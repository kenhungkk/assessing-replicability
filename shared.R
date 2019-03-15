library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(foreach)
library(doParallel)
registerDoParallel()

# constants ===================================================================
# level used in journal
alpha <- 0.05
# level used in our analysis for ucb / test
level <- 0.05
# fdp controlled
fdp.cont <- 0.1
# colors
c.red <- "#d7191c"
c.blue <- "#2c7bb6"

# functions ===================================================================
# plotting with common legend
grid.arrange.shared.legend <-
  function(...,
           ncol = length(list(...)),
           nrow = 1,
           position = c("bottom", "right")) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )
    
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
  }

# computes the one tailed naive pval
one.tail.pval <- function(type, stat, deg, n) {
  pval <- NA
  pval <- ifelse(type == "F", pf(stat, 1, deg, lower.tail = FALSE), pval)
  pval <- ifelse(type == "t", pt(stat, deg, lower.tail = FALSE), pval)
  pval <- ifelse(type == "z", pnorm(stat, lower.tail = FALSE), pval)
  pval <- ifelse(type == "r",
                 suppressWarnings(
                   # Fisher transformation
                   pnorm(atanh(stat) * sqrt(n - 3), lower.tail = FALSE)),
                 pval)
  pval
}

# pnorm for truncated normal distribution, lower.tail = TRUE
trunc.pnorm <- function(q, mu = 0, sd = 1, lo, hi) {
  A <- pnorm(lo, mu, sd) + pnorm(hi, mu, sd, lower.tail = FALSE)
  ifelse(
    q < lo, pnorm(q, mu, sd) / A,
    1 - pnorm(pmax(q, hi), mu, sd, lower.tail = FALSE) / A
  )
}

# convert t-tests for pcor to z via Fisher transformation
# eff.n should be n - p, where p is the number of variables controlled
pcor.t2z <- function(q, eff.n) {
  r <- q / sqrt(eff.n - 2 + q^2)
  atanh(r) * sqrt(eff.n - 3)
}

# Naive test as seen in paper
# z1 ~ N(mu1, 1)
# z2 ~ N(mu2, 1)
# H0: eta1 * mu1 + eta2 * mu2 = delta
# returns p-value by default
naive.test <- function(z1, z2, eta1, eta2, delta = 0) {
  d <- eta1 * z1 + eta2 * z2
  pnorm(d, delta, sqrt(eta1^2 + eta2^2))
}

# Selective z-test as derived in the paper
# z1 ~ N(mu1, 1) 1{z1 < lo or z1 > hi}
# z2 ~ N(mu2, 1)
# H0: eta1 * mu1 + eta2 * mu2 <= delta
# returns one sided p-value by default
sel.test <- function(z1, z2, eta1, eta2, lo, hi, delta = 0) {
  d <- eta1 * z1 + eta2 * z2
  m <- eta2 * z1 - eta1 * z2
  eta.sq <- eta1^2 + eta2^2
  d.lo <- (lo * eta.sq - eta2 * m) / eta1
  d.hi <- (hi * eta.sq - eta2 * m) / eta1
  trunc.pnorm(d, delta, sqrt(eta.sq), d.lo, d.hi)
}

# Selective z-test with finite sample correction for t
# t1 ~ t,df1(mu1)
# t2 ~ t,df2(mu2)
# H0: eta1 * mu1 + eta2 * mu2 <= delta
# returns one sided-pvalue by default
sel.fs.test <-
  function(t1, t2, eta1, eta2, lo, hi, df1 = Inf, df2 = Inf, delta = 0) {
    # NA stands for z-test, taking df1 and df2 as Inf
    df1[is.na(df1)] <- Inf
    df2[is.na(df2)] <- Inf
    
    # plugin estimate for true value of mu1 and mu2
    mu2.hat <- t2
    mu1.hat <- (delta - eta2 * mu2.hat) / eta1
    var1 <- 1 + 2 * mu1.hat^2 / df1
    var2 <- 1 + 2 * mu2.hat^2 / df2
    sel.test(
      t1 / sqrt(var1), t2 / sqrt(var2), sqrt(var1) * eta1, sqrt(var2) * eta2,
      lo / sqrt(var1), hi / sqrt(var1), delta
    )
  }

# Predictive interval computation
# selective predictive interval for z2
sel.pi <- function(z1, k1, k2, lo, hi, alpha) {
  z2.hi <- 8 * k2
  z2.lo <- -z2.hi
  eta1 <- 1 / k1
  eta2 <- -1 / k2
  for (i in 1:round(log2(max(z2.hi)) - log2(.Machine$double.eps) / 2)) {
    z2.mi <- (z2.hi + z2.lo) / 2
    pval <- sel.test(z1, z2.mi, eta1, eta2, lo, hi)
    z2.lo <- ifelse(pval < alpha, z2.lo, z2.mi)
    z2.hi <- ifelse(pval < alpha, z2.mi, z2.hi)
  }
  (z2.hi + z2.lo) / 2
}

# Naive predictive interval computation
# naive predictive interval for z2
naive.pi <- function(z1, k1, k2, alpha) {
  eta1 <- 1 / k1
  eta2 <- -1 / k2
  eta.sq <- eta1^2 + eta2^2
  (sqrt(eta.sq) * qnorm(alpha) - eta1 * z1) / eta2
}

# selective confidence interval for delta = theta1 - theta2
sel.ci <- function(z1, z2, k1, k2, lo, hi, alpha) {
  delta.hi <- 8
  delta.lo <- -8.5
  eta1 <- 1 / k1
  eta2 <- -1 / k2
  for (i in 1:round(log2(max(delta.hi)) - log2(.Machine$double.eps) / 2)) {
    delta.mi <- (delta.hi + delta.lo) / 2
    pval <- sel.test(z1, z2, eta1, eta2, lo, hi, delta.mi)
    delta.lo <- ifelse(pval < alpha, delta.lo, delta.mi)
    delta.hi <- ifelse(pval < alpha, delta.mi, delta.hi)
  }
  (delta.hi + delta.lo) / 2
}

# naive confidence interval for delta = theta1 - theta2
naive.ci <- function(z1, z2, k1, k2, alpha) {
  eta1 <- 1 / k1
  eta2 <- -1 / k2
  eta.sq <- eta1^2 + eta2^2
  eta1 * z1 + eta2 * z2 - qnorm(alpha) * sqrt(eta.sq)
}

# load all data ============================================================
# load main data
df.original <- read.csv("rpp_data.csv")
df <- df.original %>%
  dplyr::filter(
    # restrict scope to univariate statistically significant studies
    T_pval_USE..O. < alpha, !is.na(T_pval_USE..R.),
    T_Test.Statistic..O. %in% c("F", "t", "z", "r"), T_df1..O. %in% c(NA, 1),
    T_Test.Statistic..R. %in% c("F", "t", "z", "r"), T_df1..R. %in% c(NA, 1)
  ) %>%
  dplyr::select(
    Study.Num,
    T_Test.Statistic..O., T_df2..O., T_Test.value..O., T_N..O., T_r..O.,
    T_pval_USE..O.,
    T_Test.Statistic..R., T_df2..R., T_Test.value..R., T_N..R., T_r..R.,
    T_pval_USE..R., Direction..R.
  ) %>%
  dplyr::rename(
    study.num = Study.Num,
    type.O = T_Test.Statistic..O., deg.O = T_df2..O., stat.O = T_Test.value..O.,
    n.O = T_N..O., eff.O = T_r..O., pval.O = T_pval_USE..O.,
    type.R = T_Test.Statistic..R., deg.R = T_df2..R., stat.R = T_Test.value..R.,
    n.R = T_N..R., eff.R = T_r..R., pval.R = T_pval_USE..R.,
    dir = Direction..R.
  ) %>%
  dplyr::mutate(
    # turn all F into t
    stat.O = ifelse(type.O == "F", suppressWarnings(sqrt(stat.O)), stat.O),
    type.O = ifelse(type.O == "F", "t", as.character(type.O)),
    stat.O = abs(stat.O),
    stat.R = ifelse(type.R == "F", suppressWarnings(sqrt(stat.R)), stat.R),
    type.R = ifelse(type.R == "F", "t", as.character(type.R)),
    stat.R = ifelse(dir == "same", abs(stat.R), -abs(stat.R)),
    # recompute the number of tails
    tails.O = round(pval.O / one.tail.pval(type.O, stat.O, deg.O, n.O)),
    tails.R = round(pval.R / one.tail.pval(type.R, abs(stat.R), deg.R, n.R))
  ) %>%
  # pval adjustments for computing fdp estimates and ucbs
  dplyr::mutate(
    # adjust original pval
    adjp.O = pval.O / alpha,
    # make replication pval all one-tailed
    adjp.R =
      ifelse(
        tails.R == 1, pval.R, ifelse(dir == "same", pval.R / 2, 1 - pval.R / 2)
      )
  ) %>%
  # convert everything to z score
  dplyr::mutate(
    # convert everything to t or z
    z.O =
      ifelse(
        type.O != "r", stat.O, suppressWarnings(atanh(stat.O) * sqrt(n.O - 3))
      ),
    type.O = ifelse(type.O == "r", "z", as.character(type.O)),
    z.R =
      ifelse(
        type.R != "r", stat.R, suppressWarnings(atanh(stat.R) * sqrt(n.R - 3))
      ),
    type.R = ifelse(type.R == "r", "z", as.character(type.R)),
    z.O.lo =
      ifelse(
        type.O == "t",
        # if the test is a t-test
        ifelse(tails.O == 1, -Inf, qt(alpha / 2, df = deg.O)),
        # if the test is a z-test (including correlations)
        ifelse(tails.O == 1, -Inf, qnorm(alpha / 2))
      ),
    z.O.hi =
      ifelse(
        type.O == "t",
        # if the test is a t-test
        qt(1 - alpha / tails.O, df = deg.O),
        # if the test is a z-test (including correlations)
        qnorm(1 - alpha / tails.O)
      )
  )

# special cases, all regressions and correlations
df[df$study.num %in% c(44, 48, 93, 112, 134), c("type.O", "type.R")] <- "z"
# Study 44: Payne et al.
df[df$study.num == 44, c("z.O", "z.O.lo", "z.O.hi")] <-
  pcor.t2z(df[df$study.num == 44, c("z.O", "z.O.lo", "z.O.hi")], 71 - 2)
df[df$study.num == 44, "z.R"] <-
  pcor.t2z(df[df$study.num == 44, "z.R"], 180 - 2)
# Study 48: Cox et al.
df[df$study.num == 48, c("z.O", "z.O.lo", "z.O.hi")] <-
  pcor.t2z(df[df$study.num == 48, c("z.O", "z.O.lo", "z.O.hi")], 100 - 6)
df[df$study.num == 48, "z.R"] <-
  pcor.t2z(df[df$study.num == 48, "z.R"], 200 - 6)
# Study 93: Murray et al.
df[df$study.num == 93, c("z.O", "z.O.lo", "z.O.hi")] <-
  pcor.t2z(df[df$study.num == 93, c("z.O", "z.O.lo", "z.O.hi")], 91 - 6)
df[df$study.num == 93, "z.R"] <-
  pcor.t2z(df[df$study.num == 93, "z.R"], 76 - 6)
# Study 112: McKinstry et al.
df[df$study.num == 112, c("z.O", "z.O.lo", "z.O.hi")] <-
  pcor.t2z(df[df$study.num == 112, c("z.O", "z.O.lo", "z.O.hi")], 11)
df[df$study.num == 112, "z.R"] <-
  pcor.t2z(df[df$study.num == 112, "z.R"], 11)
# Study 134: Larsen and McKibban
df[df$study.num == 134, c("z.O", "z.O.lo", "z.O.hi")] <-
  pcor.t2z(df[df$study.num == 134, c("z.O", "z.O.lo", "z.O.hi")], 119 - 1)
df[df$study.num == 134, "z.R"] <-
  pcor.t2z(df[df$study.num == 134, "z.R"], 238 - 1)

# hardcoding k.O and k.R according to the appendix
source("effect-size-k.R")
df <- merge(df, df.k, by = "study.num")
