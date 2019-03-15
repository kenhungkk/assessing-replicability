# load data and functions if df is not already there
if (!exists("df.original") | !exists("df")) {
  source("shared.R")
}

df <- df %>%
  dplyr::filter(
    deg.O >= 30 | is.na(deg.O),
    deg.R >= 30 | is.na(deg.R)
  )
# Row 27 is an outlier for deg.R
range.degs.O <- range(df$deg.O[-27], na.rm = TRUE)
range.degs.R <- range(df$deg.R[-27], na.rm = TRUE)
degs.O <-
  round(exp(seq(log(range.degs.O[1]), log(range.degs.O[2]), length.out = 5)))
degs.R <-
  round(exp(seq(log(range.degs.R[1]), log(range.degs.R[2]), length.out = 5)))

# Show distribution of dfs ====================================================
pdf("fig/df.pdf", width = 4, height = 6)
ggplot(df) +
  geom_point(aes(x = deg.O, y = deg.R)) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  geom_point(
    data = expand.grid(x = degs.O, y = degs.R), aes(x = x, y = y),
    color = c.blue, shape = "+", size = 5
  ) +
  annotate(
    "polygon", fill = c.blue, alpha = 0.25,
    x =
      c(
        degs.O[1], degs.O[2], degs.O[3], degs.O[5], degs.O[5], degs.O[4],
        degs.O[3], degs.O[1], degs.O[1]
      ),
    y =
      c(
        degs.R[3], degs.R[4], degs.R[5], degs.R[5], degs.R[3], degs.R[2],
        degs.R[1], degs.R[1], degs.R[3]
      )
  ) +
  labs(x = expression(df[O]), y = expression(df[R])) +
  theme_minimal()
dev.off()

# Simulation ==================================================================

# one-sample t-test simulation parameters
ncps.O <- seq(0, 5, 0.1)
# size of simulation
m <- 1e5

nom.levels <- c(level, level / 46 * 4, level / 46)
filenames <- c("t-approx", "t-approx-bh", "t-approx-bonf")

for (i in 1:3) {
  nom.level <- nom.levels[i]
  filename <- filenames[i]
  
  # parallelized code for power simulation
  powers <-
    foreach(df1 = degs.O, .combine = "rbind") %:%
    foreach(df2 = degs.R, .combine = "rbind") %:%
    foreach(ncp1 = ncps.O, .combine = "rbind") %dopar% {
      cutoff <- qt(alpha / 2, df1, lower.tail = FALSE)
      # lower tail and upper tail probabilities
      lt <- pt(-cutoff, df1, ncp1)
      ut <- pt(cutoff, df1, ncp1, lower.tail = FALSE)
      p.t1 <- runif(m, max = lt + ut)
      p.t1 <- p.t1 + (p.t1 > lt) * (1 - lt - ut)
      # inverse CDF to generate t1
      t1 <- qt(p.t1, df1, ncp1)
      # ncp2 under true null of equal effect size
      ncp2 <- ncp1 / sqrt(df1 + 1) * sqrt(df2 + 1)
      t2 <- rt(m, df2, ncp2)
      # simple p-value
      pval <- sel.test(t1, t2, 1 / sqrt(df1), -1 / sqrt(df2), -cutoff, cutoff)
      pval.2s <- 2 * pmin(pval, 1 - pval)
      # p-value with finite sample correction
      pval.fs <-
        sel.fs.test(
          t1, t2, 1 / sqrt(df1), -1 / sqrt(df2), -cutoff, cutoff, df1, df2
        )
      pval.fs.2s <- 2 * pmin(pval.fs, 1 - pval.fs)
      c(df1, df2, ncp1, mean(pval.2s < nom.level), mean(pval.fs.2s < nom.level))
    }
  powers <- data.frame(powers)
  colnames(powers) <- c("deg.O", "deg.R", "ncp.O", "simple", "fs")
  powers <-
    powers %>%
    gather(key = "test", value = "power", -deg.O, -deg.R, -ncp.O) %>%
    mutate(
      df1 =
        factor(
          deg.O, labels = sapply(degs.O, function(x) paste0("df[O]=='", x, "'"))
        ),
      df2 =
        factor(
          deg.R, labels = sapply(degs.R, function(x) paste0("df[R]=='", x, "'"))
        ),
      df2 = factor(df2, levels = rev(levels(df2))),
      # facets to be greyed out
      useless =
        case_when(
          deg.O == degs.O[1] ~ deg.R %in% degs.R[c(4, 5)],
          deg.O == degs.O[2] ~ deg.R == degs.R[5],
          deg.O == degs.O[4] ~ deg.R == degs.R[1],
          deg.O == degs.O[5] ~ deg.R %in% degs.R[c(1, 2)],
          TRUE ~ FALSE
        )
    )
  
  pdf(paste0("fig/", filename, ".pdf"), width = 6, height = 6)
  ggplot(powers) +
    facet_grid(df2 ~ df1, labeller = label_parsed) +
    geom_rect(
      data = subset(powers, useless), alpha = 0.0025, fill = "black", color = NA,
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    geom_line(aes(x = ncp.O, y = power, group = test, color = test)) +
    scale_color_manual(
      values = c(c.blue, c.red),
      labels =
        c(
          "Finite-sample correction",
          expression("Standard selective"~italic(z)*"-test")
        )
    ) +
    geom_hline(yintercept = nom.level, lty = 2, color = "gray50") +
    labs(
      x = expression("Original noncentrality parameter ("*ncp[O]*")"),
      y = "Type I error rate", color = "Test"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  dev.off()
}

# Special case: Row 27 or Purdie-Vaughns ======================================
df1 <- df$deg.O[27]
df2 <- df$deg.R[27]

# parallelized code for power simulation
powers <-
  foreach(nom.level = nom.levels, .combine = "rbind") %:%
  foreach(ncp1 = ncps.O, .combine = "rbind") %dopar% {
    cutoff <- qt(alpha / 2, df1, lower.tail = FALSE)
    # lower tail and upper tail probabilities
    lt <- pt(-cutoff, df1, ncp1)
    ut <- pt(cutoff, df1, ncp1, lower.tail = FALSE)
    p.t1 <- runif(m, max = lt + ut)
    p.t1 <- p.t1 + (p.t1 > lt) * (1 - lt - ut)
    # inverse CDF to generate t1
    t1 <- qt(p.t1, df1, ncp1)
    # ncp2 under true null of equal effect size
    ncp2 <- ncp1 / sqrt(df1 + 1) * sqrt(df2 + 1)
    t2 <- rt(m, df2, ncp2)
    # simple p-value
    pval <- sel.test(t1, t2, 1 / sqrt(df1), -1 / sqrt(df2), -cutoff, cutoff)
    pval.2s <- 2 * pmin(pval, 1 - pval)
    # p-value with finite sample correction
    pval.fs <-
      sel.fs.test(
        t1, t2, 1 / sqrt(df1), -1 / sqrt(df2), -cutoff, cutoff, df1, df2
      )
    pval.fs.2s <- 2 * pmin(pval.fs, 1 - pval.fs)
    c(nom.level, ncp1, mean(pval.2s < nom.level), mean(pval.fs.2s < nom.level))
  }
powers <- data.frame(powers)
colnames(powers) <- c("nom.level", "ncp.O", "simple", "fs")
powers <-
  powers %>%
  gather(key = "test", value = "power", -nom.level, -ncp.O) %>%
  mutate(
    fac.level = factor(nom.level, labels = c("Bonferroni", "BH", "Nominal")),
    fac.level = factor(fac.level, levels = rev(levels(fac.level)))
  )

pdf(paste0("fig/t-approx-pv.pdf"), width = 9, height = 4)
ggplot(powers) +
  facet_grid(cols = vars(fac.level)) +
  geom_line(aes(x = ncp.O, y = power, group = test, color = test)) +
  scale_color_manual(
    values = c(c.blue, c.red),
    labels =
      c(
        "Finite-sample correction",
        expression("Standard selective"~italic(z)*"-test")
      )
  ) +
  geom_hline(aes(yintercept = nom.level), lty = 2, color = "gray50") +
  labs(
    x = expression("Original noncentrality parameter ("*ncp[O]*")"),
    y = "Type I error rate", color = "Test"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
dev.off()

# Effect shift reanaylsis =====================================================
df <- df %>%
  dplyr::filter(
    deg.O >= 30 | is.na(deg.O),
    deg.R >= 30 | is.na(deg.R)
  ) %>%
  dplyr::mutate(
    sel.fs.p =
      sel.fs.test(z.O, z.R, 1 / k.O, -1 / k.R, z.O.lo, z.O.hi, deg.O, deg.R),
    sel.fs.p = 2 * pmin(sel.fs.p, 1 - sel.fs.p),
    sel.fs.rej = sel.fs.p < alpha
  )

df$bh.fs.p <- p.adjust(df$sel.fs.p, method = "BH")
df$bh.fs.rej <- df$bh.fs.p < fdp.cont
df$holm.fs.p <- p.adjust(df$sel.fs.p, method = "holm")
df$holm.fs.rej <- df$holm.fs.p < 0.05

cat(
  "======== Rejections (with finite sample correction) ========",
  "\nSelective test:", sum(df$sel.fs.rej), "/", nrow(df), "=",
  mean(df$sel.fs.rej),
  "\nSelective test, adjusted with BH:",
  sum(df$bh.fs.rej), "/", nrow(df), "=", mean(df$bh.fs.rej),
  "\nSelective test, adjusted with Holm:",
  sum(df$holm.fs.rej), "/", nrow(df), "=", mean(df$holm.fs.rej),
  "\n"
)

# Effect decline reanalysis ===================================================
lambda <- 0.5
rho <- 0.2

# Plot the p-value
pval <-
  1 -
  with(
    df,
    sel.fs.test(
      z.O, z.R, 1 / k.O, -1 / (k.R * (1 - rho)), -Inf, z.O.hi, deg.O, deg.R
    )
  )

# Underestimate and lower confidence bound
R <- nrow(df)
B <- sum(pval > lambda)
V.O <- B / (1 - lambda)
V.hi.O <- sum(pbinom(B, 1:R, 1 - lambda) >= level)
cat(
  "Underestimate of effect size decline: ", 1 - V.O / R,
  " (LCB: ", 1 - V.hi.O / R, ")",
  sep = ""
)

# Overestimate and upper confidence bound
pval <- 1 - pval
R <- nrow(df)
B <- sum(pval > lambda)
V.O <- B / (1 - lambda)
V.hi.O <- sum(pbinom(B, 1:R, 1 - lambda) >= level)
cat(
  "Overestimate of effect size decline: ", min(V.O / R, 1),
  " (UCB: ", min(V.hi.O / R, 1), ")",
  sep = ""
)
