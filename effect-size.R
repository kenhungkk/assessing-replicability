# load data and functions if df is not already there
if (!exists("df.original") | !exists("df")) {
  source("shared.R")
}

# Testing (naive, selective, BH) ==============================================
df <- df %>%
  dplyr::filter(
    deg.O >= 30 | is.na(deg.O),
    deg.R >= 30 | is.na(deg.R)
  ) %>%
  dplyr::mutate(
    naive.p = naive.test(z.O, z.R, 1 / k.O, -1 / k.R),
    naive.p = 2 * pmin(naive.p, 1 - naive.p),
    # sel.1s.p is pval for H0: theta.O > theta.R
    sel.p = sel.test(z.O, z.R, 1 / k.O, -1 / k.R, z.O.lo, z.O.hi),
    # make it two sided
    sel.p = 2 * pmin(sel.p, 1 - sel.p),
    naive.rej = naive.p < alpha,
    sel.rej = sel.p < alpha
  )
df$bh.p <- p.adjust(df$sel.p, method = "BH")
df$bh.rej <- df$bh.p < fdp.cont
df$holm.p <- p.adjust(df$sel.p, method = "holm")
df$holm.rej <- df$holm.p < 0.05
cat(
  "======== Rejections ========",
  "\nNaive test:", sum(df$naive.rej), "/", nrow(df), "=", mean(df$naive.rej),
  "\nSelective test:", sum(df$sel.rej), "/", nrow(df), "=", mean(df$sel.rej),
  "\nSelective test, adjusted with BH:",
  sum(df$bh.rej), "/", nrow(df), "=", mean(df$bh.rej),
  "\nSelective test, adjusted with Holm:",
  sum(df$holm.rej), "/", nrow(df), "=", mean(df$holm.rej),
  "\n"
)

# Predictive interval for theta.hat.R ==========================================
df <- df %>%
  dplyr::mutate(
    z.R.lo = sel.pi(z.O, k.O, k.R, z.O.lo, z.O.hi, 1 - alpha / 2),
    z.R.hi = sel.pi(z.O, k.O, k.R, z.O.lo, z.O.hi, alpha / 2),
    z.R.lo.naive = naive.pi(z.O, k.O, k.R, 1 - alpha / 2),
    z.R.hi.naive = naive.pi(z.O, k.O, k.R, alpha / 2),
    eff.R.lo = z.R.lo / k.R,
    eff.R.hi = z.R.hi / k.R,
    eff.R.lo.naive = z.R.lo.naive / k.R,
    eff.R.hi.naive = z.R.hi.naive / k.R,
    eff.O = z.O / k.O,
    eff.R = z.R / k.R
  )

# Plot effect size estimate, replication vs. original
plot.pi.adjusted <-
  ggplot(df %>% filter(eff.O < 2.5)) +
  geom_abline(slope = 1, intercept = 0, size = 0.2, color = "gray50") +
  geom_vline(xintercept = 0, size = 0.2, color = "gray50") +
  geom_hline(yintercept = 0, size = 0.2, color = "gray50") +
  geom_point(aes(x = eff.O, y = eff.R, color = sel.rej)) +
  geom_errorbar(
    aes(x = eff.O, ymin = eff.R.lo, ymax = eff.R.hi, color = sel.rej),
    width = 0.05
  ) +
  scale_colour_manual(values = c("black", c.red)) +
  geom_point(
    data = df %>% filter(sel.rej), aes(x = eff.O, y = eff.R), color = c.red
  ) +
  geom_errorbar(
    data = df %>% filter(sel.rej),
    aes(x = eff.O, ymin = eff.R.lo, ymax = eff.R.hi), color = c.red,
    width = 0.05
  ) +
  coord_cartesian(ylim = c(-2, max(df$eff.R.hi[df$eff.O < 2.5]))) +
  scale_y_continuous(minor_breaks = seq(-2, 4, 0.5)) +
  labs(
    x = expression(hat(theta)[O]), y = expression(hat(theta)[R]),
    color = "Rejected", title = "Adjusted P.I."
  ) +
  theme_minimal()
plot.pi.unadjusted <-
  ggplot(df %>% filter(eff.O < 2.5)) +
  geom_abline(slope = 1, intercept = 0, size = 0.2, color = "gray50") +
  geom_vline(xintercept = 0, size = 0.2, color = "gray50") +
  geom_hline(yintercept = 0, size = 0.2, color = "gray50") +
  geom_point(aes(x = eff.O, y = eff.R, color = naive.rej)) +
  geom_errorbar(
    aes(
      x = eff.O, ymin = eff.R.lo.naive, ymax = eff.R.hi.naive, color = naive.rej
    ),
    width = 0.05
  ) +
  scale_colour_manual(values = c("black", c.red)) +
  geom_point(
    data = df %>% filter(naive.rej), aes(x = eff.O, y = eff.R), color = c.red
  ) +
  geom_errorbar(
    data = df %>% filter(naive.rej),
    aes(x = eff.O, ymin = eff.R.lo.naive, ymax = eff.R.hi.naive), color = c.red,
    width = 0.05
  ) +
  coord_cartesian(ylim = c(-2, max(df$eff.R.hi[df$eff.O < 2.5]))) +
  scale_y_continuous(minor_breaks = seq(-2, 4, 0.5)) +
  labs(
    x = expression(hat(theta)[O]), y = expression(hat(theta)[R]),
    color = "Rejected", title = "Unadjusted P.I."
  ) +
  theme_minimal()

pdf("fig/pi.pdf", width = 7.5, height = 7.5, onefile = FALSE)
grid.arrange.shared.legend(
  plot.pi.adjusted, plot.pi.unadjusted, nrow = 2, ncol = 1
)
dev.off()

# Confidence interval for theta.O - theta.R ===================================
df <- df %>%
  dplyr::mutate(
    delta.lo = sel.ci(z.O, z.R, k.O, k.R, z.O.lo, z.O.hi, 1 - alpha / 2),
    delta.hi = sel.ci(z.O, z.R, k.O, k.R, z.O.lo, z.O.hi, alpha / 2),
    delta.lo.naive = naive.ci(z.O, z.R, k.O, k.R, 1 - alpha / 2),
    delta.hi.naive = naive.ci(z.O, z.R, k.O, k.R, alpha / 2),
    delta.hat = eff.O - eff.R
  )

# Plot diff in effect size estimate, delta, and CI for it
plot.ci.adjusted <-
  ggplot(df) +
  geom_abline(slope = 1, intercept = 0, size = 0.2, color = "gray50") +
  geom_vline(xintercept = 0, size = 0.2, color = "gray50") +
  geom_hline(yintercept = 0, size = 0.2, color = "gray50") +
  geom_point(aes(x = delta.hat, y = delta.hat, color = sel.rej)) +
  geom_errorbar(
    aes(x = delta.hat, ymin = delta.lo, ymax = delta.hi, color = sel.rej),
    width = 0.05
  ) +
  scale_colour_manual(values = c("black", c.red)) +
  geom_point(
    data = df %>% filter(sel.rej), aes(x = delta.hat, y = delta.hat),
    color = c.red
  ) +
  geom_errorbar(
    data = df %>% filter(sel.rej),
    aes(x = delta.hat, ymin = delta.lo, ymax = delta.hi), color = c.red,
    width = 0.05
  ) +
  coord_cartesian(ylim = c(-3, max(df$delta.hi))) +
  scale_y_continuous(minor_breaks = seq(-3, 4, 0.5)) +
  labs(
    x = expression(hat(theta)[O] - hat(theta)[R]),
    y = expression(theta[O] - theta[R]), color = "Rejected",
    title = "Adjusted C.I."
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
plot.ci.unadjusted <-
  ggplot(df) +
  geom_abline(slope = 1, intercept = 0, size = 0.2, color = "gray50") +
  geom_vline(xintercept = 0, size = 0.2, color = "gray50") +
  geom_hline(yintercept = 0, size = 0.2, color = "gray50") +
  geom_point(aes(x = delta.hat, y = delta.hat, color = naive.rej)) +
  geom_errorbar(
    aes(
      x = delta.hat, ymin = delta.lo.naive, ymax = delta.hi.naive,
      color = naive.rej
    ),
    width = 0.05
  ) +
  scale_colour_manual(values = c("black", c.red)) +
  geom_point(
    data = df %>% filter(naive.rej), aes(x = delta.hat, y = delta.hat),
    color = c.red
  ) +
  geom_errorbar(
    data = df %>% filter(naive.rej),
    aes(x = delta.hat, ymin = delta.lo.naive, ymax = delta.hi.naive),
    color = c.red, width = 0.05
  ) +
  coord_cartesian(ylim = c(-3, max(df$delta.hi))) +
  scale_y_continuous(minor_breaks = seq(-3, 4, 0.5)) +
  labs(
    x = expression(hat(theta)[O] - hat(theta)[R]),
    y = expression(theta[O] - theta[R]), color = "Rejected",
    title = "Unadjusted C.I."
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

pdf("fig/ci.pdf", width = 7.5, height = 7.5, onefile = FALSE)
grid.arrange.shared.legend(
  plot.ci.adjusted, plot.ci.unadjusted, nrow = 2, ncol = 1
)
dev.off()

# Effect decline analysis ======================================================
lambda <- 0.5

# Plot the p-value
pval <- 1 - with(df, sel.test(z.O, z.R, 1 / k.O, -1 / k.R, -Inf, z.O.hi))
pdf("fig/effect-decline.pdf", width = 6, height = 4)
ggplot() +
  geom_vline(xintercept = c(0, 1), col = "gray50", size = 0.2) +
  stat_bin(
    aes(x = pval), breaks = seq(0, 1, 0.05), fill = NA, col = "black"
  ) +
  labs(x = expression(italic("p")*"-value"), y = "Frequency") +
  geom_vline(xintercept = lambda, lty = 2, col = c.red) +
  annotate(
    "text", x = lambda, y = 10, parse = TRUE, col = c.red, hjust = -0.1,
    label = expression(lambda * "=0.5")
  ) +
  geom_hline(
    yintercept = sum(pval > lambda) / lambda * 0.05,
    lty = 2, col = c.red
  ) +
  scale_y_continuous(breaks = seq(0, 12, 4), minor_breaks = seq(0, 16, 2)) +
  theme_minimal()
dev.off()

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

# test the null H^D: theta.O * (1 - rho) <= theta.R
rhos <- seq(-1, 1, 0.01)
df.decline <- foreach(rho = rhos, .combine = "rbind") %do% {
  if (rho == 1) {
    # same as testing theta.R >= 0, which gives p-value complement to testing
    # theta.R <= 0
    pval <- 1 - df$adjp.R
  } else {
    # use selective Gaussian test
    # k.R changes to k.R * (1 - rho)
    pval <-
      1 -
      with(
        df, sel.test(z.O, z.R, 1 / k.O, -1 / (k.R * (1 - rho)), -Inf, z.O.hi)
      )
  }
  # compute effect decline, the lower ones first
  R <- nrow(df)
  B <- sum(pval > lambda)
  V.O <- B / (1 - lambda)
  V.hi.O <- sum(pbinom(B, 1:R, 1 - lambda) >= level)
  lo.est <- max(1 - V.O / R, 0)
  lcb <- max(1 - V.hi.O / R, 0)
  
  # compute effect decline, the upper ones now
  pval <- 1 - pval
  B <- sum(pval > lambda)
  V.O <- B / (1 - lambda)
  V.hi.O <- sum(pbinom(B, 1:R, 1 - lambda) >= level)
  hi.est <- min(V.O / R, 1)
  ucb <- min(V.hi.O / R, 1)
  
  # returns FDP estimate and lcb
  c(rho, lo.est, hi.est, lcb, ucb)
}
colnames(df.decline) <- c("rho", "underest", "overest", "lcb", "ucb")
df.decline <- data.frame(df.decline)

pdf("fig/effect-decline-range.pdf", width = 6, height = 4)
ggplot(df.decline) +
  geom_ribbon(
    aes(x = rho, ymin = lcb, ymax = ucb), fill = "black", alpha = 0.25,
    color = NA
  ) +
  geom_line(aes(x = rho, y = overest), color = "black") +
  geom_line(aes(x = rho, y = underest), color = "black") +
  geom_hline(yintercept = c(0, 1), color = "gray50", size = 0.2) +
  geom_vline(xintercept = 0, color = "gray50", size = 0.2) +
  labs(
    x = expression(paste(rho, " (Relative effect size decline)")),
    y = "Proportion"
  ) +
  theme_minimal()
dev.off()
