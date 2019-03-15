library(dplyr)
library(foreach)
library(ggplot2)
library(tidyr)

# colors
c.red <- "#d7191c"
c.blue <- "#2c7bb6"

n <- 1000000
thetas <- seq(0, 2, 0.05)
alpha <- 0.05
z.alpha.2 <- qnorm(alpha / 2, lower.tail = FALSE)

# Naive assessment of directional claims ======================================

neg.prob <- pnorm(-z.alpha.2, thetas)
pos.prob <- pnorm(z.alpha.2, thetas, lower.tail = FALSE)
naive.no.rep <- 1 - (neg.prob^2 + pos.prob^2) / (neg.prob + pos.prob)
fdp.sim <- neg.prob / (neg.prob + pos.prob)
df <-
  data.frame(theta = thetas, naive.no.rep = naive.no.rep, fdp = fdp.sim) %>%
  gather(key = "key", value = "fraction", -theta)

pdf("fig/naive-same-dir.pdf", width = 4.5, height = 3)
ggplot(df) +
  geom_hline(yintercept = 0, size = 0.2, color = "gray50") +
  geom_vline(xintercept = 0, size = 0.2, color = "gray50") +
  geom_line(aes(x = theta, y = fraction, group = key, color = key)) +
  scale_color_manual(
    values = c(c.blue, c.red),
    labels = c("original claim is false", "original claim not confirmed")
  ) +
  geom_hline(yintercept = 0.64, lty = 2, color = c.red) +
  annotate(
    "text", color = c.red, x = 0, y = 0.64, hjust = -0.1, vjust = -0.5,
    label = "RP:P rate = 64%"
  ) +
  labs(x = expression(theta), y = "Rate (lower is better)", parse = TRUE) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_y_continuous(labels = scales::percent) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
dev.off()

# Simulation for naive CI coverage ============================================

set.seed(20180715)

x.rand <- rnorm(n)

coverage <- foreach(theta = thetas, .combine = "c") %do% {
  x <- x.rand + theta
  x <- x[abs(x) > z.alpha.2]
  p.cov <- pnorm(x + z.alpha.2, theta) - pnorm(x - z.alpha.2, theta)
  mean(p.cov)
}
df <- data.frame(theta = thetas, coverage = coverage)

pdf("fig/naive-ci-sim.pdf", width = 4.5, height = 3)
ggplot(df) +
  geom_hline(yintercept = 0, size = 0.2, color = "gray50") +
  geom_vline(xintercept = 0, size = 0.2, color = "gray50") +
  geom_hline(yintercept = 0.53, lty = 2, color = c.red) +
  annotate(
    "text", label = "RP:P non-coverage = 53%", color = c.red, x = 2, y = 0.53,
    hjust = 1, vjust = -0.5
  ) +
  geom_line(aes(x = theta, y = 1 - coverage), color = c.red) +
  labs(
    title =
      "Non-coverage of original point estimate by
      replication 95% C.I.",
    x = expression(theta), y = "Non-coverage (lower is better)"
  ) +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal()
dev.off()

# Effect decline simulation ===================================================

set.seed(20181018)

x1.rand <- rnorm(n)
x2.rand <- rnorm(n)

decline.factor <- foreach(theta = thetas, .combine = "c") %do% {
  x1 <- x1.rand + theta
  x2 <- x2.rand + theta
  sig <- abs(x1) > z.alpha.2
  x1 <- x1[sig]
  x2 <- x2[sig]
  
  # decline in effect size
  # since all original effect size > 0, if x1 < 0 a decline means x2 increases
  mean(x2 / x1 < 1)
}
df <- data.frame(theta = thetas, decline = decline.factor)

pdf("fig/naive-decline.pdf", width = 4.5, height = 3)
ggplot(df) +
  geom_hline(yintercept = 1, size = 0.2, color = "gray50") +
  geom_vline(xintercept = 0, size = 0.2, color = "gray50") +
  geom_hline(yintercept = 0.828, lty = 2, color = c.red) +
  annotate(
    "text", label = "RP:P fraction = 83%", color = c.red, x = 0, y = 0.828,
    hjust = -0.1, vjust = -0.5
  ) +
  geom_line(aes(x = theta, y = decline), color = c.red) +
  labs(
    title = "Fraction of effect sizes declined",
    x = expression(theta), y = "Fraction (lower is better)"
  ) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal()
dev.off()

# Simulation of certain values of theta =======================================
n <- 1000

# for naive same direction significant claims, look at theta = 1
set.seed(20181213)
theta <- 1
theta.O.hat <- rnorm(n, theta)
theta.R.hat <- rnorm(n, theta)
df <-
  data.frame(theta.O.hat = theta.O.hat, theta.R.hat = theta.R.hat) %>%
  mutate(selected = abs(theta.O.hat) > z.alpha.2)

pdf("fig/naive-same-dir-theta1.pdf", width = 3, height = 3)
ggplot(df) +
  annotate(
    "rect", xmin = -Inf, xmax = -z.alpha.2, ymin = -z.alpha.2, ymax = Inf,
    fill = c.red, alpha = 0.25
  ) +
  annotate(
    "rect", xmin = z.alpha.2, xmax = Inf, ymin = -Inf, ymax = z.alpha.2,
    fill = c.red, alpha = 0.25
  ) +
  annotate(
    "rect", xmin = -Inf, xmax = -z.alpha.2, ymin = -Inf, ymax = Inf,
    fill = c.blue, alpha = 0.25
  ) +
  annotate(
    "rect", xmin = -z.alpha.2, xmax = z.alpha.2, ymin = -Inf, ymax = Inf,
    fill = "gray50", alpha = 0.1
  ) +
  geom_vline(xintercept = c(-z.alpha.2, z.alpha.2), lty = 1, color = "black") +
  geom_hline(yintercept = c(-z.alpha.2, z.alpha.2), lty = 2, color = "gray50") +
  geom_vline(xintercept = 0, lty = 1, size = .2, color = "gray50") +
  geom_hline(yintercept = 0, lty = 1, size = .2, color = "gray50") +
  geom_point(aes(x = theta.O.hat, y = theta.R.hat, alpha = selected)) +
  scale_alpha_manual(values = c(0.05, 1)) +
  labs(
    x = "original effect size estimate", y = "replication effect size estimate"
  ) +
  annotate("point", x = theta, y = theta, size = 6, shape = "+") +
  annotate("segment", x = -6, xend = theta, y = theta, yend = theta, 
           lty = 2, col = "gray50"
  ) +
  annotate("segment", x = theta, xend = theta, y = -6, yend = theta, 
           lty = 2, col = "gray50"
  ) +
  coord_fixed(xlim = c(-4.5, 4.5), ylim = c(-4.5, 4.5)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  scale_x_continuous(
    breaks = c(-z.alpha.2, 0, theta, z.alpha.2),
    labels =
      c(
        expression(-z[alpha / 2]), "0", expression(theta),
        expression(z[alpha / 2])
      )
  ) +
  scale_y_continuous(
    breaks = c(-z.alpha.2, 0, theta, z.alpha.2),
    labels =
      c(
        expression(-z[alpha / 2]), "0", expression(theta),
        expression(z[alpha / 2])
      )
  )
dev.off()

# for naive confidence interval coverage, look at theta = 0.3
set.seed(20181213)
theta <- 0.5
theta.O.hat <- rnorm(n, theta)
theta.R.hat <- rnorm(n, theta)
df <-
  data.frame(theta.O.hat = theta.O.hat, theta.R.hat = theta.R.hat) %>%
  mutate(selected = abs(theta.O.hat) > z.alpha.2)

pdf("fig/naive-ci-theta05.pdf", width = 3, height = 3)
ggplot(df) +
  annotate(
    "rect", xmin = -z.alpha.2, xmax = z.alpha.2, ymin = -Inf, ymax = Inf,
    fill = "gray50", alpha = 0.1
  ) +
  annotate(
    "polygon", fill = c.red, alpha = 0.25, 
    x = c(z.alpha.2, z.alpha.2, 5 + z.alpha.2, 5),
    y = c(-5, 0, 5, -5)
  ) +
  annotate(
    "polygon", fill = c.red, alpha = 0.25, 
    x = c(-z.alpha.2, -z.alpha.2, -5 - z.alpha.2, -5),
    y = c(5, 0, -5, 5)
  ) +
  annotate(
    "polygon", fill = c.red, alpha = 0.25, 
    x = c(-z.alpha.2, -z.alpha.2, -5 + z.alpha.2, -5),
    y = c(-5, -2 * z.alpha.2, -5, -5)
  ) +
  annotate(
    "polygon", fill = c.red, alpha = 0.25, 
    x = c(z.alpha.2, z.alpha.2, 5 - z.alpha.2, 5),
    y = c(5, 2 * z.alpha.2, 5, 5)
  ) +
  geom_vline(xintercept = c(-z.alpha.2, z.alpha.2), lty = 1, color = "black") +
  geom_vline(xintercept = 0, lty = 1, size = 0.2, color = "gray50") +
  geom_hline(yintercept = 0, lty = 1, size = 0.2, color = "gray50") +
  annotate(
    "segment", x = -6, xend = theta, y = theta, yend = theta,
    lty = 2, col = "gray50"
  ) +
  annotate(
    "segment", x = theta, xend = theta, y = -6, yend = theta,
    lty = 2, col = "gray50"
  ) +
  geom_abline(slope = 1, intercept = 0, lty = 1, size = 0.2, color = "gray50") +
  geom_abline(slope = 1, intercept = -z.alpha.2, lty = 2, color = "gray50") +
  geom_abline(slope = 1, intercept = z.alpha.2, lty = 2, color = "gray50") +
  geom_point(aes(x = theta.O.hat, y = theta.R.hat, alpha = selected)) +
  scale_alpha_manual(values = c(0.05, 1)) +
  labs(
    x = "original effect size estimate", y = "replication effect size estimate"
  ) +
  annotate("point", x = theta, y = theta, size = 6, shape = "+") +
  coord_fixed(xlim = c(-4.5, 4.5), ylim = c(-4.5, 4.5)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    legend.position = "none"
  ) +
  scale_x_continuous(
    breaks = c(-z.alpha.2, 0, theta, z.alpha.2),
    labels =
      c(
        expression(-z[alpha / 2]), "0", expression(theta),
        expression(z[alpha / 2])
      )
  ) +
  scale_y_continuous(
    breaks = c(-z.alpha.2, 0, theta, z.alpha.2),
    labels =
      c(
        expression(-z[alpha / 2]), "0", expression(theta),
        expression(z[alpha / 2])
      )
  )
dev.off()

# for naive same direction significant claims, look at theta = 1
set.seed(20181213)
theta <- 1
theta.O.hat <- rnorm(n, theta)
theta.R.hat <- rnorm(n, theta)
df <-
  data.frame(theta.O.hat = theta.O.hat, theta.R.hat = theta.R.hat) %>%
  mutate(selected = abs(theta.O.hat) > z.alpha.2)

pdf("fig/naive-decline-theta1.pdf", width = 3, height = 3)
ggplot(df) +
  annotate(
    "rect", xmin = -z.alpha.2, xmax = z.alpha.2, ymin = -Inf, ymax = Inf,
    fill = "gray50", alpha = 0.1
  ) +
  annotate(
    "polygon", fill = c.red, alpha = 0.25, 
    x = c(z.alpha.2, z.alpha.2, 5, 5),
    y = c(-5, z.alpha.2, 5, -5)
  ) +
  annotate(
    "polygon", fill = c.red, alpha = 0.25, 
    x = c(-z.alpha.2, -z.alpha.2, -5, -5),
    y = c(5, -z.alpha.2, -5, 5)
  ) +
  geom_vline( xintercept = c(-z.alpha.2, z.alpha.2), lty = 1, color = "black") +
  geom_vline(xintercept = 0, lty = 1, size = 0.2, color = "gray50") +
  geom_hline(yintercept = 0, lty = 1, size = 0.2, color = "gray50") +
  annotate(
    "segment", x = -6, xend = theta, y = theta, yend = theta, 
    lty = 2, col = "gray50"
  ) +
  annotate(
    "segment", x = theta, xend = theta, y = -6, yend = theta,
    lty = 2, col = "gray50"
  ) +
  geom_abline(slope = 1, intercept = 0, lty = 2, color = "gray50") +
  geom_point(aes(x = theta.O.hat, y = theta.R.hat, alpha = selected)) +
  scale_alpha_manual(values = c(0.05, 1)) +
  labs(
    x = "original effect size estimate", y = "replication effect size estimate"
  ) +
  annotate("point", x = theta, y = theta, size = 6, shape = "+") +
  coord_fixed(xlim = c(-4.5, 4.5), ylim = c(-4.5, 4.5)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    legend.position = "none"
  ) +
  scale_x_continuous(
    breaks = c(-z.alpha.2, 0, theta, z.alpha.2),
    labels =
      c(
        expression(-z[alpha / 2]), "0", expression(theta),
        expression(z[alpha / 2])
      )
  ) +
  scale_y_continuous(
    breaks = c(-z.alpha.2, 0, theta, z.alpha.2),
    labels =
      c(
        expression(-z[alpha / 2]), "0", expression(theta),
        expression(z[alpha / 2])
      )
  )
dev.off()
