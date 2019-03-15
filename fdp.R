# load data and functions if df is not already there
if (!exists("df.original") | !exists("df")) {
  source("shared.R")
}

# Storey"s method using [0.5, 1.0] to estimate ================================
lambda <- 0.5
alphas <- c(0.001, 0.005, 0.01, 0.05)

# Estimates for V based on original, replication or combined p-values, also UCBs

# When alpha < 0.05, refer to method in Appendix B
B <- sum(df$adjp.O > lambda)
# a is alpha / 0.05
fdrs <- foreach(a = alphas[-length(alphas)] / alpha, .combine = "rbind") %do% {
  beta <- (1 - lambda) / (1 - lambda + a)
  R <- sum(df$adjp.O < a)
  N <- R + B
  # Estimate
  gamma.hat <- (1 / beta - 1) * B / R
  # UCB
  Q <- sum(pbinom(B, 1:N, beta) >= level)
  gamma.ucb <- (Q - B) / R
  c(gamma.hat, gamma.ucb)
}

# When alpha = 0.05 refer to method in the methodology section
R <- nrow(df)
# for original
# estimate for V
V.O <- B / (1 - lambda)
# ucb for V
V.hi.O <- sum(pbinom(B, 1:R, 1 - lambda) >= level)
# adding special row to fdrs
fdrs <- rbind(fdrs, c(V.O / R, V.hi.O / R))

# Computes the FDR based on replication and combined p-vals
more.fdrs <- foreach(a = alphas / alpha, .combine = "rbind") %do% {
  # shrink data set according to new alpha
  df.tmp <-
    df %>%
    filter(adjp.O < a) %>%
    mutate(adjp.O = adjp.O / a)
  R <- nrow(df.tmp)
  
  # for replication
  B.R <- sum(df.tmp$adjp.R > lambda)
  # estimate for V
  V.R <- B.R / (1 - lambda)
  # ucb for V
  V.hi.R <- sum(pbinom(B.R, 1:R, 1 - lambda) >= level)
  
  c(V.R / R, V.hi.R / R)
}

fdrs <- cbind(fdrs, more.fdrs)
rownames(fdrs) <- alphas
colnames(fdrs) <- c("Ori", "Ori UCB", "Rep", "Rep UCB")
# Table 2, last row gives Table 1
fdrs

# p-value distributions =======================================================
pdf("fig/fdp-original.pdf", width = 4, height = 4)
ggplot(df) +
  geom_vline(xintercept = c(0, 1), col = "gray50", size = 0.2) +
  stat_bin(
    aes(x = adjp.O), breaks = seq(0, 1, 0.05), fill = NA, col = "black"
  ) +
  labs(
    x = expression("Adjusted original" ~ italic("p") * "-values"),
    y = "Frequency"
  ) +
  geom_vline(xintercept = lambda, lty = 2, col = c.red) +
  annotate(
    "text", x = lambda, y = 20, parse = TRUE, col = c.red, hjust = -0.1,
    label = expression(lambda * "=0.5")
  ) +
  geom_hline(
    yintercept = sum(df$adjp.O > lambda) / lambda * 0.05, lty = 2, col = c.red
  ) +
  theme_minimal()
dev.off()

pdf("fig/fdp-replication.pdf", width = 4, height = 4)
ggplot(df) +
  geom_vline(xintercept = c(0, 1), col = "gray50", size = 0.2) +
  stat_bin(
    aes(x = adjp.R), breaks = seq(0, 1, 0.05), fill = NA, col = "black"
  ) +
  labs(
    x = expression("Replication" ~ italic("p") * "-values"),
    y = "Frequency"
  ) +
  geom_vline(xintercept = lambda, lty = 2, col = c.red) +
  annotate(
    "text", x = lambda, y = 20, parse = TRUE, col = c.red, hjust = -0.1,
    label = expression(lambda * "=0.5")
  ) +
  geom_hline(
    yintercept = sum(df$adjp.R > lambda) / lambda * 0.05, lty = 2, col = c.red
  ) +
  theme_minimal()
dev.off()
