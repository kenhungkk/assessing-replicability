# load data and functions if df is not already there
if (!exists("df.original") | !exists("df")) {
  source("shared.R")
}

df.pval <- na.omit(df.original[, c("T_pval_USE..O.", "T_pval_USE..R.")])
colnames(df.pval) <- c("Original", "Replication")
df.pval <- df.pval %>% gather(key = "experiment", value = "value")

pdf("fig/ecdf.pdf", width = 6, height = 4)
ggplot(df.pval) +
  geom_vline(xintercept = c(0, 1), col = "gray50", size = 0.2) +
  geom_hline(yintercept = c(0, 1), col = "gray50", size = 0.2) +
  stat_ecdf(
    aes(x = value, group = experiment, color = experiment), na.rm = TRUE
  ) +
  labs(
    x = expression(italic("p") * "-value"),
    y = "Cumulative frequency",
    color = "Experiment"
  ) +
  geom_vline(xintercept = alpha, col = c.red, lty = 2) +
  annotate(
    "text", x = 0.05, y = 0.15, parse = TRUE, col = c.red, hjust = -0.1,
    label = expression(italic("p") * "=0.05")
  ) +
  scale_color_manual(values = c(c.red, c.blue)) +
  coord_fixed() + 
  theme_minimal()
dev.off()
