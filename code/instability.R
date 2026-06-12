library(mvtnorm)
library(patchwork)
library(ggtern)
library(cowplot)

set.seed(251111)

K <- 3
mu0 <- rep(0, K)  # equal means
Sigma0 <- matrix(c(
  1, -0.99, 0,
  -0.99, 1,  0,
  0, 0, 0.01
), nrow = 3, byrow = TRUE)

n <- 500  
S <- 10000


mu_samples <- rmvnorm(S, mu0, Sigma0/n)

winners <- integer(S)
for (s in 1:S) {
  x <- mu_samples[s,]
  winners[s] <- which.min(x)
}

w <- tabulate(winners, nbins = K) / S
round(w, 3)



Bplot <- 1000
draws <- rmvnorm(Bplot, mean = mu0, sigma = Sigma0 / n)
df <- data.frame(mu1 = draws[,1], mu2 = draws[,2], mu3 = draws[,3])
df$min12 <- pmin(df$mu1, df$mu2)


lims <- range(c(df$mu1, df$mu2, df$mu3))

p1 <- ggplot(df, aes(mu1, mu2)) +
  geom_point(alpha=0.8, size = 0.001) +
  coord_cartesian(xlim = lims, ylim = lims) +
  labs(
    x = expression(mu[1]),
    y = expression(mu[2]),
    # title = bquote(mu[2]~"vs"~~mu[1])
  ) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_line(color = scales::alpha("grey70", 0.3), linewidth = 0.2),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 15),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5))
  )



p2 <- ggplot(df, aes(mu1, mu3)) +
  geom_point(alpha=0.8, size = 0.001) +
  coord_cartesian(xlim = lims, ylim = lims) +
  labs(
    x = expression(mu[1]),
    y = expression(mu[3]),
    # title = bquote(mu[3]~"vs"~~mu[1])
  ) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_line(color = scales::alpha("grey70", 0.3), linewidth = 0.2),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 15),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5))
  )

p3 <- ggplot(df, aes(min12, mu3)) +
  geom_point(alpha=0.8, size = 0.001) +
  coord_cartesian(xlim = lims, ylim = lims) +
  labs(
    x = expression(paste("min(", mu[1], ", ", mu[2], ")")),
    y = expression(mu[3]),
    # title = expression(mu[3]~"vs"~~paste("min(", mu[1], ", ", mu[2], ")"))
  ) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_line(color = scales::alpha("grey70", 0.3), linewidth = 0.2),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 15),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5))
  )

p_instability <- p1 + p2 + p3

# p_instability <- cowplot::plot_grid(
#   p1, p2, p3,
#   ncol = 3, align = "h", axis = "tb",
#   labels = c("(A)", "(B)", "(C)"),
#   label_size = 22,
#   label_fontface = "bold",
#   label_x = 0.02,
#   label_y = 0.98,
#   hjust = 0,
#   vjust = 1
# )


# save for Figure 2
# fig_dir <- here::here("output", "figures")
# ggsave(file.path(fig_dir, "instability_fig2.png"), p_instability, width = 15, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "instability_fig2.pdf"), p_instability, width = 15, height = 5, device = grDevices::cairo_pdf)






