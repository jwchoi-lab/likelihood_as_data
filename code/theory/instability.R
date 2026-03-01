library(mvtnorm)
library(patchwork)
library(ggtern)

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
p_instability

ggsave("Figures/251111_instability_example.png", p_instability, width = 15, height = 5, dpi = 300, device = "png")




df$winner <- factor(ifelse(df$mu3 < df$min12, 3, ifelse(df$mu1 < df$mu2, 1, 2)),
                    levels = c(1,2,3), labels = c("Model 1","Model 2","Model 3"))
p3_thresh <- ggplot(df, aes(min12, mu3)) +
  geom_point(aes(color = winner), alpha = 0.35, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_cartesian(xlim = lims, ylim = lims) +
  labs(
    x = expression(paste("min(", mu[1], ", ", mu[2], ")")),
    y = expression(mu[3])
  ) +
  scale_color_manual(values = c("#1f77b4","#1f77b4","#d62728")) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_line(color = scales::alpha("grey70", 0.3), linewidth = 0.2),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 15),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5)),
    legend.position = "bottom"
  )

p3_thresh








# Ternary plot

library(mvtnorm)
library(ggtern)
library(dplyr)


set.seed(251111)

# difference-covariance matrix 
diff_cov <- function(Sigma, k) {
  idx <- setdiff(seq_len(nrow(Sigma)), k)
  C <- matrix(0, nrow = length(idx), ncol = nrow(Sigma))
  for (r in seq_along(idx)) {
    e <- rep(0, nrow(Sigma))
    e[idx[r]] <- 1
    e[k] <- -1
    C[r, ] <- e
  }
  C %*% Sigma %*% t(C)
}

Lambdas <- lapply(1:K, function(k) diff_cov(Sigma0, k))

# compute pi(y) using MVN CDFs (K=3 -> bivariate CDF)
pi_of_y <- function(y) {
  out <- numeric(K)
  for (k in 1:K) {
    idx <- setdiff(1:K, k)
    out[k] <- as.numeric(pmvnorm(upper = y[idx], mean = c(0, 0), sigma = Lambdas[[k]]))
  }
  out
}

# Simulate Y ~ N(0, Sigma0) 
B <- 1000  # number of points for the ternary plot
Ydraws <- rmvnorm(B, mean = rep(0, K), sigma = Sigma0)

pi_mat <- t(apply(Ydraws, 1, pi_of_y))
row_sums <- rowSums(pi_mat)
pi_simplex <- pi_mat / row_sums # normalize to sum to 1

# hard-min winner (argmin of y)
winners <- apply(Ydraws, 1, function(v) which.min(v))
lab <- factor(winners, levels = 1:K, labels = paste0("Model ", 1:K))

df_tern <- data.frame(p1 = pi_simplex[,1], p2 = pi_simplex[,2], p3 = pi_simplex[,3],
                      winner = lab)


prop.table(table(lab))


p_tern <- ggtern(df_tern, aes(x = p1, y = p2, z = p3, color = winner)) +
  geom_point(alpha = 0.5, size = 0.8) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 14)
  ) +
  labs(
    T = expression(pi[1]), L = expression(pi[2]), R = expression(pi[3]),
    color = "Winner"
  )

print(p_tern)




