# 01_shapley_gmm.R
# Shapley galaxy Gaussian mixtures (Figures 1 and 3)

library(mclust)
library(matrixStats)
library(tidyverse)
library(cowplot)
library(here)

source(here::here("code","functions.R"))


# For exact replication of the figures in the paper, we provide precomputed results in
# output/shapley_gmm_fitted.Rdata:
# load(here::here("output", "shapley_gmm_fitted.Rdata"))


# Otherwise, the scripts below also allow recomputing the results;
# these may differ slightly due to random initialization in the Gaussian mixture fits, but yield the same qualitative conclusions.


# Permuted Shapley galaxy velocities data
x <- read.csv(here::here("data", "processed", "julia_run", "x_perm.csv"))$x_perm

# Sample sizes used in the Shapley example
sample_sizes <- c(40, 120, 400, 1200, 4000)

Kmax <- 15 # maximum K

# NIW prior hyperparameters for LaD posterior
mu0 <- rep(0, Kmax)
lambda0 <- 0.01
nu0 <- Kmax + 2
Psi0 <- diag(Kmax)
n_post_samples <- 1000


# -------------------------------------------------------------------
# Fit Gaussian mixtures, build LaD matrix Z, AIC/BIC, LaD posterior
# -------------------------------------------------------------------
set.seed(9910)  # for reproducibility of random initializations

results_Z <- vector("list", length(sample_sizes))
results_fits <- vector("list", length(sample_sizes))
results_AIC <- vector("list", length(sample_sizes))
results_BIC <- vector("list", length(sample_sizes))
results_mu <- vector("list", length(sample_sizes))
names(results_Z) <- names(results_fits) <- names(results_AIC) <- names(results_BIC) <- names(results_mu) <- as.character(sample_sizes)

for (n in sample_sizes) {
  cat("Fitting GMM for Sample size:", n, "\n")
  x_n <- x[1:n]

  fits <- vector("list", Kmax)
  Z <- matrix(NA, nrow = n, ncol = Kmax)  # raw per-point NLL
  Z_bc <- matrix(NA, nrow = n, ncol = Kmax)  # bias-corrected per-point NLL
  
  for (K in 1:Kmax) {
    best_fit <- NULL
    best_ll <- -Inf
    n_restarts <- 50
    
    for (r in 1:n_restarts) {
      # Random hierarchical clustering initialization
      init <- list(hcPairs = hcRandomPairs(matrix(x_n, ncol = 1)))
      
      fit_try <- Mclust(x_n, G = K, modelNames = "V", prior = priorControl(dof=10), initialization = init)
      
      if (!is.null(fit_try$loglik) && fit_try$loglik > best_ll) {
        best_fit <- fit_try
        best_ll <- fit_try$loglik
      }
    }
    
    fit <- best_fit
    fits[[K]] <- fit
    
    pro <- pmax(as.numeric(fit$parameters$pro), .Machine$double.xmin) 
    ms <- as.numeric(fit$parameters$mean)      
    s2 <- pmax(as.numeric(fit$parameters$variance$sigmasq), .Machine$double.eps)
   
    logW <- matrix(NA, nrow = n, ncol = K)
    for (j in seq_len(K)) {
      logW[, j] <- log(pro[j]) + dnorm(x_n, mean = ms[j], sd = sqrt(s2[j]), log = TRUE)
    }
    Z[, K] <- -matrixStats::rowLogSumExps(logW)  # per-point NLL under the K-component mixture
    
    # Bias correction dk/(2n) with dk = 3K - 1 
    dk <- 3*K-1
    Z_bc[, K] <- Z[, K] + dk/(2*n)
  }
  
  results_Z[[as.character(n)]] <- Z_bc
  results_fits[[as.character(n)]] <- fits
  
  # AIC/BIC using the same log-likelihood definition
  logLik <- -colSums(Z)
  dK <- 3*(1:Kmax) - 1 # parameter count 
  AIC <- 2*dK - 2*logLik
  BIC <- log(n)*dK - 2*logLik
  
  results_AIC[[as.character(n)]] <- AIC
  results_BIC[[as.character(n)]] <- BIC
  
  # LaD posterior draws
  mu_samples <- niw_posterior(Z_bc, mu0, lambda0, Psi0, nu0, n_post_samples)$mu
  results_mu[[as.character(n)]] <- mu_samples
}

# We onnly use K from 1 to 10 for LaD result in the figures
results_mu <- lapply(results_mu, function(M) M[, 1:10])
results_Z <- lapply(results_Z,  function(Z) Z[, 1:10])


# save results for later reuse 
save(results_mu, results_AIC, results_BIC, results_Z, results_fits, x,
     file = here::here("output", "shapley_gmm_fitted.Rdata"))






# -------------------------------------------------------------------
# Figure 1
# -------------------------------------------------------------------

## (a) Histogram of Shapley galaxy velocities
p_left <- ggplot(data.frame(x = x), aes(x)) +
  geom_histogram(bins = 150, fill = "grey50", color = "grey50") +
  labs(x = "Velocity (1000 km/s)", y = NULL, title = "Shapley galaxy data") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5)),
  )

## (b) Boxplots of MFM posterior draws of K
mfm_dir <- here::here("data", "processed", "julia_run", "MFM")

read_k_draws <- function(n) {
  f <- file.path(mfm_dir, sprintf("mfm_n%05d_k_draws.csv", n))
  readr::read_csv(f, show_col_types = FALSE) %>%
    transmute(n = as.numeric(n), K = as.integer(K))
}

k_draws_df <- bind_rows(lapply(sample_sizes, read_k_draws))

p_mid <- ggplot(k_draws_df, aes(x = factor(n), y = K)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, fill = "white", color = "black") + # outliers removed
  scale_y_continuous(breaks = c(5, 10, 15, 20, 25, 30), limits = c(1, 30)) +
  geom_vline(xintercept = k_draws_df$n, color = "grey70", linewidth = 0.2, alpha=0.3) +
  geom_hline(yintercept = seq(5, 30, 5), color = "grey70", linewidth = 0.2, alpha=0.3) +
  labs(x = "n", y = "k", title = "Bayesian posterior") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5)),
  )


## (c) AIC/BIC winners
winners <- tibble(
  sample_size = sample_sizes,
  AIC = vapply(results_AIC, which.min, integer(1)),
  BIC = vapply(results_BIC, which.min, integer(1))
) %>%
  pivot_longer(cols = c(AIC, BIC), names_to = "criterion", values_to = "K") %>%
  mutate(sample_size_fac = factor(sample_size, levels = sample_sizes))


p_right <- ggplot(winners, aes(x = sample_size_fac, y = K,
                               group = criterion, color = criterion, linetype = criterion)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_y_continuous(breaks = c(3, 6, 9, 12, 15), limits = c(1, Kmax)) +
  geom_vline(xintercept = winners$sample_size_fac, color = "grey70", linewidth = 0.2, alpha=0.3) +
  geom_hline(yintercept = seq(3, 15, 3), color = "grey70", linewidth = 0.2, alpha=0.3) +
  scale_color_manual(values = c("AIC" = "red1", "BIC" = "blue1")) +
  scale_linetype_manual(values = c("AIC" = "dashed", "BIC" = "solid"), guide="none") +
  labs(x = "n", y = "k", color = NULL, title = "Information Criteria") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 15),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = c(0.12, 0.8),    
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.key = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5)),
  )



fig1 <- gridExtra::grid.arrange(p_left, p_mid, p_right, ncol = 3)


# save for Figure 1
# fig_dir <- here::here("output", "figures")
# ggsave(file.path(fig_dir, "shapley_fig1.png"), fig1, width = 15, height = 5, dpi = 300)





# -------------------------------------------------------------------
# Figure 3: four-row GMM plot
# -------------------------------------------------------------------

ns <- c(40, 400, 4000) # subset of sample_sizes used in Figure 3

## Top row: MFM posterior on K with AIC/BIC overlay
read_k_pmf <- function(n) {
  f <- file.path(mfm_dir, sprintf("mfm_n%05d_k_posterior.csv", n))
  readr::read_csv(f, show_col_types = FALSE) %>%
    transmute(k = as.integer(k), prob = prob, n = n)
}

posterior_dfs <- bind_rows(lapply(ns, read_k_pmf)) %>%
  mutate(n = factor(n, levels = ns))

winners_long <- winners %>% transmute(n = sample_size, criterion, K) %>% 
  filter(n %in% ns)

w_AIC <- winners_long %>% filter(criterion == "AIC", n %in% ns)
w_BIC <- winners_long %>% filter(criterion == "BIC", n %in% ns)


p_row1_base <- ggplot(posterior_dfs, aes(k, prob)) +
  geom_col(fill = "lightskyblue2", color = "steelblue", width = 0.85) +
  geom_segment(data = w_AIC, aes(x = K, xend = K, y = 0, yend = 0.25),
               color = "red1", linetype = "dashed", linewidth = 1) +
  geom_segment(data = w_BIC, aes(x = K, xend = K, y = 0, yend = 0.25),
               color = "blue1", linetype = "dotdash", linewidth = 1) +
  geom_vline(xintercept = seq(1, 30, 5), color = "grey70", linewidth = 0.2, alpha=0.3) +
  geom_hline(yintercept = seq(0, 0.25, 0.05), color = "grey70", linewidth = 0.2, alpha=0.3) +
  scale_y_continuous(limits = c(0, 0.25),
                     breaks  = c(0, 0.05, 0.1, 0.15, 0.2, 0.25),
                     labels = c("0", "0.05", "0.10", "0.15", "0.2", "0.25")) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25, 30), limits = c(1, 30)) +
  facet_wrap(~ n, nrow = 1, labeller = labeller(n = function(v) paste0("n = ", v))) +
  labs(x = "k", y = expression(pi(k~"|"~x)))+
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 20),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 17),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5))
  ) 

legend_pos <- tibble::tibble(
  n = factor(ns, levels = ns),
  x0 = 20, # left edge for legend glyphs
  x1 = 24, # text anchor to the right
  y_b = 0.22, # Bayes posterior row 
  y_a = 0.2, # AIC row
  y_c = 0.18 # BIC row
)


p_row1 <- p_row1_base +
  # Bayes posterior swatch
  geom_rect(data = legend_pos,
            aes(xmin = x0, xmax = x0 + 3, ymin = y_b - 0.002, ymax = y_b + 0.002),
            inherit.aes = FALSE, fill = "lightskyblue2", color = NA) +
  geom_text(data = legend_pos,
            aes(x = x1, y = y_b, label = "posterior"),
            inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 5) +
  # AIC line swatch + label
  geom_segment(data = legend_pos,
               aes(x = x0, xend = x0 + 3, y = y_a, yend = y_a),
               inherit.aes = FALSE, color = "red1", linetype = "dashed", linewidth = 0.9) +
  geom_text(data = legend_pos,
            aes(x = x1, y = y_a, label = "AIC"),
            inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 5) +
  # BIC line swatch + label
  geom_segment(data = legend_pos,
               aes(x = x0, xend = x0 + 3, y = y_c, yend = y_c),
               inherit.aes = FALSE, color = "blue1", linetype = "dotdash", linewidth = 0.9) +
  geom_text(data = legend_pos,
            aes(x = x1, y = y_c, label = "BIC"),
            inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 5)


## Second row: posterior draws of LaD means mu_k
K_pick <- c(2,3,4)
mu_samples_4000 <- results_mu$`4000`
mu_means_4000 <- colMeans(mu_samples_4000)
names(mu_means_4000) <- as.character(seq_len(length(mu_means_4000)))

levels_common <- mu_means_4000[as.character(K_pick)]

tolerance_common <- tidyr::crossing(
  n = ns,
  Kline = K_pick
) %>%
  mutate(level = as.numeric(levels_common[as.character(Kline)]),
         Klab  = paste0("K=", Kline))

df_draws <- do.call(
  rbind,
  lapply(names(results_mu), function(nm) {
    M <- results_mu[[nm]]
    K <- ncol(M)
    tibble(n  = as.integer(nm),
           K  = rep(1:K, each = nrow(M)),
           mu = as.vector(M))
  })
) %>% filter(n %in% ns)

p_row2a <- ggplot(df_draws, aes(x = factor(K), y = mu)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3, color = "black") +
  geom_vline(xintercept = seq(1, 10, 1), color = "grey70", linewidth = 0.2, alpha=0.3) +
  geom_hline(yintercept = seq(2.25, 4, 0.25), color = "grey70", linewidth = 0.2, alpha=0.3) +
  facet_wrap(~ n, nrow = 1, labeller = labeller(n = function(v) paste0("n = ", v))) +
  labs(x = "k", y = expression(mu[k])) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 20),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 17),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5))
  ) 


## Third row: SLC scores for fixed deltas
# print(as.numeric(mu_means_4000[as.character(K_pick)] - min(mu_means_4000)))
# print(as.numeric(mu_means_4000[as.character(K_pick-1)] - min(mu_means_4000)))
# mu_samples_4000[,K_pick] - rep(min(mu_samples_4000), 3) 

delta_vals <- c(0.30, 0.12, 0.06)
names(delta_vals) <- paste0("delta_", K_pick)
colors_delta <- setNames(c("orange", "purple1", "green3"), names(delta_vals))
alpha_n <- function(n) n^(0.45)  
complexities <- 1:10

res_fixed <- imap(results_mu, function(M, nm) {
  n_val <- as.numeric(nm)
  a_n <- alpha_n(n_val)
  Kseq <- seq_len(ncol(M))
  
  imap_dfr(delta_vals, function(delta_val, lbl) {
    w_hat <- selection_probabilities(M, complexities = Kseq, delta = delta_val, alpha_n = a_n)
    tibble(n = n_val, K = Kseq, delta = delta_val, delta_label = lbl, w_hat = as.numeric(w_hat))
  })
}) %>%
  bind_rows() %>%
  filter(n %in% ns) %>%
  mutate(n = factor(n, levels = ns),
         delta_label = factor(delta_label, levels = names(delta_vals)))


expr_to_chr <- function(e) paste(deparse(e), collapse = "")

ann_df <- res_fixed %>%
  dplyr::distinct(n, delta_label, delta) %>%
  dplyr::mutate(
    label = sprintf('delta=="%s"', formatC(delta, format = "f", digits = 2))
  )

p_row3 <- ggplot(res_fixed, aes(x = K, y = w_hat, fill = delta_label)) +
  geom_col(color = "black", width = 0.8, show.legend = FALSE) +
  facet_grid(rows = vars(delta_label),
             cols = vars(n),
             labeller = labeller(n = function(v) paste0("n = ", v)),
             scales = "fixed") +
  geom_vline(xintercept = seq(1, 10, 1), color = "grey70", linewidth = 0.2, alpha=0.3) +
  geom_hline(yintercept = seq(0, 1, 0.25), color = "grey70", linewidth = 0.2, alpha=0.3) +
  geom_text(
    data = ann_df,
    aes(label = label, color = delta_label),
    x = Inf, y = 0.93, hjust = 1.5, vjust = 0.8, 
    size = 6, inherit.aes = FALSE, parse = TRUE
  ) +
  scale_fill_manual(values = colors_delta) +
  scale_color_manual(values = colors_delta, guide = "none") +
  scale_x_continuous(breaks = seq_len(max(res_fixed$K))) +
  scale_y_continuous(limits = c(0, 1),
                     breaks  = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.50", "0.75", "1")) +
  labs(x = "k", y = expression(hat(w)[delta](k))) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),                
    strip.text = element_text(size = 20),
    strip.background.x= element_rect(fill = NA, colour = "black", linewidth = 0.8),
    strip.background.y= element_blank(),
    strip.text.y = element_blank(),  
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 17),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5))
  )



## Fourth row: SLC scores as a function of tau in [0,1]
tau_grid <- seq(0, 1, length.out = 101)
bin_width  <- min(diff(tau_grid))

res <- list();

for (nm in names(results_mu)) {
  n_val <- as.numeric(nm)
  x_n <- x[1:n_val]
  mu_samples <- results_mu[[nm]]   # S x K matrix
  K <- ncol(mu_samples)
  
  # mu_star (best for each n)
  mu_bar <- colMeans(mu_samples)
  mu_star <- min(mu_bar)
  
  # mu_noise (for each n) - Uniform noise baseline on [min(x), max(x)] 
  a <- min(x_n); b <- max(x_n)
  mu_noise <- -mean(dunif(x_n, min=a, max=b, log=TRUE))
  
  denom <- mu_noise - mu_star
  delta_grid <- denom * tau_grid
  
  a_n <- n_val^(0.45) # no meaning for this example
  
  # selection probabilities
  res_n <- purrr::map2_dfr(tau_grid, delta_grid, ~{
    xleft  <- pmax(.x - bin_width/2, 0)
    xright <- pmin(.x + bin_width/2, 1)
    tibble(
      n = n_val, tau = .x, K=seq_len(K), complexities = complexities, 
      sel_prob = as.numeric( selection_probabilities(mu_samples, complexities, delta = .y, alpha_n = a_n) ),
      xleft = xleft, xright = xright, xmid = (xleft + xright)/2
    )
  })
  
  res[[nm]] <- res_n
}

res <- bind_rows(res)

n_levels <- sort(unique(res$n))

res <- res %>% 
  filter(n %in% ns) %>%
  mutate(n_f = factor(n, levels = n_levels))

heat_cols <- c("white", "orange1", "orange2", "red1", "red3")

p_row4 <- ggplot() +
  geom_rect(data = res,
            aes(xmin = xleft, xmax = xright,
                ymin = K - 0.5, ymax = K + 0.5,
                fill = sel_prob)) +
  geom_vline(xintercept = seq(0, 1, 0.1), color = "grey70", linewidth = 0.2, alpha=0.3) +
  geom_hline(yintercept = 1:10, color = "grey70", linewidth = 0.2, alpha=0.3) +
  scale_fill_gradientn(
    name = expression(hat(w)[delta](k)),
    limits = c(0, 1),
    colours = heat_cols,
    breaks = c(0, 0.5, 1),  
    labels = c("0", "0.5", "1"),
    oob = scales::squish
  ) +
  labs(y="k", x = expression(hat(tau))) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.2),
    labels = function(x) ifelse(x %in% c(0, 1), as.character(x), sprintf("%.1f", x)),
    limits = c(0, 1)
  ) +
  scale_y_continuous(breaks = complexities) +
  facet_wrap(~ n_f, nrow = 1, dir="v",
             labeller = labeller(n_f = function(v) paste0("n = ", v))) + 
  guides(
    fill = guide_colorbar(
      barwidth = unit(5, "cm"),
      barheight = unit(0.7, "cm"),
      title.position = "left",
      title.hjust = 0.7,
      title.vjust = 0.9
    )) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 20),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 17),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5)),
    legend.position = "bottom",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 20)
  )  


# Assemble the 4-row figure (Figure 3)
fig.gmm <- plot_grid(
  p_row1, p_row2a, p_row3, p_row4,
  ncol = 1, align = "v", axis = "l",
  rel_heights = c(1, 1, 1.2, 1.3) 
)
fig.gmm


# fig_dir <- here::here("output", "figures")
# ggsave(file.path(fig_dir, "shapley_fig3.png"), fig.gmm, width = 15, height = 20, dpi = 300)



