# 02_sparseMVN.R
# Sparse multivariate normal example (Section 5.2)

options(pillar.sigfig = 3)

library(MCMCpack)
library(mvtnorm)
library(tidyverse)
library(matrixStats)
library(latex2exp)
library(cowplot)
library(grid)
library(RColorBrewer)

source(here::here("code", "functions.R"))


# For exact replication of the figures in the paper, we provide precomputed results in
# output/sparseMVN_fitted.Rdata:
# load(here::here("output", "sparseMVN_fitted.Rdata"))
# After load this, you can jump directly to codes to create figures below. 

# Otherwise, the scripts below also allow recomputing the results from the raw data;
# these may differ slightly up to random sampling, but yield the same qualitative conclusions.


# ---------------------------
# Generate data from MVN 
# ---------------------------
generate_data_mvn <- function(n) {
  true_mean <- c(1, 1, 0.5, 0.5, 0.4, 0)  # True mean parameters
  true_cov <- diag(6) # 6 x 6 identity covariance matrix
  x <- rmvnorm(n, mean = true_mean, sigma = true_cov)
  return(x)
}


# ---------------------------
# LaD: compute Z (pointwise neg logliks under each model)
# ---------------------------
fit_models <- function(x, models) {
  n <- nrow(x)
  d <- ncol(x)
  K <- length(models)
  Z <- matrix(0, nrow = n, ncol = K)
  
  beta_hats <- colMeans(x) # Use the sample means as plug-in estimates for the free parameters.
  
  for (k in seq_len(K)) {
    # For model k, free parameters are those in models[[k]].
    model_coeffs <- ifelse(1:d %in% models[[k]], beta_hats, 0)
    # Compute residuals for each observation.
    resids <- x - matrix(model_coeffs, n, d, byrow = TRUE)
    # Negative log-likelihood for each observation under MVN with identity covariance
    Z[, k] <- 0.5*rowSums(resids^2) + (d/2) * log(2*pi)
  }
  return(Z)
}


# -------------------------------------------------------------------
# KL distances (Table 1)
# -------------------------------------------------------------------

KL_for_model <- function(mu_T, S) 0.5 * sum(mu_T[-S]^2)

# Example with your seven models:
mu_T <- c(1, 1, 0.5, 0.5, 0.4, 0)
D1 <- KL_for_model(mu_T, c(1,4))
D2 <- KL_for_model(mu_T, c(1,2))
D3 <- KL_for_model(mu_T, c(1,2,5))
D4 <- KL_for_model(mu_T, c(1,2,4))
D5 <- KL_for_model(mu_T, c(1,2,3))
D6 <- KL_for_model(mu_T, c(1,2,3,4,5))
D7 <- KL_for_model(mu_T, 1:6)
KL_true <- c(D1, D2, D3, D4, D5, D6, D7)
cat("True KL divergences (Table 1):", round(KL_true, 3), "\n")



#############################
# Simulation Settings
############################# 
set.seed(23512)

# Candidate models: each model is defined by the indices (columns) whose parameters are free.
models <- list(
  c(1,4),          # Model 1.
  c(1,2),          # Model 2.
  c(1,2,5),        # Model 3: Only the first three parameters free.
  c(1,2,4),        # Model 4.
  c(1,2,3),        # Model 5.
  c(1,2,3,4,5),    # Model 6: Parameter 6 fixed to 0.
  c(1,2,3,4,5,6)   # Model 7: All parameters free.
)

K <- length(models)

df_k <- vapply(models, length, integer(1)) # Number of free means per model: d_k = |J_k|

# Complexity levels for the candidate models (lower number means simpler).
# complexities <- c(1, 1, 2, 2, 2, 3, 4)
complexities <- df_k # should be identical.

# Simulation parameters
sample_sizes <- c(50, 500, 5000)
n_max <- max(sample_sizes)  # Maximum sample size used to generate sample 
n_sim <- 50 # Number of simulation replicates.
n_post_samples <- 1000      # Number of NIW posterior samples per simulation replicate.

delta_values <- c(0.75, 0.255, 0.05)

# NIW Prior Parameters
mu0 <- rep(0, K) # Prior mean (vector of zeros).
lambda0 <- 0.01 # Prior scaling for the mean.
nu0 <- K + 2 # Degrees of freedom (must be > K-1).
Psi0 <- 1 * diag(K)  # Prior inverse scale matrix.

# NIG priors
a0 <- rep((nu0 - K + 1)/2, K) 
b0 <- rep(diag(Psi0)/2, each = 1) 

# Coarsening parameters
alphas <- c(10, 100, Inf)




#############################
# Run simulation over 50 replicated datasets and sample sizes 
#############################
store_x <- list()
store_Z <- list()

lad_rows   <- list()
aicbic_rows <- list()
coarse_rows <- list()


for (sim in 1:n_sim) {
  # cat("Simulation: ", sim, "\n")
  x_full <- generate_data_mvn(n = n_max) 
  
  for (n in sample_sizes) {
    # cat("Sample size: ", n, "\n")
    x <- x_full[1:n, ]
    Z <- fit_models(x, models)

    # Bias-corrected LaD values: Z_{ik} + d_k/(2n)
    Z_bc  <- t(t(Z) + df_k/(2*n)) 
    
    # Save x and Z
    df_x <- as_tibble(x) %>%
      mutate(i = row_number()) %>%
      pivot_longer(-i, names_to = "dim", values_to = "xij") %>%
      mutate(sim = sim, sample_size = n, .before = 1)
    
    store_x[[length(store_x) + 1]] <- df_x
    
    df_Z <- as_tibble(Z_bc) %>%
      mutate(i = row_number()) %>%
      pivot_longer(-i, names_to = "model_idx", values_to = "Zik") %>%
      mutate(k = as.integer(sub("V","", model_idx))) %>%
      dplyr::select(i, k, Zik) %>%
      mutate(sim = sim, sample_size = n, .before = 1)
    
    store_Z[[length(store_Z) + 1]] <- df_Z
    
    
    # LaD posteriors
    mu_niw <- niw_posterior(Z_bc, mu0, lambda0, Psi0, nu0, n_post_samples)$mu 
    mu_nig <- nig_posterior(Z_bc, mu0, lambda0, a0, b0, n_post_samples)$mu
    
    # NIW
    df_niw <- as_tibble(mu_niw) %>%
      mutate(draw = row_number()) %>%
      pivot_longer(-draw, names_to = "k", values_to = "mu_sample") %>%
      mutate(k = as.integer(sub("V","", k)),
             complexity = complexities[k],
             method = "LaD_NIW",
             sim = sim, n = n, .before = 1) %>%
      dplyr::select(sim, n, k, complexity, method, draw, mu_sample)
    
    # NIG
    df_nig <- as_tibble(mu_nig) %>%
      mutate(draw = row_number()) %>%
      pivot_longer(-draw, names_to = "k", values_to = "mu_sample") %>%
      mutate(k = as.integer(sub("V","", k)),
             complexity = complexities[k],
             method = "LaD_NIG",
             sim = sim, n = n, .before = 1) %>%
      dplyr::select(sim, n, k, complexity, method, draw, mu_sample)
    
    lad_rows[[length(lad_rows)+1]] <- bind_rows(df_niw, df_nig)
    
    
    # AIC and BIC
    logLik <- -colSums(Z)
    AIC <- 2 * df_k - 2 * logLik              
    BIC <- log(n) * df_k - 2 * logLik       
    
    k_star_AIC <- which.min(AIC)
    k_star_BIC <- which.min(BIC)
    
    df_aic <- tibble(
      sim = sim, n = n, k = 1:K,
      complexity = complexities,
      method = "AIC",
      value = AIC,
      chosen_k = k_star_AIC
    )
    df_bic <- tibble(
      sim = sim, n = n, k = 1:K,
      complexity = complexities,
      method = "BIC",
      value = BIC,
      chosen_k = k_star_BIC
    )
    aicbic_rows[[length(aicbic_rows)+1]] <- bind_rows(df_aic, df_bic)
    
    
    # Coarsening (uniform prior)
    df_coarse <- purrr::map_dfr(
      alphas,
      ~coarsen_post_probs(
        x, models, alpha = .x,
        kappa0 = 0.01, m0 = rep(0, ncol(x)),
        prior_wts = NULL  # uniform prior
      )
    ) %>%
      dplyr::mutate(sim = sim, n = n, .before = 1) %>%
      dplyr::mutate(complexity = complexities[k], .after = k) %>%
      dplyr::select(sim, n, k, complexity, alpha, zeta_n, logml, prob)
    
    coarse_rows[[length(coarse_rows)+1]] <- df_coarse
    
  }
  
  if (sim %% 10 == 0) cat("Completed simulation replicate:", sim, "\n")
}

df_x_all <- bind_rows(store_x)
df_Z_all <- bind_rows(store_Z)

df_LaD <- bind_rows(lad_rows)           
df_AICBIC <- bind_rows(aicbic_rows)    
df_Coarsen <- bind_rows(coarse_rows)    



# tau values
df_d <- df_x_all %>%
  group_by(sim, sample_size) %>%
  summarise(d = n_distinct(dim), .groups = "drop")

# mu_noise(sim,n) = (d/2)log(2 \pi) + 0.5 * mean_i [ sum_j x_{ij}^2 ]
df_mu_noise <- df_x_all %>%
  group_by(sim, sample_size, i) %>%
  summarise(row_sqsum = sum(xij^2), .groups = "drop") %>%
  group_by(sim, sample_size) %>%
  summarise(mean_row_sqsum = mean(row_sqsum), .groups = "drop") %>%
  left_join(df_d, by = c("sim","sample_size")) %>%
  mutate(mu_noise = (d/2) * log(2*pi) + 0.5 * mean_row_sqsum) %>%
  dplyr::select(sim, sample_size, mu_noise)

# mu_star(sim,n) = min_k mean_i Z_{ik}
df_mu_star <- df_Z_all %>%
  group_by(sim, sample_size, k) %>%
  summarise(mu_k = mean(Zik), .groups = "drop") %>%
  group_by(sim, sample_size) %>%
  summarise(mu_star = min(mu_k), .groups = "drop")

# tau dataframe
df_tau <- df_mu_noise %>%
  left_join(df_mu_star, by = c("sim","sample_size")) %>%
  mutate(denom = mu_noise - mu_star) %>%
  tidyr::crossing(tibble(delta = delta_values)) %>%
  mutate(tau = delta / denom) %>%
  rename(n = sample_size)





#############################
# Compute scores and compare with other methods 
#############################

# Metrics
brier_loss <- function(p, A) sum((p-A)^2)
# gini_loss <- function(p, A) sum(A*(1-p) + (1-A)*p)
# log_loss <- function(p, A) sum(-A*log(p) - (1-A)*log(1-p))
# log_loss <- function(p, A) {
#   t1 <- ifelse(A == 1, -log(p), 0)
#   t2 <- ifelse(A == 0, -log(1 - p), 0)
#   sum(t1 + t2)
# }

truth_sets <- list(large = c(2), middle = c(4,5), small = c(6))
delta_to_truth <- tibble::tibble(delta = delta_values, set_name = names(truth_sets))

A_vec_for_delta <- function(delta, K) {
  set_name <- delta_to_truth$set_name[which.min(abs(delta_to_truth$delta - delta))]
  winners  <- truth_sets[[set_name]]
  as.integer(seq_len(K) %in% winners)
}

score_rows <- list()

for (delta in delta_values) {
  A_vec <- A_vec_for_delta(delta, K)
  
  for (n in sample_sizes) {
    alpha_n_soft <- n^(0.45) # LaD-soft / LaD-diag
    
    for (sim in 1:n_sim) {
      tau_val <- df_tau %>% filter(sim == !!sim, n == !!n, delta == !!delta) %>% pull(tau)
      
      # LaD
      for (lad_core in c("LaD_NIW", "LaD_NIG")) {
        df_this <- df_LaD %>%
          filter(sim == !!sim, n == !!n, method == lad_core) %>%
          arrange(draw, k)
        
        if (nrow(df_this) > 0) {
          # S x K mu-samples matrix
          mat_mu <- df_this %>%
            dplyr::select(draw, k, mu_sample) %>%
            pivot_wider(names_from = k, values_from = mu_sample) %>%
            dplyr::select(as.character(seq_len(K))) %>%
            as.matrix()
          
          # soft-minimum 
          what_soft <- selection_probabilities(mat_mu, complexities, delta, alpha_n = alpha_n_soft, mode = "soft")
          score_rows[[length(score_rows)+1]] <- tibble::tibble(
            sim = sim, n = n, delta = delta, tau = tau_val,
            method = if (lad_core == "LaD_NIW") "LaD_FC_soft" else "LaD_DC_soft",
            k = seq_len(K), score = what_soft,
            S_brier = brier_loss(what_soft, A_vec),
            # S_log = log_loss(what_soft, A_vec),
            # S_gini = gini_loss(what_soft, A_vec)
          )
          
          # hard-minimum; only for FC 
          if (lad_core == "LaD_NIW") {
            what_hard <- selection_probabilities(mat_mu, complexities, delta, alpha_n = NA, mode = "hard")
            score_rows[[length(score_rows)+1]] <- tibble::tibble(
              sim = sim, n = n, delta = delta, tau = tau_val,
              method = "LaD_FC_hard",
              k = seq_len(K), score = what_hard,
              S_brier = brier_loss(what_hard, A_vec),
              # S_log   = log_loss(what_hard, A_vec),
              # S_gini  = gini_loss(what_hard, A_vec)
            )
          }
        }
      }
      
      # Coarsening
      df_co <- df_Coarsen %>% dplyr::filter(sim == !!sim, n == !!n)
      if (nrow(df_co) > 0) {
        alphas_here <- sort(unique(df_co$alpha))
        for (a in alphas_here) {
          grp <- df_co %>% dplyr::filter(alpha == a) %>% dplyr::arrange(k)
          grp <- tibble::tibble(k = 1:K) %>% dplyr::left_join(grp, by = "k") %>% dplyr::arrange(k)
          p <- grp$prob
          lbl <- if (is.infinite(a)) "Bayes" else paste0("c-post ($\\alpha = ", a, "$)")
          score_rows[[length(score_rows)+1]] <- tibble::tibble(
            sim = sim, n = n, delta = delta, tau = tau_val,
            method = lbl,
            k = seq_len(K), score = p,
            S_brier = brier_loss(p, A_vec),
            # S_log   = log_loss(p, A_vec),
            # S_gini  = gini_loss(p, A_vec)
          )
        }
      }
      
      # AIC/BIC
      df_aic_slice <- df_AICBIC %>% filter(sim == !!sim, n == !!n, method == "AIC")
      df_bic_slice <- df_AICBIC %>% filter(sim == !!sim, n == !!n, method == "BIC")
      
      if (nrow(df_aic_slice) > 0) {
        k_star <- unique(df_aic_slice$chosen_k)
        p_aic <- as.numeric(seq_len(K) == k_star)
        score_rows[[length(score_rows)+1]] <- tibble::tibble(
          sim = sim, n = n, delta = delta, tau = tau_val,
          method = "AIC",
          k = seq_len(K), score = p_aic,
          S_brier = brier_loss(p_aic, A_vec),
          # S_log = log_loss(p_aic, A_vec),
          # S_gini = gini_loss(p_aic, A_vec)
        )
      }
      if (nrow(df_bic_slice) > 0) {
        k_star <- unique(df_bic_slice$chosen_k)
        p_bic <- as.numeric(seq_len(K) == k_star)
        score_rows[[length(score_rows)+1]] <- tibble::tibble(
          sim = sim, n = n, delta = delta, tau = tau_val,
          method = "BIC",
          k = seq_len(K), score = p_bic,
          S_brier = brier_loss(p_bic, A_vec),
          # S_log = log_loss(p_bic, A_vec),
          # S_gini = gini_loss(p_bic, A_vec)
        )
      }
    } # sim
  } #n
} # delta

df_scores <- bind_rows(score_rows)


# cat("Save results...\n")
# save(df_x_all, df_Z_all, df_LaD, df_AICBIC, df_Coarsen, df_tau, df_scores, file = "output/sparseMVN_fitted.Rdata")
# cat("Results saved successfully.\n")



################################################
# Figure 4
################################################
# df_tau %>% dplyr::filter(sim==25)

delta_values <- c(0.75, 0.255, 0.05)
K <- 7

# Summarizing LaD plot
labels_delta <- setNames(
  paste0("delta==", round(delta_values, 2)),
  delta_values
)

plot_dat <- df_scores %>%
  filter(n %in% c(50, 500, 5000)) %>% 
  filter(method == "LaD_FC_soft") %>%
  dplyr::select(sim, n, delta, tau, k, score) %>%
  mutate(
    sample_size = factor(n, levels = sort(unique(n))),
    model = factor(k, labels=c(1:K)), 
    delta = factor(delta, levels = delta_values, labels = labels_delta)
  )

plot_dat_bar <- plot_dat %>%
  filter(sim == 2) # One representative data set for Figure 4

summary_plot <- ggplot(plot_dat_bar, aes(x = model, y = score, fill = delta, color=delta)) +
  geom_col(color = "black", width = 0.8, show.legend = FALSE) +
  geom_vline(xintercept = seq(1, 7, 1), color = "grey70", linewidth = 0.2, alpha=0.3) +
  geom_hline(yintercept = seq(0, 1, 0.25), color = "grey70", linewidth = 0.2, alpha=0.3) +
  scale_fill_manual(values = c("orange", "purple1", "green3")) +
  scale_color_manual(values = c("orange", "purple1", "green3"), guide = "none") +
  scale_y_continuous(limits = c(0, 1),
                     breaks  = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.50", "0.75", "1")) +
  facet_grid(delta ~ sample_size, 
             # switch = "y",
             labeller = labeller(sample_size = function(x) paste0("n = ", x), 
                                 delta = label_parsed)
  ) +
  labs(x = "k", y = expression(hat(w)[delta](k))) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5))
  )


# Mu plot
df_LaD_NIW <- df_LaD %>%
  dplyr::filter(method == "LaD_NIW", n %in% c(50, 500, 5000), sim == 25) 

K_max <- 7

# boxplots of mu_k 
ygrid_mu <- pretty(range(df_LaD_NIW$mu_sample, na.rm = TRUE), n = 6)

p_mu <- ggplot(df_LaD_NIW, aes(x = factor(k), y = mu_sample)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3, color = "black") +
  geom_vline(xintercept = seq(1, K_max, 1), color = "grey70", linewidth = 0.2, alpha = 0.3) +
  geom_hline(yintercept = ygrid_mu,        color = "grey70", linewidth = 0.2, alpha = 0.3) +
  scale_y_continuous(limits = range(ygrid_mu), breaks = ygrid_mu) +
  facet_wrap(~ n, nrow = 1, labeller = labeller(n = function(v) paste0("n = ", v))) +
  labs(x = "k", y = expression(mu[k])) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 15),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5))
  )


fig.mvn <- plot_grid(
  summary_plot, p_mu, 
  ncol = 1, align = "v", axis = "l",
  rel_heights = c(2, 1) 
)
fig.mvn


# save for Figure 4
# fig_dir <- here::here("output", "figures")
# ggsave(file.path(fig_dir, "sparseMVN_fig4.png"), fig.mvn, width = 12, height = 12, dpi = 300)




# -------------------------------------------------------------------
# Figure 5 -- Metrics: mean Brier loss with s.e. across methods
# -------------------------------------------------------------------
df_scores_dataset <- df_scores %>%
  dplyr::filter(n %in% c(50, 500, 5000)) %>%
  dplyr::group_by(sim, n, delta, method) %>%
  dplyr::summarise(
    S_brier = dplyr::first(S_brier)
  )

df_scores_summary <- df_scores_dataset %>%
  dplyr::group_by(n, delta, method) %>%
  dplyr::summarise(
    M = dplyr::n(),
    E_brier = round(mean(S_brier), 3),
    SE_brier = round(sd(S_brier)/sqrt(M), 3)
  )


method_levels <- c(
  "LaD-soft", "LaD-hard", "LaD-diag",
  "c-post ($\\alpha = 10$)", "c-post ($\\alpha = 100$)", "Bayes",
  "AIC", "BIC"
)


dark2_colors <- brewer.pal(8, "Dark2")
color_vals <- setNames(dark2_colors, method_levels)


shape_vals <- c(
  "LaD-soft"                 = 1,  # open circle
  "LaD-hard"                 = 2,  # open triangle up
  "LaD-diag"                 = 0,  # open square
  "c-post ($\\alpha = 10$)"  = 5,  # open diamond
  "c-post ($\\alpha = 100$)" = 6,  # open triangle down
  "Bayes"                    = 11,  # open circle with X
  "AIC"                      = 4,  # x
  "BIC"                      = 7   # open rec with X
)

plot_df <- df_scores_summary %>%
  dplyr::mutate(
    method_chr = as.character(method),
    method_plot = dplyr::case_when(
      method_chr == "LaD_FC_soft" ~ "LaD-soft",
      method_chr == "LaD_FC_hard" ~ "LaD-hard",
      method_chr == "LaD_DC_soft" ~ "LaD-diag",
      method_chr == "Bayes" ~ "Bayes",
      grepl("^c-post", method_chr) ~ method_chr,   # e.g., "c-post ($\\alpha = 10$)"
      TRUE ~ method_chr
    ),
    method_plot = factor(method_plot, levels = method_levels),
    delta = factor(delta, levels = delta_values)
  ) %>%
  dplyr::select(-method_chr)


legend_labels <- latex2exp::TeX(method_levels)

labels_delta <- setNames(
  paste0("delta==", round(delta_values, 2)),
  delta_values
)

hgrid_df <- plot_df %>% 
  dplyr::group_by(delta) %>% 
  dplyr::summarise(ymin = min(E_brier - SE_brier, na.rm = TRUE), 
                  ymax = max(E_brier + SE_brier, na.rm = TRUE)) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(y = list(pretty(c(ymin, ymax), n = 5))) %>% 
  tidyr::unnest(y) %>% 
  dplyr::select(delta, y)

pd <- position_dodge(width = 0.7)

vlines_at <- c(1.5, 2.5)

p_metric <- ggplot(plot_df, aes(x = factor(n), y = E_brier, color = method_plot, shape = method_plot)) +
  geom_vline(xintercept = vlines_at, linetype = "dashed", color = "grey80", linewidth = 0.7) +
  geom_hline(data = hgrid_df, aes(yintercept = y), color = "grey80", linewidth = 0.3, alpha = 0.3) +
  geom_point(position = pd, size = 2.6) +
  geom_errorbar(aes(ymin = E_brier - SE_brier, ymax = E_brier + SE_brier),
                position = pd, width = 0.18, linewidth = 0.5) +
  facet_wrap(
    ~ delta, nrow = 1, scales = "free_y",
    labeller = labeller(delta = as_labeller(labels_delta, default = label_parsed))
  ) +
  scale_color_manual(
    values = color_vals,
    limits = method_levels,
    labels = legend_labels,
    name   = "Method"
  ) +
  scale_shape_manual(
    values = shape_vals,
    limits = method_levels,
    labels = legend_labels,
    name   = "Method"
  ) +
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3, linetype = 0, linewidth = 0)),
    shape = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3, linetype = 0, linewidth = 0))
  ) +
  labs(
    x = "n",
    y = "Mean Brier loss with s.e."
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 15),
    axis.text  = element_text(size = 13),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

p_metric


# save for Figure 5
# fig_dir <- here::here("output", "figures")
# ggsave(file.path(fig_dir, "sparseMVN_fig5.png"), p_metric, width = 12, height = 5, dpi = 300)







