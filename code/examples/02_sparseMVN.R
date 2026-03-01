.libPaths("~/rlibs")
options(pillar.sigfig = 3)

library(MCMCpack)
library(mvtnorm)
library(tidyverse)
library(matrixStats)
library(latex2exp)
library(cowplot)
library(grid)

# ---------------------------
# Generate data from MVN 
# ---------------------------
generate_data_mvn <- function(n) {
  # true_mean <- c(1, 1, 0.3, 0.3, 0.25, 0)  # True mean parameters
  true_mean <- c(1, 1, 0.5, 0.5, 0.4, 0)  # True mean parameters
  true_cov <- diag(6)                      # 6 x 6 identity covariance matrix
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
    # Under the assumption of a multivariate normal with identity covariance, the 
    # (un-normalized) negative log-likelihood for each observation is given by:
    Z[, k] <- 0.5*rowSums(resids^2) + (d/2) * log(2*pi)
  }
  
  return(Z)
}


# ---------------------------
# NIW posterior for LaD mean
# ---------------------------

niw_posterior <- function(Z, mu0, lambda0, Psi0, nu0, n_samples) {
  # Function to draw posterior samples of mu and Sigma using NIW conjugacy
  
  n <- nrow(Z)
  K <- ncol(Z)
  
  # Compute sufficient statistics
  Z_bar <- colMeans(Z)
  Psi_obs <- t(Z - matrix(Z_bar, n, K, byrow = TRUE)) %*% (Z - matrix(Z_bar, n, K, byrow = TRUE))
  
  # Update posterior parameters
  lambda_n <- lambda0 + n
  mu_n <- (lambda0*mu0 + n*Z_bar)/lambda_n
  nu_n <- nu0 + n
  Psi_n <- Psi0 + Psi_obs + (lambda0*n/lambda_n)*(matrix(Z_bar-mu0, ncol=1) %*% matrix(Z_bar-mu0, nrow=1))
  
  mu_samples <- matrix(0, nrow = n_samples, ncol = K)
  Sigma_samples <- array(0, dim = c(K, K, n_samples))
  
  # Draw posterior samples
  for (i in 1:n_samples) {
    # Sample Sigma from the Inverse-Wishart distribution
    # Sigma <- solve(rwish(nu_n, solve(Psi_n)))
    Sigma <- MCMCpack::riwish(nu_n, Psi_n)
    
    # Sample mu from the conditional Normal distribution given Sigma
    mu <- mvtnorm::rmvnorm(1, mu_n, Sigma/lambda_n)
    # mu <- mvrnorm(1, mu_n, Sigma / lambda_n)
    
    mu_samples[i, ] <- mu
    Sigma_samples[, , i] <- Sigma
  }

  return(list(mu = mu_samples, Sigma = Sigma_samples))
}

# ---------------------------
# NIG 
# ---------------------------
nig_posterior <- function(Z, mu0, lambda0, a0, b0, n_samples) {
  # Independent NIG posterior per column (diagonal Covariance Sigma)
  n <- nrow(Z)
  K <- ncol(Z)
  
  mu_samples <- matrix(NA, nrow = n_samples, ncol = K)
  sigma2_samples <- matrix(NA, nrow = n_samples, ncol = K)
  
  for (k in seq_len(K)) {
    xk <- Z[, k]
    xbar <- mean(xk)
    sse <- sum((xk - xbar)^2)
    
    lambda_n <- lambda0 + n
    mu_n <- (lambda0*mu0[k] + n*xbar) / lambda_n
    a_n <- a0[k] + n/2
    b_n <- b0[k] + 0.5*sse + 0.5*(lambda0*n / lambda_n) * (xbar - mu0[k])^2
    
    # Draws
    sig2_draws <- MCMCpack::rinvgamma(n_samples, shape = a_n, scale = b_n) # IG(a,b)
    mu_draws <- rnorm(n_samples, mean = mu_n, sd = sqrt(sig2_draws / lambda_n))
    
    mu_samples[, k] <- mu_draws
    sigma2_samples[, k] <- sig2_draws
  }
  list(mu = mu_samples, sigma2 = sigma2_samples)
}

# ---------------------------
# LaD selection-probabilities (given posterior draws of LaD means)
# ---------------------------
selection_probabilities <- function(mu_samples, complexities, delta, alpha,
                                    mode = c("soft", "hard"), tol = 1e-12) {
  # mu_samples: S x K matrix of posterior draws for mu_k (smaller is better)
  # complexities: length-K integer vector of complexity classes
  # delta: tolerance for "surviving" from the global best
  # alpha: temperature (used only in mode="soft"; ignored in mode="hard")
  # mode: "soft" = original exp(-alpha * (mu - class_min)) weights
  #       "hard" = hard minimum inside class (uniform over ties)
  # tol:  tie tolerance inside a class for hard rule
  
  mode <- match.arg(mode)
  S <- nrow(mu_samples); K <- ncol(mu_samples)
  
  class_selection_counts <- numeric(K)  # factor 1: class selection frequency
  weight_sums <- numeric(K)             # factor 2: within-class weights
  
  for (s in seq_len(S)) {
    mu_draw <- mu_samples[s, ]
    
    # Step 1: choose class among models within delta of global best
    overall_best <- min(mu_draw)
    surviving <- which(mu_draw <= overall_best + delta)
    chosen_class <- min(complexities[surviving])
    
    selected_models <- which(complexities == chosen_class)
    class_selection_counts[selected_models] <- class_selection_counts[selected_models] + 1
    
    # Step 2: within-class weighting
    for (k in seq_len(K)) {
      class_k <- complexities[k]
      idx <- which(complexities == class_k)
      mu_best_in_class <- min(mu_draw[idx])
      
      if (mode == "soft") {
        # soft minimum: temperature-based weight
        weight <- exp(-alpha * (mu_draw[k] - mu_best_in_class))
      } else {
        # hard minimum: split uniformly among class minimizers
        winners <- idx[abs(mu_draw[idx] - mu_best_in_class) <= tol]
        weight <- if (k %in% winners) 1 / length(winners) else 0
      }
      weight_sums[k] <- weight_sums[k] + weight
    }
  }
  
  # Monte Carlo averages of the two factors; product is the two-factor score
  p_class <- class_selection_counts / S
  avg_weight <- weight_sums / S
  p_class * avg_weight
}


# ---------------------------
# Coarsening via power posterior 
# ---------------------------
coarsen_post_probs <- function(x, models, alpha,
                               kappa0 = 0.01,      # prior precision on free means (coarsening prior)
                               prior_wts = NULL,   # length-K prior over models; renormalized inside
                               m0 = NULL) {        # prior mean vector (length d)
  # Coarsened (power) posterior model probabilities for MVN mean, Sigma = I
  
  n <- nrow(x); d <- ncol(x); K <- length(models)
  zeta_n <- if (is.infinite(alpha)) 1 else alpha / (alpha + n)
  xbar <- colMeans(x)
  if (is.null(m0)) m0 <- rep(0, d)
  
  # coefficients in (★)
  a <- 0.5 * (log(kappa0) - log(kappa0 + zeta_n * n))
  b <- 0.5 * (kappa0 * zeta_n * n) / (kappa0 + zeta_n * n)
  c <- 0.5 * zeta_n * n
  
  logml <- numeric(K)
  for (k in seq_len(K)) {
    J <- models[[k]]
    if (length(J) == 0L) { logml[k] <- 0; next }
    diffJ <- xbar[J] - m0[J]
    logml[k] <- sum(a - b * diffJ^2 + c * xbar[J]^2)
  }
  
  if (is.null(prior_wts)) prior_wts <- rep(1 / K, K)
  prior_wts <- prior_wts / sum(prior_wts)
  
  w <- logml + log(prior_wts)
  w <- w - max(w) 
  post <- exp(w); post <- post / sum(post)
  
  tibble::tibble(
    k      = seq_len(K),
    prob   = post,
    logml  = logml,
    alpha  = alpha,
    zeta_n = zeta_n
  )
}


#############################
# Simulation Settings
#############################
# seed <- 1243412
# seed <- 25102021
seed <- 251102
set.seed(seed)

# Candidate models: each model is defined by the indices (columns) whose parameters are free.
models <- list(
  c(1,4),          # Model 1.
  c(1,2),          # Model 2.
  c(1,2,5),        # Model 3: Only the first three parameters free.
  c(1,2,3),        # Model 4.
  c(1,2,4),        # Model 5.
  c(1,2,3,4,5),    # Model 6: Parameter 6 fixed to 0.
  c(1,2,3,4,5,6)   # Model 7: All parameters free.
)

K <- length(models)

df_k <- vapply(models, length, integer(1))

# Complexity levels for the candidate models (lower number means simpler).
complexities <- c(1, 1, 2, 2, 2, 3, 4)

# Simulation parameters
# sample_sizes <- c(50, 500, 5000)
sample_sizes <- c(50, 100, 500, 1000, 5000, 10000)
n_max <- max(sample_sizes)  # Maximum sample size used to generate sample 
n_sim <- 50      # Number of simulation replicates.
n_post_samples <- 1000      # Number of NIW posterior samples per simulation replicate.

# gamma_values <- c(1, 0.082, 0.01) # 0.091/1.12125
# gamma_values <- c(1, 0.3, 0.01) # 0.5 / 1.70125
# gamma_values <- c(1, 0.083, 0.085, 0.09, 0.05, 0.045, 0.01)

# noise: 0.5*(1^2 + 1^2 + 0.5^2 + 0.5^2 + 0.25^2) = 1.28125
# delta_values <- 1.12125 * gamma_values
# delta_values <- c(1, 0.092, 0.02)
# delta_values <- c(1, 0.26, 0.02)
delta_values <- c(0.75, 0.255, 0.05)
# delta_values <- 1.28125 * gamma_values
# delta_values <- c(0.01, 0.091, 0.58)  # Tolerance levels based on ture KL divergence.

# NIW Prior Parameters
mu0 <- rep(0, K) # Prior mean (vector of zeros).
lambda0 <- 0.01 # Prior scaling for the mean.
nu0 <- K + 2 # Degrees of freedom (must be > K-1).
Psi0 <- 1 * diag(K)  # Prior inverse scale matrix.

# NIG priors
# mu0 <- rep(0, K)
# lambda0 <- 0.01
a0 <- rep((nu0 - K + 1)/2, K) 
b0 <- rep(diag(Psi0)/2, each = 1) 

# Coarsening parameters
alphas <- c(10, 100, Inf)


#############################
# Simulation Phase A: run simulation over replicated datasets and sample sizes 
#############################

store_x <- list()
store_Z <- list()

lad_rows   <- list()
aicbic_rows <- list()
coarse_rows <- list()

# y <- generate_data_mvn(n=100000)
# z <- fit_models(y, models)
# true_cov <- cov(z)

for (sim in 1:n_sim) {
  x_full <- generate_data_mvn(n = n_max) 
  
  for (n in sample_sizes) {
    # cat("Sample size: ", n, "\n")
    x <- x_full[1:n, ]
    Z <- fit_models(x, models)

    Z  <- t(t(Z) + df_k/(2*n)) 
    
    # Save x and Z
    df_x <- as_tibble(x) |>
      mutate(i = row_number()) |>
      pivot_longer(-i, names_to = "dim", values_to = "xij") |>
      mutate(sim = sim, sample_size = n, .before = 1)
    
    store_x[[length(store_x) + 1]] <- df_x
    
    df_Z <- as_tibble(Z) |>
      mutate(i = row_number()) |>
      pivot_longer(-i, names_to = "model_idx", values_to = "Zik") |>
      mutate(k = as.integer(sub("V","", model_idx))) |>
      select(i, k, Zik) |>
      mutate(sim = sim, sample_size = n, .before = 1)
    
    store_Z[[length(store_Z) + 1]] <- df_Z
    
    
    # LaD posteriors
    mu_niw <- niw_posterior(Z, mu0, lambda0, Psi0, nu0, n_post_samples)$mu # NIW bayesian inference
    # Sigma_samples <- niw_posterior(Z, mu0, lambda0, Psi0, nu0, n_post_samples)$Sigma # NIW bayesian inference
    # mu_samples <- rmvnorm(n_post_samples, colMeans(Z), cov(Z) / n) # using normal distribution
    # mu_samples <- rmvnorm(n_post_samples, colMeans(Z), true_cov / n) # using true mu, cov.
    mu_nig <- nig_posterior(Z, mu0, lambda0, a0, b0, n_post_samples)$mu
    
    # NIW
    df_niw <- as_tibble(mu_niw) |>
      mutate(draw = row_number()) |>
      pivot_longer(-draw, names_to = "k", values_to = "mu_sample") |>
      mutate(k = as.integer(sub("V","", k)),
             complexity = complexities[k],
             method = "LaD_NIW",
             sim = sim, n = n, .before = 1) |>
      select(sim, n, k, complexity, method, draw, mu_sample)
    
    # NIG
    df_nig <- as_tibble(mu_nig) |>
      mutate(draw = row_number()) |>
      pivot_longer(-draw, names_to = "k", values_to = "mu_sample") |>
      mutate(k = as.integer(sub("V","", k)),
             complexity = complexities[k],
             method = "LaD_NIG",
             sim = sim, n = n, .before = 1) |>
      select(sim, n, k, complexity, method, draw, mu_sample)
    
    lad_rows[[length(lad_rows)+1]] <- bind_rows(df_niw, df_nig)
    
    
    # AIC and BIC
    # df_k <- vapply(models, length, integer(1))   # number of free means per model
    
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
        prior_wts = NULL                 # uniform prior
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



# gamma values
df_d <- df_x_all %>%
  group_by(sim, sample_size) %>%
  summarise(d = n_distinct(dim), .groups = "drop")

# mu_noise(sim,n) = (d/2)log(2π) + 0.5 * mean_i [ sum_j x_{ij}^2 ]
df_mu_noise <- df_x_all %>%
  group_by(sim, sample_size, i) %>%
  summarise(row_sqsum = sum(xij^2), .groups = "drop") %>%
  group_by(sim, sample_size) %>%
  summarise(mean_row_sqsum = mean(row_sqsum), .groups = "drop") %>%
  left_join(df_d, by = c("sim","sample_size")) %>%
  mutate(mu_noise = (d/2) * log(2*pi) + 0.5 * mean_row_sqsum) %>%
  select(sim, sample_size, mu_noise)

# mu_star(sim,n) = min_k mean_i Z_{ik}
df_mu_star <- df_Z_all %>%
  group_by(sim, sample_size, k) %>%
  summarise(mu_k = mean(Zik), .groups = "drop") %>%
  group_by(sim, sample_size) %>%
  summarise(mu_star = min(mu_k), .groups = "drop")

# gamma dataframe
df_gamma <- df_mu_noise %>%
  left_join(df_mu_star, by = c("sim","sample_size")) %>%
  mutate(denom = mu_noise - mu_star) %>%
  tidyr::crossing(tibble(delta = delta_values)) %>%
  mutate(gamma = delta / denom) %>%
  rename(n = sample_size)





#############################
# Simulation Phase B: compute scores and compare with other methods 
#############################

# Metrics
brier_loss <- function(p, A) sum((p-A)^2)
gini_loss <- function(p, A) sum(A*(1-p) + (1-A)*p)
# log_loss <- function(p, A) sum(-A*log(p) - (1-A)*log(1-p))
log_loss <- function(p, A) {
  t1 <- ifelse(A == 1, -log(p), 0)
  t2 <- ifelse(A == 0, -log(1 - p), 0)
  sum(t1 + t2)
}

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
      gamma_val <- df_gamma %>% filter(sim == !!sim, n == !!n, delta == !!delta) %>% pull(gamma)
      
      # LaD
      for (lad_core in c("LaD_NIW", "LaD_NIG")) {
        df_this <- df_LaD %>%
          filter(sim == !!sim, n == !!n, method == lad_core) %>%
          arrange(draw, k)
        
        if (nrow(df_this) > 0) {
          # S x K mu-samples matrix
          mat_mu <- df_this %>%
            select(draw, k, mu_sample) %>%
            pivot_wider(names_from = k, values_from = mu_sample) %>%
            select(as.character(seq_len(K))) %>%
            as.matrix()
          
          # soft-minimum 
          what_soft <- selection_probabilities(mat_mu, complexities, delta, alpha = alpha_n_soft, mode = "soft")
          score_rows[[length(score_rows)+1]] <- tibble::tibble(
            sim = sim, n = n, delta = delta, gamma = gamma_val,
            method = if (lad_core == "LaD_NIW") "LaD_FC_soft" else "LaD_DC_soft",
            k = seq_len(K), score = what_soft,
            S_brier = brier_loss(what_soft, A_vec),
            S_log = log_loss(what_soft, A_vec),
            S_gini = gini_loss(what_soft, A_vec)
          )
          
          # hard-minimum; only for FC 
          if (lad_core == "LaD_NIW") {
            what_hard <- selection_probabilities(mat_mu, complexities, delta, alpha = NA, mode = "hard")
            score_rows[[length(score_rows)+1]] <- tibble::tibble(
              sim = sim, n = n, delta = delta, gamma = gamma_val,
              method = "LaD_FC_hard",
              k = seq_len(K), score = what_hard,
              S_brier = brier_loss(what_hard, A_vec),
              S_log   = log_loss(what_hard, A_vec),
              S_gini  = gini_loss(what_hard, A_vec)
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
            sim = sim, n = n, delta = delta, gamma = gamma_val,
            method = lbl,
            k = seq_len(K), score = p,
            S_brier = brier_loss(p, A_vec),
            S_log   = log_loss(p, A_vec),
            S_gini  = gini_loss(p, A_vec)
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
          sim = sim, n = n, delta = delta, gamma = gamma_val,
          method = "AIC",
          k = seq_len(K), score = p_aic,
          S_brier = brier_loss(p_aic, A_vec),
          S_log = log_loss(p_aic, A_vec),
          S_gini = gini_loss(p_aic, A_vec)
        )
      }
      if (nrow(df_bic_slice) > 0) {
        k_star <- unique(df_bic_slice$chosen_k)
        p_bic <- as.numeric(seq_len(K) == k_star)
        score_rows[[length(score_rows)+1]] <- tibble::tibble(
          sim = sim, n = n, delta = delta, gamma = gamma_val,
          method = "BIC",
          k = seq_len(K), score = p_bic,
          S_brier = brier_loss(p_bic, A_vec),
          S_log = log_loss(p_bic, A_vec),
          S_gini = gini_loss(p_bic, A_vec)
        )
      }
    } # sim
  } #n
} # delta

df_scores <- bind_rows(score_rows)




cat("Save results...\n")
# save(df_x_all, df_Z_all, df_LaD, df_AICBIC, df_Coarsen, df_gamma, df_scores, file = paste("Output/251105_mvn_sim_seed_", seed, ".Rdata", sep=""))
cat("Results saved successfully.\n")

load("Output/251105_mvn_sim_seed_251102.Rdata")


df_gamma %>% dplyr::filter(sim==25)




# Summarizing LaD plot
labels_delta <- setNames(
  # paste0("delta==", signif(delta_values, 2), "*';'~gamma==", signif(gamma_values, 2)),
  paste0("delta==", round(delta_values, 2)),
  delta_values
)

plot_dat <- df_scores %>%
  filter(n %in% c(50, 500, 5000)) %>% 
  filter(method == "LaD_FC_soft") %>%
  select(sim, n, delta, gamma, k, score) %>%
  mutate(
    sample_size = factor(n, levels = sort(unique(n))),
    model = factor(k, labels=c(1:K)), 
    delta = factor(delta, levels = delta_values, labels = labels_delta)
  )

# highlight_tbl <- plot_dat %>%
#   distinct(sample_size, delta, sim) %>%
#   group_by(sample_size, delta) %>%
#   slice_sample(n = 5) %>%
#   ungroup() %>%
#   mutate(highlight = TRUE)


# highlight_tbl <- plot_dat %>%
#   distinct(sample_size, delta, sim) %>%
#   # filter(sim %in% c(1, 2, 4, 5, 10)) %>%
#   filter(sim %in% c(25)) %>% # 25, 45
#   mutate(highlight = TRUE)
# 
# plot_dat2 <- plot_dat %>%
#   left_join(highlight_tbl, by = c("sample_size","delta","sim")) %>%
#   mutate(highlight = !is.na(highlight))


plot_dat_bar <- plot_dat %>%
  filter(sim == 25)

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

summary_plot


# library(grid)
# 
# g <- ggplotGrob(summary_plot)
# 
# strip_idx <- grep("^strip-l", g$layout$name) # all left strips
# old_strip_cols <- unique(g$layout$l[strip_idx]) # their current column(s)
# strip_col_width <- g$widths[old_strip_cols[1]] # Use the existing strip column width to size the new column
# g <- gtable::gtable_add_cols(g, strip_col_width, pos = 0) # Insert a new column at the extreme left (before everything)
# new_col <- 1  # the new leftmost column index
# 
# # Move strips into the new leftmost column
# g$layout$l[strip_idx] <- new_col
# g$layout$r[strip_idx] <- new_col
# 
# old_strip_cols_after <- old_strip_cols + 1
# g$widths[old_strip_cols_after] <- unit(0, "pt")
# 
# grid.newpage()
# grid.draw(g)



# ggsave("Figures/251019_mvn_sim.png", summary_plot, width = 10, height = 8, dpi = 300, device = "png")




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

p_mu


# Boxplots for mu_{nk} - mu_n^*
df_centered <- df_LaD_NIW %>%
  dplyr::group_by(sim, n, draw) %>%
  dplyr::mutate(mu_star = min(mu_sample)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(mu_diff = mu_sample - mu_star)

delta_lines <- tidyr::expand_grid(n = sample_sizes, delta = delta_values) %>%
  dplyr::mutate(delta_lab = factor(signif(delta, 2), levels = signif(delta_values, 2)))

ymax_center <- max(0,
                   max(df_centered$mu_diff, na.rm = TRUE),
                   max(delta_values, na.rm = TRUE))
ygrid_center <- pretty(c(0, ymax_center), n = 5)

p_qoi <- ggplot(df_centered, aes(x = factor(k), y = mu_diff)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3, color = "black", fill = "white") +
  geom_hline(data = delta_lines,
             aes(yintercept = delta, color = delta_lab),
             linewidth = 0.6, show.legend = FALSE) +
  geom_vline(xintercept = seq(1, K_max, 1), color = "grey70", linewidth = 0.2, alpha = 0.3) +
  geom_hline(yintercept = ygrid_center, color = "grey70", linewidth = 0.2, alpha = 0.3) +
  scale_color_manual(values = c("orange", "purple1", "green3")) +
  scale_y_continuous(limits = range(ygrid_center), breaks = ygrid_center) +
  facet_wrap(~ n, nrow = 1, labeller = labeller(n = function(v) paste0("n = ", v))) +
  labs(x = "k", y = expression(mu[k] - mu^"*")) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 15),
    axis.text  = element_text(size = 13),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5))
  )

p_qoi



fig.mvn <- plot_grid(
  summary_plot, p_mu, 
  ncol = 1, align = "v", axis = "l",
  rel_heights = c(2, 1) 
)
fig.mvn

# ggsave("Figures/251019_p_qoi.png", p_qoi, width = 10, height = 5, dpi = 300)
ggsave("Figures/251111_mvn_result.png", fig.mvn, width = 12, height = 12, dpi = 300)





# Metrics
df_scores_dataset <- df_scores %>%
  dplyr::filter(n %in% c(50, 500, 5000)) %>%
  dplyr::group_by(sim, n, delta, method) %>%
  dplyr::summarise(
    S_brier = dplyr::first(S_brier),
    S_log = dplyr::first(S_log),
    S_gini = dplyr::first(S_gini),
    .groups = "drop"
  )

df_scores_summary <- df_scores_dataset %>%
  dplyr::group_by(n, delta, method) %>%
  dplyr::summarise(
    M = dplyr::n(),
    E_brier = round(mean(S_brier), 3),
    SE_brier = round(sd(S_brier)/sqrt(M), 3),
    E_log = round(mean(S_log), 3),
    SE_log = round(sd(S_log)/sqrt(M), 3),
    E_gini = round(mean(S_gini), 3),
    SE_gini = round(sd(S_gini)/sqrt(M), 3),
    .groups = "drop"
  )


method_levels <- c(
  "LaD-soft", "LaD-hard", "LaD-diag",
  "c-post ($\\alpha = 10$)", "c-post ($\\alpha = 100$)", "Bayes",
  "AIC", "BIC"
)

color_vals <- c(
  "LaD-soft" = "red2",
  "LaD-hard" = "magenta",
  "LaD-diag" = "darkorange1",
  
  "c-post ($\\alpha = 10$)"  = "blue",
  "c-post ($\\alpha = 100$)" = "skyblue2",
  "Bayes" = "cyan",
  
  "AIC" = "seagreen3",
  "BIC" = "green2"
)

# shape_vals <- c(
#   "LaD-soft"              = 1,  # filled circle
#   "LaD-hard"              = 19,  # filled triangle
#   "LaD-diag"              = 21,  # filled square
#   "c-post ($\\alpha = 10$)"  = 2,  
#   "c-post ($\\alpha = 100$)" = 17,  
#   "Bayes"                    = 24,  
#   "AIC"                      = 5,   # filled star
#   "BIC"                      = 23    # filled square cross
# )


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

# hgrid_df <- plot_df %>%
#   group_by(delta) %>%
#   summarise(ymax = max(E_brier + SE_brier, na.rm = TRUE), .groups = "drop") %>%
#   mutate(ymax_half = ceiling(ymax * 2) / 2,                 # round up to next 0.5
#          y = purrr::map(ymax_half, ~ seq(0, .x, by = 0.5))) %>%
#   unnest(y) %>%
#   select(delta, y)

hgrid_df <- plot_df %>% 
  dplyr::group_by(delta) %>% 
  dplyr::summarise(ymin = min(E_brier - SE_brier, na.rm = TRUE), 
                  ymax = max(E_brier + SE_brier, na.rm = TRUE), 
                  .groups = "drop" ) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(y = list(pretty(c(ymin, ymax), n = 5))) %>% 
  tidyr::unnest(y) %>% 
  dplyr::select(delta, y)

pd <- position_dodge(width = 0.7)

vlines_at <- c(1.5, 2.5)

p_metric <- ggplot(plot_df, aes(x = factor(n), y = E_brier, color = method_plot)) +
  geom_vline(xintercept = vlines_at, linetype = "dashed", color = "grey80", linewidth = 0.7) +
  geom_hline(data = hgrid_df, aes(yintercept = y), color = "grey80", linewidth = 0.3, alpha = 0.3) +
  geom_point(position = pd, size = 2.6) +
  geom_errorbar(aes(ymin = E_brier - SE_brier, ymax = E_brier + SE_brier),
                position = pd, width = 0.18, linewidth = 0.5) +
  facet_wrap(
    ~ delta, nrow = 1, scales = "free_y",
    labeller = labeller(delta = as_labeller(labels_delta, default = label_parsed))
  ) +
  scale_color_manual(values = color_vals, limits = method_levels, labels = legend_labels, name = "Method",
                     guide = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3))) +
  # scale_y_continuous(
  #   limits = c(0, NA),            
  #   breaks = scales::breaks_width(0.5),            
  #   labels = function(b) sprintf("%.1f", b)
  # ) +
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


ggsave("Figures/251111_mvn_metrics.png", p_metric, width = 12, height = 5, dpi = 300)




















# S=1000
# mu_NIW  <- niw_posterior(Z, mu0, lambda0, Psi0, nu0, S)$mu
# mu_DIAG <- nig_posterior(Z, mu0, lambda0, a0, b0, S)$mu
#   
# # Correlation and variance of the tie (models 3 vs 4)
# cor_34_NIW  <- cor(mu_NIW[,3],  mu_NIW[,4])
# cor_34_DIAG <- cor(mu_DIAG[,3], mu_DIAG[,4])
#   
# vardiff_NIW  <- var(mu_NIW[,3]  - mu_NIW[,4])
# vardiff_DIAG <- var(mu_DIAG[,3] - mu_DIAG[,4])
#   
# # Global-min distribution and # of survivors under the delta-robustness
# gmin_NIW  <- apply(mu_NIW,  1, min)
# gmin_DIAG <- apply(mu_DIAG, 1, min)
# 
# nsurv <- function(M) rowSums(M <= apply(M, 1, min) + delta)
# ns_NIW  <- nsurv(mu_NIW)
# ns_DIAG <- nsurv(mu_DIAG)
# 
# # Conditional concentration (when the class is chosen and both survive):
# softmin_concentration <- function(M, class_ids, delta, alpha){
#   S <- nrow(M); 
#   K <- ncol(M)
#   min_all <- apply(M,1,min)
#   survive <- M <= (min_all + delta)
#   chosen_class <- apply(survive,1,function(v) {
#     # class of any surviving model with smallest complexity
#     surv_idx <- which(v)
#     if(length(surv_idx)==0) return(NA_integer_)
#     min(complexities[surv_idx])
#   })
#   # rows where desired class is chosen AND both 3 and 4 survive
#   rows <- which(chosen_class == min(complexities[class_ids]) &
#                   apply(survive[,class_ids,drop=FALSE],1,all))
#   if(length(rows)==0) return(NA_real_)
#   # class-relative soft-min weights for models 3 and 4
#   cl <- min(complexities[class_ids])
#   cl_idx <- which(complexities==cl)
#   mu_cl_min <- apply(M[rows, cl_idx, drop=FALSE],1,min)
#   w <- exp(-alpha*(M[rows, class_ids] - mu_cl_min))
#   # concentration: average of max weight share between the tied pair
#   mean(apply(w,1,function(x) max(x)/sum(x)))
# }
# 
# alpha <- n^(0.45)
# conc_NIW  <- softmin_concentration(mu_NIW,  c(3,4), delta, alpha)
# conc_DIAG <- softmin_concentration(mu_DIAG, c(3,4), delta, alpha)
# c(conc_NIW=conc_NIW, conc_DIAG=conc_DIAG)
# 
# conc_overall <- function(M, alpha) {
#   mean( 1 / (1 + exp(-alpha * abs(M[,3] - M[,4]))) )
# }
# conc_overall(mu_NIW,  alpha)
# conc_overall(mu_DIAG, alpha)




# Test alphas
results_list <- list()
alpha_funcs <- list(
  "alpha_n^(0.45)" = function(n) { n^(0.45) },
  "alpha_10*log(n)" = function(n) { 10 * log(n) },
  "alpha_3*sqrt(n)" = function(n) { 3*sqrt(n) },
  "alpha_100*n" = function(n) { 100*n }
)

alpha_levels <- c("alpha_n^(0.45)", "alpha_10*log(n)", "alpha_3*sqrt(n)", "alpha_100*n")

sample_sizes <- c(3, 10, 20, 50, 500, 5000, 50000, 500000, 5000000)

for (sim in 1:3) {
  y_full <- generate_data_mvn(n = 5000000) 
  for (n in sample_sizes) {
    Y <- y_full[1:n, ]
    Z <- fit_models(Y, models)
    mu_samples <- niw_posterior(Z, mu0, lambda0, Psi0, nu0, n_post_samples)$mu  # NIW inference
    for (alpha_name in names(alpha_funcs)) {
      alpha_val <- alpha_funcs[[alpha_name]](n)
      sel_probs <- selection_probabilities(mu_samples, complexities, delta = 0.091, alpha_val)
      
      temp_df <- data.frame(
        sim = sim,
        sample_size = n,
        model = paste0("M", 1:K),
        sel_prob = sel_probs,
        alpha = alpha_name,
        stringsAsFactors = FALSE
      )
      results_list[[length(results_list) + 1]] <- temp_df
    }
  }
}

results_df <- do.call(rbind, results_list) %>%
  filter(model %in% paste0("M", 2:6)) %>%
  mutate(alpha = factor(alpha, levels = alpha_levels))

results_df_sim1 <- filter(results_df, sim == 1)
results_df_others <- filter(results_df, sim != 1)

ggplot() +
  geom_line(data = results_df_others,
            aes(x = sample_size, y = sel_prob, group = interaction(sim, model), color = model),
            linewidth = 0.5, alpha = 0.2) +
  # Thick lines for simulation 1
  geom_line(data = results_df_sim1,
            aes(x = sample_size, y = sel_prob, group = model, color = model),
            linewidth = 1, alpha = 1) +
  facet_wrap(~ alpha, ncol = 2) +
  scale_x_continuous(
    trans = "log10",
    breaks = sample_sizes,
    labels = sample_sizes
  ) +
  theme_bw() +
  labs(
    x = "Sample Size",
    y = "Selection Probability",
    color = "Model"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# Pairs Plot
library(GGally)
library(ggplot2)
ss <- 50
# delta_values <- c(0.04, 0.11, 0.6)
delta_val <- unique(simulation_results_df$delta)[2]

wide_df <- simulation_results_df %>%
  filter(sample_size == ss, delta == delta_val) %>%
  dplyr::select(sim, model, sel_prob) %>%
  pivot_wider(names_from = model, values_from = sel_prob) %>%
  arrange(sim)

my_hist <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) +
    geom_histogram(
      breaks = seq(0, 1, by = 0.1),  
      fill = "steelblue",
      color = "black",
      closed = "right"
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0,1), labels = scales::label_number(accuracy = 1)) +
    theme_minimal()
}

my_scatter <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.5, size = 1, color = "skyblue2") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0,1), labels = scales::label_number(accuracy = 1)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0,1), labels = scales::label_number(accuracy = 1)) +
    theme_minimal()
}

ggpairs(
  wide_df %>% dplyr::select(-sim), 
  diag = list(continuous = wrap(my_hist)),
  lower = list(continuous = wrap(my_scatter)),
  upper = NULL
) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = "Selection Probability",
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "none",
    strip.text = element_text(size = 18)
  )



simulation_results_df %>%
  filter(sample_size==50 & delta == 0.091 & model=="Model_4" & sel_prob < 0.01) %>%
  select(sel_prob)

simulation_results_df %>%
  filter(sample_size==200 & delta == 0.091 & model=="Model_4" & sel_prob < 0.01) %>%
  select(sel_prob)

simulation_results_df %>%
  filter(sample_size==200 & delta == 0.091 & model=="Model_6" & sel_prob > 0.9) %>%
  select(sel_prob)

simulation_results_df %>%
  filter(sample_size==20000 & delta == 0.091 & model=="Model_5" & sel_prob > 0.01) %>%
  select(sel_prob)







ggplot(plot_dat, aes(x = model, y = sel_prob, fill = sample_size)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7, position = position_dodge(width = 0.8)) +
  facet_grid(sample_size ~ delta , scales = "free_y", labeller = custom_labeller) +
  labs(
    x = "Model",
    y = "Selection Probability",
    fill = "Sample Size (n)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "none",
    strip.text = element_text(size = 16)
  )



summary_df <- simulation_results_df %>%
  group_by(sample_size, delta, model) %>%
  summarise(mean_sel_prob = mean(sel_prob),
            se_sel_prob = sd(sel_prob) / sqrt(n()),
            .groups = "drop") %>%
  mutate(delta = factor(delta))

ggplot(summary_df, aes(x = sample_size, y = mean_sel_prob, color = model)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_sel_prob - se_sel_prob, ymax = mean_sel_prob + se_sel_prob), 
                width = 0.1) +
  scale_x_log10() +  # log scale for sample size
  facet_wrap(~ delta) +  # facets by delta; adjust ncol as needed
  labs(x = "Sample Size (log scale)",
       y = "Mean Selection Probability",
       title = "Model Selection Probability by Sample Size and Delta",
       color = "Model") +
  theme_minimal()




KL_for_model <- function(mu_T, S) 0.5 * sum(mu_T[-S]^2)

# Example with your seven models:
mu_T <- c(1, 1, 0.5, 0.5, 0.4, 0)
D1 <- KL_for_model(mu_T, c(1,4))
D2 <- KL_for_model(mu_T, c(1,2))
D3 <- KL_for_model(mu_T, c(1,2,5))
D4 <- KL_for_model(mu_T, c(1,2,3))
D5 <- KL_for_model(mu_T, c(1,2,4))
D6 <- KL_for_model(mu_T, c(1,2,3,4,5))
D7 <- KL_for_model(mu_T, 1:6)
c(D1,D2,D3,D4,D5,D6,D7)

