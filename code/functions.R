# functions.R
# Shared functions for "Robust Model Selection using Likelihood as Data"
#
# Source this file at the top of each analysis script:
#   source(here::here("code", "functions.R"))
#
# Dependencies: MCMCpack, mvtnorm, tibble


#' Draw NIW posterior samples for the LaD mean vector
#'
#' Fits a Normal-Inverse-Wishart conjugate model to the LaD matrix Z and
#' returns posterior draws of mu and Sigma. Implements Step 3 of Algorithm 1 (Section 3.1).
#'
#' @param Z : n x K matrix of per-observation loss values (raw, MLE plug-in, or bias-corrected)
#' @param mu0 : length-K prior mean vector; recommended default is rep(0, K)
#' @param lambda0 : prior precision scaling (positive scalar); recommended default is 0.01
#' @param Psi0 : K x K prior inverse-scale matrix; recommended default is diag(K)
#' @param nu0 : prior degrees of freedom, must exceed K - 1; recommended default is K + 2
#' @param n_samples : number of posterior draws to return
#'
#' @return mu : n_samples x K matrix of posterior draws of mu
#' @return Sigma : K x K x n_samples array of posterior draws of Sigma

niw_posterior <- function(Z, mu0, lambda0, Psi0, nu0, n_samples) {
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
  
  mu_samples <- matrix(NA, nrow = n_samples, ncol = K)
  Sigma_samples <- array(NA, dim = c(K, K, n_samples))
  
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


#' Draw NIG posterior samples for the LaD mean vector (diagonal covariance)
#'
#' Fits K independent Normal-Inverse-Gamma models, one per column of Z.
#' This is the diagonal-covariance variant of the LaD model (Section S5),
#' used for the "LaD-diag" comparison in the sparse MVN simulation (Section 5.2.1).
#'
#' To match the NIW prior in niw_posterior(), set
#' a0 = (nu0 - K + 1) / 2 and b0 = diag(Psi0) / 2.
#'
#' @param Z : n x K matrix of per-observation loss values
#' @param mu0 : length-K prior mean vector
#' @param lambda0 : prior precision scaling (positive scalar)
#' @param a0 : length-K vector of IG shape hyperparameters
#' @param b0 : length-K vector of IG scale hyperparameters
#' @param n_samples : number of posterior draws to return
#'
#' @return mu : n_samples x K matrix of posterior draws of mu
#' @return sigma2 : n_samples x K matrix of posterior draws of the variances

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



#' Compute SLC scores for LaD model selection
#'
#' Implements Algorithm 2 to compute the Smooth LaD Criterion (SLC) score
#' w_delta(k) for each of the K candidate models (Definition 2, Section 2.2).
#' Each score is the product of a between-class selection probability and a
#' within-class soft-minimum weight. Setting mode = "hard" produces the
#' unstable alternative from Theorem 4.4(b), included for comparison only.
#'
#' @param mu_samples : S x K matrix of posterior draws of mu, e.g. the mu element from niw_posterior()
#' @param complexities : integer vector of length K giving the complexity c(k) of each model (smaller = simpler)
#' @param delta : non-negative tolerance scalar; a model is delta-optimal if its mean is within delta of the minimum (Definition 1)
#' @param alpha : temperature parameter alpha_n > 0; should satisfy alpha_n -> Inf and alpha_n = o(sqrt(n)); recommended n^0.45; ignored when mode = "hard"
#' @param mode : "soft" (default, proposed method) or "hard" (unstable alternative for comparison)
#' @param tol : tie tolerance used in mode = "hard", default is 1e-12
#'
#' @return numeric vector of length K with the SLC score for each model

selection_probabilities <- function(mu_samples, complexities, delta, alpha_n,
                                    mode = c("soft", "hard"), tol = 1e-12) {
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
        weight <- exp(-alpha_n * (mu_draw[k] - mu_best_in_class))
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


#' Coarsened posterior model probabilities for sparse MVN mean models
#'
#' Computes model probabilities under the coarsened (power) posterior of
#' Miller & Dunson (2019) for Gaussian models with identity covariance and
#' sparse mean vectors. Used as a comparison method in the sparse MVN
#' simulation (Section 5.2.1, labeled "c-posterior"). As alpha -> Inf,
#' zeta_n -> 1 and this recovers the standard Bayesian posterior.
#'
#' @param x : n x d data matrix
#' @param models : list of length K; each element is an integer vector of free coordinate indices for that model
#' @param alpha : coarsening parameter (positive scalar or Inf for standard Bayes posterior)
#' @param kappa0 : prior precision on free mean coordinates (coarsening prior), default is 0.01
#' @param prior_wts : length-K prior weights over models; uniform if NULL (default)
#' @param m0 : length-d prior mean vector; zero vector if NULL (default)
#'
#' @return tibble with columns: k (model index), prob (posterior probability),
#'   logml (log marginal power likelihood), alpha, zeta_n

coarsen_post_probs <- function(x, models, alpha,
                               kappa0 = 0.01,     
                               prior_wts = NULL,   
                               m0 = NULL) {        
  # Coarsened (power) posterior model probabilities for MVN mean, Sigma = I
  
  n <- nrow(x); d <- ncol(x); K <- length(models)
  zeta_n <- if (is.infinite(alpha)) 1 else alpha / (alpha + n)
  xbar <- colMeans(x)
  if (is.null(m0)) m0 <- rep(0, d)
  
  # coefficients 
  a <- 0.5 * (log(kappa0) - log(kappa0 + zeta_n * n))
  b <- 0.5 * (kappa0 * zeta_n * n) / (kappa0 + zeta_n * n)
  c <- 0.5 * zeta_n * n
  # Equation 138 is same as this mathematically.
  
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
    k = seq_len(K),
    prob = post,
    logml = logml,
    alpha = alpha,
    zeta_n = zeta_n
  )
}