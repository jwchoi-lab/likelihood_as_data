library(matrixStats) 
library(dplyr)
library(purrr)
library(tidyr)
library(gridExtra)
library(ggplot2)

# ----------------------------
# Parse STRUCTURE output
# ----------------------------
parse_structure_results <- function(f, K) {
  lines <- readLines(f, warn = FALSE)
  txt   <- paste(lines, collapse = " ")
  
  # alpha: try "alpha_j", else single "alpha"
  alpha <- NULL
  m1 <- gregexpr("Mean value of alpha_[0-9]+ *= *-?[0-9]+\\.?[0-9]*", txt)
  hits1 <- regmatches(txt, m1)[[1]]
  if (length(hits1) == K) {
    alpha <- as.numeric(sub(".*= *", "", hits1))
  } else {
    m2 <- regexpr("Mean value of alpha *= *-?[0-9]+\\.?[0-9]*", txt)
    if (m2 != -1) {
      a <- as.numeric(sub(".*= *", "", regmatches(txt, m2)))
      alpha <- rep(a, K)
    } else {
      alpha <- rep(NA_real_, K)
    }
  }
  
  # locate locus headers
  locus_starts <- grep("^Locus [0-9]+", lines)
  L <- length(locus_starts)
  
  # per-locus: allele labels (as printed) and per-k frequency vectors
  labels <- vector("list", L)
  theta_by_k <- vector("list", K)
  for (k in seq_len(K)) theta_by_k[[k]] <- vector("list", L)
  
  for (l in seq_len(L)) {
    start <- locus_starts[l] + 2  # skip "X alleles" & "missing data"
    end   <- if (l < L) locus_starts[l + 1] - 1 else length(lines)
    block <- lines[start:end]
    block <- block[nzchar(trimws(block))]
    
    parsed <- lapply(block, function(ln) {
      ln <- gsub("\\([^)]*\\)", "", ln)
      tok <- strsplit(ln, "\\s+")[[1]]
      tok <- tok[nzchar(tok)]
      allele <- suppressWarnings(as.integer(tok[1]))
      if (is.na(allele)) return(NULL)
      freqs <- suppressWarnings(as.numeric(tok[-1]))
      if (any(is.na(freqs))) return(NULL)
      if (K == 1 && length(freqs) >= 1) freqs <- freqs[1]
      if (K > 1 && length(freqs) != K) return(NULL)
      list(allele = allele, freqs = freqs)
    })
    parsed <- Filter(Negate(is.null), parsed)
    if (!length(parsed)) stop("No allele lines parsed at locus ", l)
    
    labels[[l]] <- vapply(parsed, `[[`, integer(1), "allele")
    for (k in seq_len(K)) {
      theta_by_k[[k]][[l]] <- vapply(parsed, function(r) {
        if (K == 1) r$freqs else r$freqs[k]
      }, numeric(1))
    }
  }
  
  # optional summary numbers
  get_num <- function(pattern) {
    m <- regexpr(pattern, txt)
    if (m == -1) return(NA_real_)
    as.numeric(sub(".*= *", "", regmatches(txt, m)))
  }
  
  list(
    alpha = alpha,            # numeric length K
    theta = theta_by_k,       # list length K; each is list length L (numeric vectors)
    labels = labels,          # list length L (allele codes per locus)
    estimated_ln_prob  = get_num("Estimated Ln Prob of Data *= *-?[0-9]+\\.?[0-9]*"),
    mean_ln_likelihood = get_num("Mean value of ln likelihood *= *-?[0-9]+\\.?[0-9]*"),
    var_ln_likelihood  = get_num("Variance of ln likelihood *= *[0-9]+\\.?[0-9]*")
  )
}

find_locus_map <- function(X, labels) {
  L  <- ncol(X) / 2
  Sx <- lapply(seq_len(L), function(d) unique(na.omit(c(X[,2*d-1], X[,2*d]))))
  Sl <- labels
  
  # Jaccard similarity between dataset-locus d and parsed-locus l
  C <- matrix(0, nrow = L, ncol = L)
  for (d in seq_len(L)) for (l in seq_len(L)) {
    a <- Sx[[d]]; b <- Sl[[l]]
    if (length(a) == 0L && length(b) == 0L) { C[d,l] <- 1; next }
    if (length(a) == 0L || length(b) == 0L) { C[d,l] <- 0; next }
    C[d,l] <- length(intersect(a,b))/length(union(a,b))
  }
  
  # greedy assignment (fine for small L). If 'clue' is installed you can swap to LSAP.
  map <- rep(NA_integer_, L); used <- rep(FALSE, L)
  order_d <- order(apply(C, 1, max), decreasing = TRUE)
  for (d in order_d) {
    s <- C[d, ]; s[used] <- -Inf
    l <- which.max(s)
    if (!is.finite(s[l])) l <- which(!used)[1]
    map[d] <- l; used[l] <- TRUE
  }
  map
}


compute_point_nll <- function(X, alpha, theta, labels, S = 1000, eps = 1e-12, auto_map = TRUE) {
  n <- nrow(X); P <- ncol(X); L <- P/2; K <- length(theta)
  
  # optional locus permutation to match X
  loc_map <- if (auto_map) find_locus_map(X, labels) else seq_len(L)
  labels <- labels[loc_map]
  theta <- lapply(theta, function(thk) thk[loc_map])
  
  idx_of <- function(l, v) match(v, labels[[l]])  # NA if allele unseen at locus l
  
  # nll <- numeric(n)
  nll <- rep(0, n)
  
  if (K == 1) {
    th1 <- theta[[1]]  # list length L of numeric vectors (probs aligned to labels[[l]])
    for (i in seq_len(n)) {
      xi <- as.integer(X[i, ])
      obs <- which(!is.na(xi))
      ll <- 0.0
      for (j in obs) {
        l <- ceiling(j/2)
        v <- xi[j]
        idx <- idx_of(l, v)
        p <- if (is.na(idx)) eps else max(th1[[l]][idx], eps)
        ll <- ll + log(p)
      }
      nll[i] <- -ll
    }
    return(data.frame(ind = seq_len(n), nll = nll))
  }
  
  # K > 1: Monte Carlo Dirichlet integral
  alpha <- as.numeric(alpha)
  for (i in seq_len(n)) {
    xi  <- as.integer(X[i, ])
    obs <- which(!is.na(xi))
    
    U <- matrix(rgamma(S * K, shape = alpha, rate = 1), nrow = S, ncol = K, byrow = TRUE)
    W <- U / rowSums(U)
    logW <- log(W)
    
    logA <- matrix(0, nrow = S, ncol = P)
    for (j in obs) {
      l   <- ceiling(j/2)
      v   <- xi[j]
      idx <- idx_of(l, v)
      logTheta <- vapply(seq_len(K), function(k) {
        if (is.na(idx)) log(eps) else log( max(theta[[k]][[l]][idx], eps) )
      }, numeric(1))
      logA[, j] <- rowLogSumExps(logW + matrix(logTheta, nrow = S, ncol = K, byrow = TRUE))
    }
    
    logp_samps <- rowSums(logA, na.rm = TRUE)
    nll[i] <- -(logSumExp(logp_samps) - log(S))
  }
  
  data.frame(ind = seq_len(n), nll = nll)
}

# Function to draw posterior samples of mu and Sigma using NIW conjugacy
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


# Selection probaiblity function
selection_probabilities <- function(mu_samples, complexities, delta, alpha_n) {
  # mu_samples: S x K matrix of posterior draws for mu_k (e.g., negative log-likelihood)
  # complexities: integer vector of length K specifying the complexity class for each model
  # delta: numeric tolerance for "surviving" from the global best
  # alpha_n: temperature parameter (positive)
  
  S <- nrow(mu_samples)  # number of posterior draws
  K <- ncol(mu_samples)  # number of models
  
  class_selection_counts <- numeric(K) # Class selection indicator counts
  weight_sums <- numeric(K) # Sum of performance weights
  
  for (s in seq_len(S)) {
    mu_draw <- mu_samples[s, ]
    
    # Step 1: Determine the chosen complexity class for this draw
    overall_best <- min(mu_draw)
    surviving <- which(mu_draw <= overall_best + delta)
    surviving_classes <- complexities[surviving]
    chosen_class <- min(surviving_classes)
    
    selected_models <- which(complexities == chosen_class)
    class_selection_counts[selected_models] <- class_selection_counts[selected_models] + 1
    
    # Step 2: Compute the weight.
    for (k in seq_len(K)) {
      class_k <- complexities[k]
      models_in_class <- which(complexities == class_k)
      mu_best_in_class <- min(mu_draw[models_in_class])
      weight <- exp(-alpha_n * (mu_draw[k] - mu_best_in_class))
      weight_sums[k] <- weight_sums[k] + weight
    }
  }
  
  # Compute the Monte Carlo estimates by averaging over S draws:
  p_class <- class_selection_counts / S
  avg_weight <- weight_sums / S
  
  # Overall model selection probability: product of the two components.
  final_probs <- p_class * avg_weight
  return(final_probs)
}

seed <- 251111
set.seed(seed)

# data
# dat <- read.table("Data/impala.txt")
dat <- read.table("Data/brooktrout.txt")
X <- dat
# X_mat <- as.matrix(dat)

grid <- expand.grid(
  K = 1:10,
  rep = 1:20
) %>%
  as_tibble() %>%
  mutate(
    file = sprintf("Structure/250914_output/results_K%d_rep%d_f", K, rep)
  )

results <- grid %>%
  mutate(
    parsed = map2(file, K, parse_structure_results),
    # mean_ll = map_dbl(parsed, "mean_ln_likelihood")
    mean_ll = map_dbl(parsed, "estimated_ln_prob")
  ) %>%
  dplyr::select(K, rep, mean_ll)

best_reps <- results %>%
  group_by(K) %>%
  slice_max(mean_ll, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(file = sprintf("Structure/250914_output/results_K%d_rep%d_f", K, rep))

final_rep_nll <- best_reps %>%
  mutate(parsed = map2(file, K, parse_structure_results)) %>%
  mutate(nlldf  = map(parsed, ~ compute_point_nll(X, .x$alpha, .x$theta, labels = .x$labels))) %>%
  dplyr::select(K, rep, nlldf) %>%
  unnest(nlldf)

final_nll <- final_rep_nll %>% dplyr::select(-rep)



## Apply LaD
nll_matrix <- final_nll %>%
  # filter(K %in% c(1:10)) %>%
  pivot_wider(names_from = K, values_from = nll) %>%
  arrange(ind)

Z <- as.matrix(nll_matrix[,-1]) # remove index column 


K <- length(colnames(Z))
n <-nrow(Z)
complexities <- c(1,2,3,4,5,6,7,8,9,10)

L <- ncol(X) / 2
V_l <- sapply(seq_len(L), function(l) {
  x1 <- X[, 2*l - 1]; x2 <- X[, 2*l]
  length(unique(na.omit(c(x1, x2))))
})

Ks <- seq_len(ncol(Z))

d_k <- Ks * sum(V_l - 1) + Ks

Z  <- t(t(Z) + d_k/(2*n))  # bias-corrected. 



# Noise model
M <- (!is.na(X[, seq(1, 2*L, 2)])) + (!is.na(X[, seq(2, 2*L, 2)]))
nll_noise <- as.vector(M %*% log(V_l))
mu_noise_uniform <- mean(nll_noise)
print(mu_noise_uniform )

# mu_noise <- as.numeric(colMeans(Z)[1])
mu_k <- as.numeric(colMeans(Z))
print(mu_k)

# NIW Prior Parameters
mu0 <- rep(0, K)  # Prior mean (vector of zeros).
lambda0 <- 0.01 # Prior scaling for the mean.
nu0 <- K + 2  # Degrees of freedom (must be > K-1).
Psi0 <- 1 * diag(K)  # Prior inverse scale matrix.

n_post_samples <- 1000
mu_samples <- niw_posterior(Z, mu0, lambda0, Psi0, nu0, n_post_samples)$mu # NIW bayesian inference




# Step plot
denom <- (mu_noise_uniform - min(mu_k)) # total recoverable info
tau_grid <- seq(0, 1, length.out = 101)
delta_grid <- denom * tau_grid
bin_width  <- min(diff(tau_grid))

a_n <- n^(0.45) # no meaning for this example

res <- purrr::map2_dfr(tau_grid, delta_grid, ~{
  xleft  <- pmax(.x - bin_width/2, 0)
  xright <- pmin(.x + bin_width/2, 1)
  tibble(
    n = n, tau = .x, K=seq_len(K), complexities = complexities, 
    sel_prob = as.numeric( selection_probabilities(mu_samples, complexities, delta = .y, alpha = a_n) ),
    xleft = xleft, xright = xright, xmid = (xleft + xright)/2
  )
})


heat_cols <- c("white", "orange1", "orange2", "red1", "red3")

p_heat <- ggplot() +
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
  guides(
    fill = guide_colorbar(
      barwidth = unit(5, "cm"),
      barheight = unit(0.5, "cm"),
      title.position = "left",
      title.hjust = 0.5,
      title.vjust = 0.9
    )) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 15),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5)),
    legend.position = "bottom",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.justification = c(0.5, 0),
    legend.box.margin = margin(t = -10),
    plot.margin = margin(5, 5, 5, 5) 
  )  

p_heat

# ggsave("Figures/250915_brooktrout_lad_pp.png", p_heat, width = 6, height = 4, dpi = 300, device = "png")




# Mu plot


S <- nrow(mu_samples); K <- ncol(mu_samples)
colnames(mu_samples) <- paste0("k", seq_len(K))

df_mu <- as.data.frame(mu_samples)
df_mu$draw <- seq_len(S)

df_mu_long <- df_mu %>%
  pivot_longer(
    cols = -draw,
    names_to = "k_lab",
    values_to = "mu_sample"
  ) %>%
  mutate(
    # numeric k index matching the column order of mu_samples
    k = as.integer(factor(k_lab, levels = colnames(mu_samples)))
  )

K_max <- K
ygrid_mu <- pretty(range(df_mu_long$mu_sample, na.rm = TRUE), n = 6)

p_mu <- ggplot(df_mu_long, aes(x = factor(k), y = mu_sample)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3, color = "black", fill = "white") +
  geom_vline(xintercept = seq(1, K_max, by = 1), color = "grey70", linewidth = 0.2, alpha = 0.3) +
  geom_hline(yintercept = ygrid_mu, color = "grey70", linewidth = 0.2, alpha = 0.3) +
  scale_y_continuous(limits = range(ygrid_mu), breaks = ygrid_mu) +
  labs(x = "k", y = expression(mu[k])) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 15),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5))
  )
p_mu



library(cowplot)
fig.ad <- plot_grid(
  p_mu, p_heat,
  ncol = 2, align = "v", axis = "l",
  rel_widths = c(1, 1.2) 
)
fig.ad






# Evanno's method

# Panel A: L(K) mean over replicates 
LK_stats <- results %>%
  group_by(K) %>%
  summarise(
    LK_mean = mean(mean_ll, na.rm = TRUE),
    LK_sd   = sd(mean_ll, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(K)

# Per-replicate finite differences
by_rep <- results %>%
  group_by(rep) %>%
  arrange(K, .by_group = TRUE) %>%
  mutate(
    Lprime = mean_ll - lag(mean_ll),  # L'(K)
    L2 = lead(mean_ll) - 2*mean_ll + lag(mean_ll) # L''(K)
  ) %>%
  ungroup()

# Panel B: L'(K) mean across replicates
Lprime_stats <- by_rep %>%
  group_by(K) %>%
  summarise(
    Lp_mean = mean(Lprime, na.rm = TRUE),
    Lp_sd = sd(Lprime, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(K)

# Panel C: |L''(K)| mean across replicates
L2_stats <- by_rep %>%
  group_by(K) %>%
  summarise(
    L2_mean = mean(abs(L2), na.rm = TRUE),
    L2_sd   = sd(abs(L2),   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(K)

# Panel D: delta K = mean(|L''(K)|) / sd(L(K)) 
evanno_df <- LK_stats %>%
  left_join(L2_stats, by = "K") %>%
  mutate(
    DeltaK = ifelse(LK_sd > 0, L2_mean / LK_sd, NA)     # undefined if sd = 0
  )

ygrid_LK <- pretty(range(LK_stats$LK_mean, na.rm = TRUE), n = 6)
ygrid_Lprime <- pretty(range(Lprime_stats$Lp_mean, na.rm = TRUE), n = 6)
ygrid_L2 <- pretty(range(L2_stats$L2_mean, na.rm = TRUE), n = 6)
ygrid_DeltaK <- pretty(range(evanno_df$DeltaK, na.rm = TRUE), n = 6)


pA <- ggplot(LK_stats, aes(K, LK_mean)) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin = LK_mean - LK_sd, ymax = LK_mean + LK_sd), width = 0.2) +
  geom_vline(xintercept = seq(1, K_max, by = 1), color = "grey70", linewidth = 0.2, alpha = 0.3) +
  geom_hline(yintercept = ygrid_LK, color = "grey70", linewidth = 0.2, alpha = 0.3) +
  labs(title = expression("Mean " * L(k)), x = "k", y = "") +
  scale_x_continuous(breaks = LK_stats$K) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size=18, hjust = 0.5),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black", size = 18),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 15),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5))
  ) + theme(plot.margin = margin(10, 12, 10, 12))

pB <- ggplot(Lprime_stats, aes(K, Lp_mean)) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin = Lp_mean - Lp_sd, ymax = Lp_mean + Lp_sd), width = 0.2) +
  geom_vline(xintercept = seq(1, K_max, by = 1), color = "grey70", linewidth = 0.2, alpha = 0.3) +
  geom_hline(yintercept = ygrid_Lprime, color = "grey70", linewidth = 0.2, alpha = 0.3) +
  labs(title = expression(L*"\u2032"(k)), x = "k", y = "") +
  scale_x_continuous(breaks = Lprime_stats$K) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size=18, hjust = 0.5),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black", size = 18),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 15),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5))
  ) + theme(plot.margin = margin(10, 12, 10, 12))

pC <- ggplot(L2_stats, aes(K, L2_mean)) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin = L2_mean - L2_sd, ymax = L2_mean + L2_sd), width = 0.2) +
  geom_vline(xintercept = seq(1, K_max, by = 1), color = "grey70", linewidth = 0.2, alpha = 0.3) +
  geom_hline(yintercept = ygrid_L2, color = "grey70", linewidth = 0.2, alpha = 0.3) +
  labs(title = expression("|" * L*"\u2033"(k) * "|"), x = "k", y = "") +
  scale_x_continuous(breaks = L2_stats$K) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size=15, hjust = 0.5),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black", size = 18),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5))
  ) + theme(plot.margin = margin(10, 12, 10, 12))

pD <- ggplot(evanno_df, aes(K, DeltaK)) +
  geom_point() + geom_line() +
  labs(title = expression(Delta * k), x = "k", y = "") +
  geom_vline(xintercept = seq(1, K_max, by = 1), color = "grey70", linewidth = 0.2, alpha = 0.3) +
  geom_hline(yintercept = ygrid_DeltaK, color = "grey70", linewidth = 0.2, alpha = 0.3) +
  scale_x_continuous(breaks = evanno_df$K) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size=18, hjust = 0.5),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black", size = 18),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 15),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5))
  ) + theme(plot.margin = margin(10, 12, 10, 12)) 

evanno_plot <- grid.arrange(pA, pB, pD, ncol = 3)




ggsave("Figures/260225_brooktrout_result.png", fig.ad, width = 15, height = 5, dpi = 300)

ggsave("Figures/251111_brooktrout_evanno.png", evanno_plot, width = 15, height = 5, dpi = 300, device = "png")


save.image("admixture_251111.RData") 


