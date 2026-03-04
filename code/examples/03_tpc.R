# 03_tpc.R
# Thermal performance curves in microbial ecology (Section 5.3)

library(tidyverse)
library(ggplot2)
library(nls.multstart)
library(minpack.lm)
library(purrr)
library(cowplot)
library(here)

source(here::here("code", "functions.R"))


# Read data
dat.all <- read_csv(
  here::here("data", "raw", "tpc", "thermal_performance_datasets.csv")
)

# dat.all %>%
#   filter(trait %in% c("Growth rate", "growth rate")) %>%
#   count(taxon, sort = TRUE) %>%
#   slice_max(n, n = 10) 

pm_dat <- dat.all %>%
  filter(taxon == "Prorocentrum minimum",
         trait == "growth rate") %>%
  filter(str_detect(citation, "Loma")==FALSE) %>%
  dplyr::select(temperature, trait_value) %>%
  rename(temp = temperature) 

# ggplot(pm_dat, aes(x = temp, y = trait_value)) +
#   geom_point() +
#   geom_smooth(se = FALSE, method = "loess") +
#   labs(title = "pm minor: Population Growth vs Temperature",
#        x = "Temperature (°C)",
#        y = "Population Growth rate")




# For exact replication of the figures in the paper, we provide precomputed results in
# output/tpc_fitted.Rdata:
# load(here::here("output", "tpc_fitted.Rdata"))

# Otherwise, the scripts below also allow recomputing the results from the raw data;
# these may differ slightly up to random sampling, but yield the same qualitative conclusions.




set.seed(251111)


# Models
fit_Eubank_3_pars <- function(dataset) {
  a_start <- 300
  b_start <- 50
  T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
  function_to_be_fitted <- function(a, T_pk, b, temp){
    result <- a / ((temp - T_pk)^2 + b)
    if ( any(result <= 0) ) {
      return(rep(1e10, length(temp)))
    } else {
      return(result)
    }
  }
  fit <- NULL
  try(
    fit <- nls_multstart(
      trait_value ~ function_to_be_fitted(a, T_pk, b, temp=temp),
      data=dataset, iter=1000,
      start_lower=c(a=0.5*a_start, T_pk=0.5*T_pk_start, b=0.5*b_start),
      start_upper=c(a=1.5*a_start, T_pk=1.5*T_pk_start, b=1.5*b_start),
      supp_errors='Y', convergence_count=FALSE,
      control = nls.lm.control(
        ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
        maxfev = 100000
      ),
      lower=c(0, 0, 0), upper=c(Inf, 150, Inf)
    )
  )
  return(fit)
}

fit_Gaussian_3_pars <- function(dataset) {
  B_pk_start <- max(dataset$trait_value)
  T_pk_start <- max(dataset$temp[dataset$trait_value==B_pk_start])
  a_start    <- 90
  function_to_be_fitted <- function(B_pk, T_pk, a, temp) {
    B_pk * exp(-0.5 * (abs(temp - T_pk)/a)^2)
  }
  fit <- NULL
  try(
    fit <- nls_multstart(
      trait_value ~ function_to_be_fitted(B_pk, T_pk, a, temp=temp),
      data=dataset, iter=1000,
      start_lower = c(B_pk=0.5*B_pk_start, T_pk=0.5*T_pk_start, a=0.5*a_start),
      start_upper = c(B_pk=1.5*B_pk_start, T_pk=1.5*T_pk_start, a=1.5*a_start),
      supp_errors = 'Y', convergence_count = FALSE,
      control = nls.lm.control(ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, maxfev = 100000),
      lower = c(0,0,0), upper = c(Inf,150,Inf)
    )
  )
  return(fit)
}

fit_Mitchell_Angilletta_3_pars <- function(dataset) {
  a_start <- 1000
  b_start <- max(dataset$temp) - min(dataset$temp)
  T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
  function_to_be_fitted <- function(a, b, T_pk, temp) {
    result <- ( 1 / (2 * b) ) * ( 1 + cos( ( (temp - T_pk) / b ) * pi ) ) * a
    # Make sure that all values are positive and that the resulting  curve is not multimodal.
    if ( any(result <= 0 ) || any( ((temp - T_pk)/b) > 2 ) || any( ((temp - T_pk)/b) < -2 ) ) {
      return(rep(1e10, length(temp)))
    } else {
      return(result)
    }
  }
  fit <- NULL
  try(
    fit <- nls_multstart(
      trait_value ~ function_to_be_fitted(a, b, T_pk, temp=temp),
      data=dataset, iter=1000,
      start_lower = c(a=0.5*a_start, b=0.5*b_start, T_pk=0.5*T_pk_start),
      start_upper = c(a=1.5*a_start, b=1.5*b_start, T_pk=1.5*T_pk_start),
      supp_errors = 'Y', convergence_count = FALSE,
      control = nls.lm.control(ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, maxfev = 100000),
      lower = c(0,0,0), upper = c(Inf,150,150)
    )
  )
  return(fit)
}

fit_Briere_I_3_pars <- function(dataset) {
  a_start <- 1
  Tmin_start <- min(dataset$temp) - 1
  Tmax_start <- max(dataset$temp) + 1
  function_to_be_fitted <- function(a, Tmin, Tmax, temp){
    result <- a * temp * (temp - Tmin) * sqrt(Tmax - temp)
    if ( any(result <= 0) || any(is.nan(result)) ) {
      return(rep(1e10, length(temp)))
    } else {
      return(result)
    }
  }
  fit <- NULL
  try(
    fit <- nls_multstart(
      trait_value ~ function_to_be_fitted(a, Tmin, Tmax, temp=temp),
      data=dataset, iter=1000,
      start_lower=c(a=0.5*a_start, Tmin=0.5*Tmin_start, Tmax=0.5*Tmax_start),
      start_upper=c(a=1.5*a_start, Tmin=1.5*Tmin_start, Tmax=1.5*Tmax_start),
      supp_errors='Y', convergence_count=FALSE,
      control = nls.lm.control(
        ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
        maxfev = 100000
      ),
      lower=c(0, -20, 0), upper=c(Inf, 150, 150)
    )
  )
  return(fit)
}


fit_Weibull_4_pars <- function(dataset) {
  B_pk_start <- max(dataset$trait_value)
  T_pk_start <- max(dataset$temp[dataset$trait_value == B_pk_start])
  d_start <- 90
  e_start <- 2
  
  function_to_be_fitted <- function(B_pk, T_pk, d, e, temp) {
    result <- B_pk *
      (( (e - 1)/e )^((1 - e)/e)) *
      (((temp - T_pk)/d + ((e - 1)/e)^(1/e))^(e - 1)) *
      exp(-((temp - T_pk)/d + ((e - 1)/e)^(1/e))^e + ( (e - 1)/e ))
    if (any(is.na(result)) || any(result <= 0)) {
      return(rep(1e10, length(temp)))
    } else {
      return(result)
    }
  }
  
  fit <- NULL
  try(
    fit <- nls_multstart(
      trait_value ~ function_to_be_fitted(B_pk, T_pk, d, e, temp=temp),
      data = dataset,
      iter = 1000,
      start_lower = c(
        B_pk = 0.5 * B_pk_start, T_pk = 0.5 * T_pk_start,
        d = 0.5 * d_start, e = 0.5 * e_start
      ),
      start_upper = c(
        B_pk = 1.5 * B_pk_start, T_pk = 1.5 * T_pk_start,
        d = 1.5 * d_start, e = 1.5 * e_start
      ),
      supp_errors = 'Y',
      convergence_count = FALSE,
      control = nls.lm.control(
        ftol = .Machine$double.eps, ptol = .Machine$double.eps,
        maxiter = 1024, maxfev = 100000
      ),
      lower = c(0, 0, -Inf, -Inf),
      upper = c(Inf, 150, Inf, Inf)
    )
  )
  return(fit)
}


fit_modified_Gaussian_4_pars <- function(dataset) {
  B_pk_start <- max(dataset$trait_value)
  T_pk_start <- dataset$temp[which.max(dataset$trait_value)]
  a_start    <- 90
  b_start    <- 2
  
  function_to_be_fitted <- function(B_pk, T_pk, a, b, temp) {
    result <- B_pk * exp(-0.5 * (abs(temp - T_pk)/a)^b)
    if (any(is.na(result)) || any(result <= 0)) {
      return(rep(1e10, length(temp)))
    } else {
      return(result)
    }
  }
  
  fit <- NULL
  try(
    fit <- nls_multstart(
      trait_value ~ function_to_be_fitted(B_pk, T_pk, a, b, temp=temp),
      data = dataset,
      iter = 1000,
      start_lower = c(B_pk = 0.5 * B_pk_start, T_pk = 0.5 * T_pk_start,
                      a = 0.5 * a_start, b = 0.5 * b_start),
      start_upper = c(B_pk = 1.5 * B_pk_start, T_pk = 1.5 * T_pk_start,
                      a = 1.5 * a_start, b = 1.5 * b_start),
      lower = c(0, 0, 0, 0),
      upper = c(Inf, 150, Inf, Inf),
      supp_errors = 'Y',
      convergence_count = FALSE,
      control = nls.lm.control(
        ftol = .Machine$double.eps,
        ptol = .Machine$double.eps,
        maxiter = 1024,
        maxfev = 100000
      )
    )
  )
  return(fit)
}



fit_extended_Briere_5_pars <- function(dataset) {
  a_start      <- 1
  T_min_start  <- min(dataset$temp) - 1
  T_max_start  <- max(dataset$temp) + 1
  b_start      <- 1
  d_start      <- 1.5
  
  function_to_be_fitted <- function(a, T_min, T_max, b, d, temp) {
    if (T_min >= T_max || any(temp < T_min) || any(temp > T_max)) {
      return(rep(1e10, length(temp)))
    } else {
      result <- a * temp * ((temp - T_min)^b) * ((T_max - temp)^d)
      if (any(is.na(result)) || any(result <= 0)) {
        return(rep(1e10, length(temp)))
      } else {
        return(result)
      }
    }
  }
  
  fit <- NULL
  try(
    fit <- nls_multstart(
      trait_value ~ function_to_be_fitted(a, T_min, T_max, b, d, temp = temp),
      data = dataset,
      iter = 1000,
      start_lower = c(
        a = 0.5 * a_start, T_min = 0.5 * T_min_start,
        T_max = 0.5 * T_max_start, b = 0.5 * b_start, d = 0.5 * d_start
      ),
      start_upper = c(
        a = 1.5 * a_start, T_min = 1.5 * T_min_start,
        T_max = 1.5 * T_max_start, b = 1.5 * b_start, d = 1.5 * d_start
      ),
      lower = c(0, -20, 0, 0, 0),
      upper = c(Inf, 150, 150, Inf, Inf),
      control = nls.lm.control(
        ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
        maxfev = 100000
      ),
      supp_errors = 'Y',
      convergence_count = FALSE
    )
  )
  
  return(fit)
}



fit_poly5_6_pars <- function(dataset) {
  a_start <- min(dataset$trait_value)
  b_start <- 0.1
  c_start <- 0.003
  d_start <- 2e-7
  e_start <- 1e-10
  g_start <- 1e-10
  
  f <- function(a, b, c, d, e, g, temp) {
    res <- a + b*temp + c*temp^2 + d*temp^3 + e*temp^4 + g*temp^5
    if(any(is.na(res)) || any(res <= 0)) rep(1e10, length(res)) else res
  }
  
  fit <- NULL
  
  try(
    nls_multstart(
      trait_value ~ f(a, b, c, d, e, g, temp=temp),
      data = dataset, iter = 1000,
      start_lower = c(a=0.5*a_start, b=0.5*b_start, c=0.5*c_start,
                      d=0.5*d_start, e=0.5*e_start, g=0.5*g_start),
      start_upper = c(a=1.5*a_start, b=1.5*b_start, c=1.5*c_start,
                      d=1.5*d_start, e=1.5*e_start, g=1.5*g_start),
      lower = rep(-Inf,6), upper = rep(Inf,6),
      supp_errors='Y', convergence_count=FALSE,
      control = nls.lm.control(ftol = .Machine$double.eps,
                               ptol = .Machine$double.eps,
                               maxiter = 1024, maxfev = 100000)
    )
  )
}


fit_Sharpe_Schoolfield_7_pars <- function(dataset) {
  R <- 1.987
  Tref <- 273.15
  
  B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])
  T_pk <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
  H_A_start <- 15000
  H_L_start <- 70000
  H_H_start <- 50000
  T_L50_start <- min(dataset$temp) + 273.15
  T_H50_start <- T_pk + 273.15
  a_start <- 1
  
  f <- function(B_0, H_A, H_L, T_L50, H_H, T_H50, a, temp) {
    TK <- temp + 273.15
    num <- B_0 * TK / Tref * exp(H_A / R * (1 / Tref - 1 / TK))
    den <- a + exp(H_L / R * (1 / T_L50 - 1 / TK)) + exp(H_H / R * (1 / T_H50 - 1 / TK))
    res <- num / den
    if (any(is.na(res)) || any(res <= 0)) rep(1e10, length(res)) else res
  }
  
  fit <- NULL
  try(
    fit <- nls_multstart(
      trait_value ~ f(B_0, H_A, H_L, T_L50, H_H, T_H50, a, temp = temp),
      data = dataset,
      iter = 1000,
      start_lower = c(
        B_0 = 0.5 * B_0_start, H_A = 0.5 * H_A_start, H_L = 0.5 * H_L_start,
        T_L50 = 0.5 * T_L50_start, H_H = 0.5 * H_H_start, T_H50 = 0.5 * T_H50_start,
        a = 0.5 * a_start
      ),
      start_upper = c(
        B_0 = 1.5 * B_0_start, H_A = 1.5 * H_A_start, H_L = 1.5 * H_L_start,
        T_L50 = 1.5 * T_L50_start, H_H = 1.5 * H_H_start, T_H50 = 1.5 * T_H50_start,
        a = 1.5 * a_start
      ),
      lower = c(0, 0, 0, 273.15, 0, 273.15, 0),
      upper = c(Inf, Inf, Inf, 273.15 + 150, Inf, 273.15 + 150, Inf),
      supp_errors = "Y",
      convergence_count = FALSE,
      control = nls.lm.control(
        ftol = .Machine$double.eps,
        ptol = .Machine$double.eps,
        maxiter = 1024,
        maxfev = 100000
      )
    )
  )
  return(fit)
}



# Model fitting
fit_list <- list(
  eubank   = fit_Eubank_3_pars(pm_dat),
  gaussian = fit_Gaussian_3_pars(pm_dat),
  mitchell = fit_Mitchell_Angilletta_3_pars(pm_dat),
  briere1  = fit_Briere_I_3_pars(pm_dat),
  weibull  = fit_Weibull_4_pars(pm_dat),
  m_gaussian = fit_modified_Gaussian_4_pars(pm_dat),
  e_briere = fit_extended_Briere_5_pars(pm_dat),
  poly5 = fit_poly5_6_pars(pm_dat),
  schoolfield = fit_Sharpe_Schoolfield_7_pars(pm_dat)
)

names(fit_list) <- c(
  "Eubank","Gaussian","Mitchell–Angilletta",
  "Brière I","Weibull","Modified Gaussian",
  "Extended Brière","Poly-5","Sharpe–Schoolfield"
)


# save(fit_list, file = paste("output/tpc_fitted.Rdata"))


model_labels <- c(
  "1. Eubank",
  "2. Gaussian",
  "3. Mitchell–Angilletta",
  "4. Brière I",
  "5. Weibull",
  "6. Modified Gaussian",
  "7. Extended Brière",
  "8. Poly-5",
  "9. Sharpe–Schoolfield"
)


resid_df <- purrr::imap_dfr(fit_list, function(fit, model) {
  if (is.null(fit)) return(NULL)
  tibble(
    model = as.factor(model),
    temp = pm_dat$temp,
    observed = pm_dat$trait_value,
    fitted = predict(fit, newdata = pm_dat),
    residual = observed - fitted
  )
})

curve_df <- resid_df %>%
  dplyr::select(model, temp, fitted) %>%
  distinct()


# Convert s^-1 → divisions per day: k = (μ_s × 86400) / ln(2)
pm_dat_plot   <- pm_dat  %>% mutate(trait_value_divd = trait_value * 86400 / log(2))
curve_df_plot <- curve_df %>% mutate(fitted_divd = fitted * 86400 / log(2))
curve_df_plot$model <- factor(curve_df_plot$model,
                              levels = names(fit_list),
                              labels = model_labels)


tpc_fit <- ggplot(pm_dat_plot, aes(x = temp, y = trait_value_divd)) +
  geom_point(color = "black", size=0.5) +
  geom_line(data = curve_df_plot, aes(x = temp, y = fitted_divd), 
            color = "steelblue3", size = 1) +   
  facet_wrap(~ model, scales = "free") +
  labs(
    x = "Temperature (°C)",
    y = expression("Growth rate"~(div~d^{-1}))
  ) +
  theme_bw(base_size = 20) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black", size = 17),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 18),
    legend.position = "none"
  )




# Apply to our method
compute_LaD <- function(dataset, fit_list) {
  n <- nrow(dataset)
  K <- length(fit_list)
  Z <- matrix(NA, nrow=n, ncol=K)
  colnames(Z) <- names(fit_list)
  for (k in seq_along(fit_list)) {
    model_name <- names(fit_list)[k]
    fit <- fit_list[[k]]
    if (!is.null(fit)) {
      preds <- predict(fit, newdata=dataset)
      resids <- dataset$trait_value - preds
      sigma2 <- mean(resids^2, na.rm=TRUE)
      Z[,k] <- 0.5*log(2*pi*sigma2) + (resids^2)/(2*sigma2)
    }
  }
  return(Z)
}
Z <- compute_LaD(pm_dat, fit_list)


complexities <- c(3,3,3,3,4,4,5,6,7)
df_k <- complexities + 1
n <- nrow(pm_dat)

Z_bc <- t(t(Z) + df_k/(2*n)) # Bias-adjusted


# Noise (baseline) model
compute_noise_model_Z <- function(y) {
  mu_hat <- mean(y)
  sigma2_hat <- mean((y - mu_hat)^2)
  Z_noise <- 0.5 * log(2 * pi * sigma2_hat) + ((y - mu_hat)^2) / (2 * sigma2_hat)
  return(Z_noise)
}

Z_noise <- compute_noise_model_Z(pm_dat$trait_value)
mu_noise <- mean(Z_noise) # should be strictly worse (greater) than mu_ks

# Proxy (best possible) model
proxy_fit <- pm_dat %>%
  group_by(temp) %>%
  mutate(mu_hat = mean(trait_value)) %>%
  ungroup()

sigma2_hat <- with(proxy_fit, mean((trait_value - mu_hat)^2))

proxy_fit <- proxy_fit %>%
  mutate(g_proxy = 0.5 * (log(2*pi*sigma2_hat) + (trait_value - mu_hat)^2 / sigma2_hat))

mu_flex <- mean(proxy_fit$g_proxy)


K <- length(colnames(Z_bc))
# NIW Prior Parameters
mu0 <- rep(0, K)       # Prior mean (vector of zeros).
lambda0 <- 0.01        # Prior scaling for the mean.
nu0 <- K + 2           # Degrees of freedom (must be > K-1).
Psi0 <- 1 * diag(K)  # Prior inverse scale matrix.

n_post_samples <- 1000
mu_samples <- niw_posterior(Z_bc, mu0, lambda0, Psi0, nu0, n_post_samples)$mu # NIW bayesian inference



# Choosing tau and delta
tau <- c(0.25, 0.1, 0.05, 0.001)  
delta_values <- tau * (mu_noise - mu_flex)

all_results <- list()
for (i in seq_along(delta_values)) {
  delta <- delta_values[i]
  tau_i <- tau[i]
  sel_probs <- selection_probabilities(mu_samples, complexities, delta, alpha_n = n^(0.45))
  all_results[[i]] <- data.frame(
    model = colnames(Z_bc),
    complexity = complexities,
    sel_prob = sel_probs,
    delta = delta,
    tau = tau_i
  )
}

final_results <- do.call(rbind, all_results)  %>%
  mutate(
    model = factor(model, levels = names(fit_list)),
    k = as.integer(model)
  ) 

delta_palette <- c("orange", "purple1", "green3", "steelblue3")
plot_results <- final_results %>%
  mutate(
    model = factor(model, levels = colnames(Z_bc), labels = colnames(Z_bc)),
    facet_label = paste0("paste(delta==", round(delta, 3), ")")
  ) %>%
  mutate(delta_f = factor(signif(delta, 3), levels = unique(signif(delta, 3))))

plot_results$facet_label <- factor(
  plot_results$facet_label,
  levels = unique(plot_results$facet_label),
  ordered = TRUE
)


tpc_plot <- ggplot(plot_results, aes(x = factor(k), y = sel_prob, fill = delta_f)) +
  geom_col(width = 0.8, color = "black", linewidth = 0.2) +
  facet_wrap(~ facet_label, ncol = 1, labeller = label_parsed) +
  scale_fill_manual(values = delta_palette, guide = "none") +
  labs(x = "k", y = expression(hat(w)[delta](k))) +
  scale_y_continuous(limits = c(0, 1),
                     breaks  = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.50", "0.75", "1")) +
  theme_bw(base_size = 20) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black", size = 18),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 18),
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )



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


df_centered <- df_mu_long %>%
  group_by(draw) %>%
  mutate(mu_flex = min(mu_sample, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(mu_diff = mu_sample - mu_flex)


delta_df <- NULL
if (exists("delta_values") && length(delta_values) > 0) {
  delta_df <- tibble(
    delta = delta_values,
    delta_lab = factor(signif(delta_values, 2), levels = signif(delta_values, 2))
  )
}

ymax_center <- max(0, max(df_centered$mu_diff, na.rm = TRUE))
ygrid_center <- pretty(c(0, ymax_center), n = 5)

p_qoi <- ggplot(df_centered, aes(x = factor(k), y = mu_diff)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3, color = "black", fill = "white") +
  geom_hline(data = delta_df, aes(yintercept = delta, color = delta_lab),
             linewidth = 0.6, show.legend = FALSE) +
  geom_vline(xintercept = seq(1, K_max, by = 1), color = "grey70", linewidth = 0.2, alpha = 0.3) +
  geom_hline(yintercept = ygrid_center, color = "grey70", linewidth = 0.2, alpha = 0.3) +
  scale_y_continuous(limits = range(ygrid_center), breaks = ygrid_center) +
  scale_color_manual(values = c("orange", "purple1", "green3", "steelblue3")) +
  labs(x = "k", y = expression(mu[k] - mu^"*")) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 18),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5))
  )


spacer <- ggplot() + theme_void()
fig.tpc.1 <- plot_grid(
  tpc_fit,
  plot_grid(p_qoi, spacer, ncol = 2, rel_widths = c(1, 0.2)),  # shrink p_qoi width
  ncol = 1, align = "v", axis = "l",
  rel_heights = c(1, 0.55),
  labels = c("(A)", "(B)"),
  label_size = 22, label_fontface = "bold",
  label_x = -0.01, label_y = 1
)


fig.tpc <- plot_grid(
  fig.tpc.1, spacer, tpc_plot,
  ncol = 3, rel_widths = c(3, 0.1, 1.5), 
  labels = c("", "", "(C)"),
  label_size = 22, label_fontface = "bold",
  label_x = 0.01, label_y = 1
)


# save for Figure 6
# fig_dir <- here::here("output", "figures")
# ggsave(file.path(fig_dir, "tpc_fig6.png"), fig.tpc, width = 15, height = 10, dpi = 300)




