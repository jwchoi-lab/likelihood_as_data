# 04_structure_admixture.R
# Population structure admixture models (Section 5.4)

library(matrixStats) 
library(dplyr)
library(purrr)
library(tidyr)
library(gridExtra)
library(ggplot2)
library(cowplot)

source(here::here("code", "functions.R"))
source(here::here("code", "structure_utils.R"))


# Read Data
dat <- read.table(
  here::here("data", "raw", "admixture", "brooktrout.txt")
)

X <- dat


# For exact replication of the figures in the paper, we provide precomputed results in
# output/admixture_fitted.Rdata:
# load(here::here("output", "admixture_fitted.Rdata"))

# Otherwise, the scripts below also allow recomputing the results;
# these may differ slightly up to random sampling, but yield the same qualitative conclusions.




set.seed(251111)

grid <- expand.grid(
  K = 1:10,
  rep = 1:20
) %>%
  as_tibble() %>%
  mutate(
    file = here::here("data", "processed", "structure_run",
                      sprintf("results_K%d_rep%d_f", K, rep))
  )

results <- grid %>%
  mutate(
    parsed = map2(file, K, parse_structure_results),
    mean_ll = map_dbl(parsed, "estimated_ln_prob")
  ) %>%
  dplyr::select(K, rep, mean_ll)

best_reps <- results %>%
  group_by(K) %>%
  slice_max(mean_ll, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(file = here::here("data", "processed", "structure_run",
                           sprintf("results_K%d_rep%d_f", K, rep))
         )

final_rep_nll <- best_reps %>%
  mutate(parsed = map2(file, K, parse_structure_results)) %>%
  mutate(nlldf  = map(parsed, ~ compute_point_nll(X, .x$alpha, .x$theta, labels = .x$labels))) %>%
  dplyr::select(K, rep, nlldf) %>%
  unnest(nlldf)

final_nll <- final_rep_nll %>% dplyr::select(-rep)



# save(final_nll, file = paste("output/admixture_fitted.Rdata"))



## Apply LaD
nll_matrix <- final_nll %>%
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
Z_bc  <- t(t(Z) + d_k/(2*n))  # bias-corrected. 



# Noise model
M <- (!is.na(X[, seq(1, 2*L, 2)])) + (!is.na(X[, seq(2, 2*L, 2)]))
nll_noise <- as.vector(M %*% log(V_l))
mu_noise_uniform <- mean(nll_noise)
print(mu_noise_uniform )

# mu_noise <- as.numeric(colMeans(Z_bc)[1])
mu_k <- as.numeric(colMeans(Z_bc))
print(mu_k)

# NIW Prior Parameters
mu0 <- rep(0, K)  # Prior mean (vector of zeros).
lambda0 <- 0.01 # Prior scaling for the mean.
nu0 <- K + 2  # Degrees of freedom (must be > K-1).
Psi0 <- 1 * diag(K)  # Prior inverse scale matrix.

n_post_samples <- 1000
mu_samples <- niw_posterior(Z_bc, mu0, lambda0, Psi0, nu0, n_post_samples)$mu 




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
    sel_prob = as.numeric( selection_probabilities(mu_samples, complexities, delta = .y, alpha_n = a_n) ),
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

fig.ad <- plot_grid(
  p_mu, p_heat,
  ncol = 2, align = "v", axis = "l",
  rel_widths = c(1, 1.2) 
)





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


# save for Figure S1 and S2
# fig_dir <- here::here("output", "figures")
# ggsave(file.path(fig_dir, "admixture_figS1.png"), fig.ad, width = 15, height = 5, dpi = 300)
# ggsave(file.path(fig_dir, "admixture_figS2.png"), evanno_plot, width = 15, height = 5, dpi = 300)


