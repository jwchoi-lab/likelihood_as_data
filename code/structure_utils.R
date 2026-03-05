## ------------------------------------------------------------------
## Utilities for the STRUCTURE admixture example
##
## These helpers are used in the population structure example to:
##  - read a STRUCTURE output file for a given K,
##  - align the locus ordering in the STRUCTURE output with the genotype matrix X, and
##  - compute per-individual negative log-likelihoods under the fitted admixture model, 
##    integrating over ancestry proportions and latent population assignments (see Supplement S3–S4.4).
## ------------------------------------------------------------------


## ------------------------------------------------------------------
## Parse STRUCTURE output
##
## For a single STRUCTURE run and a fixed number of populations K, this function pulls out alpha_j, theta_{k,l}, and the allele labels used at each locus.
##
## It also records a few summary diagnostics reported by STRUCTURE:
##  - Estimated Ln Prob of Data,
##  - Mean value of ln likelihood,
##  - Variance of ln likelihood.
## ------------------------------------------------------------------
parse_structure_results <- function(f, K) {
  # Read the file once;
  lines <- readLines(f, warn = FALSE)
  txt <- paste(lines, collapse = " ")

  # alpha: try alpha_j, else single alpha
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
      alpha <- rep(NA, K)
    }
  }
  
  # locate locus headers
  locus_starts <- grep("^Locus [0-9]+", lines)
  L <- length(locus_starts)
  
  # per-locus: allele labels and per-k frequency vectors
  labels <- vector("list", L)
  theta_by_k <- vector("list", K)
  for (k in seq_len(K)) theta_by_k[[k]] <- vector("list", L)
  
  for (l in seq_len(L)) {
    start <- locus_starts[l] + 2  # skip "X alleles" & "missing data"
    end <- if (l < L) locus_starts[l + 1] - 1 else length(lines)
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
  
  get_num <- function(pattern) {
    m <- regexpr(pattern, txt)
    if (m == -1) return(NA_real_)
    as.numeric(sub(".*= *", "", regmatches(txt, m)))
  }
  
  list(
    alpha = alpha, 
    theta = theta_by_k, 
    labels = labels, 
    estimated_ln_prob  = get_num("Estimated Ln Prob of Data *= *-?[0-9]+\\.?[0-9]*"),
    mean_ln_likelihood = get_num("Mean value of ln likelihood *= *-?[0-9]+\\.?[0-9]*"),
    var_ln_likelihood  = get_num("Variance of ln likelihood *= *[0-9]+\\.?[0-9]*")
  )
}


## ------------------------------------------------------------------
## Match locus ordering between X and STRUCTURE output
##
## The raw genotype matrix X stores alleles as:
##   (locus 1 copy 1, locus 1 copy 2, locus 2 copy 1, locus 2 copy 2, ...)
##
## The STRUCTURE output may list loci in a different order. This function constructs a permutation that best aligns the two orderings, 
## based on the overlap of observed allele labels at each locus (Jaccard similarity).
##
## Returns: an integer vector 'map' of length L such that map[d] = l means "locus d in X corresponds to locus l in the STRUCTURE output".
## ------------------------------------------------------------------

find_locus_map <- function(X, labels) {
  L  <- ncol(X) / 2
  Sx <- lapply(seq_len(L), function(d) unique(na.omit(c(X[,2*d-1], X[,2*d]))))
  Sl <- labels
  
  C <- matrix(0, nrow = L, ncol = L)
  for (d in seq_len(L)) for (l in seq_len(L)) {
    a <- Sx[[d]]; b <- Sl[[l]]
    if (length(a) == 0L && length(b) == 0L) { C[d,l] <- 1; next }
    if (length(a) == 0L || length(b) == 0L) { C[d,l] <- 0; next }
    C[d,l] <- length(intersect(a,b)) / length(union(a,b))
  }
  
  map <- rep(NA, L); used <- rep(FALSE, L)
  order_d <- order(apply(C, 1, max), decreasing = TRUE)
  for (d in order_d) {
    s <- C[d, ]; s[used] <- -Inf
    l <- which.max(s)
    if (!is.finite(s[l])) l <- which(!used)[1]
    map[d] <- l; used[l] <- TRUE
  }
  map
}



## ------------------------------------------------------------------
## Per-individual negative log-likelihoods for the admixture model
##
## Given:
##  - X: n x (2L) genotype matrix; columns are allele copies (locus 1 copy 1, locus 1 copy 2, ..., locus L copy 2),
##  - alpha: length-K Dirichlet parameters (one per population),
##  - theta: list of length K; theta[[k]][[l]] is a vector of allele frequencies for population k at locus l,
##  - labels: list of length L; labels[[l]] is the vector of allele codes corresponding to theta[[k]][[l]],
##
## this function returns an n x 2 data.frame with:
##  - ind: individual index (row of X),
##  - nll: Monte Carlo estimate of -log p(x_i | alpha, theta),
## integrating out both the ancestry proportions q_i ~ Dirichlet(alpha) and the latent population assignments as in Supplement S4.4.
##
## For K = 1, the admixture model reduces to a single population and we can compute the likelihood exactly. 
## For K > 1, we approximate the Dirichlet integral with S Monte Carlo samples of q_i.
## ------------------------------------------------------------------
compute_point_nll <- function(X, alpha, theta, labels, S = 1000, eps = 1e-12, auto_map = TRUE) {
  n <- nrow(X); P <- ncol(X); L <- P/2; K <- length(theta)
  
  # locus permutation to match X
  loc_map <- if (auto_map) find_locus_map(X, labels) else seq_len(L)
  labels <- labels[loc_map]
  theta <- lapply(theta, function(thk) thk[loc_map])
  
  idx_of <- function(l, v) match(v, labels[[l]])  # NA if allele unseen at locus l
  
  nll <- rep(0, n)
  
  if (K == 1) {
    th1 <- theta[[1]]  
    for (i in seq_len(n)) {
      xi <- as.integer(X[i, ])
      obs <- which(!is.na(xi))
      ll <- 0
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
    xi <- as.integer(X[i, ])
    obs <- which(!is.na(xi))
    
    U <- matrix(rgamma(S * K, shape = alpha, rate = 1), nrow = S, ncol = K, byrow = TRUE)
    W <- U / rowSums(U)
    logW <- log(W)
    
    logA <- matrix(0, nrow = S, ncol = P)
    for (j in obs) {
      l <- ceiling(j/2)
      v <- xi[j]
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
