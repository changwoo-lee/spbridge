
# Build all within-cluster pairs once (optionally subsample per cluster)
.make_pairs <- function(id, pair_frac = 1) {
  id <- as.factor(id)
  idx_by_cluster <- split(seq_along(id), id)
  Ps <- vector("list", length(idx_by_cluster))
  m <- 0L
  for (g in seq_along(idx_by_cluster)) {
    idx <- idx_by_cluster[[g]]
    ni <- length(idx)
    if (ni < 2L) next
    P <- utils::combn(idx, 2L)
    if (pair_frac < 1) {
      keep <- sample(ncol(P), max(1L, floor(pair_frac * ncol(P))))
      P <- P[, keep, drop = FALSE]
    }
    m <- m + 1L
    Ps[[m]] <- P
  }
  Ps <- Ps[seq_len(m)]
  if (!length(Ps)) return(matrix(integer(), nrow = 2L))
  do.call(cbind, Ps)
}


.pair_probs_vec <- function(mu1, mu2, phi, rel.tol) {
  vapply(
    seq_along(mu1),
    function(i) pair_probs_logistic_bridge(mu1[i], mu2[i], phi, rel.tol),
    FUN.VALUE = numeric(4L)
  )
}

expit <- function(x) 1/(1+exp(-x))

composite_ll_phi_integrate_fast <- function(y, id, mu_hat, rel.tol = 1e-8, pair_frac = 1) {
  P <- .make_pairs(id, pair_frac)
  if (ncol(P) == 0L) {
    return(function(phi) if (is.finite(phi) && phi > 0 && phi < 1) 0.0 else -Inf)
  }

  # Precompute
  j <- P[1L, ]; k <- P[2L, ]
  mu1 <- mu_hat[j]; mu2 <- mu_hat[k]

  # Map outcomes (00,01,10,11) -> (1,2,3,4)
  y <- as.integer(y)
  key <- 1L + 2L * y[j] + y[k]

  function(phi) {
    if (!(is.finite(phi) && phi > 0 && phi < 1)) return(-Inf)
    probs_mat <- .pair_probs_vec(mu1, mu2, phi, rel.tol) # 4 x M
    if (!all(is.finite(probs_mat))) return(-Inf)
    # Pick the probability for each pair's observed (y_j,y_k)
    sel <- probs_mat[cbind(key, seq_along(key))]
    if (any(sel <= 0 | !is.finite(sel))) return(-Inf)
    sum(log(sel))
  }
}

# \int g(b) dbridge(b; phi) db
int_bridge <- function(gfun, phi, rel.tol = 1e-8) {
  out <- try(integrate(function(b) gfun(b) * dbridge(b, phi),
                       lower = -Inf, upper = Inf, rel.tol = rel.tol),
             silent = TRUE)
  if (inherits(out, "try-error") || !is.finite(out$value)) return(NA_real_)
  out$value
}


pair_probs_logistic_bridge <- function(mu_j, mu_k, phi, rel.tol = 1e-8) {
  if (!(is.finite(phi) && phi > 0 && phi < 1)) return(rep(NA_real_, 4))
  eta_j <- mu_j / phi
  eta_k <- mu_k / phi

  # Joint success probability p11
  g <- function(b) expit(eta_j + b) * expit(eta_k + b)
  p11 <- int_bridge(g, phi, rel.tol = rel.tol)

  # Marginal means
  m_j <- expit(mu_j)
  m_k <- expit(mu_k)
  # Remaining joint probs
  p10 <- m_j - p11
  p01 <- m_k - p11
  p00 <- 1 - m_j - m_k + p11

  eps <- 1e-12
  p <- pmax(pmin(c(p00, p01, p10, p11), 1 - eps), eps)
  p <- p / sum(p)
  p
}

estimate_phi_pairCL <- function(y, id, mu_hat, rel.tol = 1e-8, pair_frac = 1) {
  ll <- composite_ll_phi_integrate_fast(y, id, mu_hat, rel.tol, pair_frac)
  eps <- 1e-8
  obj <- function(phi) -ll(phi)
  opt <- optimize(f = obj, interval = c(eps, 1 - eps), tol = 1e-8)
  list(phi = opt$minimum, value = -opt$objective, conv = 0L)
}
