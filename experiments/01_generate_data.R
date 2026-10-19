#' Simulate a target population with a nested trial
#'
#' @param ssize Expected trial size
#' @param xtau pCATE model: "simple" (linear) or "complex" (adds sin(5 X_tau))
#' @param xmu Control outcome model: "simple" (linear) or "nonlinear"
#' @param delta Number of selection covariates not in X_tau, |X_S \ X_tau|
#' @param N Size of the target population
#' @param ntau Number of treatment effect modifiers X_tau
#' @param nsel Number of selection covariates X_S
#' @param nmu Number of outcome covariates (X_S plus prognostic Z)
#' @param seed Seed set at the start
#' @param subsample Draw an external random subsample of size nsubsample?
#'
#' @return Data frame of the target population with trial indicator S
generate_data <- function(
  ssize,
  xtau,
  xmu,
  delta = 1,
  N = 20000,
  ntau = 2,
  nsel = 2,
  nmu = 5,
  seed = NULL,
  subsample = FALSE,
  nsubsample = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  stopifnot(xtau %in% c("simple", "complex"))
  stopifnot(xmu %in% c("simple", "nonlinear"))
  stopifnot(is.numeric(ssize))
  stopifnot(is.numeric(ntau) && ntau >= 1)
  stopifnot(is.numeric(nsel) && nsel >= 1)
  stopifnot(is.numeric(nmu) && nmu >= 1)
  stopifnot(is.numeric(delta) && delta >= 0 && delta <= nsel && nsel - delta <= ntau)
  stopifnot(ntau <= nsel || nsel <= nmu)
  stopifnot(is.logical(subsample) && length(subsample) == 1)
  if (subsample) {
    stopifnot(!is.null(nsubsample))
    stopifnot(is.numeric(nsubsample), nsubsample >= 1, nsubsample <= N)
  }

  X_sel <- matrix(runif(N * nsel, -1, 1), nrow = N, ncol = nsel)

  n_extra_outcome <- max(0, nmu - nsel)
  X_prognostic <- matrix(
    runif(N * n_extra_outcome, -1, 1),
    nrow = N,
    ncol = n_extra_outcome
  )

  X_outcome <- cbind(X_sel, X_prognostic)

  X_tau <- matrix(runif(N * ntau, -1, 1), nrow = N, ncol = ntau)
  n_shared <- nsel - delta
  if (n_shared > 0) {
    X_tau[, seq_len(n_shared)] <- X_sel[, seq_len(n_shared)]
  }
  sel_only <- if (delta > 0) (n_shared + 1):nsel else integer(0)

  pcate <- if (xtau == "simple") {
    1 + rowSums(X_tau) / sqrt(ntau)
  } else {
    1 + rowSums(X_tau + sin(5 * X_tau)) / sqrt(ntau)
  }

  gamma_sel <- 3
  kappa_tail <- 3

  # Scaled by sqrt(nsel), not sqrt(delta), so the sCATE-pCATE bias grows
  # linearly in delta
  sel_effect <- if (delta > 0) {
    gamma_sel * rowSums(X_sel[, sel_only, drop = FALSE]) / sqrt(nsel)
  } else {
    0
  }

  # Effect concentrated where trial units are sparse; mean 0 and sd kappa_tail
  # from a fixed reference draw, and independent of X_tau so the pCATE is unchanged
  tail_effect <- if (delta > 0) {
    tail_f <- function(w) -pmax(-rowSums(w) / sqrt(nsel), 0)^2
    ref <- withr::with_seed(1, tail_f(matrix(runif(1e6 * delta, -1, 1), ncol = delta)))
    kappa_tail * (tail_f(X_sel[, sel_only, drop = FALSE]) - mean(ref)) / sd(ref)
  } else {
    0
  }

  effect_indiv <- pcate + sel_effect + tail_effect

  mu0 <- if (xmu == "nonlinear") {
    mu0_nonlinear(X_outcome, nmu)
  } else {
    rowSums(X_outcome[, seq_len(nmu), drop = FALSE]) / sqrt(nmu)
  }
  Y0 <- mu0
  Y1 <- Y0 + effect_indiv

  marg <- ssize / N

  sel_df <- data.frame(rep(1, N), sel_score = rowSums(X_sel) / sqrt(nsel))
  intercept <- determine_intercept(1, marg, sel_df)

  lp_sel <- as.vector(as.matrix(sel_df) %*% c(intercept, 1))
  selection_prob <- expit(lp_sel)

  S <- rbinom(N, 1, selection_prob)

  A <- rep(0, N)
  if (sum(S) > 0) {
    A[S == 1] <- rbinom(sum(S), 1, 0.5)
  }

  e <- rnorm(N, 0, 1)

  subsample_ind <- rep(0L, N)
  if (subsample) {
    idx <- sample.int(N, nsubsample, replace = FALSE)
    subsample_ind[idx] <- 1
  }

  df_prognostic <- as.data.frame(X_prognostic)
  names(df_prognostic) <- sprintf("Z%d", seq_len(n_extra_outcome))
  df_tau <- as.data.frame(X_tau)
  names(df_tau) <- paste0("T", seq_len(ncol(df_tau)))
  df_selection <- as.data.frame(X_sel)
  names(df_selection) <- paste0("X", seq_len(ncol(df_selection)))

  data.frame(
    df_prognostic,
    df_tau,
    df_selection,
    tau = pcate,
    ite = effect_indiv,
    S = S,
    subsample = subsample_ind,
    Y0 = Y0,
    Y1 = Y1,
    Y = (1 - A) * Y0 + A * Y1 + e,
    A = A,
    e = e,
    selection_prob = selection_prob,
    check.names = FALSE
  )
}


# High-frequency main effects plus all pairwise interactions, each block O(1)
mu0_nonlinear <- function(X, nmu) {
  X <- X[, seq_len(nmu), drop = FALSE]
  main <- rowSums(X + sin(4 * pi * X)) / sqrt(nmu)
  if (nmu < 2) {
    return(main)
  }
  pairs <- utils::combn(nmu, 2)
  inter <- vapply(
    seq_len(ncol(pairs)),
    function(k) sin(2 * pi * X[, pairs[1, k]] * X[, pairs[2, k]]),
    numeric(nrow(X))
  )
  main + rowSums(as.matrix(inter)) / sqrt(ncol(pairs))
}


expit <- function(x) {
  exp(x) / (1 + exp(x))
}

# Bisection for the intercept b0 with mean(expit(b0 + beta * x)) = marg, adapted
# from https://github.com/serobertson/BalancingInterceptSolver. Coefficients are
# pasted into a formula string as in the original, which rounds them to 15
# significant digits; kept to reproduce the simulated data exactly.
determine_intercept <- function(
  beta_vec,
  marg = 0.5,
  external_dataset,
  tolerance = .001,
  lower_bound = -10,
  upper_bound = 10
) {
  mean_prob <- function(coefs) {
    linearformula <- paste(coefs, colnames(external_dataset), sep = "*")
    lp <- with(external_dataset, eval(parse(text = paste(linearformula, collapse = " + "))))
    mean(expit(lp))
  }

  beta_vec <- c((upper_bound + lower_bound) / 2, beta_vec)
  mean.oc <- mean_prob(beta_vec)
  count <- 1

  while (
    abs(mean.oc - marg) > tolerance &&
      upper_bound > lower_bound &&
      count <= 99
  ) {
    if (mean.oc > marg) {
      upper_bound <- beta_vec[1]
    } else if (mean.oc < marg) {
      lower_bound <- beta_vec[1]
    }
    beta_vec[1] <- (upper_bound + lower_bound) / 2
    mean.oc <- mean_prob(beta_vec)
    count <- count + 1
  }

  if (abs(mean.oc - marg) > tolerance) {
    warning("Algorithm did not converge within tolerance.")
  }
  beta_vec[1]
}
