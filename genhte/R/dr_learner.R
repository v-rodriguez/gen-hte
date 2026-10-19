#' DR-Learner for generalising HTEs
#'
#' @param x a data frame of predictors of treatment response
#' @param y a response vector
#' @param s a vector of indicators of study participation
#' @param a a vector indicator of treatment assignment
#' @param x_s (optional) a data frame of predictors of selection
#' @param x_a (optional) a data frame of predictors of treatment
#' @param x_tau (optional) vector of covariates for the second step regression
#' @param x_g (optional) vector of covariates for g_a step (defaults to x_tau if not provided)
#' @param x_test (optional) a data frame of predictors for testing
#' @param reweight should the data be re-weighted using IPSW?
#' @param trim should the propensity scores be trimmed?
#' @param trim_perc (only when trim = TRUE) percentile used for trimming
#' @param s.weights sample weights for study selection propensities
#' @param subsample (optional) vector indicating external subsample membership for improved propensity estimation
#' @param pool_selection (only when subsample is supplied) fit the selection model on the pooled trial and subsample rows, inverting the pooled regression back onto the omega scale. FALSE fits on the subsample rows alone, which needs no inversion but sees far fewer S=1 rows.
#' @param use_g_step (logical) if TRUE and subsample provided, compute intermediate g_a for marginalization over unmeasured vars
#' @param seed sets the seed at the start
#' @param K number of folds used for cross-fitting (default K=5)
#' @param out_method Outcome regression learner: "rf" (default), "xgboost", "superlearner", "xgb_cv" (shallow xgboost, rounds chosen by CV), "xgb_stack" (super learner over xgboost settings) or "xgb_outcome" (deliberately under-fitted xgboost).
#' @param prop_method Treatment propensity learner: "xgboost" (default), "rf", "rf_cv" (leaf size chosen out-of-bag) or "superlearner".
#' @param select_method Selection propensity learner, with the same options as prop_method.
#' @param final_method Second-stage regression learner for the CATE and the potential-outcome regressions: "rf" (default), "rf_cv", "xgboost", "superlearner", "xgb_cv", "xgb_cv_tuned" (xgb_cv with min_child_weight chosen by CV) or "xgb_stack".
#'
#' @return fitted values of the CATE
#' @export
dr_learner <- function(
  x,
  y,
  s,
  a,
  x_s = NULL,
  x_a = NULL,
  x_tau = NULL,
  x_g = NULL,
  x_test = NULL,
  reweight = FALSE,
  trim = FALSE,
  trim_perc = 0.9,
  s.weights = NULL,
  subsample = NULL,
  pool_selection = TRUE,
  use_g_step = FALSE,
  K = 5,
  seed = NULL,
  out_method = "rf",
  prop_method = "xgboost",
  select_method = "xgboost",
  final_method = "rf"
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- length(y)

  if (is.null(x_a)) {
    x_a <- x
  }
  if (is.null(x_s)) {
    x_s <- x
  }
  if (is.null(x_tau)) {
    x_tau <- x
  }
  if (is.null(x_g)) {
    x_g <- x_tau
  }
  if (!check_lengths(x[, 1], y, a, s, x_s[, 1], x_a[, 1])) {
    stop("x, y, a and s must have the same length")
  }
  if (!is.null(x_test) && !setequal(names(x_tau), names(x_test))) {
    stop("variables in x_test must have the same names as variables in x_tau")
  }
  check_subsample(subsample, s)

  if (select_method == "xgboost") {
    x_s <- as.matrix(x_s)
  }

  # P(S'=1), needed to invert the pooled selection regression back onto the omega
  # scale. A scalar summary of the sampling design, so it is taken over all rows
  # rather than per fold. NULL whenever the pooled fit is not in play, which is
  # also how study_prop is told no inversion is wanted.
  pool_selection <- pool_selection && !is.null(subsample)
  p_subsample <- if (pool_selection) mean(subsample) else NULL

  folds <- groupdata2::fold(data.frame(x, s = s), k = K, cat_col = "s")$.folds

  mu1_list <- rep(NA_real_, n)
  mu0_list <- rep(NA_real_, n)
  pi_list <- rep(NA_real_, n)
  ex_list <- rep(NA_real_, n)
  phi_list <- rep(NA_real_, n)
  phi0_list <- rep(NA_real_, n)
  phi1_list <- rep(NA_real_, n)

  for (k in 1:K) {
    d1 <- folds != k
    d2 <- folds == k

    pi_mod <- treat_prop(x_a[d1, ], a[d1], s[d1], method = prop_method)

    A <- a[d2]
    S <- s[d2]
    Y <- y[d2]

    mu_0 <- outcome_train(subset(x, d1), y[d1], a[d1], s[d1], 0, method = out_method)
    mu_1 <- outcome_train(subset(x, d1), y[d1], a[d1], s[d1], 1, method = out_method)

    mu_0_hat <- get_predictions(mu_0, subset(x, d2))
    mu_1_hat <- get_predictions(mu_1, subset(x, d2))

    mu1_list[d2] <- mu_1_hat
    mu0_list[d2] <- mu_0_hat

    # The pseudo-outcomes are base + weight * residual. base is the plug-in
    # mu_a_hat, or g_a_hat once the g-step has marginalised it; with the g-step
    # the residual correction applies to trial rows only.
    base_0 <- mu_0_hat
    base_1 <- mu_1_hat
    rows <- rep(TRUE, length(Y))

    if (reweight) {
      sel_rows <- selection_rows(d1, s, subsample, pool_selection)

      ex_mod <- if (is.null(s.weights)) {
        study_prop(
          subset(x_s, sel_rows),
          s[sel_rows],
          method = select_method,
          p_subsample = p_subsample
        )
      }

      pi_hat <- get_predictions(pi_mod, subset(x_a, d2))

      if (is.null(s.weights)) {
        ex_hat <- get_predictions(ex_mod, subset(x_s, d2))
        if (trim) {
          ex_hat <- WeightIt::trim(ex_hat, at = trim_perc, lower = TRUE)
        }
      } else {
        ex_hat <- s.weights[d2]
      }

      ex_list[d2] <- ex_hat

      trial_d1_global <- d1 & s == 1
      if (use_g_step && any(trial_d1_global)) {
        # Train g_a on the d1 trial rows, evaluate on d2 for proper cross-fitting
        trial_d1_local <- s[d1] == 1
        x_g_d1 <- subset(x_g, d1)

        mu_0_hat_d1 <- get_predictions(mu_0, subset(x, d1))
        mu_1_hat_d1 <- get_predictions(mu_1, subset(x, d1))

        if (is.null(s.weights)) {
          ex_hat_d1 <- get_predictions(ex_mod, subset(x_s, trial_d1_global))
          if (trim) {
            ex_hat_d1 <- truncate_vector(ex_hat_d1, trim_perc)
          }
        } else {
          ex_hat_d1 <- s.weights[trial_d1_global]
        }

        # mu_a is already fit unweighted here, so only the g_a step is gated.
        g_step_weights <- if (needs_ipsw(x_s, x_g)) 1 / ex_hat_d1 else NULL

        g_0_mod <- final_regression(
          phi = mu_0_hat_d1[trial_d1_local],
          x_tau = x_g_d1[trial_d1_local, , drop = FALSE],
          method = final_method,
          weights = g_step_weights
        )
        g_1_mod <- final_regression(
          phi = mu_1_hat_d1[trial_d1_local],
          x_tau = x_g_d1[trial_d1_local, , drop = FALSE],
          method = final_method,
          weights = g_step_weights
        )

        base_0 <- get_predictions(g_0_mod, subset(x_g, d2))
        base_1 <- get_predictions(g_1_mod, subset(x_g, d2))
        rows <- S == 1
      }
    } else {
      pi_hat <- get_predictions(pi_mod, subset(x_a, d2))
    }

    # The I(S=1) factor keeps outcomes of non-trial units out of the
    # pseudo-outcome: they contribute the plug-in only.
    wts <- if (reweight) dr_weights(A, pi_hat, S, ex_hat) else dr_weights(A, pi_hat)
    w <- wts$w
    w_1 <- wts$w_1
    w_0 <- wts$w_0

    pi_list[d2] <- pi_hat

    if (any(is.na(w))) {
      warning(
        "There are NA's in the weights. Consider reducing trimming percentile."
      )
    }

    residual <- Y - A * mu_1_hat - (1 - A) * mu_0_hat
    phi_0 <- base_0
    phi_1 <- base_1
    phi <- base_1 - base_0
    phi_0[rows] <- w_0[rows] * (Y[rows] - mu_0_hat[rows]) + base_0[rows]
    phi_1[rows] <- w_1[rows] * (Y[rows] - mu_1_hat[rows]) + base_1[rows]
    phi[rows] <- w[rows] * residual[rows] + (base_1[rows] - base_0[rows])

    phi_list[d2] <- phi
    phi0_list[d2] <- phi_0
    phi1_list[d2] <- phi_1
  }

  # Keep this fit/predict order: ranger's predict draws from the R RNG.
  dr_mod <- final_regression(phi = phi_list, x_tau = x_tau, method = final_method)
  tau_hat <- get_predictions(dr_mod, x_tau)
  mu_0_mod <- final_regression(phi = phi0_list, x_tau = x_tau, method = final_method)
  mu_1_mod <- final_regression(phi = phi1_list, x_tau = x_tau, method = final_method)

  train <- list(
    id = seq_len(n),
    tau_hat = tau_hat,
    mu_0_hat = get_predictions(mu_0_mod, x_tau),
    mu_1_hat = get_predictions(mu_1_mod, x_tau),
    phi = phi_list,
    mu_0 = mu0_list,
    mu_1 = mu1_list,
    pi = pi_list,
    ex = ex_list
  )

  if (is.null(x_test)) {
    return(train)
  }

  list(
    train = train,
    test = list(
      id = seq_len(nrow(x_test)),
      tau_hat = get_predictions(dr_mod, x_test),
      mu_0_hat = get_predictions(mu_0_mod, x_test),
      mu_1_hat = get_predictions(mu_1_mod, x_test)
    )
  )
}
