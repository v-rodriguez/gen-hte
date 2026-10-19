#' T-Learner for generalising HTEs
#'
#' @param x a data frame of predictors for training
#' @param y a response vector
#' @param s a vector of indicators of study participation
#' @param a a vector indicator of treatment assignment
#' @param x_s (optional) a data frame of predictors of selection
#' @param x_a (optional) a data frame of predictors of treatment
#' @param x_test a data frame of predictors for testing
#' @param reweight should the data be re-weighted using IPSW?
#' @param trim should the propensity scores be trimmed?
#' @param trim_perc (only when trim = TRUE) percentile used for trimming
#' @param s.weights sample weights for study selection propensities
#' @param subsample (optional) vector indicating external subsample membership for improved propensity estimation
#' @param pool_selection (only when subsample is supplied) fit the selection model on the pooled trial and subsample rows, inverting the pooled regression back onto the omega scale. FALSE fits on the subsample rows alone, which needs no inversion but sees far fewer S=1 rows.
#' @param use_g_step (logical) if TRUE and subsample provided, compute intermediate g_a for marginalization
#' @param x_tau (optional) a data frame of treatment effect modifiers (when use_g_step=TRUE)
#' @param x_g (optional) variables for g_a regression (marginalizes over variables NOT in x_g). Defaults to x_tau if not provided.
#' @param seed sets the seed at the start
#' @param K number of folds used for cross-fitting (default K=5)
#' @param out_method Outcome regression learner: "rf" (default), "xgboost", "superlearner", "xgb_cv" (shallow xgboost, rounds chosen by CV), "xgb_stack" (super learner over xgboost settings) or "xgb_outcome" (deliberately under-fitted xgboost).
#' @param prop_method Treatment propensity learner: "xgboost" (default), "rf", "rf_cv" (leaf size chosen out-of-bag) or "superlearner".
#' @param select_method Selection propensity learner, with the same options as prop_method.
#' @param final_method G-step regression learner: "rf" (default), "rf_cv", "xgboost", "superlearner", "xgb_cv", "xgb_cv_tuned" (xgb_cv with min_child_weight chosen by CV) or "xgb_stack".
#'
#' @return fitted values of the CATE
#' @export
t_learner <- function(
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
  if (!check_lengths(x[, 1], y, a, s, x_s[, 1], x_a[, 1])) {
    stop("x, y, a, and s must have the same length")
  }

  if (!is.null(x_test)) {
    if (!is.null(x_tau)) {
      if (!all(names(x_tau) %in% names(x_test))) {
        stop("x_test must contain all variables in x_tau")
      }
    } else if (!setequal(names(x), names(x_test))) {
      stop("Variables in x_test must have the same names as variables in x")
    }
  }
  check_subsample(subsample, s)

  if (use_g_step) {
    if (is.null(x_tau)) {
      stop("use_g_step requires x_tau to be provided")
    }
    if (!reweight) {
      stop("use_g_step requires reweight=TRUE")
    }
  }

  if (is.null(x_g)) {
    x_g <- x_tau
  }

  if (select_method == "xgboost") {
    x_s <- as.matrix(x_s)
  }

  # Both regressions below are conditional expectations, so the IPSW weights are
  # redundant wherever x_s is already among the regressors (see needs_ipsw). x_g
  # is NULL when use_g_step is FALSE and x_tau is unset, in which case
  # weight_g_step is never read.
  weight_outcome <- needs_ipsw(x_s, x)
  weight_g_step <- needs_ipsw(x_s, x_g)

  # P(S'=1), needed to invert the pooled selection regression back onto the omega
  # scale. A scalar summary of the sampling design, so it is taken over all rows
  # rather than per fold. NULL whenever the pooled fit is not in play, which is
  # also how study_prop is told no inversion is wanted.
  pool_selection <- pool_selection && !is.null(subsample)
  p_subsample <- if (pool_selection) mean(subsample) else NULL

  folds <- groupdata2::fold(data.frame(x, s), k = K, cat_col = "s")$.folds

  mu1_list <- rep(NA_real_, n)
  mu0_list <- rep(NA_real_, n)

  mu_0_models <- list()
  mu_1_models <- list()
  g_0_models <- list()
  g_1_models <- list()

  for (k in 1:K) {
    d1 <- folds != k
    d2 <- folds == k

    if (reweight) {
      pi_mod <- treat_prop(subset(x_a, d1), a[d1], s[d1], method = prop_method)
      pi_hat <- get_predictions(pi_mod, subset(x_a, d1 & s == 1))

      if (is.null(s.weights)) {
        sel_rows <- selection_rows(d1, s, subsample, pool_selection)
        ex_mod <- study_prop(
          subset(x_s, sel_rows),
          s[sel_rows],
          method = select_method,
          p_subsample = p_subsample
        )
        ex_hat <- get_predictions(ex_mod, subset(x_s, d1 & s == 1))
        if (trim) {
          ex_hat <- truncate_vector(ex_hat, trim_perc)
        }
      } else {
        # Supplied weights are taken as given and left untrimmed.
        ex_hat <- s.weights[d1 & s == 1]
      }

      A <- a[d1 & s == 1]
      S <- s[d1]

      # Stabilised weights
      pi_a <- mean(A)
      p_s <- mean(S)
      w <- (A * pi_a * p_s) /
        (pi_hat * ex_hat) +
        ((1 - A) * (1 - pi_a) * p_s) / (ex_hat * (1 - pi_hat))

      if (any(is.na(w))) {
        warning(
          "There are NA's in the weights. These points have been removed. Consider reducing trimming percentile."
        )
      }

      w_outcome <- if (weight_outcome) w
    } else {
      w_outcome <- NULL
    }

    mu_0_models[[k]] <- outcome_train(
      subset(x, d1), y[d1], a[d1], s[d1], 0, w = w_outcome, method = out_method
    )
    mu_1_models[[k]] <- outcome_train(
      subset(x, d1), y[d1], a[d1], s[d1], 1, w = w_outcome, method = out_method
    )

    mu_0_hat <- get_predictions(mu_0_models[[k]], subset(x, d2))
    mu_1_hat <- get_predictions(mu_1_models[[k]], subset(x, d2))

    trial_d1_global <- d1 & s == 1
    if (use_g_step && any(trial_d1_global)) {
      trial_d1_local <- s[d1] == 1
      x_g_d1 <- subset(x_g, d1)

      mu_0_hat_d1 <- get_predictions(mu_0_models[[k]], subset(x, d1))
      mu_1_hat_d1 <- get_predictions(mu_1_models[[k]], subset(x, d1))

      if (is.null(s.weights)) {
        ex_hat_d1 <- get_predictions(ex_mod, subset(x_s, trial_d1_global))
        if (trim) {
          ex_hat_d1 <- truncate_vector(ex_hat_d1, trim_perc)
        }
      } else {
        ex_hat_d1 <- s.weights[trial_d1_global]
      }

      g_step_weights <- if (weight_g_step) 1 / ex_hat_d1 else NULL

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

      x_g_d2 <- subset(x_g, d2)
      mu_0_hat <- get_predictions(g_0_mod, x_g_d2)
      mu_1_hat <- get_predictions(g_1_mod, x_g_d2)

      if (identical(names(x_g), names(x_tau))) {
        g_0_models[[k]] <- g_0_mod
        g_1_models[[k]] <- g_1_mod
      } else {
        # Two-stage marginalisation x_g -> x_tau. x_g is observed for all units,
        # so this regression is fit unweighted on all of d1.
        x_tau_d1 <- subset(x_tau, d1)
        g_0_hat_d1 <- get_predictions(g_0_mod, x_g_d1)
        g_1_hat_d1 <- get_predictions(g_1_mod, x_g_d1)

        g_0_models[[k]] <- final_regression(
          phi = g_0_hat_d1, x_tau = x_tau_d1, method = final_method
        )
        g_1_models[[k]] <- final_regression(
          phi = g_1_hat_d1, x_tau = x_tau_d1, method = final_method
        )

        x_tau_d2 <- subset(x_tau, d2)
        mu_0_hat <- get_predictions(g_0_models[[k]], x_tau_d2)
        mu_1_hat <- get_predictions(g_1_models[[k]], x_tau_d2)
      }
    }

    mu0_list[d2] <- mu_0_hat
    mu1_list[d2] <- mu_1_hat
  }

  train <- list(
    id = seq_len(n),
    tau_hat = mu1_list - mu0_list,
    mu_0_hat = mu0_list,
    mu_1_hat = mu1_list,
    mu_0_models = mu_0_models,
    mu_1_models = mu_1_models
  )

  if (is.null(x_test)) {
    return(train)
  }

  # Ensemble test predictions across the K fold models
  if (use_g_step && length(g_0_models) > 0) {
    models_0 <- g_0_models
    models_1 <- g_1_models
    x_test <- x_test[, names(x_tau), drop = FALSE]
  } else {
    models_0 <- mu_0_models
    models_1 <- mu_1_models
    x_test <- x_test[, names(x), drop = FALSE]
  }

  mu_0_hat.test <- 0
  mu_1_hat.test <- 0
  tau_test <- 0
  for (k in 1:K) {
    mu_0_hat.test_k <- get_predictions(models_0[[k]], x_test)
    mu_1_hat.test_k <- get_predictions(models_1[[k]], x_test)
    mu_0_hat.test <- mu_0_hat.test + mu_0_hat.test_k
    mu_1_hat.test <- mu_1_hat.test + mu_1_hat.test_k
    tau_test <- tau_test + (mu_1_hat.test_k - mu_0_hat.test_k)
  }

  list(
    train = train,
    test = list(
      id = seq_len(NROW(x_test)),
      tau_hat = tau_test / K,
      mu_0_hat = mu_0_hat.test / K,
      mu_1_hat = mu_1_hat.test / K
    )
  )
}
