#' Compute causal prediction intervals using split conformal inference for generalization
#'
#' @param x matrix of outcome predictors
#' @param y vector of outcomes
#' @param s indicator of selection
#' @param a indicator of treatment (binary)
#' @param x_s matrix of selection predictors
#' @param x_a matrix of treatment predictors
#' @param x_tau matrix of treatment effect modifiers
#' @param train_ind indicator of training data
#' @param reweight (TRUE/FALSE) should the models use IPSW?
#' @param trim (TRUE/FALSE) should the propensity scores be trimmed?
#' @param trim_perc (only when trim = TRUE) percentile used for trimming
#' @param s.weights vector of given selection weights
#' @param seed sets the initial seed
#' @param method T-Learner or DR-learner
#' @param alpha nominal confidence level
#' @param K number of folds used for cross-fitting (default K=5)
#' @param range optional range to truncate the intervals
#' @param out_method method for outcome model
#' @param prop_method method for propensity model
#' @param select_method method for selection of nuisance models
#' @param final_method method for final regression model (DR-Learner only)
#'
#' @return a list of 2-d matrices with the prediction intervals for the
#'  CATE and potential outcomes Y0 and Y1
#' @export
split_causal_conformal <- function(
  x,
  y,
  s,
  a,
  x_s = NULL,
  x_a = NULL,
  x_tau = NULL,
  train_ind = NULL,
  reweight = FALSE,
  trim = TRUE,
  trim_perc = 0.9,
  s.weights = NULL,
  K = 5,
  seed = NULL,
  method = c("DR", "TL"),
  alpha = 0.05,
  range = NULL,
  out_method = "rf",
  prop_method = "rf",
  select_method = "xgboost",
  final_method = "rf"
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(x_a)) {
    x_a <- x
  }
  if (is.null(x_s)) {
    x_s <- x
  }
  if (is.null(x_tau)) {
    x_tau <- x
  }
  if (!check_lengths(x[, 1], y, a, s, x_s[, 1], x_a[, 1])) {
    stop("x, y, a and s must have the same length")
  }
  method <- match.arg(method)

  split <- rep(3, length(y))
  if (!is.null(train_ind)) {
    split[train_ind] <- groupdata2::partition(
      subset(data.frame(x, s), train_ind),
      p = c(0.8),
      cat_col = 's',
      list_out = FALSE
    )$.partitions
  } else {
    split <- groupdata2::partition(
      data.frame(x, s),
      p = c(0.8),
      cat_col = 's',
      list_out = FALSE
    )$.partitions
  }

  train <- (split == 1)
  cal <- (split == 2)

  learner <- if (method == "DR") dr_learner else t_learner
  pred <- learner(
    x = subset(x, train),
    y = y[train],
    s = s[train],
    a = a[train],
    x_tau = subset(x_tau, train),
    x_s = subset(x_s, train),
    x_a = subset(x_a, train),
    x_test = x_tau,
    reweight = reweight,
    s.weights = s.weights[train],
    K = K,
    trim = trim,
    trim_perc = trim_perc,
    out_method = out_method,
    prop_method = prop_method,
    select_method = select_method,
    final_method = final_method
  )

  cal_arm_1 <- cal & s == 1 & a == 1
  cal_arm_0 <- cal & s == 1 & a == 0

  # non-conformity scores (absolute residuals)
  scores_0 <- abs(pred$test$mu_0_hat[cal_arm_0] - y[cal_arm_0])
  scores_1 <- abs(pred$test$mu_1_hat[cal_arm_1] - y[cal_arm_1])

  arm_weights <- function(arm, cal_arm) {
    w <- calculate_weights(
      subset(x_a, train),
      subset(x_s, train),
      s[train],
      a[train],
      subset(x_a, cal_arm),
      subset(x_s, cal_arm),
      x_a,
      x_s,
      arm = arm,
      reweight = reweight,
      prop_method = prop_method,
      select_method = select_method,
      trim = trim,
      trim_perc = trim_perc
    )
    normalise_weights(w$cal, w$test)
  }
  weights_0 <- arm_weights(0, cal_arm_0)
  weights_1 <- arm_weights(1, cal_arm_1)

  # quantile of the weighted score distribution, per test point
  q_0 <- apply(weights_0, 2, function(w) {
    weighted_quantile(c(scores_0, Inf), prob = 1 - alpha, w = w)
  })
  q_1 <- apply(weights_1, 2, function(w) {
    weighted_quantile(c(scores_1, Inf), prob = 1 - alpha, w = w)
  })
  c_0 <- cbind(pred$test$mu_0_hat - q_0, pred$test$mu_0_hat + q_0)
  c_1 <- cbind(pred$test$mu_1_hat - q_1, pred$test$mu_1_hat + q_1)

  if (!is.null(range)) {
    c_0 <- pmin(pmax(c_0, range[1]), range[2])
    c_1 <- pmin(pmax(c_1, range[1]), range[2])
  }

  c_tau <- cbind(
    pmin(c_1[, 1] - c_0[, 2], c_1[, 2] - c_0[, 1]),
    pmax(c_1[, 1] - c_0[, 2], c_1[, 2] - c_0[, 1])
  )

  list(
    c_tau = c_tau,
    c_0 = c_0,
    c_1 = c_1,
    tau_hat = pred$test$tau_hat,
    mu_0_hat = pred$test$mu_0_hat,
    mu_1_hat = pred$test$mu_1_hat,
    q_0 = q_0,
    q_1 = q_1,
    max_score_0 = max(scores_0),
    max_score_1 = max(scores_1)
  )
}


#' Calculate IPSW weights
#'
#' @param x_a.train matrix of treatment predictors for training
#' @param x_s.train matrix of selection predictors for training
#' @param s.train vector of selection indicators for training
#' @param a.train vector of treatment indicators for training
#' @param x_a.cal matrix of treatment predictors for calibration
#' @param x_s.cal matrix of selection predictors for calibration
#' @param x_a.test matrix of treatment predictors for testing
#' @param x_s.test matrix of selection predictors for testing
#' @param arm indicator of treatment (1) or control (0) arm
#' @param reweight (TRUE/FALSE) should the weights carry the selection factor
#'   p_s / omega-hat? When FALSE no selection model is fit and the factor is
#'   dropped, which is what the trial-only arms want: their training rows are all
#'   S=1, so there is no omega to estimate.
#' @param prop_method learner for the treatment propensity.
#' @param select_method learner for the selection model.
#' @param trim (TRUE/FALSE) should the fitted propensities be trimmed?
#' @param trim_perc (only when trim = TRUE) percentile used for trimming
#'
#' @return list of vectors of weights
#' @export
calculate_weights <- function(
  x_a.train,
  x_s.train,
  s.train,
  a.train,
  x_a.cal,
  x_s.cal,
  x_a.test,
  x_s.test,
  arm = c(0, 1),
  reweight = TRUE,
  prop_method = "rf",
  select_method = "xgboost",
  trim = TRUE,
  trim_perc = 0.9
) {
  # Trimming clips the fitted propensities, with one pair of thresholds taken
  # from the training predictions and applied to all three sets. Per-set
  # quantiles would cap the calibration rows (trial members, omega near 1)
  # differently from the test rows, which span the whole population.
  clip_from <- function(reference) {
    if (!trim) {
      return(identity)
    }
    lo <- stats::quantile(reference, 1 - trim_perc, na.rm = TRUE, names = FALSE)
    hi <- stats::quantile(reference, trim_perc, na.rm = TRUE, names = FALSE)
    function(p) pmin(pmax(p, lo), hi)
  }

  pi_mod <- treat_prop(x_a.train, a.train, s.train, method = prop_method)
  pi_hat.train <- get_predictions(pi_mod, x_a.train)
  clip_pi <- clip_from(pi_hat.train)

  pi_hat.train <- clip_pi(pi_hat.train)
  pi_hat.cal <- clip_pi(get_predictions(pi_mod, x_a.cal))
  pi_hat.test <- clip_pi(get_predictions(pi_mod, x_a.test))

  # The selection factor is skipped when reweight is FALSE: fit on a constant S
  # it would be a constant that cancels in normalise_weights anyway.
  if (reweight) {
    ex_mod <- study_prop(x_s.train, s.train, method = select_method)
    p_s <- mean(s.train)

    ex_hat.train <- get_predictions(ex_mod, x_s.train)
    clip_ex <- clip_from(ex_hat.train)

    sel.train <- p_s / clip_ex(ex_hat.train)
    sel.cal <- p_s / clip_ex(get_predictions(ex_mod, x_s.cal))
    sel.test <- p_s / clip_ex(get_predictions(ex_mod, x_s.test))
  } else {
    sel.train <- sel.cal <- sel.test <- 1
  }

  p_a <- mean(a.train)
  num <- if (arm == 1) p_a else 1 - p_a
  den <- function(p) if (arm == 1) p else 1 - p

  list(
    train = sel.train * (num / den(pi_hat.train)),
    cal = sel.cal * (num / den(pi_hat.cal)),
    test = sel.test * (num / den(pi_hat.test))
  )
}


#' Normalise weights
#'
#' @param w.cal calibration weights
#' @param w.test test weights
#'
#' @return matrix with one column per test point: the calibration weights and
#'   that test point's weight, normalised to sum to one
#' @export
normalise_weights <- function(w.cal, w.test) {
  vapply(
    w.test,
    function(w) c(w.cal, w) / sum(w.cal, w),
    numeric(length(w.cal) + 1)
  )
}


#' Compute quantile of weighted scores
#' function from https://github.com/ryantibs/conformal/
#'
#' @param v scores
#' @param prob quantile level
#' @param w weights
#' @param sorted have the scores and weights been sorted?
#'
#' @return scalar of the quantile at `prob`
#' @export
weighted_quantile = function(v, prob, w = NULL, sorted = FALSE) {
  if (is.null(w)) {
    w = rep(1, length(v))
  }
  if (!sorted) {
    o = order(v)
    v = v[o]
    w = w[o]
  }
  i = which(cumsum(w / sum(w)) >= prob)
  if (length(i) == 0) {
    # Can happen with infinite weights
    return(Inf)
  } else {
    return(v[min(i)])
  }
}


#' Inexact causal conformal inference for generalisation based on Lei & Candés 2021
#'
#' @param x_test Optional covariates for independent predictions.
#' @param range Optional limits for truncating prediction intervals.
#' @param x matrix of outcome predictors
#' @param y vector of outcomes
#' @param s indicator of selection
#' @param a indicator of treatment (binary)
#' @param x_s matrix of selection predictors
#' @param x_a matrix of treatment predictors
#' @param x_tau matrix of treatment effect modifiers
#' @param reweight should the data be re-weighted using IPSW?
#' @param trim should the propensity scores be trimmed?
#' @param trim_perc (only when trim = TRUE) percentile used for trimming
#' @param s.weights vector of selection weights
#' @param seed sets the initial seed
#' @param method T-Learner or DR-learner
#' @param alpha nominal confidence level
#' @param gamma level for quantile regression
#' @param train_perc percentage of data used for training the meta-learner
#' @param K number of folds used for cross-fitting (default K=5)
#' @param out_method method for outcome model
#' @param prop_method method for propensity model
#' @param select_method method for selection of nuisance models
#' @param final_method method for final regression model (DR-Learner only)
#'
#' @return a list of 2-d matrices with the prediction intervals for the inexact
#'  CATE, naive CATE and potential outcomes Y0 and Y1
#' @export
inexact_conformal <- function(
  x,
  y,
  s,
  a,
  x_s = NULL,
  x_a = NULL,
  x_tau = NULL,
  x_test = NULL,
  reweight = FALSE,
  trim = TRUE,
  trim_perc = 0.9,
  s.weights = NULL,
  K = 5,
  seed = NULL,
  method = c("DR", "TL"),
  alpha = 0.05,
  gamma = 0.3,
  train_perc = 0.5,
  range = NULL,
  out_method = "rf",
  prop_method = "rf",
  select_method = "xgboost",
  final_method = "rf"
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(x_a)) {
    x_a <- x
  }
  if (is.null(x_s)) {
    x_s <- x
  }
  if (is.null(x_tau)) {
    x_tau <- x
  }
  if (!check_lengths(x[, 1], y, a, s, x_s[, 1], x_a[, 1])) {
    stop("x, y, a and s must have the same length")
  }
  if (!is.null(x_test) && !setequal(names(x_tau), names(x_test))) {
    stop("variables in x_test must have the same names as variables in x_tau")
  }
  method <- match.arg(method)

  split <- groupdata2::partition(
    data.frame(x, s),
    p = train_perc,
    cat_col = 's',
    list_out = FALSE
  )$.partitions
  d1 <- split == 1
  d2 <- split == 2

  conformal <- split_causal_conformal(
    x = data.frame(x),
    y = y,
    s = s,
    a = a,
    x_s = x_s,
    x_a = data.frame(x_a),
    x_tau = x_tau,
    reweight = reweight,
    train_ind = d1,
    trim = trim,
    trim_perc = trim_perc,
    s.weights = s.weights,
    seed = seed,
    method = method,
    K = K,
    alpha = alpha,
    range = range,
    out_method = out_method,
    prop_method = prop_method,
    select_method = select_method,
    final_method = final_method
  )

  A <- a[d2]
  Y <- y[d2]

  # Infinite quantiles are capped at the arm's largest calibration score so
  # high-weight trial units still inform the quantile regression
  capped_band <- function(mu_hat, q, max_score) {
    q <- ifelse(is.infinite(q), max_score, q)
    band <- cbind(mu_hat - q, mu_hat + q)
    if (!is.null(range)) band <- pmin(pmax(band, range[1]), range[2])
    band
  }
  c_0 <- capped_band(conformal$mu_0_hat, conformal$q_0, conformal$max_score_0)
  c_1 <- capped_band(conformal$mu_1_hat, conformal$q_1, conformal$max_score_1)

  treated_intervals <- cbind(
    (Y - c_0[d2, 2]), # lower: Y1 - Y0_upper
    (Y - c_0[d2, 1]) # upper: Y1 - Y0_lower
  )
  control_intervals <- cbind(
    (c_1[d2, 1] - Y), # lower: Y1_lower - Y0
    (c_1[d2, 2] - Y) # upper: Y1_upper - Y0
  )

  intervals <- matrix(NA, nrow = length(A), ncol = 2)
  intervals[A == 1, ] <- treated_intervals[A == 1, ]
  intervals[A == 0, ] <- control_intervals[A == 0, ]

  # Only trial rows have an observed Y, so only they can calibrate
  trial <- s[d2] == 1
  finite_idx <- trial &
    !is.infinite(intervals[, 1]) &
    !is.infinite(intervals[, 2])

  q_used <- ifelse(A == 1, conformal$q_0[d2], conformal$q_1[d2])
  message(
    sum(trial & is.infinite(q_used)), " of ", sum(trial),
    " trial intervals in the calibration fold capped at the largest calibration score; ",
    sum(finite_idx), " intervals used for quantile regression."
  )

  out <- list(
    c_0 = conformal$c_0,
    c_1 = conformal$c_1,
    c_tau_naive = conformal$c_tau,
    tau_hat = conformal$tau_hat,
    mu_0_hat = conformal$mu_0_hat,
    mu_1_hat = conformal$mu_1_hat
  )

  if (sum(finite_idx) < 10) {
    warning(
      "Too few finite intervals for training quantile regression. Returning naive conformal intervals."
    )
    return(out)
  }

  x_q <- x_tau[d2, , drop = FALSE][finite_idx, , drop = FALSE]
  c_quant_l <- quantregForest::quantregForest(x = x_q, y = intervals[finite_idx, 1])
  c_quant_r <- quantregForest::quantregForest(x = x_q, y = intervals[finite_idx, 2])

  predict_band <- function(newdata) {
    unname(cbind(
      predict(c_quant_l, newdata = newdata, what = gamma),
      predict(c_quant_r, newdata = newdata, what = 1 - gamma)
    ))
  }

  out <- c(list(c_tau = predict_band(x_tau)), out)
  if (!is.null(x_test)) {
    out$c_tau_test <- predict_band(x_test)
  }
  out
}
