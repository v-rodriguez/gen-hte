# Default SuperLearner library for outcome and second-stage regressions
SL_LIB <- c("SL.glm", "SL.xgboost", "SL.ranger", "SL.mean")

# Fit a gradient-boosted model.
#
# Wraps xgb.train rather than xgboost(): from xgboost 2.0 onwards xgboost()
# renames its arguments and returns a booster whose predict method rejects
# xgb.DMatrix input. xgb.train keeps the same signature across 1.7.x and 3.x and
# fits the same model.
fit_xgb <- function(x, y, objective, nrounds = 20, weight = NULL) {
  dtrain <- xgboost::xgb.DMatrix(
    data = as.matrix(x),
    label = y,
    weight = weight
  )

  xgboost::xgb.train(
    params = list(objective = objective),
    data = dtrain,
    nrounds = nrounds,
    verbose = 0
  )
}

# Forest with min.node.size chosen by out-of-bag log-loss (probability) or MSE.
# Large leaves keep 1/omega and 1/pi from blowing up on a few overfit rows.
RF_CV_NODE <- c(20, 50, 100, 200)
RF_CV_MTRY <- 1 / 3
RF_CV_NODE_FINAL <- c(100, 200, 400, 800)

fit_rf_cv <- function(x, y, weights = NULL, probability = TRUE, nodes = RF_CV_NODE) {
  if (probability) y <- as.factor(y)
  fit <- function(node, mtry) {
    ranger::ranger(
      y = y, x = x, probability = probability, num.trees = 200, case.weights = weights,
      min.node.size = node, mtry = max(1, floor(mtry * NCOL(x))), num.threads = 1
    )
  }
  grid <- expand.grid(node = nodes, mtry = RF_CV_MTRY)
  fits <- Map(fit, grid$node, grid$mtry)
  oob_loss <- vapply(fits, function(f) {
    if (!probability) return(f$prediction.error)
    p <- f$predictions[cbind(seq_along(y), as.integer(y))]
    -mean(log(pmax(p, 1e-6)), na.rm = TRUE)
  }, numeric(1))
  fits[[which.min(oob_loss)]]
}

# Treatment propensity, fit on trial members only.
treat_prop <- function(
  x,
  a,
  s,
  method = "rf",
  SL.lib = c("SL.glm", "SL.xgboost", "SL.mean", "SL.gam")
) {
  x <- subset(x, s == 1)
  a <- a[s == 1]
  if (method == "rf") {
    ranger::ranger(y = as.factor(a), x = x, probability = TRUE)
  } else if (method == "rf_cv") {
    fit_rf_cv(x, a)
  } else if (method == "superlearner") {
    SuperLearner::SuperLearner(
      Y = a,
      X = x,
      family = binomial(),
      SL.library = SL.lib
    )
  } else if (method == "xgboost") {
    fit_xgb(x, a, "binary:logistic")
  } else {
    stop("method not recognized")
  }
}

# Study participation propensity.
#
# Rows from the external subsample alone estimate P(S=1 | x_s, S'=1), which is
# omega when S' is independent of S given x_s, but see few S=1 cases. Rows pooled
# over R = {S=1 or S'=1} bring the whole trial to bear but estimate
# q = omega / (omega + (1 - omega)p). Supply p_subsample = P(S'=1) with pooled
# rows and the returned "study_prop_pooled" object inverts q back to omega on
# predict, via omega = qp / (1 - q(1 - p)).
study_prop <- function(
  x,
  s,
  method = "xgboost",
  SL.library = SL_LIB,
  p_subsample = NULL
) {
  fit <- train_selection_model(x, s, method, SL.library)

  if (is.null(p_subsample)) {
    return(fit)
  }

  if (!is.numeric(p_subsample) || length(p_subsample) != 1) {
    stop("p_subsample must be a single number")
  }
  if (!is.finite(p_subsample) || p_subsample <= 0 || p_subsample > 1) {
    stop("p_subsample must lie in (0, 1]; got ", p_subsample)
  }

  structure(
    list(fit = fit, p_subsample = p_subsample),
    class = "study_prop_pooled"
  )
}

train_selection_model <- function(x, s, method, SL.library) {
  # A constant response silently yields a model that predicts that constant
  # everywhere, which is indistinguishable from a fitted propensity downstream.
  if (length(unique(s)) < 2) {
    stop(
      "selection model requires both S=0 and S=1 in the training data; ",
      "got only S=",
      unique(s),
      " across ",
      length(s),
      " rows"
    )
  }

  if (method == "xgboost") {
    fit_xgb(x, s, "binary:logistic")
  } else if (method == "superlearner") {
    SuperLearner::SuperLearner(
      Y = s,
      X = as.data.frame(x),
      family = binomial(),
      SL.library = SL.library
    )
  } else if (method == "rf") {
    ranger::ranger(y = as.factor(s), x = as.data.frame(x), probability = TRUE)
  } else if (method == "rf_cv") {
    fit_rf_cv(as.data.frame(x), s)
  } else {
    stop("method not recognized")
  }
}

# Outcome model for arm ind_A, fit on trial members only. w are optional case
# weights (T-learner IPSW); rows with NA weights are dropped.
outcome_train <- function(x, y, a, s, ind_A, w = NULL, method = "rf") {
  x <- subset(x, s == 1)
  y <- y[s == 1]
  a <- a[s == 1]

  probability <- length(unique(y)) == 2

  if (any(is.na(w))) {
    x <- subset(x, !is.na(w))
    y <- y[!is.na(w)]
    a <- a[!is.na(w)]
    w <- w[!is.na(w)]
  }

  arm <- a == ind_A
  x <- subset(x, arm)
  y <- y[arm]
  w <- w[arm]

  if (method == "rf") {
    ranger::ranger(x = x, y = y, case.weights = w, probability = probability)
  } else if (method == "superlearner") {
    SuperLearner::SuperLearner(
      Y = y,
      X = x,
      family = if (probability) binomial() else gaussian(),
      SL.library = SL_LIB,
      obsWeights = w
    )
  } else if (method == "xgboost") {
    fit_xgb(
      x,
      y,
      if (probability) "binary:logistic" else "reg:squarederror",
      weight = w
    )
  } else if (method == "xgb_stack" && !probability) {
    fit_xgb_stack(x, y, w)
  } else if (method == "xgb_cv" && !probability) {
    fit_xgb_cv(x, y, w)
  } else if (method == "xgb_outcome" && !probability) {
    # deliberately underfit; base_score pinned so xgboost versions agree
    xgboost::xgb.train(
      params = list(
        objective = "reg:squarederror", eta = 0.1, max_depth = 2,
        base_score = stats::weighted.mean(y, if (is.null(w)) rep(1, length(y)) else w),
        nthread = 1
      ),
      data = xgboost::xgb.DMatrix(as.matrix(x), label = y, weight = w),
      nrounds = 5,
      verbose = 0
    )
  } else {
    stop("method not recognized")
  }
}

final_regression <- function(phi, x_tau, method = "rf", weights = NULL) {
  if (method == "rf") {
    ranger::ranger(y = phi, x = x_tau, case.weights = weights)
  } else if (method == "superlearner") {
    SuperLearner::SuperLearner(
      Y = phi,
      X = as.data.frame(x_tau),
      family = gaussian(),
      SL.library = SL_LIB,
      obsWeights = weights
    )
  } else if (method == "xgboost") {
    fit_xgb(x_tau, phi, "reg:squarederror", weight = weights)
  } else if (method == "xgb_cv") {
    fit_xgb_cv(x_tau, phi, weights)
  } else if (method == "rf_cv") {
    fit_rf_cv(x_tau, phi, weights, probability = FALSE, nodes = RF_CV_NODE_FINAL)
  } else if (method == "xgb_cv_tuned") {
    fit_xgb_cv(x_tau, phi, weights, mcw = c(10, 30, 100))
  } else if (method == "xgb_stack") {
    fit_xgb_stack(x_tau, phi, weights)
  } else {
    stop("method not recognized")
  }
}

# Convex NNLS combination of cross-validated xgboost fits over depth x rounds,
# plus the mean
XGB_STACK_DEPTHS <- c(1, 2, 4, 6)
XGB_STACK_ROUNDS <- c(50, 200, 600)

fit_xgb_stack <- function(x, y, weights = NULL, V = 5) {
  x <- as.matrix(x)
  w <- if (is.null(weights)) rep(1, length(y)) else weights
  fit_depth <- function(rows, depth) {
    xgboost::xgb.train(
      params = list(
        objective = "reg:squarederror", eta = 0.03, max_depth = depth,
        min_child_weight = 50, subsample = 0.8, nthread = 1
      ),
      data = xgboost::xgb.DMatrix(x[rows, , drop = FALSE], label = y[rows], weight = w[rows]),
      nrounds = max(XGB_STACK_ROUNDS),
      verbose = 0
    )
  }
  library_preds <- function(fits, mean_y, newx) {
    cbind(do.call(cbind, lapply(fits, xgb_stack_rounds, newx = newx)), mean_y)
  }

  fold <- sample(rep(seq_len(V), length.out = length(y)))
  z <- matrix(NA_real_, length(y), length(XGB_STACK_DEPTHS) * length(XGB_STACK_ROUNDS) + 1)
  for (v in seq_len(V)) {
    tr <- fold != v
    fits_v <- lapply(XGB_STACK_DEPTHS, fit_depth, rows = tr)
    z[!tr, ] <- library_preds(fits_v, stats::weighted.mean(y[tr], w[tr]), x[!tr, , drop = FALSE])
  }
  coef <- nnls::nnls(z * sqrt(w), y * sqrt(w))$x
  coef <- if (sum(coef) > 0) coef / sum(coef) else c(rep(0, ncol(z) - 1), 1)

  structure(
    list(
      fits = lapply(XGB_STACK_DEPTHS, fit_depth, rows = rep(TRUE, length(y))),
      mean_y = stats::weighted.mean(y, w),
      coef = coef,
      library_preds = library_preds
    ),
    class = "xgb_stack"
  )
}

# xgboost >= 2 treats the end of iterationrange as inclusive, 1.7.x as exclusive
xgb_stack_rounds <- function(fit, newx) {
  inclusive <- utils::packageVersion("xgboost") >= "2.0.0"
  vapply(XGB_STACK_ROUNDS, function(r) {
    predict(fit, newx, iterationrange = c(1, if (inclusive) r else r + 1))
  }, numeric(NROW(newx)))
}

# Shallow, regularised xgboost with rounds chosen by 5-fold CV; with several
# min_child_weight values (mcw), the one with the lowest CV error is kept
fit_xgb_cv <- function(x, y, weights = NULL, mcw = 100) {
  # min_child_weight is in weight units
  if (!is.null(weights)) weights <- weights / mean(weights)
  dtrain <- xgboost::xgb.DMatrix(as.matrix(x), label = y, weight = weights)
  cv_fit <- function(m) {
    params <- list(
      objective = "reg:squarederror", eta = 0.03, max_depth = 2,
      # capped so small fits can still split
      min_child_weight = min(m, NROW(x) / 4), subsample = 0.8, lambda = 10, nthread = 1,
      base_score = stats::weighted.mean(y, if (is.null(weights)) rep(1, length(y)) else weights)
    )
    cv <- xgboost::xgb.cv(
      params, dtrain,
      nrounds = 1000, nfold = 5, early_stopping_rounds = 30, verbose = 0,
      # stratified folds fail on a constant y
      stratified = FALSE
    )
    # xgboost 1.7.x stores the best round on the object, >= 2 under early_stop
    best <- if (!is.null(cv$best_iteration)) cv$best_iteration else cv$early_stop$best_iteration
    best <- max(1, best)
    list(params = params, best = best, err = cv$evaluation_log$test_rmse_mean[best])
  }
  fits <- lapply(mcw, cv_fit)
  f <- fits[[which.min(vapply(fits, `[[`, numeric(1), "err"))]]
  xgboost::xgb.train(f$params, dtrain, nrounds = f$best, verbose = 0)
}
