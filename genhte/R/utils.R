#' Should a regression carry IPSW weights?
#'
#' omega is a function of x_s, so weighting a regression whose regressors already
#' contain x_s leaves its estimand unchanged and only costs effective sample size.
#' Weight only when x_s is absent from the regressors, where the weights perform
#' the marginalisation over the target law of x_s. Unrecognisable names fall back
#' to weighting: over-weighting costs variance, under-weighting costs consistency.
#'
#' @param x_s selection covariates
#' @param regressors the covariates the regression is conditioning on
#' @noRd
needs_ipsw <- function(x_s, regressors) {
  xs_names <- colnames(x_s)
  xr_names <- colnames(regressors)
  is.null(xs_names) || is.null(xr_names) || !all(xs_names %in% xr_names)
}

#' Rows the selection model is fit on
#'
#' Without an external subsample, omega is fit on the whole training half of the
#' fold. With one, the choice is between the pooled rows R = \{S=1 or S'=1\}, which
#' bring the entire trial to bear but make the regression q rather than omega, and
#' the subsample rows alone, whose only S=1 rows are the subsample members who
#' also joined the trial. Pooling is the default; see \code{study_prop} for the
#' inversion it requires.
#'
#' @param d1 logical, the training half of the fold
#' @param s study participation indicator
#' @param subsample external subsample indicator, or NULL when there is none
#' @param pool fit on the pooled rows rather than the subsample alone
#' @noRd
selection_rows <- function(d1, s, subsample, pool) {
  if (is.null(subsample)) {
    return(d1)
  }
  if (pool) d1 & (s == 1 | subsample == 1) else d1 & subsample == 1
}

check_lengths <- function(...) length(unique(lengths(list(...)))) == 1

check_subsample <- function(subsample, s) {
  if (is.null(subsample)) {
    return(invisible())
  }
  if (length(subsample) != length(s)) {
    stop("subsample must have the same length as other vectors")
  }
  if (!all(subsample %in% c(0, 1))) {
    stop("subsample must be a binary (0/1) vector")
  }
  if (!any(s == 1 | subsample == 1)) {
    stop("At least some observations must be in trial or subsample")
  }
}

get_predictions <- function(object, newdata) {
  # A pooled selection fit predicts q = P(S=1 | x_s, S=1 or S'=1), not omega.
  # Invert here so every call site downstream sees the omega scale. The map is
  # increasing in q, so trimming the result still trims the same observations.
  if (inherits(object, "study_prop_pooled")) {
    q <- pmin(pmax(get_predictions(object$fit, newdata), 0), 1)
    p <- object$p_subsample
    return(q * p / pmax(1 - q * (1 - p), .Machine$double.eps))
  }

  if (inherits(object, "SuperLearner")) {
    preds <- predict(object, newdata)$pred
  } else if (inherits(object, "xgb_stack")) {
    z <- object$library_preds(object$fits, object$mean_y, as.matrix(newdata))
    preds <- drop(z %*% object$coef)
  } else if (inherits(object, "xgb.Booster")) {
    preds <- predict(object, as.matrix(newdata))
  } else {
    preds <- predict(object, newdata)$predictions
  }

  if (is.matrix(preds) & NCOL(preds) == 2) {
    preds <- preds[, 2]
  }
  preds
}

# Hajek-stabilised residual weights S I(A=a) / (omega pi_a), each arm rescaled
# to mean 1
dr_weights <- function(A, pi_hat, S = 1, ex_hat = 1) {
  w_1 <- S * A / (pi_hat * ex_hat)
  w_0 <- S * (1 - A) / ((1 - pi_hat) * ex_hat)
  w_1 <- w_1 / mean(w_1)
  w_0 <- w_0 / mean(w_0)
  list(w = w_1 - w_0, w_1 = w_1, w_0 = w_0)
}

# Clip x at its (1 - alpha) and alpha quantiles
truncate_vector <- function(x, alpha = 0.9) {
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1 (exclusive)")
  }
  upper <- quantile(x, alpha, na.rm = TRUE, names = FALSE)
  lower <- quantile(x, 1 - alpha, na.rm = TRUE, names = FALSE)
  pmin(pmax(x, lower), upper)
}
