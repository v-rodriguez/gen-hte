#' Calculate Performance Metrics
#'
#' Computes training and validation metrics for model evaluation
#'
#' @param train_pred Predicted treatment effects on training set
#' @param train_true True treatment effects on training set
#' @param test_pred Predicted treatment effects on test set
#' @param test_true True treatment effects on test set
#' @return List of performance metrics (MSE, RMSE, bias, SE, etc.)
calculate_metrics <- function(train_pred, train_true, test_pred, test_true) {
  train_error <- train_pred - train_true
  test_error <- test_pred - test_true

  list(
    # Training metrics
    train_mse = mean(train_error^2),
    train_rmse = sqrt(mean(train_error^2)),
    train_bias = mean(train_error),
    train_abs_bias = mean(abs(train_error)),
    train_se = sd(train_error),
    # Validation metrics
    mse = mean(test_error^2),
    rmse = sqrt(mean(test_error^2)),
    bias = mean(test_error),
    abs_bias = mean(abs(test_error)),
    median_mse = median(test_error^2),
    median_bias = median(test_error),
    median_abs_bias = median(abs(test_error)),
    se = sd(test_error)
  )
}

#' Calculate Conformal Inference Metrics
#'
#' Computes coverage and interval width metrics for conformal prediction intervals
#'
#' @param intervals n x 2 matrix with lower and upper bounds
#' @param true_tau Vector of true treatment effects
#' @return List with coverage rate, average width, and percent infinite
calculate_conformal_metrics <- function(intervals, true_tau) {
  # Check coverage: true value within interval
  coverage <- mean(
    true_tau >= intervals[, 1] & true_tau <= intervals[, 2],
    na.rm = TRUE
  )
  widths <- intervals[, 2] - intervals[, 1]
  avg_width <- mean(widths[is.finite(widths)], na.rm = TRUE)
  perc_infinite <- mean(!is.finite(widths), na.rm = TRUE)

  list(
    coverage = coverage,
    avg_width = avg_width,
    perc_infinite = perc_infinite
  )
}

#' Run one configuration of a simulation experiment
#'
#' Parses --config_index, runs one_rep(i, config) for every replication and
#' writes one CSV to <script_dir>/results.
#'
#' @param configs Data frame of configurations (one row per config)
#' @param one_rep Function (rep index, config row) returning a data frame
#' @param script_dir Directory of the calling runner script
run_experiment <- function(configs, one_rep, script_dir) {
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--config_index", type = "integer", default = NULL,
      help = "Configuration index to run", metavar = "integer")
  )))
  config_index <- opt$config_index
  if (is.null(config_index)) stop("--config_index is required")
  if (config_index < 1L || config_index > nrow(configs)) {
    stop("Invalid config_index: ", config_index,
      " (valid range: 1-", nrow(configs), ")")
  }
  config <- configs[config_index, ]
  reps <- config$nreps

  output_dir <- file.path(script_dir, "results")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Processing configuration:", config$config_id, "\n")

  if (parallelly::supportsMulticore() && reps > 1L) {
    future::plan(future::multicore,
      workers = min(reps, max(1L, parallelly::availableCores() - 1L)))
  } else {
    future::plan(future::sequential)
  }
  results <- furrr::future_map(
    seq_len(reps),
    function(i) one_rep(i, config),
    .options = furrr::furrr_options(seed = TRUE)
  ) |> dplyr::bind_rows()
  future::plan(future::sequential)

  output_file <- file.path(output_dir, paste0(config$config_id, ".csv"))
  write.csv(results, output_file, row.names = FALSE)
  cat("Saved:", output_file, "\n")
}

#' Stack per-model metric lists into one results data frame
#'
#' @param config Config row (supplies config_id)
#' @param i Replication index
#' @param metrics Named list (by model code) of metric lists
#' @param labels Named character vector mapping model code to method label
results_frame <- function(config, i, metrics, labels) {
  dplyr::bind_rows(Map(
    function(model, m) {
      data.frame(
        config_id = config$config_id,
        rep = i,
        model = model,
        method = labels[[model]],
        m
      )
    },
    names(metrics),
    metrics
  ))
}
