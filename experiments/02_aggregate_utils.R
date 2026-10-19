# Shared validation and aggregation for the three local simulation runners.

aggregate_experiment <- function(script_dir, experiment, configs, metadata, expected_models) {
  result_dir <- file.path(script_dir, "results")
  result_files <- list.files(result_dir, pattern = "\\.csv$", full.names = TRUE)
  if (!length(result_files)) stop("No result CSV files found in ", result_dir)

  required <- c(
    "config_id", "rep", "model", "method", "mse", "bias", "abs_bias",
    "rmse", "median_mse", "median_bias", "median_abs_bias"
  )
  if (experiment == 1L) {
    required <- c(required, "coverage_split_95", "avg_width_split_95",
      "perc_infinite_split", "coverage_inexact_95", "avg_width_inexact_95",
      "perc_infinite_inexact", "coverage_split_95_ite", "coverage_inexact_95_ite")
  }
  parts <- lapply(result_files, function(path) {
    rows <- read.csv(path, stringsAsFactors = FALSE)
    missing_columns <- setdiff(required, names(rows))
    if (!nrow(rows)) stop("Empty result file: ", path)
    if (length(missing_columns)) {
      stop("Result file ", path, " is missing columns: ",
        paste(missing_columns, collapse = ", "))
    }
    rows
  })
  all_results <- dplyr::bind_rows(parts)
  keys <- all_results[c("config_id", "rep", "model")]
  if (anyNA(keys) || any(!nzchar(all_results$config_id)) ||
      any(!nzchar(all_results$model))) stop("Result keys must not be missing or empty")
  if (!is.numeric(all_results$rep) || any(!is.finite(all_results$rep)) ||
      any(all_results$rep != floor(all_results$rep))) {
    stop("Replication IDs must be finite integers")
  }
  if (anyDuplicated(keys)) stop("Duplicate config_id/rep/model rows in result files")

  unknown <- setdiff(unique(all_results$config_id), configs$config_id)
  if (length(unknown)) stop("Unknown config_id values: ", paste(unknown, collapse = ", "))

  unexpected_models <- setdiff(unique(all_results$model), expected_models)
  if (length(unexpected_models)) {
    stop("Unknown model values: ", paste(unexpected_models, collapse = ", "))
  }
  # With no duplicate keys, a config is complete exactly when every
  # (rep, model) pair in 1..nreps x expected_models is present.
  validation <- data.frame(
    config_id = configs$config_id,
    expected_reps = configs$nreps,
    n_rows = vapply(configs$config_id, function(id) sum(all_results$config_id == id), integer(1), USE.NAMES = FALSE)
  )
  validation$is_complete <- with(
    validation,
    n_rows == expected_reps * length(expected_models)
  ) & vapply(seq_len(nrow(validation)), function(i) {
    reps <- all_results$rep[all_results$config_id == validation$config_id[i]]
    all(reps %in% seq_len(validation$expected_reps[i]))
  }, logical(1))

  if (any(!validation$is_complete)) {
    print(validation[!validation$is_complete, ])
    stop("Experiment ", experiment, " results are incomplete; final CSV was not written")
  }

  idx <- match(all_results$config_id, configs$config_id)
  for (column in metadata) all_results[[column]] <- configs[[column]][idx]

  final_dir <- file.path(script_dir, "final")
  dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(
    all_results,
    file.path(final_dir, paste0("exp", experiment, "_results.csv")),
    row.names = FALSE
  )
  write.csv(validation, file.path(final_dir, "validation_summary.csv"), row.names = FALSE)
  cat("Saved complete Experiment", experiment, "results to", final_dir, "\n")
}
