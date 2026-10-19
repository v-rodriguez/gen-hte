all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", all_args[grepl("^--file=", all_args)])
if (length(file_arg) != 1L) stop("Run this file with Rscript")
script_dir <- dirname(normalizePath(file_arg, winslash = "/"))
source(file.path(script_dir, "..", "00_config.R"))
source(file.path(script_dir, "..", "02_aggregate_utils.R"))

aggregate_experiment(
  script_dir,
  experiment = 3,
  configs = get_exp3_configs(),
  metadata = c("scenario", "ssize", "nsubsample"),
  expected_models = c("tl", "dr", "tl_sub", "dr_sub", "tl_sub_marg", "dr_sub_marg")
)
