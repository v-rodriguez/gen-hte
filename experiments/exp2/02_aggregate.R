all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", all_args[grepl("^--file=", all_args)])
if (length(file_arg) != 1L) stop("Run this file with Rscript")
script_dir <- dirname(normalizePath(file_arg, winslash = "/"))
source(file.path(script_dir, "..", "00_config.R"))
source(file.path(script_dir, "..", "02_aggregate_utils.R"))

aggregate_experiment(
  script_dir,
  experiment = 2,
  configs = get_exp2_configs(),
  metadata = c("dimension", "complexity", "nmu", "ssize"),
  expected_models = c("tl", "tl_ipw", "tl_ipw_z", "dr", "dr_ipw", "dr_ipw_z")
)
