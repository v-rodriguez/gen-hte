# Usage: Rscript experiments/exp3/01_exp3.R --config_index N

all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", all_args[grepl("^--file=", all_args)])
if (length(file_arg) != 1L) stop("Run this file with Rscript")
script_dir <- dirname(normalizePath(file_arg, winslash = "/"))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/")

renv::load(repo_root)
library(dplyr)

devtools::load_all(file.path(repo_root, "genhte"), quiet = TRUE)
source(file.path(repo_root, "experiments", "00_config.R"))
source(file.path(repo_root, "experiments", "01_generate_data.R"))
source(file.path(repo_root, "experiments", "utils.R"))

one_rep <- function(i, config) {
  cat("Running config:", config$config_id, "Rep:", i, "\n")
  seed <- config$base_seed + i

  # Training population with an external subsample, and a validation population
  # at a fixed seed
  draw <- function(s, subsample) {
    generate_data(
      ssize = config$ssize,
      xtau = config$xtau,
      xmu = config$xmu,
      delta = config$delta,
      N = config$N,
      ntau = config$ntau,
      nsel = config$nsel,
      nmu = config$nmu,
      seed = s,
      subsample = subsample,
      nsubsample = config$nsubsample
    )
  }
  data <- draw(seed, TRUE)
  data_val <- draw(20, FALSE)

  vars <- colnames(data)[grepl("T|X", colnames(data))]
  t_vars <- colnames(data)[grepl("T", colnames(data))]
  s_vars <- colnames(data)[grepl("X", colnames(data))]
  trial_data <- data |> filter(S == 1)

  common_params <- list(
    trim = TRUE,
    trim_perc = 0.9,
    K = 5,
    seed = seed,
    out_method = "xgb_cv",
    prop_method = "xgboost",
    select_method = "xgboost",
    final_method = "xgb_cv"
  )
  on_trial <- c(common_params, list(y = trial_data$Y, s = trial_data$S, a = trial_data$A))
  on_full <- c(
    common_params,
    list(y = data$Y, s = data$S, a = data$A, subsample = data$subsample)
  )

  score <- function(fit, d) {
    calculate_metrics(fit$train$tau_hat, d$tau, fit$test$tau_hat, data_val$tau)
  }

  dr_mod <- do.call(dr_learner, c(list(
    x = data.frame(trial_data[vars]),
    x_s = trial_data[s_vars],
    x_tau = as.matrix(trial_data[t_vars]),
    x_test = as.matrix(data_val[t_vars]),
    reweight = FALSE
  ), on_trial))

  gdr_mod <- do.call(dr_learner, c(list(
    x = data.frame(data[vars]),
    x_s = data[s_vars],
    x_tau = as.matrix(data[t_vars]),
    x_test = as.matrix(data_val[t_vars]),
    reweight = TRUE
  ), on_full))

  gdr_marg_mod <- do.call(dr_learner, c(list(
    x = data.frame(data[vars]),
    x_s = data[s_vars],
    x_tau = as.matrix(data[t_vars]),
    x_test = as.matrix(data_val[t_vars]),
    reweight = TRUE,
    use_g_step = TRUE
  ), on_full))

  t_mod <- do.call(t_learner, c(list(
    x = data.frame(trial_data[t_vars]),
    x_test = data.frame(data_val[t_vars]),
    reweight = FALSE
  ), on_trial))

  gt_mod <- do.call(t_learner, c(list(
    x = data.frame(data[t_vars]),
    x_s = data[s_vars],
    x_test = data.frame(data_val[t_vars]),
    reweight = TRUE
  ), on_full))

  gt_marg_mod <- do.call(t_learner, c(list(
    x = data.frame(data[vars]),
    x_s = data[s_vars],
    x_test = data.frame(data_val[t_vars]),
    use_g_step = TRUE,
    x_tau = data.frame(data[t_vars]),
    reweight = TRUE
  ), on_full))

  results_frame(
    config,
    i,
    list(
      dr = score(dr_mod, trial_data),
      dr_sub = score(gdr_mod, data),
      dr_sub_marg = score(gdr_marg_mod, data),
      tl = score(t_mod, trial_data),
      tl_sub = score(gt_mod, data),
      tl_sub_marg = score(gt_marg_mod, data)
    ),
    c(
      dr = "DR-Learner",
      dr_sub = "gDR (Subsample)",
      dr_sub_marg = "gDR (Restricted)",
      tl = "T-Learner",
      tl_sub = "gT (Subsample)",
      tl_sub_marg = "gT (Restricted)"
    )
  )
}

run_experiment(get_exp3_configs(), one_rep, script_dir)
