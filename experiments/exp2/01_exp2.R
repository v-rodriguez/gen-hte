# Usage: Rscript experiments/exp2/01_exp2.R --config_index N

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

  # Training population, and a validation population at a fixed seed
  draw <- function(s) {
    generate_data(
      ssize = config$ssize,
      xtau = config$xtau,
      xmu = config$xmu,
      delta = config$delta,
      N = config$N,
      ntau = config$ntau,
      nsel = config$nsel,
      nmu = config$nmu,
      seed = s
    )
  }
  data <- draw(seed)
  data_val <- draw(20)

  # X_sel: selection variables; X_tau: treatment effect modifiers;
  # Z: prognostic variables observed only in the trial (S=1)
  x_sel_vars <- grep("^X", colnames(data), value = TRUE)[1:config$nsel]
  x_tau_vars <- grep("^T", colnames(data), value = TRUE)
  z_prog_vars <- grep("^Z", colnames(data), value = TRUE)
  all_vars <- c(x_sel_vars, x_tau_vars, z_prog_vars)
  no_z_vars <- c(x_sel_vars, x_tau_vars)

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
  on_full <- c(common_params, list(y = data$Y, s = data$S, a = data$A))

  score <- function(fit, d) {
    calculate_metrics(fit$train$tau_hat, d$tau, fit$test$tau_hat, data_val$tau)
  }

  # DR (baseline): (X_sel, X_tau) only, trial data
  dr_mod <- do.call(dr_learner, c(list(
    x = data.frame(trial_data[no_z_vars]),
    x_s = trial_data[x_sel_vars],
    x_tau = as.matrix(trial_data[x_tau_vars]),
    x_test = as.matrix(data_val[x_tau_vars]),
    reweight = FALSE
  ), on_trial))

  # gDR: (X_sel, X_tau), IPSW reweighting
  gdr_mod <- do.call(dr_learner, c(list(
    x = data.frame(data[no_z_vars]),
    x_s = data[x_sel_vars],
    x_tau = as.matrix(data[x_tau_vars]),
    x_test = as.matrix(data_val[x_tau_vars]),
    reweight = TRUE
  ), on_full))

  # gDR (Enr.): fit mu_a(X_sel, X_tau, Z), marginalise Z via g_a(X_sel, X_tau),
  # then regress on X_tau
  gdr_marg_mod <- do.call(dr_learner, c(list(
    x = data.frame(data[all_vars]),
    x_s = data[x_sel_vars],
    x_tau = as.matrix(data[x_tau_vars]),
    x_test = as.matrix(data_val[x_tau_vars]),
    x_g = data[no_z_vars],
    reweight = TRUE,
    use_g_step = TRUE
  ), on_full))

  # T (baseline), trial data
  t_mod <- do.call(t_learner, c(list(
    x = data.frame(trial_data[x_tau_vars]),
    x_test = data.frame(data_val[x_tau_vars]),
    reweight = FALSE
  ), on_trial))

  # gT: IPSW reweighting
  gt_mod <- do.call(t_learner, c(list(
    x = data.frame(data[x_tau_vars]),
    x_s = data[x_sel_vars],
    x_test = data.frame(data_val[x_tau_vars]),
    reweight = TRUE
  ), on_full))

  # gT (Enr.): as gDR (Enr.)
  gt_marg_mod <- do.call(t_learner, c(list(
    x = data.frame(data[all_vars]),
    x_s = data[x_sel_vars],
    x_tau = data.frame(data[x_tau_vars]),
    x_test = data.frame(data_val[x_tau_vars]),
    x_g = data[no_z_vars],
    reweight = TRUE,
    use_g_step = TRUE
  ), on_full))

  results_frame(
    config,
    i,
    list(
      dr = score(dr_mod, trial_data),
      dr_ipw = score(gdr_mod, data),
      dr_ipw_z = score(gdr_marg_mod, data),
      tl = score(t_mod, trial_data),
      tl_ipw = score(gt_mod, data),
      tl_ipw_z = score(gt_marg_mod, data)
    ),
    c(
      dr = "DR-Learner",
      dr_ipw = "gDR-Learner",
      dr_ipw_z = "gDR (Enr.)",
      tl = "T-Learner",
      tl_ipw = "gT-Learner",
      tl_ipw_z = "gT (Enr.)"
    )
  )
}

run_experiment(get_exp2_configs(), one_rep, script_dir)
