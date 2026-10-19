# Usage: Rscript experiments/exp1/01_exp1.R --config_index N

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
  # Independent noise per arm for ITE coverage, on its own seed so the main RNG
  # stream is unchanged
  e_cf <- withr::with_seed(21, rnorm(NROW(data_val)))
  ite_indep <- data_val$ite +
    ifelse(data_val$A == 1, data_val$e - e_cf, e_cf - data_val$e)

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

  # DR fits the outcome model on all covariates and regresses on X_tau; T uses
  # X_tau throughout
  run_model <- function(method, reweight) {
    dr <- method == "DR"
    d <- if (reweight) data else trial_data
    x_vars <- if (dr) vars else t_vars
    in_train <- c(rep(TRUE, NROW(d)), rep(FALSE, NROW(data_val)))
    on_d <- c(common_params, list(y = d$Y, s = d$S, a = d$A))

    point <- if (dr) {
      do.call(dr_learner, c(list(
        x = data.frame(d[vars]),
        x_s = d[s_vars],
        x_tau = as.matrix(d[t_vars]),
        x_test = as.matrix(data_val[t_vars]),
        reweight = reweight
      ), on_d))
    } else {
      do.call(t_learner, c(list(
        x = data.frame(d[t_vars]),
        x_s = d[s_vars],
        x_test = data.frame(data_val[t_vars]),
        reweight = reweight
      ), on_d))
    }

    split_conf <- do.call(split_causal_conformal, c(list(
      x = bind_rows(data.frame(d[x_vars]), data.frame(data_val[x_vars])),
      y = c(d$Y, data_val$Y),
      s = c(d$S, data_val$S),
      a = c(d$A, data_val$A),
      x_s = bind_rows(d[s_vars], data_val[s_vars]),
      x_tau = if (dr) as.matrix(rbind(d[t_vars], data_val[t_vars])),
      train_ind = in_train,
      reweight = reweight,
      method = method
    ), common_params))

    inexact_conf <- do.call(inexact_conformal, c(list(
      x = data.frame(d[x_vars]),
      y = d$Y,
      s = d$S,
      a = d$A,
      x_s = d[s_vars],
      x_tau = if (dr) d[t_vars],
      x_test = data.frame(data_val[t_vars]),
      reweight = reweight,
      method = method
    ), common_params))

    split_metrics <- calculate_conformal_metrics(
      split_conf$c_tau[!in_train, , drop = FALSE],
      data_val$tau
    )
    inexact_metrics <- calculate_conformal_metrics(
      inexact_conf$c_tau_test,
      data_val$tau
    )
    split_ite <- calculate_conformal_metrics(
      split_conf$c_tau[!in_train, , drop = FALSE],
      ite_indep
    )
    inexact_ite <- calculate_conformal_metrics(
      inexact_conf$c_tau_test,
      ite_indep
    )

    c(
      calculate_metrics(point$train$tau_hat, d$tau, point$test$tau_hat, data_val$tau),
      list(
        coverage_split_95 = split_metrics$coverage,
        avg_width_split_95 = split_metrics$avg_width,
        perc_infinite_split = split_metrics$perc_infinite,
        coverage_inexact_95 = inexact_metrics$coverage,
        avg_width_inexact_95 = inexact_metrics$avg_width,
        perc_infinite_inexact = inexact_metrics$perc_infinite,
        coverage_split_95_ite = split_ite$coverage,
        coverage_inexact_95_ite = inexact_ite$coverage
      )
    )
  }

  results_frame(
    config,
    i,
    list(
      dr = run_model("DR", FALSE),
      dr_ipw = run_model("DR", TRUE),
      tl = run_model("TL", FALSE),
      tl_ipw = run_model("TL", TRUE)
    ),
    c(
      dr = "DR-Learner",
      dr_ipw = "gDR-Learner",
      tl = "T-Learner",
      tl_ipw = "gT-Learner"
    )
  )
}

run_experiment(get_exp1_configs(), one_rep, script_dir)
