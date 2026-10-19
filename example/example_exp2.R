
all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", all_args[grepl("^--file=", all_args)])
if (length(file_arg) != 1L) stop("Run this file with Rscript")
script_dir <- dirname(normalizePath(file_arg, winslash = "/"))
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/")
renv::load(repo_root)

library(dplyr)
library(ggplot2)

devtools::load_all(file.path(repo_root, "genhte"), quiet = TRUE)

source(file.path(repo_root, "experiments", "00_config.R"))
source(file.path(repo_root, "experiments", "01_generate_data.R"))
source(file.path(repo_root, "experiments", "utils.R"))

reps <- 50
configs <- get_exp2_configs()
config <- configs[configs$config_id == "low_complex_nmu6_ssize1500", , drop = FALSE]
stopifnot(nrow(config) == 1L)

cat("Configuration:", config$config_id, "\n")
cat("Replications:", reps, "\n\n")


one_rep <- function(i) {
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
      seed = s,
      subsample = FALSE
    )
  }
  data <- draw(seed)
  data_val <- draw(20)

  x_sel_vars <- grep("^X", colnames(data), value = TRUE)[1:config$nsel]
  x_tau_vars <- grep("^T", colnames(data), value = TRUE)
  z_prog_vars <- grep("^Z", colnames(data), value = TRUE)

  all_vars <- c(x_sel_vars, x_tau_vars, z_prog_vars)
  no_z_vars <- c(x_sel_vars, x_tau_vars)

  trial_data <- data |> filter(S == 1)

  common <- list(
    trim = TRUE,
    trim_perc = 0.9,
    K = 5,
    seed = seed,
    out_method = "xgb_cv",
    prop_method = "xgboost",
    select_method = "xgboost",
    final_method = "xgb_cv"
  )
  on_trial <- c(common, list(y = trial_data$Y, s = trial_data$S, a = trial_data$A))
  on_full <- c(common, list(y = data$Y, s = data$S, a = data$A))

  rmse <- function(fit) {
    sqrt(mean((fit$test$tau_hat - data_val$tau)^2))
  }

  # DR-Learner: trial data only, sCATE
  dr <- do.call(dr_learner, c(list(
    x = data.frame(trial_data[no_z_vars]),
    x_s = trial_data[x_sel_vars],
    x_tau = as.matrix(trial_data[x_tau_vars]),
    x_test = as.matrix(data_val[x_tau_vars]),
    reweight = FALSE
  ), on_trial))

  # gDR-Learner: IPSW reweighting to the target population
  gdr <- do.call(dr_learner, c(list(
    x = data.frame(data[no_z_vars]),
    x_s = data[x_sel_vars],
    x_tau = as.matrix(data[x_tau_vars]),
    x_test = as.matrix(data_val[x_tau_vars]),
    reweight = TRUE
  ), on_full))

  # gDR (Enr.): fit mu_a on (X, T, Z), marginalize over X_P, then regress on X_tau
  gdr_z <- do.call(dr_learner, c(list(
    x = data.frame(data[all_vars]),
    x_s = data[x_sel_vars],
    x_tau = as.matrix(data[x_tau_vars]),
    x_test = as.matrix(data_val[x_tau_vars]),
    x_g = data[no_z_vars],
    reweight = TRUE,
    use_g_step = TRUE
  ), on_full))

  tl <- do.call(t_learner, c(list(
    x = data.frame(trial_data[x_tau_vars]),
    x_test = data.frame(data_val[x_tau_vars]),
    reweight = FALSE
  ), on_trial))

  gt <- do.call(t_learner, c(list(
    x = data.frame(data[x_tau_vars]),
    x_s = data[x_sel_vars],
    x_test = data.frame(data_val[x_tau_vars]),
    reweight = TRUE
  ), on_full))

  gt_z <- do.call(t_learner, c(list(
    x = data.frame(data[all_vars]),
    x_s = data[x_sel_vars],
    x_tau = data.frame(data[x_tau_vars]),
    x_test = data.frame(data_val[x_tau_vars]),
    x_g = data[no_z_vars],
    reweight = TRUE,
    use_g_step = TRUE
  ), on_full))

  fits <- list(
    "DR-Learner" = dr,
    "gDR-Learner" = gdr,
    "gDR (Enr.)" = gdr_z,
    "T-Learner" = tl,
    "gT-Learner" = gt,
    "gT (Enr.)" = gt_z
  )

  data.frame(
    rep = i,
    method = names(fits),
    rmse = vapply(fits, rmse, numeric(1)),
    row.names = NULL
  )
}

results <- bind_rows(lapply(seq_len(reps), function(i) {
  cat("rep", i, "of", reps, "\n")
  one_rep(i)
}))

method_levels <- c(
  "T-Learner",
  "gT-Learner",
  "gT (Enr.)",
  "DR-Learner",
  "gDR-Learner",
  "gDR (Enr.)"
)
results$method <- factor(results$method, levels = method_levels)

summary_table <- results |>
  group_by(method) |>
  summarise(
    mean_rmse = mean(rmse),
    se = sd(rmse) / sqrt(n()),
    median_rmse = median(rmse),
    .groups = "drop"
  )

cat("\nRMSE against the true population CATE, over", reps, "replications:\n\n")
print(as.data.frame(summary_table), row.names = FALSE, digits = 3)

p <- ggplot(results, aes(x = method, y = rmse, fill = method)) +
  geom_boxplot(alpha = 0.8, show.legend = FALSE) +
  labs(
    title = "Population CATE estimation error",
    subtitle = paste0(config$config_id, ", ", reps, " replications"),
    x = NULL,
    y = "RMSE"
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(script_dir, "example_results.png"), p, width = 7, height = 4.5, dpi = 150)
