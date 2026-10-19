path <- this.path::this.path()
setwd(dirname(path))
load("data.RData")
source("01_analysis_functions.R")

devtools::load_all(paste0(dirname(getwd()), "/genhte"))

library(dplyr)
library(mice)
library(caret)
library(ggplot2)
library(tidyr)
library(patchwork)
library(furrr)

plan(future.callr::callr, workers = availableCores() - 1)

# Which variables do we want in x_tau: "ejec_frac" (ejection fraction only),
# "sex" or "age" (ejection fraction plus that variable)
cli_args <- commandArgs(trailingOnly = TRUE)
tau_opt <- if (length(cli_args) > 0) cli_args[1] else "ejec_frac"
stopifnot(tau_opt %in% c("ejec_frac", "sex", "age"))

M <- imputed$m
results <- data.frame()

set.seed(88)
for (m in 1:M) {
  d_imp <- mice::complete(imputed, action = m)

  dmv <- caret::dummyVars(" ~ .", data = d_imp, fullRank = TRUE)
  d_imp <- data.frame(predict(dmv, newdata = d_imp))

  cor <- cor(d_imp[, colSums(is.na(d_imp)) == 0])
  cor <- cor[, colSums(is.na(cor)) == 0]

  cor_cut <- findCorrelation(cor, cutoff = 0.9, exact = TRUE)

  df <- d_imp[, -cor_cut] %>%
    dplyr::select(-"vessel_disease")
  names(df)[(names(df)) == "sex.Male"] <- "sex"

  x <- df %>%
    dplyr::select(contains(var_names))
  x_a <- df %>%
    dplyr::select(
      age,
      contains("angina"),
      ejec_frac,
      sys_bp,
      prox_lad_obstr,
      lv_wall_score
    )

  reps <- 50
  x_grid <- seq(30, 80, by = 0.1)
  N_test = length(x_grid)

  if (tau_opt == "ejec_frac") {
    x_tau <- df %>%
      mutate(hist_of_mi = prev_mi_hist.Single | prev_mi_hist.Multiple) %>%
      dplyr::select(ejec_frac, hist_of_mi)
    x_test <- bind_cols(
      ejec_frac = rep(x_grid, 2),
      hist_of_mi = c(rep(FALSE, N_test), rep(TRUE, N_test))
    )
  } else {
    x_tau <- df %>%
      mutate(hist_of_mi = prev_mi_hist.Single | prev_mi_hist.Multiple) %>%
      dplyr::select(all_of(c("ejec_frac", "hist_of_mi", tau_opt)))

    if (tau_opt == "age") {
      x_tau$age <- ifelse(x_tau$age <= median(x_tau$age), 0, 1)
    }
    x_test <- bind_cols(
      ejec_frac = rep(x_grid, 4),
      hist_of_mi = c(rep(FALSE, N_test * 2), rep(TRUE, N_test * 2)),
      tau_opt = rep(c(0, 1), N_test * 2)
    )
    names(x_test)[length(names(x_test))] <- tau_opt
  }

  # Both learners fit mu_a on the full covariates and regress onto x_tau with
  # final_method: the pseudo-outcomes for gDR, the g-step for gT
  args_list <- list(
    x = x,
    x_tau = x_tau,
    y = df$Y,
    s = df$S,
    a = df$A,
    x_s = x,
    x_a = x_a,
    x_test = x_test,
    reweight = TRUE,
    trim = FALSE,
    prop_method = "rf_cv",
    out_method = "rf",
    select_method = "rf_cv",
    final_method = "rf_cv"
  )

  tl_est <- repeat_function_parallel(
    func = t_learner,
    args = c(args_list, use_g_step = TRUE),
    n = reps
  )

  dr_est <- repeat_function_parallel(dr_learner, args_list, reps)

  print(m)
  results <- rbind(
    results,
    cbind(
      `gT-Learner` = tl_est$tau_hat,
      `gDR-Learner` = dr_est$tau_hat,
      x_test,
      m
    )
  )
}

colors <- c(`gDR-Learner` = "#0072B2", `gT-Learner` = "#D55E00")

collect_res <- results %>%
  group_by_at(vars("ejec_frac", "hist_of_mi", all_of(tau_opt))) %>%
  summarise(across(c(`gT-Learner`, `gDR-Learner`), mean))

collect_res_piv <- collect_res %>%
  pivot_longer(
    cols = c(`gT-Learner`, `gDR-Learner`),
    names_to = "Method",
    values_to = "tau_hat"
  )

if (tau_opt == "ejec_frac") {
  collect_res_sum <- collect_res %>%
    group_by(hist_of_mi)
} else {
  collect_res_sum <- collect_res %>%
    group_by_at(vars("hist_of_mi", all_of(tau_opt)))
}


collect_res_sum <- collect_res_sum %>%
  summarise(across(c(`gT-Learner`, `gDR-Learner`), mean)) %>%
  pivot_longer(
    cols = c(`gT-Learner`, `gDR-Learner`),
    names_to = "Method",
    values_to = "mean"
  )


mi_labels <- c(`TRUE` = "History of MI", `FALSE` = "No history of MI")

seed <- 240

args_list$x_test <- NULL
args_list$x <- x_tau
args_list$train_perc <- 0.9
args_list$range <- c(0, 1)
args_list$gamma <- 0.25


set.seed(seed)
tl_ci <- data.frame(do.call(
  inexact_conformal,
  args = append(args_list, c(method = "TL"))
))
tl_ci <- bind_cols(tl_ci, x_tau) %>%
  mutate(Method = "gT-Learner")

args_list$x <- x

set.seed(seed)
dr_ci <- data.frame(do.call(
  inexact_conformal,
  args = append(args_list, c(method = "DR"))
))
dr_ci <- bind_cols(dr_ci, x_tau) %>%
  mutate(Method = "gDR-Learner")

ci_res <- bind_rows(tl_ci, dr_ci)

# replace infinite values
ci_res <- ci_res %>%
  group_by(Method) %>%
  mutate(across(
    c(c_tau.1, c_tau.2),
    ~ case_when(
      . == -Inf ~ min(.[is.finite(.)]),
      . == Inf ~ max(.[is.finite(.)]),
      TRUE ~ .
    )
  ))


p2 <- ggplot(data = collect_res_piv) +
  geom_hline(aes(yintercept = 0)) +
  geom_line(aes(y = tau_hat, x = ejec_frac, color = Method), linewidth = 0.8) +
  geom_hline(
    data = collect_res_sum,
    aes(yintercept = mean, color = Method),
    lty = "dashed"
  ) +
  facet_grid(
    cols = vars(hist_of_mi),
    labeller = labeller(hist_of_mi = mi_labels)
  ) +
  theme_minimal(base_size = 22, base_family = "serif") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  xlim(c(30, 80)) +
  xlab("Ejection Fraction") +
  ylab("pCATE") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

if (tau_opt != "ejec_frac") {
  labels_names <- if (tau_opt == "sex") {
    c(`1` = "Male", `0` = "Female")
  } else {
    c(`0` = "Age <= 50", `1` = "Age >50")
  }
  p2 <- p2 +
    facet_grid(
      cols = vars(hist_of_mi),
      rows = vars(.data[[tau_opt]]),
      labeller = labeller(hist_of_mi = mi_labels, !!tau_opt := labels_names)
    )
}

# p1 adds the inexact conformal bands
p1 <- p2 +
  geom_ribbon(
    aes(x = ejec_frac, ymin = c_tau.1, ymax = c_tau.2, fill = Method),
    data = ci_res,
    alpha = 0.2
  )

file_name <- paste0(
  "plots/",
  Sys.Date(),
  "/CATE comparison_inexact_constrained",
  " ",
  tau_opt,
  ".pdf"
)
file_name_2 <- paste0("plots/CATE points", " ", tau_opt, ".pdf")

if (tau_opt == "ejec_frac") {
  ggsave(
    p1,
    file = file_name,
    create.dir = TRUE,
    width = 12,
    height = 7
  )
}

ggsave(
  p2,
  file = file_name_2,
  create.dir = TRUE,
  width = 12,
  height = ifelse(tau_opt == "ejec_frac", 7, 8)
)
