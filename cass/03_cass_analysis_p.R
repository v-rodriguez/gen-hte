path <- this.path::this.path()
setwd(dirname(path))

library(dplyr)
library(fda)
library(tidyr)
library(ggplot2)
library(patchwork)
library(mice)

source("01_analysis_functions.R")
load("data.RData")


df_list <- mice::complete(imputed_param, action = "all")
M <- imputed_param$m

# Which variables do we want in x_tau: 'ejec_frac' (ejection fraction only),
# 'sex' or 'age' (ejection fraction plus that variable)
cli_args <- commandArgs(trailingOnly = TRUE)
x_tau <- if (length(cli_args) > 0) cli_args[1] else 'ejec_frac'
stopifnot(x_tau %in% c('ejec_frac', 'sex', 'age'))

if (x_tau == "ejec_frac") {
  X_grid <- seq(30, 80, 0.1)
  dim <- 3
} else {
  X_grid <- rep(seq(30, 80, 0.1), 2)
  ind <- c(rep(0, 0.5 * length(X_grid)), rep(1, 0.5 * length(X_grid)))
  dim <- 4
}

B <- 200
alpha <- 0.05

subsets <- c(0, 1)

results <- array(NA, dim = c(2, length(X_grid), dim))
standard_errors <- array(NA, dim = c(2, length(X_grid), 2))

conf_set <- list()
point_conf <- list()
unif_conf <- list()

set.seed(88)
for (i in 1:2) {
  results0 <- array(NA, dim = c(M, length(X_grid), dim))
  standard_errors0 <- array(NA, dim = c(M, length(X_grid), 2))

  dr_tvec_quant <- c()
  t_tvec_quant <- c()

  for (m in seq_along(df_list)) {
    print(m)
    df_all <- df_list[[m]] %>%
      mutate(across(prev_mi_hist, ~ ifelse(. == "No", 0, 1))) %>%
      dplyr::select(all_of(c(var_names, "S", "A", "Y")))

    # subset data by history of MI
    df <- df_all %>%
      filter(prev_mi_hist == subsets[i])

    # nuisance models: selection, outcome per arm, treatment
    omega_hat <- predict(train_omega(df), df, type = "response")
    mu0hat <- predict(train_outcome(df, 0), df, type = "response")
    mu1hat <- predict(train_outcome(df, 1), df, type = "response")
    pi_hat <- predict(train_pi(df), df, type = "response")

    S <- df$S
    A <- df$A
    Y <- df$Y

    # construct pseudo-outcomes
    pseud <- (S / omega_hat) *
      ((A - pi_hat) / (pi_hat * (1 - pi_hat))) *
      (Y - A * mu1hat - (1 - A) * mu0hat) +
      mu1hat -
      mu0hat

    # train second stage regression
    dr_fit.mod <- dr_fit_mod(df, pseud, x_tau)
    dr_fit <- dr_fit.mod$mod
    regressors <- dr_fit.mod$mod.mat

    # Grid points
    grid_df <- data.frame(ejec_frac = X_grid)
    if (x_tau == "ejec_frac") {
      regressors.grid <- tau_design(grid_df, x_tau, df$ejec_frac)
    } else {
      regressors.grid <- tau_design(grid_df, x_tau, df$ejec_frac, ind)
      df_ind <- tau_indicator(df, x_tau, df$age)
    }

    # Fitted values using gDR-Learner
    dr_cate <- regressors.grid %*% coef(dr_fit)

    ## Training t-learner

    w1 <- (S / omega_hat) * (1 / (pi_hat))
    w1 <- w1[S == 1 & A == 1]
    w0 <- (S / omega_hat) * (1 / (1 - pi_hat))
    w0 <- w0[S == 1 & A == 0]

    # train weighted outcome models
    mu1_w <- train_outcome(df, 1, w1, x_tau)
    mu0_w <- train_outcome(df, 0, w0, x_tau)

    # Fitted values using gT-Learner
    tl_cate <- plogis(regressors.grid %*% coef(mu1_w)) -
      plogis(regressors.grid %*% coef(mu0_w))

    # Calculate standard errors for gDR-Learner
    dr_cov <- sandwich_est(regressors, pseud, dr_fit)
    dr_se <- get_se(dr_cov, regressors.grid)

    t_df <- bind_cols(Y = df$Y, regressors, S = df$S, A = df$A, )

    subgroup <- switch(x_tau, sex = "sexMale", age = "age_ind")

    # Calculate standard errors for gT-Learner
    t_boot <- t_boot_se(
      t_df,
      dim(regressors)[2] - 1,
      w1,
      w0,
      regressors.grid,
      subgroup = subgroup,
      B = B
    )

    # Get weighted bootstrap fitted values
    dr_boot.vals <- dr_wboot(
      regressors,
      pseud,
      dim(regressors)[2] - 1,
      regressors.grid,
      B = B
    )

    # Calculate max t-statistics
    dr_tvec <- max_tstat(dr_boot.vals, dr_cate[, 1], dr_se)
    t_tvec <- max_tstat(t_boot$boot.vals, tl_cate[, 1], t_boot$cate_se)

    dr_tvec_quant <- c(dr_tvec_quant, quantile(dr_tvec, 1 - alpha, names = FALSE))
    t_tvec_quant <- c(t_tvec_quant, quantile(t_tvec, 1 - alpha, names = FALSE))

    results0[m, , ] <- if (x_tau == "ejec_frac") {
      cbind(dr_cate[, 1], tl_cate[, 1], X_grid)
    } else {
      cbind(dr_cate[, 1], tl_cate[, 1], X_grid, ind)
    }

    standard_errors0[m, , ] <- cbind(dr_se, t_boot$cate_se)
  }

  # conformal prediction on last imputed dataset
  dr_set <- bind_cols(
    causal_conformal_param(df, method = "DR", x_tau = x_tau),
    ejec_frac = df$ejec_frac,
    Method = "gDR-Learner"
  )

  tl_set <- bind_cols(
    causal_conformal_param(df, method = "TL", x_tau = x_tau),
    ejec_frac = df$ejec_frac,
    Method = "gT-Learner"
  )

  if (x_tau != "ejec_frac") {
    dr_set$ind <- df_ind
    tl_set$ind <- df_ind
  }

  conf_set[[i]] <- data.frame(rbind(dr_set, tl_set))

  # pooling point predictions
  results[i, , ] <- apply(results0, c(2, 3), mean)

  # pooling standard errors (Rubin's rules)
  var_within <- apply(standard_errors0^2, c(2, 3), mean)
  var_across <- apply(results0[,, 1:2], c(2, 3), var)
  total_var <- var_within + (1 + 1 / M) * var_across

  standard_errors[i, , ] <- sqrt(total_var)

  # pointwise inference
  dr.lower.point <- results[i, , 1] + qnorm(0.025) * standard_errors[i, , 1]
  dr.upper.point <- results[i, , 1] + qnorm(0.975) * standard_errors[i, , 1]

  t.lower.point <- results[i, , 2] + qnorm(0.025) * standard_errors[i, , 2]
  t.upper.point <- results[i, , 2] + qnorm(0.975) * standard_errors[i, , 2]

  # uniform confidence bands from the max-t quantile, averaged over imputations
  dr_t_quantile <- mean(dr_tvec_quant)
  t_t_quantile <- mean(t_tvec_quant)

  dr_cate.upper <- results[i, , 1] + dr_t_quantile * standard_errors[i, , 1]
  dr_cate.lower <- results[i, , 1] - dr_t_quantile * standard_errors[i, , 1]

  tl_cate.upper <- results[i, , 2] + t_t_quantile * standard_errors[i, , 2]
  tl_cate.lower <- results[i, , 2] - t_t_quantile * standard_errors[i, , 2]

  methods <- rep(c("gDR-Learner", "gT-Learner"), each = length(X_grid))

  point_conf[[i]] <- data.frame(
    lower.point = c(dr.lower.point, t.lower.point),
    upper.point = c(dr.upper.point, t.upper.point),
    Method = methods,
    ejec_frac = rep(X_grid, 2)
  )

  unif_conf[[i]] <- data.frame(
    lower.unif = c(dr_cate.lower, tl_cate.lower),
    upper.unif = c(dr_cate.upper, tl_cate.upper),
    Method = methods,
    ejec_frac = rep(X_grid, 2)
  )

  if (x_tau != "ejec_frac") {
    point_conf[[i]][, "ind"] <- rep(ind, 2)
    unif_conf[[i]][, "ind"] <- rep(ind, 2)
  }
}

# create plots
plots <- list()

colors <- c(`gDR-Learner` = "#0072B2", `gT-Learner` = "#D55E00")
titles <- c("No History of MI", "History of MI")

for (i in 1:2) {
  res <- as.data.frame(results[i, , ])
  names(res) <- c("gDR-Learner", "gT-Learner", "ejec_frac", if (x_tau != "ejec_frac") "ind")

  res <- res %>%
    pivot_longer(
      cols = c(`gDR-Learner`, `gT-Learner`),
      names_to = "Method",
      values_to = "pCATE"
    )

  param_res <- full_join(res, full_join(point_conf[[i]], unif_conf[[i]]))

  # Faceting (sex/age) needs param_res as the plot data
  plots[[i]] <- ggplot(if (x_tau != "ejec_frac") param_res) +
    geom_hline(aes(yintercept = 0)) +
    geom_line(aes(x = ejec_frac, y = pCATE, colour = Method), data = param_res, linewidth = 1) +
    geom_line(aes(x = ejec_frac, y = lower.point, colour = Method), data = param_res, lty = "dashed") +
    geom_line(aes(x = ejec_frac, y = upper.point, colour = Method), data = param_res, lty = "dashed") +
    geom_ribbon(
      aes(x = ejec_frac, ymin = lower.unif, ymax = upper.unif, fill = Method),
      alpha = 0.1,
      data = param_res %>% filter(Method == "gDR-Learner")
    ) +
    theme_minimal(base_size = 22, base_family = "serif") +
    scale_color_manual(values = colors) +
    coord_cartesian(ylim = c(-1, 1)) +
    ggtitle(titles[i]) +
    scale_fill_manual(values = colors) +
    xlab("Ejection Fraction") +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    guides(fill = "none")

  if (x_tau != "ejec_frac") {
    labels_names <- if (x_tau == "sex") {
      c(`1` = "Male", `0` = "Female")
    } else {
      c(`0` = "Age <= 50", `1` = "Age >50")
    }
    plots[[i]] <- plots[[i]] +
      facet_grid(rows = vars(ind), labeller = labeller(ind = labels_names))
  }
}


p1 <- plots[[1]] +
  plots[[2]] +
  plot_layout(axis_titles = "collect", guides = "collect")

filename <- paste0("plots/", Sys.Date(), "/CATE", x_tau, "B", B)

ggsave(
  paste0(filename, " param se.pdf"),
  p1,
  width = 12,
  height = ifelse(x_tau == "ejec_frac", 7, 8),
  dpi = 350,
  create.dir = TRUE
)
