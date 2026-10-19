# B-spline settings for the ejection-fraction (and age) terms
knot <- 2
norder <- 3

# B-spline basis of x (intercept column dropped), knots at the quantiles of ref
bs_term <- function(x, ref = x) {
  bsplineS(x, breaks = quantile(ref, c(0:knot) / knot), norder = norder)[, -1]
}


tau_formula <- function(lhs, x_tau, env = parent.frame()) {
  rhs <- switch(
    x_tau,
    ejec_frac = "bs_term(ejec_frac)",
    sex = "sex + sex * bs_term(ejec_frac)",
    "age_ind + age_ind * bs_term(ejec_frac)"
  )
  as.formula(paste(lhs, "~", rhs), env = env)
}


# Run func(args) n times in parallel and return the per-observation median of
# the test-set CATE predictions.
repeat_function_parallel <- function(func, args = list(), n) {
  results <- future_map(
    1:n,
    function(i) do.call(func, args),
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
  )
  preds <- sapply(results, function(r) r$test$tau_hat)
  data.frame(tau_hat = apply(preds, 1, median))
}


sandwich_est <- function(X, y, model, est = NULL) {
  # calculate Huber-White sandwich estimator
  N <- length(y)

  Q <- t(X) %*% X / N

  if (is.null(est)) {
    est <- X %*% coef(model)
  }
  resid <- (y - est)^2

  # Explicit loop: crossprod() sums in a different order and shifts the SEs
  U <- matrix(0, dim(X)[2], dim(X)[2])
  for (i in 1:N) {
    U <- U + X[i, ] %*% t(X[i, ]) * resid[i]
  }
  U <- U / N
  solve(Q) %*% U %*% solve(Q) / N
}


get_se <- function(hw_est, pred.grid) {
  # standard error on a grid of values
  sqrt(diag(pred.grid %*% hw_est %*% t(pred.grid)))
}

dr_wboot <- function(X, y, dim, pred.grid, B = 200) {
  # Algorithm 2 & 6, weighted bootstrap for gDR-Learner
  n <- length(y)
  coef.b <- matrix(0, B, dim + 1)
  for (b in 1:B) {
    v <- rexp(n)
    v <- v / sum(v)
    coef.b[b, ] <- coef(lm(y ~ X - 1, weights = v))
  }

  pred.grid %*% t(coef.b)
}

t_boot_se <- function(
  df,
  dim,
  ipsw1,
  ipsw0,
  pred.grid,
  subgroup = NULL,
  B = 200
) {

  df$ipsw <- ifelse(df$A == 1, ipsw1, ipsw0)
  df1 <- df %>% filter(A == 1 & S == 1)
  df0 <- df %>% filter(A == 0 & S == 1)

  # Bootstrap one treatment arm, stratified by subgroup.
  # If the logistic fit fit fails to converge its sample is redrawn with sample_n.
  resample <- function(d, retry) {
    d <- d %>% group_by(across(all_of(subgroup)))
    d <- if (retry) {
      sample_n(d, nrow(d), replace = TRUE)
    } else {
      slice_sample(d, prop = 1, replace = TRUE)
    }
    dplyr::select(d, -c("A", "S"))
  }

  fit_arm <- function(sample) {
    X <- as.matrix(sample %>% dplyr::select(-c("Y", "ipsw")))
    coef(glm(sample$Y ~ X - 1, weights = sample$ipsw, family = "binomial"))
  }

  coef.b1 <- matrix(0, B, dim + 1)
  coef.b0 <- matrix(0, B, dim + 1)
  for (b in 1:B) {
    retry <- FALSE
    repeat {
      sample1 <- resample(df1, retry)
      sample0 <- resample(df0, retry)
      coef.b1[b, ] <- fit_arm(sample1)
      coef.b0[b, ] <- fit_arm(sample0)
      if (!anyNA(c(coef.b1[b, ], coef.b0[b, ]))) break
      retry <- TRUE
    }
  }

  # fitted values on grid points, and their standard errors over the B replicates
  boot.vals <- plogis(pred.grid %*% t(coef.b1)) -
    plogis(pred.grid %*% t(coef.b0))
  list(boot.vals = boot.vals, cate_se = apply(boot.vals, 1, sd))
}

max_tstat <- function(boot.vals, cate.vals, se) {
  # maximum absolute t-statistic over the grid, per bootstrap replicate
  apply(abs((boot.vals - cate.vals) / se), 2, max)
}


train_omega <- function(df) {
  # selection into the study
  glm(
    S ~ bs_term(age) + bs_term(ejec_frac) + . - age - ejec_frac - A,
    data = df,
    family = "binomial"
  )
}

train_pi <- function(df) {
  # treatment propensity, fit on trial members
  glm(
    A ~ age + angina + ejec_frac + sys_bp + prox_lad_obstr + lv_wall_score,
    data = subset(df, df$S == 1),
    family = "binomial"
  )
}

train_outcome <- function(df, g, w = NULL, x_tau) {
  # outcome model for arm g, fit on trial members

  df_s <- df %>%
    filter(S == 1)

  if (!is.null(w)) {
    # weighted outcome regression 
    df_g <- subset(df_s, df_s$A == g)
    age_ind <- ifelse(df_g$age <= median(df$age), 0, 1)
    return(glm(
      tau_formula("Y", x_tau),
      data = df_g,
      family = "binomial",
      weights = w
    ))
  }

  glm(
    Y ~ bs_term(age) + bs_term(ejec_frac) + . - age - ejec_frac - S - A,
    data = subset(df_s, A == g),
    family = "binomial"
  )
}

p_dr_learner <- function(df, contrast = TRUE, g = NULL, x_tau = "ejec_frac") {
  if (!contrast & is.null(g)) {
    stop("g cannot be NULL when not computing the contrast")
  }

  S <- df$S
  A <- df$A
  Y <- df$Y

  # Train nuisance functions
  o_mod <- train_omega(df)
  pi_mod <- train_pi(df)

  omega_hat <- predict(o_mod, df, type = "response")
  pi_hat <- predict(pi_mod, df, type = "response")

  if (contrast) {
    mu0hat <- predict(train_outcome(df, 0), df, type = "response")
    mu1hat <- predict(train_outcome(df, 1), df, type = "response")

    pseud <- (S / omega_hat) *
      ((A - pi_hat) / (pi_hat * (1 - pi_hat))) *
      (Y - A * mu1hat - (1 - A) * mu0hat) +
      mu1hat -
      mu0hat
  } else {
    mu_hat <- predict(train_outcome(df, g), df, type = "response")

    pseud <- (S / omega_hat) *
      ((g / pi_hat) + ((1 - g) / (1 - pi_hat))) *
      (Y - mu_hat) +
      mu_hat
  }

  # Second-stage regression 
  dr_fit_mod(df, pseud, x_tau)$mod
}

dr_fit_mod <- function(df, pseud, x_tau) {
  # Second-stage regression model
  age_ind <- ifelse(df$age <= median(df$age), 0, 1)
  dr_fit <- lm(tau_formula("pseud", x_tau), data = df)
  list(mod = dr_fit, mod.mat = model.matrix(dr_fit))
}

p_t_learner <- function(df, contrast = TRUE, g = NULL, x_tau = "ejec_frac") {
  if (!contrast & is.null(g)) {
    stop("g cannot be NULL when not computing the contrast")
  }

  S <- df$S
  A <- df$A

  # Train nuisance models
  omega_mod <- train_omega(df)
  pi_mod <- train_pi(df)

  omega_hat <- predict(omega_mod, df, type = "response")
  pi_hat <- predict(pi_mod, df, type = "response")

  if (contrast) {
    # ipsw weights and weighted outcome regressions per arm
    w1 <- ((S / omega_hat) * (1 / (pi_hat)))[S == 1 & A == 1]
    w0 <- ((S / omega_hat) * (1 / (1 - pi_hat)))[S == 1 & A == 0]
    return(list(
      train_outcome(df, 1, w1, x_tau),
      train_outcome(df, 0, w0, x_tau)
    ))
  }

  w <- (S / omega_hat) * ((g / pi_hat) + ((1 - g) / (1 - pi_hat)))
  train_outcome(df, g, w[S == 1 & A == g], x_tau)
}

tau_design <- function(d, x_tau, ref_ef, ind = NULL) {
  if (x_tau == "ejec_frac") {
    cbind(rep(1, nrow(d)), bs_term(d$ejec_frac, ref_ef))
  } else {
    model.matrix(~ ind + ind * bs_term(ejec_frac, ref_ef), data = d)
  }
}

# 0/1 indicator for the sex or age interaction of x_tau, age split at the
# median of ref_age
tau_indicator <- function(d, x_tau, ref_age) {
  if (x_tau == "sex") {
    as.numeric(d$sex) - 1
  } else {
    ifelse(d$age <= median(ref_age), 0, 1)
  }
}

conformal_param_g <- function(
  df,
  g = 0,
  method = c("DR", "TL"),
  x_tau = "ejec_frac",
  alpha = 0.05
) {
  method <- match.arg(method)

  # create training and calibration splits
  split <- groupdata2::partition(
    df,
    p = c(0.8),
    cat_col = "S",
    list_out = FALSE
  )$.partitions

  train <- df[split == 1, ]
  cal <- df[split == 2 & df$A == g, ]

  # test and calibration predictors
  if (x_tau == "ejec_frac") {
    cal.x <- tau_design(cal, x_tau, df$ejec_frac)
    test.x <- tau_design(df, x_tau, df$ejec_frac)
  } else {
    cal.x <- tau_design(cal, x_tau, df$ejec_frac, tau_indicator(cal, x_tau, df$age))
    test.x <- tau_design(df, x_tau, df$ejec_frac, tau_indicator(df, x_tau, df$age))
  }

  # train model
  if (method == "DR") {
    mod <- p_dr_learner(train, FALSE, g, x_tau)
  } else {
    mod <- p_t_learner(train, FALSE, g, x_tau)
  }

  # fitted values on calibration and test sets
  cal.pred <- cal.x %*% coef(mod)
  test.pred <- test.x %*% coef(mod)

  if (method == "TL") {
    cal.pred <- plogis(cal.pred)
    test.pred <- plogis(test.pred)
  }

  # train nuisance models and predict on calibration and test sets
  pi_mod <- train_pi(train)
  omega_mod <- train_omega(train)

  omega.cal <- predict(omega_mod, cal, type = "response")
  omega.test <- predict(omega_mod, df, type = "response")
  pi.cal <- predict(pi_mod, cal, type = "response")
  pi.test <- predict(pi_mod, df, type = "response")

  # ipsw weights
  w.cal <- (1 / omega.cal) * ((g / pi.cal) + (1 - g) / (1 - pi.cal))
  w.test <- (1 / omega.test) * ((g / pi.test) + (1 - g) / (1 - pi.test))

  weights <- genhte::normalise_weights(w.cal, w.test)

  # calibration scores
  cal.score <- abs(cal$Y - cal.pred)

  # prediction set from the quantile of the weighted score distribution
  q <- apply(weights, 2, function(w) {
    genhte::weighted_quantile(c(cal.score, Inf), prob = 1 - alpha, w = w)
  })
  set <- cbind(pmin(test.pred - q, test.pred + q), pmax(test.pred + q, test.pred - q))

  list(
    ejec_frac = df$ejec_frac,
    pred = test.pred,
    set = replace_infs(set),
    score = cal.score
  )
}


causal_conformal_param <- function(
  df,
  method = c("DR", "TL"),
  x_tau = "ejec_frac",
  alpha = 0.025
) {
  set.seed(120)
  conf_set0 <- conformal_param_g(df, 0, method = method, x_tau = x_tau, alpha = alpha)
  conf_set1 <- conformal_param_g(df, 1, method = method, x_tau = x_tau, alpha = alpha)
  set0 <- conf_set0$set
  set1 <- conf_set1$set

  list(
    lower = pmin(set1[, 1] - set0[, 2], set1[, 2] - set0[, 1]),
    upper = pmax(set1[, 1] - set0[, 2], set1[, 2] - set0[, 1]),
    lower0 = set0[, 1],
    upper0 = set0[, 2],
    pred0 = conf_set0$pred,
    lower1 = set1[, 1],
    upper1 = set1[, 2],
    pred1 = conf_set1$pred
  )
}


# replace -Inf with the minimum of the column and Inf with the maximum to get finite values
replace_infs <- function(mat) {
  for (j in 1:ncol(mat)) {
    col <- mat[, j]
    col[col == -Inf] <- min(col[!is.infinite(col)], na.rm = TRUE)
    col[col == Inf] <- max(col[!is.infinite(col)], na.rm = TRUE)
    mat[, j] <- col
  }
  return(mat)
}

