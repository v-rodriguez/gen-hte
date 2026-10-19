SCENARIO_COLUMNS <- c("xtau", "xmu", "ntau", "nsel", "nmu")

add_fixed_parameters <- function(configs) {
  configs$N <- 15000
  configs$nreps <- 1000
  configs$base_seed <- 9243
  configs
}

#' Setting 1: trial size sweep for all four scenarios at baseline overlap
#' (delta = 1 low-dim, 5 high-dim), plus a delta sweep for the high-dim
#' scenarios at n = 1000.
get_exp1_configs <- function() {
  scenarios <- data.frame(
    scenario = c("LD_linear", "LD_nonlinear", "HD_linear", "HD_nonlinear"),
    xtau = c("simple", "complex", "simple", "complex"),
    xmu = "nonlinear",
    ntau = c(1, 1, 5, 5),
    nsel = c(2, 2, 10, 10),
    nmu = c(2, 2, 10, 10)
  )

  size_sweep <- expand.grid(
    scenario = scenarios$scenario,
    ssize = c(500, 1000, 1500, 2000, 2500),
    stringsAsFactors = FALSE
  )
  size_sweep$delta <- ifelse(grepl("^LD", size_sweep$scenario), 1L, 5L)

  delta_sweep <- expand.grid(
    scenario = c("HD_linear", "HD_nonlinear"),
    ssize = 1000,
    delta = 6:10,
    stringsAsFactors = FALSE
  )

  configs <- merge(rbind(size_sweep, delta_sweep), scenarios, by = "scenario")
  configs <- add_fixed_parameters(configs)
  configs$config_id <- sprintf(
    "%s_delta%d_ssize%d",
    configs$scenario,
    configs$delta,
    configs$ssize
  )

  configs <- configs[order(configs$scenario, configs$delta, configs$ssize), ]
  rownames(configs) <- NULL
  configs[c("config_id", "scenario", SCENARIO_COLUMNS, "ssize", "delta", "N", "nreps", "base_seed")]
}

#' Setting 2: sweep the number of trial-only prognostic covariates
#' d_P = nmu - nsel over {0, 2, 4, 6, 8} at n = 1500, and the trial size at
#' d_P = 4, for each dimension and outcome model.
get_exp2_configs <- function() {
  scenarios <- expand.grid(
    dimension = c("low", "high"),
    complexity = c("simple", "complex"),
    stringsAsFactors = FALSE
  )
  scenarios$xtau <- "simple"
  scenarios$xmu <- ifelse(scenarios$complexity == "simple", "simple", "nonlinear")
  scenarios$ntau <- ifelse(scenarios$dimension == "low", 1, 5)
  scenarios$nsel <- ifelse(scenarios$dimension == "low", 2, 10)

  prognostic_sweep <- merge(
    scenarios[c("dimension", "complexity", "nsel")],
    data.frame(d_p = c(0, 2, 4, 6, 8))
  )
  prognostic_sweep$ssize <- 1500

  size_sweep <- merge(
    scenarios[c("dimension", "complexity", "nsel")],
    data.frame(ssize = c(500, 1000, 1500, 2000, 2500))
  )
  size_sweep$d_p <- 4

  configs <- unique(rbind(prognostic_sweep, size_sweep))
  configs$nmu <- configs$nsel + configs$d_p
  configs <- merge(
    configs[c("dimension", "complexity", "nmu", "ssize")],
    scenarios,
    by = c("dimension", "complexity")
  )
  configs$delta <- ifelse(configs$dimension == "low", 1, 5)
  configs <- add_fixed_parameters(configs)
  configs$config_id <- sprintf(
    "%s_%s_nmu%d_ssize%d",
    configs$dimension,
    configs$complexity,
    configs$nmu,
    configs$ssize
  )

  configs <- configs[
    order(configs$dimension, configs$complexity, configs$nmu, configs$ssize),
  ]
  rownames(configs) <- NULL
  configs[c(
    "config_id", "dimension", "complexity", SCENARIO_COLUMNS, "delta", "ssize",
    "N", "nreps", "base_seed"
  )]
}

#' Setting 3: sweep the external subsample size at n = 2000 and the trial size
#' at n_sub = 5000, for the two non-linear scenarios.
get_exp3_configs <- function() {
  scenarios <- data.frame(
    scenario = c("LD_nonlinear", "HD_nonlinear"),
    xtau = "complex",
    xmu = "nonlinear",
    ntau = c(2, 5),
    nsel = c(2, 10),
    nmu = c(2, 10),
    delta = c(1, 9)
  )

  subsample_sweep <- expand.grid(
    scenario = scenarios$scenario,
    ssize = 2000,
    nsubsample = c(500, 1000, 2000, 5000, 10000),
    stringsAsFactors = FALSE
  )
  size_sweep <- expand.grid(
    scenario = scenarios$scenario,
    ssize = c(500, 1000, 1500, 2000, 2500),
    nsubsample = 5000,
    stringsAsFactors = FALSE
  )

  configs <- merge(unique(rbind(subsample_sweep, size_sweep)), scenarios, by = "scenario")
  configs <- add_fixed_parameters(configs)
  configs$config_id <- sprintf(
    "%s_ssize%d_nsubsample%d",
    configs$scenario,
    configs$ssize,
    configs$nsubsample
  )

  configs <- configs[order(configs$scenario, configs$ssize, configs$nsubsample), ]
  rownames(configs) <- NULL
  configs[c(
    "config_id", "scenario", SCENARIO_COLUMNS, "ssize", "nsubsample", "delta",
    "N", "nreps", "base_seed"
  )]
}
