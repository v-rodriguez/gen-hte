# DR-Learners for Generalizing CATEs from a Randomized Controlled Trial

Code for the paper **"DR-Learners for Generalizing CATEs from a Randomized Controlled Trial"** by Vanessa Rodriguez, Brieuc Lehmann, and Karla Diaz-Ordaz.

---

## Overview

Estimating heterogeneous treatment effects from a randomized controlled trial (RCT) is challenging when the trial population differs from the target population of interest. This paper develops and evaluates causal machine learning methods for **generalizing conditional average treatment effects (CATEs)** from an RCT nested within a cohort drawn from the target population.

## Installation

### Option 1: Install `genhte` from GitHub

```r
library(devtools)
install_github(repo = "v-rodriguez/gen-hte", subdir = "genhte")
```

### Option 2: Clone and install locally

```r
# Clone the repository, then from its root directory:
library(devtools)
devtools::install_local("genhte")
```

### Reproducing the exact package environment

This project uses [`renv`](https://rstudio.github.io/renv/) to record all package versions. To restore the full environment used in the paper:

```r
renv::restore()
```

---

## Quick Start

### gDR-Learner and gT-Learner

```r
library(genhte)

# data: data frame with columns for outcomes (Y), treatment (A), study indicator (S),
#       covariates (X_S, X_tau, X_P)

# Estimate the pCATE using the gDR-Learner:
results <- dr_learner(
  x        = data[c("X1", "X2", "T1")],   # outcome model covariates
  y        = data$Y,
  s        = data$S,                        # S=1 trial, S=0 population
  a        = data$A,
  x_s      = data[c("X1", "X2")],          # X_S: selection model covariates
  x_tau    = data[c("T1")],                 # X_tau: treatment effect modifiers
  reweight = TRUE,                          # enable IPSW generalization
  trim     = TRUE,
  seed     = 123
)

# results$tau_hat: population CATE predictions on the supplied data

# The gT-Learner takes the same arguments
```

Setting `reweight = FALSE` and passing only the trial rows gives the corresponding DR-Learner and T-Learner, which target the sCATE.

### gDR (Enr.) and gT (Enr.)

When prognostic variables $X_P$ (columns `Z*`) are measured in the trial but not in the target population, add them to `x` and marginalize them out with the g-step. `x_g` is the covariate set the outcome model is marginalized onto, here $X_S \cup X_\tau$:

```r
results <- dr_learner(
  x          = data[c("X1", "X2", "T1", "Z1", "Z2")],  # includes X_P
  y          = data$Y,
  s          = data$S,
  a          = data$A,
  x_s        = data[c("X1", "X2")],
  x_tau      = data[c("T1")],
  x_g        = data[c("X1", "X2", "T1")],  # marginalize X_P out -> g_a(X_S, X_tau)
  reweight   = TRUE,
  use_g_step = TRUE,
  seed       = 123
)

# gT (Enr.) takes the same arguments
```

### gDR (Subsample), gT (Subsample), and Restricted

```r
# gDR (Subsample): pooled selection model, outcome model left unmarginalized
results <- dr_learner(
  x              = data[c("X1", "X2", "T1")],
  y              = data$Y,
  s              = data$S,
  a              = data$A,
  x_s            = data[c("X1", "X2")],
  x_tau          = data[c("T1")],
  subsample      = data$subsample,   
  pool_selection = TRUE,
  reweight       = TRUE,
  seed           = 123
)

# gDR (Restricted): additionally marginalize the outcome model onto X_tau, so
# X_S is never needed outside the pooled sample. x_g defaults to x_tau.
results <- dr_learner(
  x              = data[c("X1", "X2", "T1")],
  y              = data$Y,
  s              = data$S,
  a              = data$A,
  x_s            = data[c("X1", "X2")],
  x_tau          = data[c("T1")],
  subsample      = data$subsample,
  pool_selection = TRUE,
  use_g_step     = TRUE,             # -> g_a(X_tau)
  reweight       = TRUE,
  seed           = 123
)
```

`t_learner()` takes the same arguments for gT (Subsample) and gT (Restricted).

### Prediction Intervals

`inexact_conformal()` constructs prediction intervals for the CATE using inexact conformal inference (Lei & Candès 2021). It supports both the DR-Learner (`method = "DR"`) and T-Learner (`method = "TL"`).

**Note:** this function only works in the standard setting where the target population has all covariates needed for prediction. 

```r
intervals <- inexact_conformal(
  x        = data[c("X1", "X2", "T1")],
  y        = data$Y,
  s        = data$S,
  a        = data$A,
  x_s      = data[c("X1", "X2")],
  x_tau    = data[c("T1")],
  reweight = TRUE,
  trim     = TRUE,
  method   = "DR",   # or "TL" for T-Learner
  alpha    = 0.1,    # 90% prediction intervals
  seed     = 123
)

# intervals$c_tau:       [lower, upper] inexact CATE prediction intervals
# intervals$c_tau_naive: naive (split conformal) CATE intervals
# intervals$c_0:         potential outcome Y(0) intervals
# intervals$c_1:         potential outcome Y(1) intervals
# intervals$tau_hat:     point estimates
```

---

## Running the Example

`example/example_exp2.R` is a self-contained script for Setting 2 in the paper. It runs 50 replications of a low-dimensional simulation comparing all six methods, and produces a summary table and boxplot.

**Prerequisites:** install packages via `renv::restore()`.

```bash
Rscript example/example_exp2.R
```

## Reproducing the Simulation Study

The three settings use the current designs and data-generating processes from
the paper. Full runs use many configurations and replications and are
computationally expensive. Each invocation runs one configuration locally and
writes one ignored CSV file beneath that setting's `results/` directory. 

For Setting 1, run every configuration index returned by `get_exp1_configs()`, then validate
and combine the results:
```bash
Rscript experiments/exp1/01_exp1.R --config_index 1
Rscript experiments/exp1/02_aggregate.R
```

Settings 2 and 3 are run the same way using the scripts in `experiments/exp2/`
and `experiments/exp3/`.

After all three final CSVs exist, render
`experiments/results.qmd` from the repository root with:

```bash
quarto render experiments/results.qmd
```

## CASS Analysis

The `cass/` directory contains the analysis functions and the non-parametric
and parametric analyses used in the paper. CASS data can be requested from BioLINCC.

```bash
Rscript cass/02_cass_analysis_np.R ejec_frac
Rscript cass/03_cass_analysis_p.R ejec_frac
```

The final argument may also be `age` or `sex` for the appendix subgroup
analyses. 


## Citation

Rodriguez, V., Lehmann, B., & Diaz-Ordaz, K. (manuscript in preparation). *DR-Learners for Generalizing CATEs from a Randomized Controlled Trial.*

## Licence

This repository, including the `genhte` package and all analysis scripts, is licensed under the GNU General Public License v3.0 or later. See [LICENSE.md](LICENSE.md). The CASS data are not included and are subject to BioLINCC's data use terms.
