library(dplyr)
library(ggplot2)
library(tidyr)

theme_set(
  theme_classic(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "gray95", color = "black"),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "bottom",
      legend.key.width = unit(1.5, "cm")
    )
)

PAPER_FONT_SIZE <- 10 / 0.65
paper_text <- theme(
  text = element_text(family = "Times", size = PAPER_FONT_SIZE),
  strip.text = element_text(face = "plain", size = PAPER_FONT_SIZE),
  strip.background = element_blank()
)

y_axis_title <- function(label) {
  patchwork::wrap_elements(
    full = grid::textGrob(
      label,
      rot = 90,
      gp = grid::gpar(fontsize = PAPER_FONT_SIZE, fontfamily = "Times")
    )
  )
}

# Re-facet a one-row panel plot as one row of a stacked grid, with a right-hand
# strip (plotmath row_lab) naming what the row holds fixed
as_grid_row <- function(p, col_vars, col_labels, row_lab, x_lab) {
  p$data$row <- row_lab
  p +
    ggh4x::facet_nested(
      rows = vars(row),
      cols = vars(!!!syms(col_vars)),
      scales = "free_y",
      independent = "y",
      axes = "x",
      remove_labels = "x",
      nest_line = element_line(),
      labeller = labeller(.rows = label_parsed, .cols = as_labeller(col_labels))
    ) +
    labs(x = x_lab, y = NULL, colour = NULL, linetype = NULL, shape = NULL) +
    theme(aspect.ratio = 0.88, legend.title = element_blank())
}

method_colors <- c(
  "T-Learner" = "#E69F00",
  "gT-Learner" = "#D55E00",
  "DR-Learner" = "#56B4E9",
  "gDR-Learner" = "#0072B2",
  "gT (Enr.)" = "#D55E00",
  "gDR (Enr.)" = "#0072B2",
  "gT (Subsample)" = "#D55E00",
  "gDR (Subsample)" = "#0072B2",
  "gDR (Restricted)" = "#0072B2",
  "gT (Restricted)" = "#D55E00"
)

method_shapes <- c(
  "T-Learner" = 16,
  "gT-Learner" = 17,
  "DR-Learner" = 15,
  "gDR-Learner" = 18,
  "gT (Enr.)" = 17,
  "gDR (Enr.)" = 18,
  "gT (Subsample)" = 17,
  "gDR (Subsample)" = 18,
  "gT (Restricted)" = 17,
  "gDR (Restricted)" = 18
)

SCENARIO_LABELS <- c(
  LD_linear = "Low-Dim: Linear CATE",
  LD_nonlinear = "Low-Dim: Non-Linear CATE",
  HD_linear = "High-Dim: Linear CATE",
  HD_nonlinear = "High-Dim: Non-Linear CATE"
)
scenario_label <- function(scenario) {
  droplevels(factor(unname(SCENARIO_LABELS[scenario]), levels = SCENARIO_LABELS))
}

# Mean and standard error of each metric over replications, per group
aggregate_results <- function(
  df,
  keys,
  metrics = c(
    "mse",
    "bias",
    "abs_bias",
    "rmse",
    "median_mse",
    "median_bias",
    "median_abs_bias"
  )
) {
  df |>
    group_by(across(all_of(keys))) |>
    summarize(
      n = n(),
      across(
        all_of(metrics),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          se = ~ sd(.x, na.rm = TRUE) / sqrt(n())
        )
      ),
      .groups = "drop"
    )
}

# ------------------------------------------------------------------------------
# Setting 1
# ------------------------------------------------------------------------------

aggregate_exp1 <- function(df) {
  aggregate_results(
    df,
    c("scenario", "ssize", "delta", "method"),
    c(
      "mse",
      "bias",
      "abs_bias",
      "rmse",
      "median_mse",
      "median_bias",
      "median_abs_bias",
      "coverage_split_95",
      "avg_width_split_95",
      "coverage_inexact_95",
      "avg_width_inexact_95",
      "coverage_split_95_ite",
      "coverage_inexact_95_ite"
    )
  ) |>
    mutate(baseline = delta == ifelse(grepl("^HD", scenario), 5, 1))
}

plot_exp1_interval_width <- function(agg) {
  plot_data <- agg |>
    filter(baseline) |>
    select(scenario, ssize, method, avg_width_split_95_mean, avg_width_inexact_95_mean) |>
    pivot_longer(
      c(avg_width_split_95_mean, avg_width_inexact_95_mean),
      names_to = "conformal_label",
      values_to = "avg_width"
    ) |>
    mutate(
      conformal_label = ifelse(grepl("split", conformal_label), "Naive", "Inexact")
    )

  ggplot(
    plot_data,
    aes(
      x = ssize,
      y = avg_width,
      color = method,
      linetype = conformal_label,
      group = interaction(method, conformal_label)
    )
  ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.2) +
    scale_color_manual(values = method_colors) +
    scale_linetype_manual(values = c("Naive" = "solid", "Inexact" = "dashed"))
}

# ite = TRUE scores the intervals against the ITE instead of the CATE
make_coverage_table <- function(agg, output_path = NULL, ite = FALSE) {
  sfx <- if (ite) "_ite" else ""
  tab <- agg |>
    filter(baseline) |>
    mutate(
      scenario_label = scenario_label(scenario),
      split = sprintf("%.3f", .data[[paste0("coverage_split_95", sfx, "_mean")]]),
      inexact = sprintf("%.3f", .data[[paste0("coverage_inexact_95", sfx, "_mean")]])
    ) |>
    arrange(scenario_label, method, ssize)

  ssizes <- sort(unique(tab$ssize))

  split_wide <- tab |>
    select(scenario_label, method, ssize, split) |>
    pivot_wider(names_from = ssize, values_from = split, names_prefix = "split_")
  inexact_wide <- tab |>
    select(scenario_label, method, ssize, inexact) |>
    pivot_wider(names_from = ssize, values_from = inexact, names_prefix = "inexact_")
  wide <- left_join(split_wide, inexact_wide, by = c("scenario_label", "method")) |>
    arrange(scenario_label, method)

  n_sizes <- length(ssizes)
  col_spec <- paste0("ll", strrep("c", 2 * n_sizes))

  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    sprintf(
      "\\caption{Conformal inference coverage%s (nominal 95\\%%) across scenarios and sample sizes with $\\Delta=1$ (low-dimensional), $\\Delta=5$ (high-dimensional).}",
      if (ite) " of the ITE" else ""
    ),
    sprintf("\\label{tab:conformal_coverage%s}", sfx),
    "\\resizebox{\\textwidth}{!}{%",
    sprintf("\\begin{tabular}{%s}", col_spec),
    "\\toprule",
    sprintf(
      " & & \\multicolumn{%d}{c}{Naive Conformal} & \\multicolumn{%d}{c}{Inexact Conformal} \\\\",
      n_sizes,
      n_sizes
    ),
    sprintf(
      "\\cmidrule(lr){3-%d} \\cmidrule(lr){%d-%d}",
      2 + n_sizes,
      3 + n_sizes,
      2 + 2 * n_sizes
    ),
    paste0(
      "Scenario & Method & ",
      paste(c(ssizes, ssizes), collapse = " & "),
      " \\\\"
    ),
    "\\midrule"
  )

  current_scenario <- ""
  for (i in seq_len(nrow(wide))) {
    row <- wide[i, ]
    scenario_str <- ""
    if (as.character(row$scenario_label) != current_scenario) {
      if (current_scenario != "") lines <- c(lines, "\\midrule")
      current_scenario <- as.character(row$scenario_label)
      scenario_str <- current_scenario
    }
    lines <- c(
      lines,
      sprintf(
        "%s & %s & %s & %s \\\\",
        scenario_str,
        row$method,
        paste(row[, paste0("split_", ssizes)], collapse = " & "),
        paste(row[, paste0("inexact_", ssizes)], collapse = " & ")
      )
    )
  }

  lines <- c(lines, "\\bottomrule", "\\end{tabular}}", "\\end{table}")
  latex_code <- paste(lines, collapse = "\n")

  if (!is.null(output_path)) writeLines(latex_code, output_path)
  cat(latex_code)
  invisible(latex_code)
}

# ------------------------------------------------------------------------------
# Setting 2
# ------------------------------------------------------------------------------

aggregate_exp2 <- function(df) {
  aggregate_results(df, c("dimension", "complexity", "nmu", "ssize", "model", "method"))
}

prepare_exp2 <- function(df) {
  df |>
    mutate(
      panel_label = factor(
        paste0(
          ifelse(dimension == "low", "Low-Dim: ", "High-Dim: "),
          ifelse(complexity == "simple", "Linear Outcome", "Non-Linear Outcome")
        ),
        levels = c(
          "Low-Dim: Linear Outcome",
          "Low-Dim: Non-Linear Outcome",
          "High-Dim: Linear Outcome",
          "High-Dim: Non-Linear Outcome"
        )
      ),
      linetype_var = case_when(
        model %in% c("tl", "dr") ~ "dashed",
        model %in% c("tl_ipw", "dr_ipw") ~ "dotted",
        model %in% c("tl_ipw_z", "dr_ipw_z") ~ "solid"
      ),
      method = factor(
        method,
        levels = c(
          "T-Learner",
          "DR-Learner",
          "gT-Learner",
          "gDR-Learner",
          "gT (Enr.)",
          "gDR (Enr.)"
        )
      )
    )
}

# Invisible points at each dimension's min and max, so the two outcome panels
# within a dimension share a y range under scales = "free"
exp2_shared_y_data <- function(plot_data, x_var, y) {
  rng <- plot_data |>
    group_by(dimension) |>
    summarise(
      .lo = min(.data[[y]], na.rm = TRUE),
      .hi = max(.data[[y]], na.rm = TRUE),
      .groups = "drop"
    )
  panels <- plot_data |>
    distinct(dimension, panel_label) |>
    left_join(rng, by = "dimension")
  bind_rows(
    transmute(panels, panel_label, .y = .lo),
    transmute(panels, panel_label, .y = .hi)
  ) |>
    mutate(.x = min(plot_data[[x_var]], na.rm = TRUE))
}

plot_exp2_panels <- function(plot_data, x_var, y) {
  linetypes <- plot_data |>
    distinct(method, linetype_var) |>
    tibble::deframe()

  ggplot(
    plot_data,
    aes(
      x = .data[[x_var]],
      y = .data[[y]],
      color = method,
      linetype = method,
      shape = method
    )
  ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.2) +
    geom_blank(
      data = exp2_shared_y_data(plot_data, x_var, y),
      aes(x = .x, y = .y),
      inherit.aes = FALSE
    ) +
    facet_wrap(~panel_label, scales = "free", nrow = 1) +
    scale_color_manual(values = method_colors) +
    scale_linetype_manual(values = linetypes) +
    scale_shape_manual(values = method_shapes)
}

# x is the number of additional prognostic covariates d_P = nmu - nsel, at n = 1500
plot_exp2_prognostic <- function(agg, y) {
  agg |>
    filter(ssize == 1500) |>
    mutate(d_p = nmu - ifelse(dimension == "low", 2, 10)) |>
    prepare_exp2() |>
    plot_exp2_panels("d_p", y) +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10))
}

plot_exp2_trial_size <- function(agg, y) {
  agg |>
    filter(nmu %in% c(6, 14)) |>
    prepare_exp2() |>
    plot_exp2_panels("ssize", y)
}

# ------------------------------------------------------------------------------
# Setting 3
# ------------------------------------------------------------------------------

aggregate_exp3 <- function(df) {
  aggregate_results(df, c("scenario", "ssize", "nsubsample", "model", "method"))
}

prepare_exp3 <- function(df) {
  df |>
    mutate(
      scenario_label = scenario_label(scenario),
      linetype_var = case_when(
        model %in% c("tl", "dr") ~ "dashed",
        model %in% c("tl_sub", "dr_sub") ~ "dotted",
        model %in% c("tl_sub_marg", "dr_sub_marg") ~ "solid"
      ),
      method = factor(
        method,
        levels = c(
          "T-Learner",
          "DR-Learner",
          "gT (Subsample)",
          "gDR (Subsample)",
          "gT (Restricted)",
          "gDR (Restricted)"
        )
      )
    )
}

plot_exp3_panels <- function(plot_data, x_var, y) {
  linetypes <- plot_data |>
    distinct(method, linetype_var) |>
    tibble::deframe()

  ggplot(
    plot_data,
    aes(
      x = .data[[x_var]],
      y = .data[[y]],
      color = method,
      linetype = method,
      shape = method
    )
  ) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.2) +
    scale_color_manual(values = method_colors) +
    scale_linetype_manual(values = linetypes) +
    scale_shape_manual(values = method_shapes)
}

plot_exp3_subsample_size <- function(agg, y) {
  agg |>
    filter(ssize == 2000) |>
    prepare_exp3() |>
    plot_exp3_panels("nsubsample", y)
}

plot_exp3_trial_size <- function(agg, y) {
  agg |>
    filter(nsubsample == 5000) |>
    prepare_exp3() |>
    plot_exp3_panels("ssize", y)
}
