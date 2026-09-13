# Fits the δ-HDR candidate models and supports outputs from the manuscript analysis.

hdr_model_colors <- c(
  linear = "#2ca02c",
  logarithmic = "#7b3294",
  quadratic = "#d7191c"
)

hdr_model_linetypes <- c(
  linear = "dashed",
  logarithmic = "solid",
  quadratic = "dotdash"
)

hdr_factor_model_family <- function(x) {
  factor(as.character(x), levels = names(hdr_model_colors))
}

hdr_scale_color <- function(name = "Candidate model") {
  ggplot2::scale_color_manual(
    name = name,
    values = hdr_model_colors,
    limits = names(hdr_model_colors),
    drop = FALSE
  )
}

hdr_scale_fill <- function(name = "Candidate model") {
  ggplot2::scale_fill_manual(
    name = name,
    values = hdr_model_colors,
    limits = names(hdr_model_colors),
    drop = FALSE
  )
}

hdr_scale_linetype <- function(name = "Candidate model") {
  ggplot2::scale_linetype_manual(
    name = name,
    values = hdr_model_linetypes,
    limits = names(hdr_model_linetypes),
    drop = FALSE
  )
}


save_png_path <- function(plot, path, width = 8, height = 5, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  suppressWarnings(
    ggplot2::ggsave(path, plot, width = width, height = height, device = "png", dpi = dpi, bg = "white")
  )
}


# Information criterion helpers for model comparisons in the manuscript.
model_aicc <- function(model) {
  n <- stats::nobs(model)
  k <- attr(stats::logLik(model), "df")
  aic <- stats::AIC(model)
  if (!is.finite(n) || !is.finite(k) || (n - k - 1) <= 0) return(NA_real_)
  aic + (2 * k * (k + 1)) / (n - k - 1)
}

model_r2 <- function(model) {
  if (inherits(model, "lm")) {
    g <- broom::glance(model)
    return(tibble(r_squared = g$r.squared, adj_r_squared = g$adj.r.squared))
  }
  tibble(r_squared = NA_real_, adj_r_squared = NA_real_)
}

safe_term_coef <- function(model, term) {
  cf <- stats::coef(model)
  if (term %in% names(cf)) as.numeric(cf[term]) else NA_real_
}

# Information criteria rank candidate curve families. These helpers separately
# record coefficient direction, magnitude, P values, and 95% confidence intervals.
hdr_tidy_model_terms <- function(model, conf_level = 0.95) {
  out <- tryCatch(
    broom::tidy(model, conf.int = TRUE, conf.level = conf_level),
    error = function(e) broom::tidy(model) %>% mutate(conf.low = NA_real_, conf.high = NA_real_)
  )
  for (nm in c("std.error", "statistic", "p.value", "conf.low", "conf.high")) {
    if (!nm %in% names(out)) out[[nm]] <- NA_real_
  }
  out
}

hdr_term_stats <- function(tidy_tbl, term_name) {
  empty <- list(estimate = NA_real_, std_error = NA_real_, statistic = NA_real_,
                p_value = NA_real_, conf_low = NA_real_, conf_high = NA_real_)
  if (length(term_name) == 0L || is.na(term_name[1]) || !nzchar(as.character(term_name[1]))) return(empty)
  row <- tidy_tbl %>% filter(.data$term == .env$term_name[1]) %>% slice(1)
  if (nrow(row) == 0) return(empty)
  list(
    estimate = hdr_scalar_real(row$estimate),
    std_error = hdr_scalar_real(row$std.error),
    statistic = hdr_scalar_real(row$statistic),
    p_value = hdr_scalar_real(row$p.value),
    conf_low = hdr_scalar_real(row$conf.low),
    conf_high = hdr_scalar_real(row$conf.high)
  )
}

# Coerce one coefficient at a time for summaries and handle columns in term tables
# with explicitly vectorized operations.
hdr_scalar_real <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  x <- tryCatch(
    suppressWarnings(as.numeric(x)),
    error = function(e) NA_real_
  )
  if (length(x) == 0L) return(NA_real_)
  as.numeric(x[1])
}

hdr_real_vector <- function(x, n = NULL) {
  x <- tryCatch(
    suppressWarnings(as.numeric(x)),
    error = function(e) NA_real_
  )
  if (is.null(n)) return(x)
  if (n == 0L) return(numeric(0))
  if (length(x) == 0L) return(rep(NA_real_, n))
  if (length(x) == 1L) return(rep(as.numeric(x[1]), n))
  as.numeric(rep_len(x, n))
}

hdr_clear_positive <- function(estimate, p_value, conf_low, conf_high, alpha = 0.05) {
  estimate <- hdr_scalar_real(estimate)
  p_value <- hdr_scalar_real(p_value)
  conf_low <- hdr_scalar_real(conf_low)
  conf_high <- hdr_scalar_real(conf_high)
  alpha <- hdr_scalar_real(alpha)
  is.finite(estimate) && is.finite(p_value) && is.finite(conf_low) &&
    is.finite(conf_high) && is.finite(alpha) &&
    estimate > 0 && p_value < alpha && conf_low > 0 && conf_high > 0
}

hdr_clear_negative <- function(estimate, p_value, conf_low, conf_high, alpha = 0.05) {
  estimate <- hdr_scalar_real(estimate)
  p_value <- hdr_scalar_real(p_value)
  conf_low <- hdr_scalar_real(conf_low)
  conf_high <- hdr_scalar_real(conf_high)
  alpha <- hdr_scalar_real(alpha)
  is.finite(estimate) && is.finite(p_value) && is.finite(conf_low) &&
    is.finite(conf_high) && is.finite(alpha) &&
    estimate < 0 && p_value < alpha && conf_low < 0 && conf_high < 0
}

hdr_term_evidence_label <- function(stats, alpha = 0.05) {
  estimate <- hdr_scalar_real(stats$estimate)
  p_value <- hdr_scalar_real(stats$p_value)
  conf_low <- hdr_scalar_real(stats$conf_low)
  conf_high <- hdr_scalar_real(stats$conf_high)

  if (hdr_clear_positive(estimate, p_value, conf_low, conf_high, alpha)) {
    return("clear_positive_estimate_p_lt_0.05_CI_excludes_zero")
  }
  if (hdr_clear_negative(estimate, p_value, conf_low, conf_high, alpha)) {
    return("clear_negative_estimate_p_lt_0.05_CI_excludes_zero")
  }
  if (is.finite(estimate) && estimate > 0) return("positive_estimate_but_uncertain")
  if (is.finite(estimate) && estimate < 0) return("negative_estimate_but_uncertain")
  if (is.finite(estimate) && estimate == 0) return("zero_estimate")
  "not_estimable"
}

hdr_term_evidence_vector <- function(estimate, p_value, conf_low, conf_high, alpha = 0.05) {
  n <- max(length(estimate), length(p_value), length(conf_low), length(conf_high), 0L)
  if (n == 0L) return(character(0))
  estimate <- hdr_real_vector(estimate, n)
  p_value <- hdr_real_vector(p_value, n)
  conf_low <- hdr_real_vector(conf_low, n)
  conf_high <- hdr_real_vector(conf_high, n)
  purrr::map_chr(seq_len(n), function(i) {
    hdr_term_evidence_label(
      list(
        estimate = estimate[i],
        p_value = p_value[i],
        conf_low = conf_low[i],
        conf_high = conf_high[i]
      ),
      alpha = alpha
    )
  })
}

hdr_shape_support_label <- function(model_family, shape_interpretation, term1, term2,
                                     vertex_inside, high_side_declines, alpha = 0.05) {
  term1_pos <- hdr_clear_positive(term1$estimate, term1$p_value, term1$conf_low, term1$conf_high, alpha)
  term1_neg <- hdr_clear_negative(term1$estimate, term1$p_value, term1$conf_low, term1$conf_high, alpha)
  term2_neg <- hdr_clear_negative(term2$estimate, term2$p_value, term2$conf_low, term2$conf_high, alpha)
  if (model_family == "linear" && shape_interpretation == "positive_monotonic" && term1_pos) {
    return(list(label = "positive_monotonic_supported_by_coefficient_CI", supported = TRUE))
  }
  if (model_family == "linear" && shape_interpretation == "negative_monotonic" && term1_neg) {
    return(list(label = "negative_monotonic_supported_by_coefficient_CI", supported = TRUE))
  }
  if (model_family == "logarithmic" && shape_interpretation == "positive_diminishing_return" && term1_pos) {
    return(list(label = "positive_diminishing_return_supported_by_coefficient_CI", supported = TRUE))
  }
  if (model_family == "quadratic" && shape_interpretation == "unimodal_in_observed_range" &&
      term1_pos && term2_neg && isTRUE(vertex_inside) && isTRUE(high_side_declines)) {
    return(list(label = "unimodal_supported_by_coefficients_CI_and_geometry", supported = TRUE))
  }
  list(label = paste0("IC_selected_", model_family, "_shape_but_coefficient_or_geometry_uncertain"),
       supported = FALSE)
}

classify_candidate_shape <- function(model, family_name, x_min = -Inf, x_max = Inf) {
  b1 <- safe_term_coef(model, "x")
  b2 <- safe_term_coef(model, "I(x^2)")
  blog <- safe_term_coef(model, "logx")

  if (family_name == "linear") {
    if (is.na(b1)) return("linear_unknown")
    if (b1 > 0) return("positive_monotonic")
    if (b1 < 0) return("negative_monotonic")
    return("flat")
  }

  if (family_name == "logarithmic") {
    if (is.na(blog)) return("logarithmic_unknown")
    if (blog > 0) return("positive_diminishing_return")
    if (blog < 0) return("negative_diminishing_return")
    return("flat")
  }

  if (family_name == "quadratic") {
    if (is.na(b1) || is.na(b2)) return("quadratic_unknown")
    vertex <- if (b2 == 0) NA_real_ else -b1 / (2 * b2)
    if (b2 < 0 && b1 > 0 && is.finite(vertex) && vertex > x_min && vertex < x_max) return("unimodal_in_observed_range")
    if (b2 < 0) return("concave_down_peak_outside_range")
    if (b2 > 0) return("convex_or_U_shaped")
    return("linear_equivalent")
  }

  "unknown"
}


# 03 MANUSCRIPT HELPERS AND ANALYSIS HIERARCHY ####
# ============================================================ #
# These reusable helpers support the five manuscript analysis sets and their
# figure, table, and supporting outputs. Display labels remain separate from
# stable machine identifiers.
#
manuscript_delta_symbol <- function() "\u03b4"
manuscript_Delta_symbol <- function() "\u0394"
manuscript_squared_symbol <- function() "\u00B2"

# Dataset label convention for outputs in the manuscript.
# Internal analysis identifiers still use Topography/Foliage labels, while
# display labels in outputs use the source study names from the manuscript.
manuscript_dataset_label_table <- tibble::tribble(
  ~dataset_key, ~dataset_short, ~dataset_phrase, ~source_study, ~response_label, ~delta_axis_label, ~environmental_mean_label, ~primary_covariate_label,
  "topography", "Allouche et al. (2012)", "Allouche et al. (2012) dataset", "Allouche et al. (2012) Catalonia breeding-bird reanalysis", "Breeding-bird species richness", paste0("Topographic ", manuscript_delta_symbol()), "Mean elevation (m)", paste0("Mean elevation + mean elevation", manuscript_squared_symbol()),
  "foliage", "MacArthur & MacArthur (1961)", "MacArthur & MacArthur (1961) dataset", "MacArthur & MacArthur (1961) forest-bird reanalysis", "Bird species diversity", paste0("Foliage-height ", manuscript_delta_symbol()), "Mean foliage height", "Mean foliage height",
  "both", "Allouche et al. (2012) and MacArthur & MacArthur (1961)", "Allouche et al. (2012) and MacArthur & MacArthur (1961) datasets", "Allouche et al. (2012); MacArthur & MacArthur (1961)", "Biodiversity response", paste0(manuscript_delta_symbol(), " heterogeneity"), "Environmental mean", "Dataset-specific mean covariates",
  "unknown", "Unknown dataset", "unknown dataset", NA_character_, "Biodiversity response", paste0(manuscript_delta_symbol(), " heterogeneity"), "Environmental mean", "Dataset-specific covariates"
)

manuscript_dataset_key <- function(dataset) {
  dataset <- as.character(dataset)
  dplyr::case_when(
    is.na(dataset) ~ "unknown",
    grepl("Catalonia|Allouche|Topography|elevation|topograph", dataset, ignore.case = TRUE) ~ "topography",
    grepl("MacArthur|Foliage|foliage|forest", dataset, ignore.case = TRUE) ~ "foliage",
    grepl("Both", dataset, ignore.case = TRUE) ~ "both",
    TRUE ~ "unknown"
  )
}

manuscript_dataset_lookup <- function(dataset, column) {
  key <- manuscript_dataset_key(dataset)
  idx <- match(key, manuscript_dataset_label_table$dataset_key)
  unknown_idx <- match("unknown", manuscript_dataset_label_table$dataset_key)
  idx[is.na(idx)] <- unknown_idx
  as.character(manuscript_dataset_label_table[[column]][idx])
}

manuscript_dataset_short <- function(dataset) manuscript_dataset_lookup(dataset, "dataset_short")
manuscript_dataset_phrase <- function(dataset) manuscript_dataset_lookup(dataset, "dataset_phrase")
manuscript_dataset_source <- function(dataset) manuscript_dataset_lookup(dataset, "source_study")
manuscript_dataset_response <- function(dataset) manuscript_dataset_lookup(dataset, "response_label")
manuscript_dataset_delta_axis <- function(dataset) manuscript_dataset_lookup(dataset, "delta_axis_label")
manuscript_dataset_mean_axis <- function(dataset) manuscript_dataset_lookup(dataset, "environmental_mean_label")
manuscript_delta_output_set_label <- function(dataset) {
  paste0(manuscript_dataset_short(dataset), " ", manuscript_delta_symbol(), " covariate sensitivity")
}

format_term_label <- function(term, dataset = NULL) {
  values <- as.character(term)
  dataset_values <- if (is.null(dataset)) rep("unknown", length(values)) else rep_len(as.character(dataset), length(values))
  purrr::map2_chr(values, dataset_values, function(value, ds) {
    key <- manuscript_dataset_key(ds)
    if (is.na(value) || !nzchar(value)) return(NA_character_)
    if (value == "x") return(manuscript_delta_symbol())
    if (value == "logx") return(paste0("log(", manuscript_delta_symbol(), ")"))
    if (value == "I(x^2)") return(paste0(manuscript_delta_symbol(), manuscript_squared_symbol()))
    if (value == "mu" && key == "topography") return("mean elevation")
    if (value == "mu" && key == "foliage") return("mean foliage height")
    if (value == "mu") return("environmental mean")
    if (value == "mu2" && key == "topography") return(paste0("mean elevation", manuscript_squared_symbol()))
    if (value == "mu2") return(paste0("environmental mean", manuscript_squared_symbol()))
    value
  })
}

format_candidate_formula_label <- function(model_family, covariates, dataset = NULL) {
  model_family_values <- as.character(model_family)
  if (length(model_family_values) == 0L) return(character(0))

  # Single-model calls pass covariates as a character vector, e.g. c("mu", "mu2").
  # Table calls pass one covariate string per row, e.g. "mu;mu2". Keep both forms.
  single_model_call <- length(model_family_values) == 1L &&
    (is.null(dataset) || length(dataset) <= 1L) &&
    length(covariates) > 1L
  covariate_values <- if (length(covariates) == 0L) {
    "none"
  } else if (single_model_call) {
    paste(as.character(covariates), collapse = ";")
  } else {
    as.character(covariates)
  }

  n <- max(length(model_family_values), length(covariate_values), ifelse(is.null(dataset), 1L, length(dataset)))
  model_family_values <- rep_len(model_family_values, n)
  covariate_values <- rep_len(covariate_values, n)
  dataset_values <- if (is.null(dataset)) rep("unknown", n) else rep_len(as.character(dataset), n)

  purrr::map_chr(seq_len(n), function(i) {
    covariate_text <- format_covariate_label(covariate_values[i], dataset_values[i])
    heterogeneity_text <- dplyr::case_when(
      model_family_values[i] == "linear" ~ manuscript_delta_symbol(),
      model_family_values[i] == "logarithmic" ~ paste0("log(", manuscript_delta_symbol(), ")"),
      model_family_values[i] == "quadratic" ~ paste0(manuscript_delta_symbol(), " + ", manuscript_delta_symbol(), manuscript_squared_symbol()),
      TRUE ~ paste0("f(", manuscript_delta_symbol(), ")")
    )
    rhs <- if (covariate_text == "none") heterogeneity_text else paste(covariate_text, heterogeneity_text, sep = " + ")
    paste(manuscript_dataset_response(dataset_values[i]), "~", rhs)
  })
}

nature_manuscript_theme <- function(base_size = 8.4) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 0.8, margin = ggplot2::margin(0, 0, 2, 0)),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.7, margin = ggplot2::margin(0, 0, 4, 0)),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 1.1, colour = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 0.4),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = base_size - 0.7),
      legend.text = ggplot2::element_text(size = base_size - 1.1),
      legend.key.width = grid::unit(1.25, "lines"),
      legend.key.height = grid::unit(0.70, "lines"),
      plot.margin = ggplot2::margin(3, 5, 3, 3)
    )
}


manuscript_add_dataset_labels <- function(tbl, dataset_col = "dataset") {
  if (!dataset_col %in% names(tbl)) return(tbl)
  dataset_values <- tbl[[dataset_col]]
  tbl %>%
    mutate(
      dataset_short = manuscript_dataset_short(dataset_values),
      dataset_phrase = manuscript_dataset_phrase(dataset_values),
      source_study = manuscript_dataset_source(dataset_values),
      response_label = manuscript_dataset_response(dataset_values),
      delta_axis_label = manuscript_dataset_delta_axis(dataset_values),
      environmental_mean_label = manuscript_dataset_mean_axis(dataset_values)
    )
}

candidate_model_covariates <- function(candidate_set) {
  covariates <- candidate_set$summary$covariates[1]
  if (is.na(covariates) || covariates == "none" || covariates == "") return(character(0))
  strsplit(covariates, ";", fixed = TRUE)[[1]]
}

candidate_model_adjusted_points <- function(candidate_set) {
  dat <- candidate_set$data
  fam <- candidate_set$summary$family[1]
  covariates <- candidate_model_covariates(candidate_set)
  dataset_label <- if (grepl("Topography", candidate_set$label, ignore.case = TRUE)) "topography" else if (grepl("Foliage", candidate_set$label, ignore.case = TRUE)) "foliage" else candidate_set$label
  ok_covariates <- length(covariates) > 0 && all(covariates %in% names(dat))
  adjusted <- dat$y
  adjustment_label <- "unadjusted response"

  if (ok_covariates) {
    cov_formula <- stats::as.formula(paste("y ~", paste(covariates, collapse = " + ")))
    cov_fit <- stats::lm(cov_formula, data = dat)
    adjusted <- residuals(cov_fit, type = "response") + mean(stats::predict(cov_fit, type = "response"), na.rm = TRUE)
    adjustment_label <- paste("response adjusted for", format_covariate_label(paste(covariates, collapse = ";"), dataset_label))
  }

  tibble::tibble(
    analysis = candidate_set$label,
    x = dat$x,
    adjusted_response = as.numeric(adjusted),
    y_raw = dat$y,
    family = fam,
    covariates = ifelse(length(covariates) == 0, "none", paste(covariates, collapse = ";")),
    adjustment_label = adjustment_label,
    n = nrow(dat)
  )
}

candidate_model_prediction_weight_data <- function(candidate_set) {
  weight_cols <- candidate_set$summary %>%
    mutate(
      delta_for_plot = delta_IC,
      weight_for_plot = Akaike_weight,
      criterion_for_plot = selection_criterion
    ) %>%
    select(analysis, model_family, delta_for_plot, weight_for_plot, criterion_for_plot, shape_interpretation)
  pred <- candidate_set$predictions %>%
    left_join(weight_cols, by = c("analysis", "model_family")) %>%
    mutate(model_family = hdr_factor_model_family(model_family))
  weights <- weight_cols %>%
    mutate(model_family = hdr_factor_model_family(model_family))
  list(predictions = pred, weights = weights)
}

plot_candidate_model_set_with_weights <- function(candidate_set, show_points = TRUE) {
  pts <- candidate_model_adjusted_points(candidate_set)
  pw <- candidate_model_prediction_weight_data(candidate_set)
  pred <- pw$predictions %>% filter(is.finite(yhat), is.finite(x))
  weights <- pw$weights %>%
    mutate(model_family_short = factor(format_candidate_curve_abbrev(model_family), levels = c("lin", "log", "quad")))

  dataset_label <- if (grepl("Topography", candidate_set$label, ignore.case = TRUE)) {
    "topography"
  } else if (grepl("Foliage", candidate_set$label, ignore.case = TRUE)) {
    "foliage"
  } else {
    candidate_set$label
  }
  criterion <- unique(weights$criterion_for_plot)[1]
  response_axis <- manuscript_dataset_response(dataset_label)
  if (length(candidate_model_covariates(candidate_set)) > 0) response_axis <- paste0(response_axis, " (adjusted)")

  p_pred <- ggplot() +
    {if (show_points) geom_point(data = pts, aes(x, adjusted_response), alpha = 0.45, size = ifelse(nrow(pts) <= 20, 1.65, 0.8)) else NULL} +
    geom_line(data = pred, aes(x, yhat, color = model_family, linetype = model_family), linewidth = 0.86) +
    hdr_scale_color("Candidate curve") +
    hdr_scale_linetype("Candidate curve") +
    labs(
      x = manuscript_dataset_delta_axis(dataset_label),
      y = response_axis,
      title = paste0(manuscript_dataset_short(dataset_label), ": ", manuscript_delta_symbol(), " candidate curves")
    ) +
    nature_manuscript_theme(base_size = 8.2) +
    theme(
      plot.title = element_text(size = 9.6, face = "bold"),
      legend.position = "bottom"
    )

  p_weight <- ggplot(weights, aes(model_family_short, weight_for_plot, fill = model_family)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = sprintf("%.2f", weight_for_plot)), vjust = -0.25, size = 2.8) +
    hdr_scale_fill("Candidate curve") +
    ylim(0, max(1, weights$weight_for_plot, na.rm = TRUE) * 1.08) +
    labs(x = NULL, y = paste0(criterion, " weight")) +
    nature_manuscript_theme(base_size = 8.2) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(14, 4, 3, 4)
    )

  (p_pred | p_weight) + patchwork::plot_layout(widths = c(3.25, 0.95))
}



manuscript_csv_safe <- function(x) {
  # Keep figure labels typographically clean, but write machine-readable CSVs
  # with plain ASCII "delta" rather than the Greek symbol. This avoids mojibake
  # in spreadsheet programs while preserving the statistical meaning.
  if (!is.data.frame(x)) return(x)
  out <- x
  names(out) <- gsub(manuscript_delta_symbol(), "delta", names(out), fixed = TRUE)
  names(out) <- gsub(manuscript_Delta_symbol(), "Delta", names(out), fixed = TRUE)
  text_cols <- vapply(out, function(z) is.character(z) || is.factor(z), logical(1))
  out[text_cols] <- lapply(out[text_cols], function(z) {
    z <- as.character(z)
    z <- gsub(manuscript_delta_symbol(), "delta", z, fixed = TRUE)
    z <- gsub(manuscript_Delta_symbol(), "Delta", z, fixed = TRUE)
    z
  })
  out
}

write_analysis_hierarchy_csv <- function(x, dir_key, filename, write_outputs = TRUE) {
  if (isTRUE(write_outputs)) {
    out_path <- file.path(analysis_hierarchy_dirs[[dir_key]], filename)
    ensure_parent_dir(out_path)
    readr::write_csv(manuscript_csv_safe(x), out_path)
  }
  x
}


fit_delta_candidate_curve_models <- function(data, response, predictor, covariates = character(0),
                                      label, family = "gaussian",
                                      predictor_role = "heterogeneity",
                                      information_criterion = c("AIC", "AICc"),
                                      coefficient_conf_level = 0.95,
                                      coefficient_alpha = 0.05) {
  # Before fitting, we specified linear, quadratic, and logarithmic candidate
  # families to represent the three hypotheses. This function fits each candidate
  # with Gaussian ordinary least squares; the Catalonia comparisons are nonspatial.
  # It ranks MacArthur models by AICc (n = 13) and Catalonia models by AIC.
  # Coefficient estimates, P values, and confidence intervals assess sign support
  # for each fitted response shape. Quadratic support also requires a vertex within
  # the observed predictor range and fitted decline toward its upper end.
  family <- "gaussian"
  information_criterion <- match.arg(information_criterion)
  required <- c(response, predictor, covariates)
  if (!all(required %in% names(data))) {
    missing <- setdiff(required, names(data))
    stop("Missing columns for ", label, ": ", paste(missing, collapse = ", "))
  }

  tmp <- data %>%
    transmute(
      y = as.numeric(.data[[response]]),
      x = as.numeric(.data[[predictor]]),
      !!!rlang::syms(covariates)
    ) %>%
    filter(is.finite(y), is.finite(x))

  if (nrow(tmp) < 5) {
    warning("Too few rows for candidate fit: ", label)
    return(NULL)
  }

  if (any(tmp$x <= 0, na.rm = TRUE)) {
    bad_n <- sum(tmp$x <= 0, na.rm = TRUE)
    stop("Logarithmic candidate requires strictly positive delta values when fitting log(delta). ",
         "Analysis ", label, " has ", bad_n, " non-positive predictor value(s).")
  }

  x_min <- min(tmp$x, na.rm = TRUE)
  x_max <- max(tmp$x, na.rm = TRUE)
  if (!is.finite(x_min) || !is.finite(x_max) || x_min == x_max) {
    warning("Predictor has no finite range for candidate fit: ", label)
    return(NULL)
  }

  tmp <- tmp %>%
    mutate(logx = log(x)) %>%
    filter(is.finite(x), is.finite(logx))

  cov_string <- if (length(covariates) == 0) "" else paste(covariates, collapse = " + ")
  plus_cov <- if (length(covariates) == 0) "" else paste0(cov_string, " + ")
  formulas <- list(
    linear = stats::as.formula(paste0("y ~ ", plus_cov, "x")),
    logarithmic = stats::as.formula(paste0("y ~ ", plus_cov, "logx")),
    quadratic = stats::as.formula(paste0("y ~ ", plus_cov, "x + I(x^2)"))
  )

  models <- purrr::map(formulas, ~ stats::lm(.x, data = tmp))

  summary_tbl <- purrr::imap_dfr(models, function(model, model_family) {
    aic <- stats::AIC(model)
    aicc <- model_aicc(model)
    ic_value <- if (information_criterion == "AIC") aic else aicc
    r2s <- model_r2(model)
    b1 <- safe_term_coef(model, "x")
    b2 <- safe_term_coef(model, "I(x^2)")
    vertex <- if (is.finite(b1) && is.finite(b2) && b2 != 0) -b1 / (2 * b2) else NA_real_
    vertex_inside <- is.finite(vertex) && vertex > x_min && vertex < x_max
    high_side_declines <- is.finite(b1) && is.finite(b2) && is.finite(vertex) &&
      x_max > vertex && ((b1 * x_max + b2 * x_max^2) < (b1 * vertex + b2 * vertex^2))
    shape_interpretation <- classify_candidate_shape(model, model_family, x_min = x_min, x_max = x_max)

    tidy_ci <- hdr_tidy_model_terms(model, conf_level = coefficient_conf_level)
    term_1 <- dplyr::case_when(
      model_family == "linear" ~ "x",
      model_family == "logarithmic" ~ "logx",
      model_family == "quadratic" ~ "x",
      TRUE ~ NA_character_
    )
    term_2 <- ifelse(model_family == "quadratic", "I(x^2)", NA_character_)
    term1 <- hdr_term_stats(tidy_ci, term_1)
    term2 <- hdr_term_stats(tidy_ci, term_2)
    term1_clear_pos <- hdr_clear_positive(term1$estimate, term1$p_value, term1$conf_low, term1$conf_high, coefficient_alpha)
    term1_clear_neg <- hdr_clear_negative(term1$estimate, term1$p_value, term1$conf_low, term1$conf_high, coefficient_alpha)
    term2_clear_pos <- hdr_clear_positive(term2$estimate, term2$p_value, term2$conf_low, term2$conf_high, coefficient_alpha)
    term2_clear_neg <- hdr_clear_negative(term2$estimate, term2$p_value, term2$conf_low, term2$conf_high, coefficient_alpha)
    effect_evidence <- if (model_family == "quadratic") {
      if (term1_clear_pos && term2_clear_neg && vertex_inside && high_side_declines) "clear_unimodal_terms_and_geometry"
      else if (is.finite(term1$estimate) && is.finite(term2$estimate) && term1$estimate > 0 && term2$estimate < 0) "unimodal_estimates_but_uncertain"
      else "quadratic_terms_uncertain_or_not_unimodal"
    } else {
      hdr_term_evidence_label(term1, alpha = coefficient_alpha)
    }
    shape_support <- hdr_shape_support_label(
      model_family, shape_interpretation, term1, term2,
      vertex_inside = vertex_inside,
      high_side_declines = high_side_declines,
      alpha = coefficient_alpha
    )

    tibble(
      analysis = label,
      predictor = predictor,
      predictor_role = predictor_role,
      model_family = model_family,
      response = response,
      family = family,
      covariates = ifelse(length(covariates) == 0, "none", paste(covariates, collapse = ";")),
      formula = deparse(formulas[[model_family]]),
      formula_readable = format_candidate_formula_label(model_family, covariates, dataset = label),
      n = stats::nobs(model),
      k = attr(stats::logLik(model), "df"),
      AIC = aic,
      AICc = aicc,
      selection_criterion = information_criterion,
      IC = ic_value,
      coefficient_conf_level = coefficient_conf_level,
      coefficient_alpha = coefficient_alpha,
      predictor_min = x_min,
      predictor_max = x_max,
      slope_linear_x = b1,
      slope_logx = safe_term_coef(model, "logx"),
      quadratic_x_2 = b2,
      quadratic_vertex_x = vertex,
      quadratic_vertex_inside_observed_range = vertex_inside,
      quadratic_high_side_declines = high_side_declines,
      heterogeneity_term_1 = term_1,
      heterogeneity_term_1_estimate = term1$estimate,
      heterogeneity_term_1_std_error = term1$std_error,
      heterogeneity_term_1_statistic = term1$statistic,
      heterogeneity_term_1_p_value = term1$p_value,
      heterogeneity_term_1_conf_low = term1$conf_low,
      heterogeneity_term_1_conf_high = term1$conf_high,
      heterogeneity_term_1_clear_positive = term1_clear_pos,
      heterogeneity_term_1_clear_negative = term1_clear_neg,
      heterogeneity_term_2 = term_2,
      heterogeneity_term_2_estimate = term2$estimate,
      heterogeneity_term_2_std_error = term2$std_error,
      heterogeneity_term_2_statistic = term2$statistic,
      heterogeneity_term_2_p_value = term2$p_value,
      heterogeneity_term_2_conf_low = term2$conf_low,
      heterogeneity_term_2_conf_high = term2$conf_high,
      heterogeneity_term_2_clear_positive = term2_clear_pos,
      heterogeneity_term_2_clear_negative = term2_clear_neg,
      heterogeneity_effect_evidence = effect_evidence,
      shape_support_with_uncertainty = shape_support$label,
      shape_supported_by_coefficients = shape_support$supported,
      shape_support_notes = ifelse(
        shape_support$supported,
        "heterogeneity coefficient evidence supports the selected curve-shape interpretation",
        "heterogeneity coefficient evidence is weak, uncertain, or inconsistent with a strong directional interpretation"
      ),
      shape_interpretation = shape_interpretation
    ) %>% bind_cols(r2s)
  }) %>%
    mutate(
      delta_IC = IC - min(IC, na.rm = TRUE),
      Akaike_weight = exp(-0.5 * delta_IC) / sum(exp(-0.5 * delta_IC), na.rm = TRUE),
      best_by_IC = delta_IC == min(delta_IC, na.rm = TRUE)
    ) %>%
    arrange(IC)

  term_tbl <- purrr::imap_dfr(models, function(model, model_family) {
    hdr_tidy_model_terms(model, conf_level = coefficient_conf_level) %>%
      mutate(
        analysis = label,
        predictor = predictor,
        model_family = model_family,
        family = family,
        information_criterion = information_criterion,
          coefficient_conf_level = coefficient_conf_level,
        coefficient_alpha = coefficient_alpha,
        term_pretty = format_term_label(term, label),
        term_role = dplyr::case_when(
          term == "x" & model_family == "linear" ~ "heterogeneity_linear_term",
          term == "logx" & model_family == "logarithmic" ~ "heterogeneity_logarithmic_term",
          term == "x" & model_family == "quadratic" ~ "heterogeneity_quadratic_linear_term",
          term == "I(x^2)" & model_family == "quadratic" ~ "heterogeneity_quadratic_squared_term",
          term == "(Intercept)" ~ "intercept",
          TRUE ~ "covariate_or_other_term"
        ),
        term_evidence = hdr_term_evidence_vector(
          .data$estimate,
          .data$p.value,
          .data$conf.low,
          .data$conf.high,
          alpha = coefficient_alpha
        ),
        .before = 1
      )
  })

  # Hold each covariate at its sample mean when constructing sensitivity
  # predictions. For mu and mu2, use mean(mu) and mean(mu2) separately.
  pred_grid <- tibble(x = seq(x_min, x_max, length.out = 200), logx = log(x))
  for (cv in covariates) {
    pred_grid[[cv]] <- mean(tmp[[cv]], na.rm = TRUE)
  }

  pred_tbl <- purrr::imap_dfr(models, function(model, model_family) {
    pred <- stats::predict(model, newdata = pred_grid)
    pred_grid %>%
      select(x) %>%
      mutate(yhat = as.numeric(pred), model_family = model_family)
  }) %>%
    mutate(analysis = label, predictor = predictor, family = family,
           selection_criterion = information_criterion)

  list(label = label, data = tmp, models = models, summary = summary_tbl,
       terms = term_tbl, predictions = pred_tbl,
       information_criterion = information_criterion,
       coefficient_conf_level = coefficient_conf_level,
       coefficient_alpha = coefficient_alpha)
}

combine_candidate_curve_model_sets <- function(candidate_sets, dir_key, file_prefix,
                                               write_outputs = TRUE) {
  keep <- candidate_sets[!vapply(candidate_sets, is.null, logical(1))]
  summary_tbl <- purrr::map_dfr(keep, "summary")
  terms_tbl <- purrr::map_dfr(keep, "terms")
  predictions_tbl <- purrr::map_dfr(keep, "predictions")
  write_analysis_hierarchy_csv(
    summary_tbl, dir_key, paste0(file_prefix, "_model_comparison.csv"),
    write_outputs = write_outputs
  )
  write_analysis_hierarchy_csv(
    terms_tbl, dir_key, paste0(file_prefix, "_terms.csv"),
    write_outputs = write_outputs
  )
  write_analysis_hierarchy_csv(
    predictions_tbl, dir_key, paste0(file_prefix, "_predictions.csv"),
    write_outputs = write_outputs
  )
  list(sets = keep, summary = summary_tbl, terms = terms_tbl, predictions = predictions_tbl)
}

figure1_select_plotted_family <- function(summary_tbl) {
  required <- c("model_family", "delta_IC")
  if (!all(required %in% names(summary_tbl))) {
    stop("Figure 1 model selection requires model_family and full-precision delta_IC.")
  }

  tie_order <- c("linear", "logarithmic", "quadratic")
  candidates <- summary_tbl %>%
    transmute(
      model_family = as.character(model_family),
      delta_IC = as.numeric(delta_IC),
      tie_index = match(model_family, .env$tie_order)
    ) %>%
    filter(model_family %in% .env$tie_order, is.finite(delta_IC))

  if (nrow(candidates) == 0L) {
    stop("Figure 1 model selection found no candidate with a finite delta_IC value.")
  }

  minimum_delta_ic <- min(candidates$delta_IC)
  candidates %>%
    filter(delta_IC == minimum_delta_ic) %>%
    arrange(tie_index) %>%
    pull(model_family) %>%
    dplyr::first()
}

summarize_best_curve_by_analysis <- function(summary_tbl) {
  summary_tbl %>%
    group_by(analysis) %>%
    arrange(delta_IC, .by_group = TRUE) %>%
    summarise(
      predictor = predictor[1],
      predictor_role = predictor_role[1],
      response = response[1],
      covariates = covariates[1],
      selection_criterion = selection_criterion[1],
      n = n[1],
      best_model_family = model_family[1],
      best_shape_interpretation = shape_interpretation[1],
      best_IC = IC[1],
      best_AIC = AIC[1],
      best_AICc = AICc[1],
      best_Akaike_weight = Akaike_weight[1],
      best_heterogeneity_term_1 = heterogeneity_term_1[1],
      best_heterogeneity_term_1_estimate = heterogeneity_term_1_estimate[1],
      best_heterogeneity_term_1_std_error = heterogeneity_term_1_std_error[1],
      best_heterogeneity_term_1_statistic = heterogeneity_term_1_statistic[1],
      best_heterogeneity_term_1_p_value = heterogeneity_term_1_p_value[1],
      best_heterogeneity_term_1_conf_low = heterogeneity_term_1_conf_low[1],
      best_heterogeneity_term_1_conf_high = heterogeneity_term_1_conf_high[1],
      best_heterogeneity_term_1_clear_positive = heterogeneity_term_1_clear_positive[1],
      best_heterogeneity_term_1_clear_negative = heterogeneity_term_1_clear_negative[1],
      best_heterogeneity_term_2 = heterogeneity_term_2[1],
      best_heterogeneity_term_2_estimate = heterogeneity_term_2_estimate[1],
      best_heterogeneity_term_2_std_error = heterogeneity_term_2_std_error[1],
      best_heterogeneity_term_2_statistic = heterogeneity_term_2_statistic[1],
      best_heterogeneity_term_2_p_value = heterogeneity_term_2_p_value[1],
      best_heterogeneity_term_2_conf_low = heterogeneity_term_2_conf_low[1],
      best_heterogeneity_term_2_conf_high = heterogeneity_term_2_conf_high[1],
      best_heterogeneity_term_2_clear_positive = heterogeneity_term_2_clear_positive[1],
      best_heterogeneity_term_2_clear_negative = heterogeneity_term_2_clear_negative[1],
      best_quadratic_vertex_x = quadratic_vertex_x[1],
      best_quadratic_vertex_inside_observed_range = quadratic_vertex_inside_observed_range[1],
      best_quadratic_high_side_declines = quadratic_high_side_declines[1],
      best_heterogeneity_effect_evidence = heterogeneity_effect_evidence[1],
      best_shape_support_with_uncertainty = shape_support_with_uncertainty[1],
      best_shape_supported_by_coefficients = shape_supported_by_coefficients[1],
      best_shape_support_notes = shape_support_notes[1],
      second_best_delta_IC = {
        d <- sort(unique(delta_IC[is.finite(delta_IC)]))
        if (length(d) < 2) NA_real_ else d[2]
      },
      plausible_models_delta_IC_lt_2 = paste(model_family[delta_IC < 2], collapse = ";"),
      .groups = "drop"
    ) %>%
    mutate(
      evidence_strength = dplyr::case_when(
        !is.finite(second_best_delta_IC) ~ "not estimable",
        second_best_delta_IC < 2 ~ paste0("ambiguous: multiple models have ", manuscript_Delta_symbol(), "IC < 2"),
        second_best_delta_IC < 6 ~ "moderate support for best model",
        second_best_delta_IC < 10 ~ "strong support for best model",
        TRUE ~ "very strong support for best model"
      ),
      recommended_interpretation = dplyr::case_when(
        best_model_family == "logarithmic" & best_shape_interpretation == "positive_diminishing_return" ~ "positive diminishing-return HDR",
        best_model_family == "linear" & best_shape_interpretation == "positive_monotonic" ~ "positive monotonic HDR",
        best_model_family == "linear" & best_shape_interpretation == "negative_monotonic" ~ "negative monotonic relationship",
        best_model_family == "quadratic" & best_shape_interpretation == "unimodal_in_observed_range" ~ "unimodal HDR candidate",
        TRUE ~ "ambiguous or unsupported shape"
      ),
      coefficient_qualified_interpretation = dplyr::case_when(
        best_model_family == "logarithmic" & best_shape_supported_by_coefficients == TRUE ~ "positive diminishing-return HDR supported by coefficient sign, p-value, and 95% CI",
        best_model_family == "logarithmic" ~ "IC-selected logarithmic candidate; positive diminishing-return interpretation requires coefficient check",
        best_model_family == "linear" & best_shape_supported_by_coefficients == TRUE & best_shape_interpretation == "positive_monotonic" ~ "positive monotonic HDR supported by coefficient sign, p-value, and 95% CI",
        best_model_family == "linear" & best_shape_supported_by_coefficients == TRUE & best_shape_interpretation == "negative_monotonic" ~ "negative monotonic relationship supported by coefficient sign, p-value, and 95% CI",
        best_model_family == "linear" ~ "IC-selected linear candidate; monotonic interpretation requires coefficient check",
        best_model_family == "quadratic" & best_shape_supported_by_coefficients == TRUE ~ "unimodal HDR supported by coefficient signs, p-values, 95% CIs, and geometry",
        best_model_family == "quadratic" ~ "IC-selected quadratic candidate; unimodal interpretation requires coefficient and geometry checks",
        TRUE ~ "coefficient evidence does not support a directional HDR interpretation"
      )
    )
}

