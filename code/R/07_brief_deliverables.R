# Builds the figures, tables, captions, and manifests for the manuscript analysis.

# 04 MANUSCRIPT DELIVERABLE FIGURES AND TABLES ####
# ============================================================ #
# These functions reuse the fitted candidate models rather than refitting them.
# When adjustment covariates are present, get_adjusted_points() fits a separate
# covariate-only model.
#
# Output location, relative to the project root:
#   results/delta_hdr_reanalysis/manuscript_tables_figures/
#     main_text/
#     extended_data/
#     supporting_analysis/
#     captions/
#     manifests/
#
# Main text item created:
#   Figure 1: positive HDRs with diminishing returns in landmark datasets.
#
# Extended Data items created:
#   Extended Data Table 1: MacArthur & MacArthur model comparisons and coefficient checks.
#   Extended Data Table 2: Allouche et al. model comparisons and coefficient checks.
#   Extended Data Fig. 1 and Extended Data Fig. 2 show model support across all analyses of the
#     two landmark datasets.

manuscript_deliverable_dirs <- list(
  main_text = file.path(manuscript_outputs_root, "main_text"),
  extended_data = file.path(manuscript_outputs_root, "extended_data"),
  supporting_analysis = file.path(manuscript_outputs_root, "supporting_analysis"),
  captions = file.path(manuscript_outputs_root, "captions"),
  manifests = file.path(manuscript_outputs_root, "manifests")
)
workflow_setup_dirs <- c(workflow_setup_dirs, unlist(manuscript_deliverable_dirs, use.names = FALSE))
brief_output_dirs <- c(brief_output_dirs, unlist(manuscript_deliverable_dirs, use.names = FALSE))

manuscript_write_csv <- function(x, dir_key, filename) {
  out_path <- file.path(manuscript_deliverable_dirs[[dir_key]], filename)
  ensure_parent_dir(out_path)
  readr::write_csv(manuscript_csv_safe(x), out_path)
  x
}

manuscript_write_display_csv <- function(x, dir_key, filename) {
  out_path <- file.path(manuscript_deliverable_dirs[[dir_key]], filename)
  ensure_parent_dir(out_path)
  readr::write_csv(x, out_path)
  x
}

manuscript_save_png <- function(plot, dir_key, filename, width = 9, height = 6, dpi = 300) {
  out_path <- file.path(manuscript_deliverable_dirs[[dir_key]], sub("\\.pdf$", ".png", filename))
  save_png_path(plot, out_path, width = width, height = height, dpi = dpi)
  out_path
}

manuscript_copy_file <- function(source_path, dir_key, filename) {
  out_path <- file.path(manuscript_deliverable_dirs[[dir_key]], filename)
  ensure_parent_dir(out_path)
  ok <- file.copy(source_path, out_path, overwrite = TRUE)
  if (!isTRUE(ok)) {
    stop("Failed to copy manuscript asset: ", source_path, " -> ", out_path, call. = FALSE)
  }
  out_path
}

manuscript_fmt_num <- function(x, digits = 3) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "fg", flag = "#"), NA_character_)
}

# Formatter with fixed decimals for manuscript Word tables. This is deliberately
# separate from manuscript_fmt_num(), which is used elsewhere in the workflow
# for compact labels with significant digits in plots, captions, and CSV summaries.
manuscript_fmt_num_fixed <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_character_, length(x))
  ok <- is.finite(x)
  if (any(ok)) {
    x_ok <- x[ok]
    x_ok[abs(x_ok) < 0.5 * 10^(-digits)] <- 0
    out[ok] <- formatC(x_ok, digits = digits, format = "f")
  }
  out
}

manuscript_round_numbers_in_text <- function(x, digits = 2) {
  # Matches ordinary numeric substrings, including scientific notation, while
  # avoiding embedded digits within words. Used only for text displayed in table cells.
  pattern <- "(?<![A-Za-z])[-+]?\\d+(?:\\.\\d+)?(?:[eE][-+]?\\d+)?(?![A-Za-z])"
  purrr::map_chr(as.character(x), function(value) {
    if (is.na(value)) return(NA_character_)
    matches <- gregexpr(pattern, value, perl = TRUE)
    found <- regmatches(value, matches)[[1]]
    if (length(found) == 0L || identical(found, "")) return(value)
    replacements <- manuscript_fmt_num_fixed(found, digits)
    regmatches(value, matches) <- list(replacements)
    value
  })
}

format_model_family_table_label <- function(x) {
  dplyr::case_when(
    as.character(x) == "linear" ~ "Linear",
    as.character(x) == "logarithmic" ~ "Logarithmic",
    as.character(x) == "quadratic" ~ "Quadratic",
    TRUE ~ trimws(gsub("\\s+candidate\\b", "", as.character(x), ignore.case = TRUE))
  )
}

format_model_family_label <- function(x) {
  dplyr::case_when(
    x == "linear" ~ "Linear candidate",
    x == "logarithmic" ~ "Logarithmic candidate",
    x == "quadratic" ~ "Quadratic candidate",
    TRUE ~ as.character(x)
  )
}

format_metric_label <- function(x, dataset = NULL) {
  values <- as.character(x)
  dataset_values <- if (is.null(dataset)) rep("unknown", length(values)) else rep_len(as.character(dataset), length(values))
  purrr::map2_chr(values, dataset_values, function(value, ds) {
    key <- manuscript_dataset_key(ds)
    if (is.na(value) || !nzchar(value)) return(NA_character_)
    if (value == "delta" && key == "topography") return(paste0("topographic ", manuscript_delta_symbol(), " (m; elevation variance / mean elevation)"))
    if (value == "delta" && key == "foliage") return(paste0("foliage-height ", manuscript_delta_symbol(), " (foliage-height variance / mean foliage height)"))
    if (value == "delta") return(paste0(manuscript_delta_symbol(), " (variance / mean)"))
    dplyr::case_when(
      value == "range" ~ "range",
      value == "variance" ~ "variance",
      value %in% c("sd", "sigma") ~ "standard deviation",
      value == "CV" ~ "coefficient of variation",
      value == "FHD3" ~ "three-layer foliage-height diversity",
      value == "gini" ~ "Gini coefficient",
      value == "FH_entropy_51" ~ "51-bin foliage-height entropy",
      value == "diff_entropy" ~ "differential entropy",
      TRUE ~ value
    )
  })
}

format_covariate_label <- function(x, dataset = NULL) {
  values <- as.character(x)
  dataset_values <- if (is.null(dataset)) rep("unknown", length(values)) else rep_len(as.character(dataset), length(values))
  purrr::map2_chr(values, dataset_values, function(value, ds) {
    key <- manuscript_dataset_key(ds)
    if (is.na(value) || value %in% c("", "none", "no_mean_covariates", "unadjusted")) return("none")
    if (value == "mean_elevation_quadratic") return("mean elevation + mean elevation\u00B2")
    if (value == "mean_elevation_linear") return("mean elevation only")
    if (value == "mean_foliage_height") return("mean foliage height")
    if (value == "mu;mu2") return(ifelse(key == "topography" || key == "unknown", "mean elevation + mean elevation\u00B2", "environmental mean + environmental mean\u00B2"))
    if (value == "mu" && key == "topography") return("mean elevation")
    if (value == "mu" && key == "foliage") return("mean foliage height")
    if (value == "mu") return("environmental mean")
    value
  })
}

format_candidate_curve_short_label <- function(x) {
  dplyr::case_when(
    as.character(x) == "linear" ~ "linear",
    as.character(x) == "logarithmic" ~ "logarithmic",
    as.character(x) == "quadratic" ~ "quadratic",
    TRUE ~ as.character(x)
  )
}

format_candidate_curve_abbrev <- function(x) {
  dplyr::case_when(
    as.character(x) == "linear" ~ "lin",
    as.character(x) == "logarithmic" ~ "log",
    as.character(x) == "quadratic" ~ "quad",
    TRUE ~ as.character(x)
  )
}

format_sensitivity_panel_label <- function(dataset, covariate_scenario, include_dataset = TRUE) {
  dataset_short <- manuscript_dataset_short(dataset)
  dataset_key <- manuscript_dataset_key(dataset)
  covariate_scenario <- as.character(covariate_scenario)

  scenario <- dplyr::case_when(
    dataset_key == "topography" & covariate_scenario %in% c("no_mean_covariates", "none") ~ "Primary\nno mean covariate",
    dataset_key == "topography" & covariate_scenario == "mean_elevation_linear" ~ "Sensitivity\nmean elevation",
    dataset_key == "topography" & covariate_scenario == "mean_elevation_quadratic" ~ paste0("Sensitivity\nmean elevation + mean elevation", manuscript_squared_symbol()),
    dataset_key == "foliage" & covariate_scenario %in% c("unadjusted", "none", "no_mean_covariates") ~ "Primary\nno mean covariate",
    dataset_key == "foliage" & covariate_scenario == "mean_foliage_height" ~ "Sensitivity\nmean foliage height",
    TRUE ~ format_covariate_label(covariate_scenario, dataset_short)
  )

  if (isTRUE(include_dataset)) paste0(dataset_short, ": ", scenario) else scenario
}

sensitivity_panel_order <- function(dataset, covariate_scenario) {
  dataset_key <- manuscript_dataset_key(dataset)
  covariate_scenario <- as.character(covariate_scenario)
  dplyr::case_when(
    dataset_key == "topography" & covariate_scenario %in% c("no_mean_covariates", "none") ~ 1,
    dataset_key == "topography" & covariate_scenario == "mean_elevation_linear" ~ 2,
    dataset_key == "topography" & covariate_scenario == "mean_elevation_quadratic" ~ 3,
    dataset_key == "foliage" & covariate_scenario %in% c("unadjusted", "none", "no_mean_covariates") ~ 1,
    dataset_key == "foliage" & covariate_scenario == "mean_foliage_height" ~ 2,
    TRUE ~ 99
  )
}

extract_best_model_row <- function(best_tbl, analysis_id) {
  out <- best_tbl %>% filter(analysis == analysis_id) %>% slice(1)
  if (nrow(out) == 0) {
    warning("No best-model row found for ", analysis_id)
    return(tibble::tibble(
      analysis = analysis_id,
      n = NA_integer_,
      best_model_family = NA_character_,
      best_Akaike_weight = NA_real_,
      second_best_delta_IC = NA_real_,
      recommended_interpretation = NA_character_,
      best_shape_interpretation = NA_character_,
      plausible_models_delta_IC_lt_2 = NA_character_
    ))
  }
  out
}

format_best_model_support <- function(row) {
  criterion <- if ("selection_criterion" %in% names(row)) row$selection_criterion[1] else "IC"
  paste0(
    format_model_family_table_label(row$best_model_family[1]),
    "; ", criterion, " weight = ", manuscript_fmt_num_fixed(row$best_Akaike_weight[1], 2),
    ifelse(
      is.finite(row$second_best_delta_IC[1]),
      paste0("; ", manuscript_Delta_symbol(), criterion, " to second = ", manuscript_fmt_num_fixed(row$second_best_delta_IC[1], 2)),
      ""
    )
  )
}

format_weight_range <- function(rows) {
  rows <- rows %>% arrange(covariate_scenario, analysis)
  weights <- rows$best_Akaike_weight
  models <- unique(format_model_family_table_label(rows$best_model_family))
  if (length(weights) == 0 || all(!is.finite(weights))) return(NA_character_)
  if (length(unique(round(weights, 6))) == 1) {
    weight_txt <- manuscript_fmt_num_fixed(weights[1], 2)
  } else {
    weight_txt <- paste0(manuscript_fmt_num_fixed(min(weights, na.rm = TRUE), 2), "-", manuscript_fmt_num_fixed(max(weights, na.rm = TRUE), 2))
  }
  paste0(paste(models, collapse = "; "), "; weights = ", weight_txt)
}

manuscript_ensure_column <- function(df, column_name, default = NA) {
  if (!column_name %in% names(df)) df[[column_name]] <- default
  df
}


assemble_delta_model_comparison_table <- function() {
  topography_registry <- results_manuscript_analysis_hierarchy$topography_delta_models$registry %>%
    select(dataset, any_of(c("dataset_short", "dataset_phrase", "source_study")), analysis, analysis_tier, heterogeneity_metric, covariate_scenario, metric_role, destination, manuscript_role, primary_or_sensitivity, dataset_display_order, analysis_display_order, rationale)
  foliage_registry <- results_manuscript_analysis_hierarchy$foliage_delta_models$registry %>%
    select(dataset, any_of(c("dataset_short", "dataset_phrase", "source_study")), analysis, analysis_tier, heterogeneity_metric, covariate_scenario, metric_role, destination, manuscript_role, primary_or_sensitivity, dataset_display_order, analysis_display_order, rationale)
  joined <- bind_rows(
    results_manuscript_analysis_hierarchy$topography_delta_models$summary %>% mutate(output_set = manuscript_delta_output_set_label("topography")),
    results_manuscript_analysis_hierarchy$foliage_delta_models$summary %>% mutate(output_set = manuscript_delta_output_set_label("foliage"))
  ) %>%
    left_join(bind_rows(topography_registry, foliage_registry), by = "analysis", suffix = c("", "_registry"))

  joined <- manuscript_ensure_column(joined, "dataset", NA_character_)
  joined <- manuscript_ensure_column(joined, "dataset_registry", NA_character_)
  joined <- manuscript_ensure_column(joined, "dataset_short", NA_character_)
  joined <- manuscript_ensure_column(joined, "dataset_short_registry", NA_character_)
  joined <- manuscript_ensure_column(joined, "dataset_phrase", NA_character_)
  joined <- manuscript_ensure_column(joined, "dataset_phrase_registry", NA_character_)
  joined <- manuscript_ensure_column(joined, "source_study", NA_character_)
  joined <- manuscript_ensure_column(joined, "source_study_registry", NA_character_)
  joined <- manuscript_ensure_column(joined, "heterogeneity_metric", NA_character_)
  joined <- manuscript_ensure_column(joined, "predictor", NA_character_)
  joined <- manuscript_ensure_column(joined, "covariate_scenario", NA_character_)
  joined <- manuscript_ensure_column(joined, "covariates", NA_character_)
  joined <- manuscript_ensure_column(joined, "analysis_tier", NA_character_)
  joined <- manuscript_ensure_column(joined, "manuscript_role", NA_character_)
  joined <- manuscript_ensure_column(joined, "primary_or_sensitivity", NA_character_)
  joined <- manuscript_ensure_column(joined, "dataset_display_order", NA_integer_)
  joined <- manuscript_ensure_column(joined, "analysis_display_order", NA_integer_)

  joined %>%
    mutate(
      dataset = dplyr::coalesce(dataset, dataset_registry, "unknown dataset"),
      dataset_short = dplyr::coalesce(dataset_short, dataset_short_registry, manuscript_dataset_short(dataset)),
      dataset_phrase = dplyr::coalesce(dataset_phrase, dataset_phrase_registry, manuscript_dataset_phrase(dataset)),
      source_study = dplyr::coalesce(source_study, source_study_registry, manuscript_dataset_source(dataset)),
      heterogeneity_metric = dplyr::coalesce(heterogeneity_metric, predictor),
      covariate_scenario = dplyr::coalesce(covariate_scenario, covariates),
      candidate_model = model_family,
      candidate_model_pretty = format_model_family_label(model_family),
      heterogeneity_metric_pretty = format_metric_label(heterogeneity_metric, dataset_short),
      covariates_pretty = format_covariate_label(covariate_scenario, dataset_short),
      response_pretty = manuscript_dataset_response(dataset_short),
      formula_display = formula_readable,
      interpretation = dplyr::case_when(
        model_family == "logarithmic" & shape_interpretation == "positive_diminishing_return" ~ "positive diminishing-return HDR",
        model_family == "linear" & shape_interpretation == "positive_monotonic" ~ "positive monotonic HDR",
        model_family == "quadratic" & shape_interpretation == "unimodal_in_observed_range" ~ "unimodal HDR candidate",
        TRUE ~ shape_interpretation
      )
    ) %>%
    filter(heterogeneity_metric == "delta") %>%
    select(
      output_set, dataset, dataset_short, dataset_phrase, source_study, analysis_tier, destination, analysis, response, response_pretty,
      manuscript_role, primary_or_sensitivity, dataset_display_order, analysis_display_order,
      heterogeneity_metric, heterogeneity_metric_pretty, metric_role, predictor_role,
      covariate_scenario, covariates, covariates_pretty, formula_display,
      candidate_model, candidate_model_pretty, formula, formula_readable, n, k, AIC, AICc,
      selection_criterion, IC, delta_IC, Akaike_weight, best_by_IC,
      predictor_min, predictor_max, slope_linear_x, slope_logx, quadratic_x_2,
      quadratic_vertex_x, quadratic_vertex_inside_observed_range,
      quadratic_high_side_declines,
      heterogeneity_term_1,
      heterogeneity_term_1_estimate,
      heterogeneity_term_1_std_error,
      heterogeneity_term_1_statistic,
      heterogeneity_term_1_p_value,
      heterogeneity_term_1_conf_low,
      heterogeneity_term_1_conf_high,
      heterogeneity_term_1_clear_positive,
      heterogeneity_term_1_clear_negative,
      heterogeneity_term_2,
      heterogeneity_term_2_estimate,
      heterogeneity_term_2_std_error,
      heterogeneity_term_2_statistic,
      heterogeneity_term_2_p_value,
      heterogeneity_term_2_conf_low,
      heterogeneity_term_2_conf_high,
      heterogeneity_term_2_clear_positive,
      heterogeneity_term_2_clear_negative,
      heterogeneity_effect_evidence,
      shape_support_with_uncertainty,
      shape_supported_by_coefficients,
      shape_support_notes,
      shape_interpretation, interpretation, r_squared, adj_r_squared, rationale
    ) %>% arrange(dataset_display_order, analysis_display_order, delta_IC)
}


assemble_delta_model_term_table <- function() {
  bind_rows(
    results_manuscript_analysis_hierarchy$topography_delta_models$terms %>% mutate(output_set = manuscript_delta_output_set_label("topography")),
    results_manuscript_analysis_hierarchy$foliage_delta_models$terms %>% mutate(output_set = manuscript_delta_output_set_label("foliage"))
  ) %>%
    filter(predictor == "delta") %>%
    mutate(
      dataset_short = manuscript_dataset_short(analysis),
      term_pretty = format_term_label(term, dataset_short),
      model_family_pretty = format_model_family_label(model_family)
    ) %>%
    arrange(output_set, analysis, model_family, term)
}

assemble_legacy_model_support_table <- function() {
  topography_best <- results_manuscript_analysis_hierarchy$topography_delta_models$best
  foliage_best <- results_manuscript_analysis_hierarchy$foliage_delta_models$best

  topography_primary <- extract_best_model_row(topography_best, "Topography_delta_no_mean_covariates_gaussian")
  topography_cov_sens <- topography_best %>%
    filter(heterogeneity_metric == "delta", covariate_scenario %in% c("mean_elevation_linear", "mean_elevation_quadratic"))
  foliage_primary <- extract_best_model_row(foliage_best, "Foliage_delta_unadjusted_gaussian")
  foliage_mean_sens <- extract_best_model_row(foliage_best, "Foliage_delta_mean_foliage_height_gaussian")

  tbl <- tibble::tibble(
    Dataset = c(rep(manuscript_dataset_short("topography"), 2), rep(manuscript_dataset_short("foliage"), 2)),
    `Dataset phrase` = c(rep(manuscript_dataset_phrase("topography"), 2), rep(manuscript_dataset_phrase("foliage"), 2)),
    Role = c(
      paste0("Primary ", manuscript_delta_symbol(), " model"),
      paste0(manuscript_delta_symbol(), " covariate sensitivity"),
      paste0("Primary ", manuscript_delta_symbol(), " model"),
      paste0(manuscript_delta_symbol(), " mean-covariate sensitivity")
    ),
    Response = c(rep(manuscript_dataset_response("topography"), 2), rep(manuscript_dataset_response("foliage"), 2)),
    `delta axis` = c(rep(manuscript_dataset_delta_axis("topography"), 2), rep(manuscript_dataset_delta_axis("foliage"), 2)),
    Covariates = c(
      "None",
      "Mean elevation; mean elevation + mean elevation\u00B2",
      "None",
      "Mean foliage height"
    ),
    n = c(
      topography_primary$n[1],
      ifelse(nrow(topography_cov_sens) > 0, topography_cov_sens$n[1], NA_integer_),
      foliage_primary$n[1],
      foliage_mean_sens$n[1]
    ),
    `Best-supported candidate` = c(
      format_best_model_support(topography_primary),
      format_weight_range(topography_cov_sens),
      format_best_model_support(foliage_primary),
      format_best_model_support(foliage_mean_sens)
    ),
    `Coefficient-qualified interpretation` = c(
      topography_primary$coefficient_qualified_interpretation[1],
      format_weight_range(topography_cov_sens),
      foliage_primary$coefficient_qualified_interpretation[1],
      foliage_mean_sens$coefficient_qualified_interpretation[1]
    ),
    Interpretation = c(
      "Primary Allouche et al. response-shape support",
      paste0("Checks whether Allouche et al. ", manuscript_delta_symbol(), " inference depends on mean-elevation covariate structure"),
      "Primary MacArthur & MacArthur response-shape support",
      paste0("Checks whether MacArthur & MacArthur ", manuscript_delta_symbol(), " inference persists after accounting for mean foliage height")
    )
  )
  manuscript_write_csv(tbl, "supporting_analysis", "support_legacy_analysis_hierarchy_and_model_support.csv")
}


write_extended_data_and_supporting_analysis_tables <- function(all_comparisons, all_terms) {
  table_1 <- build_extended_data_model_check_table(all_comparisons, "foliage")
  table_2 <- build_extended_data_model_check_table(all_comparisons, "topography")
  table_1_display <- build_extended_data_model_check_display_table(all_comparisons, "foliage")
  table_2_display <- build_extended_data_model_check_display_table(all_comparisons, "topography")

  support_comparisons <- all_comparisons %>%
    select(-dplyr::any_of(c("adj_r_squared", "adjusted_r_squared", "Adjusted R²", "Adjusted R2"))) %>%
    arrange(dataset_display_order, analysis_display_order, delta_IC)
  support_terms <- all_terms %>%
    arrange(output_set, analysis, model_family, term)
  support_registry <- results_manuscript_analysis_hierarchy$hierarchy$registry %>%
    arrange(dataset_display_order, analysis_display_order, analysis)

  manuscript_write_display_csv(table_1, "extended_data", "extended_data_table_1_macarthur_model_comparisons_and_coefficient_checks.csv")
  manuscript_write_display_csv(table_2, "extended_data", "extended_data_table_2_allouche_model_comparisons_and_coefficient_checks.csv")
  manuscript_write_csv(support_comparisons, "supporting_analysis", "support_complete_delta_model_comparisons.csv")
  manuscript_write_csv(support_terms, "supporting_analysis", "support_complete_delta_model_terms.csv")
  manuscript_write_csv(support_registry, "supporting_analysis", "support_analysis_hierarchy_registry.csv")

  list(
    extended_data_table_1 = table_1,
    extended_data_table_2 = table_2,
    extended_data_table_1_display = table_1_display,
    extended_data_table_2_display = table_2_display,
    support_complete_delta_model_comparisons = support_comparisons,
    support_complete_delta_model_terms = support_terms,
    support_analysis_hierarchy_registry = support_registry
  )
}


# 04A WORD TABLE DELIVERABLES ####
# ------------------------------------------------------------ #
# These helpers convert the machine-readable tables of model comparisons into
# Extended Data table drafts suitable for manuscript assembly.

manuscript_restore_delta_text <- function(x) {
  if (!is.character(x)) return(x)
  out <- x
  out <- gsub("\\bDelta\\b", manuscript_Delta_symbol(), out)
  out <- gsub("\\bdelta\\b", manuscript_delta_symbol(), out)
  out <- gsub(paste0("log\\(", manuscript_delta_symbol(), "\\)"), paste0("log(", manuscript_delta_symbol(), ")"), out)
  out <- gsub("delta axis", paste0(manuscript_delta_symbol(), " axis"), out, fixed = TRUE)
  out
}

manuscript_table_clean_text <- function(x) {
  if (is.numeric(x)) return(x)
  out <- as.character(x)
  out[is.na(out)] <- ""
  out <- manuscript_restore_delta_text(out)
  out <- gsub("_", " ", out)
  out <- gsub("\\s+", " ", out)
  trimws(out)
}

manuscript_fmt_p <- function(p, digits = 2) {
  p <- suppressWarnings(as.numeric(p))
  out <- rep(NA_character_, length(p))
  threshold <- 10^(-digits)
  out[is.finite(p) & p < threshold] <- paste0("<", manuscript_fmt_num_fixed(threshold, digits))
  idx <- is.finite(p) & p >= threshold
  out[idx] <- manuscript_fmt_num_fixed(p[idx], digits)
  out
}

manuscript_fmt_p_table <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  out <- rep(NA_character_, length(p))
  out[is.finite(p) & p < 0.001] <- "<0.001"
  idx <- is.finite(p) & p >= 0.001 & p < 0.01
  out[idx] <- manuscript_fmt_num_fixed(p[idx], 4)
  idx <- is.finite(p) & p >= 0.01 & p < 0.1
  out[idx] <- manuscript_fmt_num_fixed(p[idx], 3)
  idx <- is.finite(p) & p >= 0.1
  out[idx] <- manuscript_fmt_num_fixed(p[idx], 3)
  out
}

manuscript_fmt_weight_table <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_character_, length(x))
  idx <- is.finite(x)
  out[idx] <- ifelse(x[idx] < 0.01, manuscript_fmt_num_fixed(x[idx], 4), manuscript_fmt_num_fixed(x[idx], 3))
  out
}

manuscript_fmt_delta_ic_display <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_character_, length(x))
  ok <- is.finite(x)
  out[ok] <- ifelse(abs(x[ok]) < 0.005, "0", manuscript_fmt_num_fixed(x[ok], 2))
  out
}

manuscript_fmt_weight_display <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_character_, length(x))
  ok <- is.finite(x)
  out[ok] <- dplyr::case_when(
    x[ok] > 0 & x[ok] < 0.01 ~ "<0.01",
    abs(x[ok]) < 0.005 ~ "0",
    TRUE ~ manuscript_fmt_num_fixed(x[ok], 2)
  )
  out
}

manuscript_fmt_num_adaptive <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_character_, length(x))
  ok <- is.finite(x)
  if (any(ok)) {
    ax <- abs(x[ok])
    digits <- dplyr::case_when(
      ax == 0 ~ 0L,
      ax >= 0.01 ~ 2L,
      ax >= 0.001 ~ 3L,
      ax >= 0.0001 ~ 4L,
      ax >= 0.00001 ~ 5L,
      TRUE ~ 6L
    )
    out[ok] <- purrr::map2_chr(x[ok], digits, ~ if (.y == 0L) "0" else formatC(.x, digits = .y, format = "f"))
  }
  out
}

manuscript_fmt_ci <- function(estimate, conf_low, conf_high, digits = 2) {
  estimate <- suppressWarnings(as.numeric(estimate))
  conf_low <- suppressWarnings(as.numeric(conf_low))
  conf_high <- suppressWarnings(as.numeric(conf_high))
  out <- rep("", max(length(estimate), length(conf_low), length(conf_high), 1L))
  estimate <- rep_len(estimate, length(out))
  conf_low <- rep_len(conf_low, length(out))
  conf_high <- rep_len(conf_high, length(out))
  ok <- is.finite(estimate) & is.finite(conf_low) & is.finite(conf_high)
  out[ok] <- paste0(
    manuscript_fmt_num_fixed(estimate[ok], digits),
    " [", manuscript_fmt_num_fixed(conf_low[ok], digits),
    ", ", manuscript_fmt_num_fixed(conf_high[ok], digits), "]"
  )
  out
}

manuscript_fmt_ci_adaptive <- function(estimate, conf_low, conf_high) {
  estimate <- suppressWarnings(as.numeric(estimate))
  conf_low <- suppressWarnings(as.numeric(conf_low))
  conf_high <- suppressWarnings(as.numeric(conf_high))
  out <- rep("", max(length(estimate), length(conf_low), length(conf_high), 1L))
  estimate <- rep_len(estimate, length(out))
  conf_low <- rep_len(conf_low, length(out))
  conf_high <- rep_len(conf_high, length(out))
  ok <- is.finite(estimate) & is.finite(conf_low) & is.finite(conf_high)
  out[ok] <- paste0(
    manuscript_fmt_num_adaptive(estimate[ok]),
    " [", manuscript_fmt_num_adaptive(conf_low[ok]),
    ", ", manuscript_fmt_num_adaptive(conf_high[ok]), "]"
  )
  out
}

manuscript_term_label_for_table <- function(term) {
  dplyr::case_when(
    as.character(term) == "x" ~ manuscript_delta_symbol(),
    as.character(term) == "logx" ~ paste0("log(", manuscript_delta_symbol(), ")"),
    as.character(term) == "I(x^2)" ~ paste0(manuscript_delta_symbol(), manuscript_squared_symbol()),
    is.na(as.character(term)) | !nzchar(as.character(term)) ~ "",
    TRUE ~ as.character(term)
  )
}

manuscript_format_term_summary <- function(term_1, estimate_1, conf_low_1, conf_high_1, p_1,
                                           term_2 = NA_character_, estimate_2 = NA_real_,
                                           conf_low_2 = NA_real_, conf_high_2 = NA_real_,
                                           p_2 = NA_real_, digits = 2, table_p_values = FALSE) {
  n <- max(
    length(term_1), length(estimate_1), length(conf_low_1), length(conf_high_1), length(p_1),
    length(term_2), length(estimate_2), length(conf_low_2), length(conf_high_2), length(p_2), 1L
  )
  term_1 <- rep_len(term_1, n)
  estimate_1 <- rep_len(estimate_1, n)
  conf_low_1 <- rep_len(conf_low_1, n)
  conf_high_1 <- rep_len(conf_high_1, n)
  p_1 <- rep_len(p_1, n)
  term_2 <- rep_len(term_2, n)
  estimate_2 <- rep_len(estimate_2, n)
  conf_low_2 <- rep_len(conf_low_2, n)
  conf_high_2 <- rep_len(conf_high_2, n)
  p_2 <- rep_len(p_2, n)

  purrr::map_chr(seq_len(n), function(i) {
    t1 <- manuscript_term_label_for_table(term_1[i])
    ci1 <- manuscript_fmt_ci(estimate_1[i], conf_low_1[i], conf_high_1[i], digits = digits)[1]
    ptxt1 <- if (isTRUE(table_p_values)) manuscript_fmt_p_table(p_1[i])[1] else manuscript_fmt_p(p_1[i])[1]
    part1 <- if (nzchar(t1) && nzchar(ci1)) {
      paste0(t1, ": ", ci1, ", P = ", ptxt1)
    } else {
      ""
    }

    t2 <- manuscript_term_label_for_table(term_2[i])
    ci2 <- manuscript_fmt_ci(estimate_2[i], conf_low_2[i], conf_high_2[i], digits = digits)[1]
    ptxt2 <- if (isTRUE(table_p_values)) manuscript_fmt_p_table(p_2[i])[1] else manuscript_fmt_p(p_2[i])[1]
    part2 <- if (nzchar(t2) && nzchar(ci2)) {
      paste0(t2, ": ", ci2, ", P = ", ptxt2)
    } else {
      ""
    }

    parts <- c(part1, part2)
    parts <- parts[nzchar(parts)]
    if (length(parts) == 0L) "not estimated" else paste(parts, collapse = "; ")
  })
}

manuscript_format_term_summary_display <- function(term_1, estimate_1, conf_low_1, conf_high_1,
                                                   term_2 = NA_character_, estimate_2 = NA_real_,
                                                   conf_low_2 = NA_real_, conf_high_2 = NA_real_) {
  n <- max(
    length(term_1), length(estimate_1), length(conf_low_1), length(conf_high_1),
    length(term_2), length(estimate_2), length(conf_low_2), length(conf_high_2), 1L
  )
  term_1 <- rep_len(term_1, n)
  estimate_1 <- rep_len(estimate_1, n)
  conf_low_1 <- rep_len(conf_low_1, n)
  conf_high_1 <- rep_len(conf_high_1, n)
  term_2 <- rep_len(term_2, n)
  estimate_2 <- rep_len(estimate_2, n)
  conf_low_2 <- rep_len(conf_low_2, n)
  conf_high_2 <- rep_len(conf_high_2, n)

  purrr::map_chr(seq_len(n), function(i) {
    t1 <- manuscript_term_label_for_table(term_1[i])
    ci1 <- manuscript_fmt_ci_adaptive(estimate_1[i], conf_low_1[i], conf_high_1[i])[1]
    part1 <- if (nzchar(t1) && nzchar(ci1)) paste0(t1, ": ", ci1) else ""

    t2 <- manuscript_term_label_for_table(term_2[i])
    ci2 <- manuscript_fmt_ci_adaptive(estimate_2[i], conf_low_2[i], conf_high_2[i])[1]
    part2 <- if (nzchar(t2) && nzchar(ci2)) paste0(t2, ": ", ci2) else ""

    parts <- c(part1, part2)
    parts <- parts[nzchar(parts)]
    if (length(parts) == 0L) "not estimated" else paste(parts, collapse = "; ")
  })
}

manuscript_compact_model_support_text <- function(x) {
  out <- manuscript_table_clean_text(x)
  out <- gsub("\\s+candidate\\b", "", out, ignore.case = TRUE)
  out <- gsub("AICc weight =", "AICc w =", out, fixed = TRUE)
  out <- gsub("AIC weight =", "AIC w =", out, fixed = TRUE)
  out <- gsub("IC weight =", "IC w =", out, fixed = TRUE)
  out <- gsub("to second =", "to second-best =", out, fixed = TRUE)
  manuscript_round_numbers_in_text(out, 2)
}

manuscript_compact_inference_text <- function(x) {
  out <- manuscript_table_clean_text(x)
  lower <- tolower(out)
  dplyr::case_when(
    grepl("positive diminishing-return hdr supported", lower) ~
      "Positive diminishing-return HDR; coefficient sign, P value, and 95% CI support the fitted shape.",
    grepl("positive monotonic hdr supported", lower) ~
      "Positive monotonic HDR; coefficient sign, P value, and 95% CI support the fitted shape.",
    grepl("unimodal hdr supported", lower) ~
      "Unimodal HDR; coefficient signs, P values, 95% CIs, and within-range geometry support the fitted shape.",
    grepl("ic-selected logarithmic", lower) ~
      "Logarithmic curve selected by information criterion; coefficient evidence should qualify the diminishing-return interpretation.",
    grepl("ic-selected quadratic", lower) ~
      "Quadratic curve selected by information criterion; coefficient and geometry checks should qualify any unimodal interpretation.",
    grepl("ic-selected linear", lower) ~
      "Linear curve selected by information criterion; coefficient evidence should qualify any monotonic interpretation.",
    TRUE ~ out
  )
}

manuscript_compact_adjustment_text <- function(delta_axis, covariates) {
  axis <- manuscript_table_clean_text(delta_axis)
  cov <- manuscript_table_clean_text(covariates)
  cov <- ifelse(!nzchar(cov) | tolower(cov) %in% c("none", "unadjusted"), "unadjusted", paste0("adjusted for ", cov))
  paste0(axis, "; ", cov)
}

build_main_text_nature_model_support_table <- function(main_text_table) {
  required_cols <- c(
    "Dataset", "Role", "Response", "delta axis", "Covariates", "n",
    "Best-supported candidate", "Coefficient-qualified interpretation"
  )
  missing_cols <- setdiff(required_cols, names(main_text_table))
  if (length(missing_cols) > 0L) {
    stop("Missing column(s) in main text model-support table: ", paste(missing_cols, collapse = ", "))
  }

  main_text_table %>%
    filter(grepl("^Primary", .data$Role)) %>%
    mutate(
      Dataset = factor(.data$Dataset, levels = c(manuscript_dataset_short("topography"), manuscript_dataset_short("foliage"))),
      `Modelled heterogeneity and adjustment` = manuscript_compact_adjustment_text(.data[["delta axis"]], .data$Covariates),
      `Best-supported candidate curve` = manuscript_compact_model_support_text(.data[["Best-supported candidate"]]),
      `Inference` = manuscript_compact_inference_text(.data[["Coefficient-qualified interpretation"]])
    ) %>%
    arrange(.data$Dataset) %>%
    mutate(Dataset = as.character(.data$Dataset)) %>%
    select(
      Dataset,
      `Biodiversity response` = Response,
      `Modelled heterogeneity and adjustment`,
      n,
      `Best-supported candidate curve`,
      Inference
    )
}

manuscript_short_shape_check <- function(x) {
  raw <- as.character(x)
  out <- manuscript_table_clean_text(x)
  dplyr::case_when(
    raw == "positive_diminishing_return_supported_by_coefficient_CI" ~ "supported positive diminishing return",
    raw == "positive_monotonic_supported_by_coefficient_CI" ~ "supported positive monotonic",
    raw == "negative_monotonic_supported_by_coefficient_CI" ~ "supported negative monotonic",
    raw == "unimodal_supported_by_coefficients_CI_and_geometry" ~ "supported unimodal",
    grepl("IC selected", out, fixed = TRUE) | grepl("IC_selected", raw, fixed = TRUE) ~ "selected shape uncertain by coefficient/geometry checks",
    TRUE ~ out
  )
}

manuscript_extended_data_scenario_label <- function(dataset, covariate_scenario) {
  dataset_key <- manuscript_dataset_key(dataset)
  covariate_scenario <- as.character(covariate_scenario)
  dplyr::case_when(
    dataset_key == "foliage" & covariate_scenario %in% c("unadjusted", "none", "no_mean_covariates") ~ "no mean covariate",
    dataset_key == "foliage" & covariate_scenario == "mean_foliage_height" ~ "mean foliage height",
    dataset_key == "topography" & covariate_scenario %in% c("no_mean_covariates", "none") ~ "no mean covariate",
    dataset_key == "topography" & covariate_scenario == "mean_elevation_linear" ~ "mean elevation only",
    dataset_key == "topography" & covariate_scenario == "mean_elevation_quadratic" ~ paste0("mean elevation + mean elevation", manuscript_squared_symbol()),
    TRUE ~ format_covariate_label(covariate_scenario, dataset)
  )
}

manuscript_extended_data_adjustment_display_label <- function(dataset, covariate_scenario) {
  dataset_key <- manuscript_dataset_key(dataset)
  covariate_scenario <- as.character(covariate_scenario)
  dplyr::case_when(
    dataset_key == "foliage" & covariate_scenario %in% c("unadjusted", "none", "no_mean_covariates") ~ "none",
    dataset_key == "foliage" & covariate_scenario == "mean_foliage_height" ~ "mean foliage height",
    dataset_key == "topography" & covariate_scenario %in% c("no_mean_covariates", "none") ~ "none",
    dataset_key == "topography" & covariate_scenario == "mean_elevation_linear" ~ "mean elevation only",
    dataset_key == "topography" & covariate_scenario == "mean_elevation_quadratic" ~ paste0("mean elevation + mean elevation", manuscript_squared_symbol()),
    TRUE ~ format_covariate_label(covariate_scenario, dataset)
  )
}

manuscript_extended_data_role_label <- function(primary_or_sensitivity) {
  dplyr::case_when(
    as.character(primary_or_sensitivity) == "primary" ~ "Primary",
    as.character(primary_or_sensitivity) == "sensitivity" ~ "Sensitivity",
    TRUE ~ manuscript_table_clean_text(primary_or_sensitivity)
  )
}

manuscript_ensure_extended_data_order_columns <- function(all_comparisons) {
  if (!"dataset_display_order" %in% names(all_comparisons)) {
    all_comparisons$dataset_display_order <- ifelse(manuscript_dataset_key(all_comparisons$dataset_short) == "foliage", 1L, 2L)
  }
  if (!"analysis_display_order" %in% names(all_comparisons)) {
    all_comparisons$analysis_display_order <- sensitivity_panel_order(all_comparisons$dataset_short, all_comparisons$covariate_scenario)
  }
  if (!"primary_or_sensitivity" %in% names(all_comparisons)) {
    all_comparisons$primary_or_sensitivity <- ifelse(all_comparisons$analysis_display_order == 1L, "primary", "sensitivity")
  }
  all_comparisons
}

build_extended_data_model_check_table <- function(all_comparisons, dataset_key) {
  dataset_key <- manuscript_dataset_key(dataset_key)
  criterion <- if (dataset_key == "foliage") "AICc" else "AIC"
  ic_col <- criterion
  delta_col <- paste0(manuscript_Delta_symbol(), criterion)
  weight_col <- paste0(criterion, " weight")

  required_cols <- c(
    "dataset_short", "analysis", "response_pretty", "heterogeneity_metric_pretty",
    "covariate_scenario", "covariates_pretty", "candidate_model_pretty",
    "formula_readable", criterion, "delta_IC", "Akaike_weight",
    "best_by_IC", "heterogeneity_term_1", "heterogeneity_term_1_estimate",
    "heterogeneity_term_1_p_value", "heterogeneity_term_1_conf_low",
    "heterogeneity_term_1_conf_high", "heterogeneity_term_2",
    "heterogeneity_term_2_estimate", "heterogeneity_term_2_p_value",
    "heterogeneity_term_2_conf_low", "heterogeneity_term_2_conf_high",
    "shape_support_with_uncertainty", "interpretation"
  )
  missing_cols <- setdiff(required_cols, names(all_comparisons))
  if (length(missing_cols) > 0L) {
    stop("Missing column(s) for Extended Data model-check table: ", paste(missing_cols, collapse = ", "))
  }

  all_comparisons <- manuscript_ensure_extended_data_order_columns(all_comparisons)

  all_comparisons %>%
    filter(manuscript_dataset_key(.data$dataset_short) == .env$dataset_key) %>%
    mutate(
      dataset_display_order = dplyr::coalesce(.data$dataset_display_order, ifelse(.env$dataset_key == "foliage", 1L, 2L)),
      analysis_display_order = dplyr::coalesce(.data$analysis_display_order, sensitivity_panel_order(.data$dataset_short, .data$covariate_scenario)),
      primary_or_sensitivity = dplyr::coalesce(.data$primary_or_sensitivity, ifelse(sensitivity_panel_order(.data$dataset_short, .data$covariate_scenario) == 1L, "primary", "sensitivity")),
      Dataset = manuscript_dataset_short(.data$dataset_short),
      `Analysis role` = manuscript_extended_data_role_label(.data$primary_or_sensitivity),
      Scenario = manuscript_extended_data_scenario_label(.data$dataset_short, .data$covariate_scenario),
      Response = manuscript_table_clean_text(.data$response_pretty),
      `δ metric` = manuscript_table_clean_text(.data$heterogeneity_metric_pretty),
      `Mean covariate adjustment` = manuscript_table_clean_text(.data$covariates_pretty),
      `Candidate model` = format_model_family_table_label(.data$candidate_model_pretty),
      Formula = manuscript_restore_delta_text(.data$formula_readable),
      !!ic_col := manuscript_fmt_num_fixed(.data[[criterion]], 2),
      !!delta_col := manuscript_fmt_num_fixed(.data$delta_IC, 2),
      !!weight_col := manuscript_fmt_weight_table(.data$Akaike_weight),
      `Best?` = ifelse(.data$best_by_IC, "yes", ""),
      `Focal heterogeneity term(s)` = manuscript_format_term_summary(
        .data$heterogeneity_term_1,
        .data$heterogeneity_term_1_estimate,
        .data$heterogeneity_term_1_conf_low,
        .data$heterogeneity_term_1_conf_high,
        .data$heterogeneity_term_1_p_value,
        .data$heterogeneity_term_2,
        .data$heterogeneity_term_2_estimate,
        .data$heterogeneity_term_2_conf_low,
        .data$heterogeneity_term_2_conf_high,
        .data$heterogeneity_term_2_p_value,
        digits = 3,
        table_p_values = TRUE
      ),
      `Coefficient / geometry evidence` = manuscript_short_shape_check(.data$shape_support_with_uncertainty),
      Interpretation = manuscript_table_clean_text(.data$interpretation)
    ) %>%
    arrange(.data$analysis_display_order, .data$delta_IC) %>%
    select(
      Dataset,
      `Analysis role`,
      Scenario,
      Response,
      `δ metric`,
      `Mean covariate adjustment`,
      `Candidate model`,
      Formula,
      n,
      Criterion = selection_criterion,
      dplyr::all_of(ic_col),
      dplyr::all_of(delta_col),
      dplyr::all_of(weight_col),
      `Best?`,
      `Focal heterogeneity term(s)`,
      `Coefficient / geometry evidence`,
      Interpretation
    )
}

build_extended_data_model_check_display_table <- function(all_comparisons, dataset_key) {
  dataset_key <- manuscript_dataset_key(dataset_key)
  criterion <- if (dataset_key == "foliage") "AICc" else "AIC"
  ic_col <- criterion
  delta_col <- paste0(manuscript_Delta_symbol(), criterion)

  required_cols <- c(
    "dataset_short", "covariate_scenario", "candidate_model_pretty",
    criterion, "delta_IC", "Akaike_weight", "heterogeneity_term_1",
    "heterogeneity_term_1_estimate", "heterogeneity_term_1_conf_low",
    "heterogeneity_term_1_conf_high", "heterogeneity_term_2",
    "heterogeneity_term_2_estimate", "heterogeneity_term_2_conf_low",
    "heterogeneity_term_2_conf_high"
  )
  missing_cols <- setdiff(required_cols, names(all_comparisons))
  if (length(missing_cols) > 0L) {
    stop("Missing column(s) for Extended Data model-check display table: ", paste(missing_cols, collapse = ", "))
  }

  all_comparisons <- manuscript_ensure_extended_data_order_columns(all_comparisons)

  base <- all_comparisons %>%
    filter(manuscript_dataset_key(.data$dataset_short) == .env$dataset_key) %>%
    mutate(
      analysis_display_order = dplyr::coalesce(.data$analysis_display_order, sensitivity_panel_order(.data$dataset_short, .data$covariate_scenario)),
      primary_or_sensitivity = dplyr::coalesce(.data$primary_or_sensitivity, ifelse(sensitivity_panel_order(.data$dataset_short, .data$covariate_scenario) == 1L, "primary", "sensitivity")),
      .role = manuscript_extended_data_role_label(.data$primary_or_sensitivity),
      .adjustment = manuscript_extended_data_adjustment_display_label(.data$dataset_short, .data$covariate_scenario),
      Candidate = format_model_family_table_label(.data$candidate_model_pretty),
      !!ic_col := manuscript_fmt_num_fixed(.data[[criterion]], 2),
      !!delta_col := manuscript_fmt_delta_ic_display(.data$delta_IC),
      Weight = manuscript_fmt_weight_display(.data$Akaike_weight),
      `Focal heterogeneity term(s)` = manuscript_format_term_summary_display(
        .data$heterogeneity_term_1,
        .data$heterogeneity_term_1_estimate,
        .data$heterogeneity_term_1_conf_low,
        .data$heterogeneity_term_1_conf_high,
        .data$heterogeneity_term_2,
        .data$heterogeneity_term_2_estimate,
        .data$heterogeneity_term_2_conf_low,
        .data$heterogeneity_term_2_conf_high
      )
    ) %>%
    arrange(.data$analysis_display_order, .data$delta_IC) %>%
    group_by(.data$analysis_display_order) %>%
    mutate(Hierarchy = ifelse(dplyr::row_number() == 1L, .data$.role, "")) %>%
    ungroup() %>%
    select(
      .analysis_display_order = analysis_display_order,
      Hierarchy,
      Adjustment = .adjustment,
      Candidate,
      dplyr::all_of(ic_col),
      dplyr::all_of(delta_col),
      Weight,
      `Focal heterogeneity term(s)`
    )

  orders <- unique(base$.analysis_display_order)
  display_cols <- setdiff(names(base), ".analysis_display_order")
  blank <- tibble::as_tibble(setNames(rep(list(""), length(display_cols)), display_cols))
  dplyr::bind_rows(lapply(seq_along(orders), function(i) {
    block <- base %>%
      filter(.data$.analysis_display_order == orders[[i]]) %>%
      select(dplyr::all_of(display_cols))
    if (i < length(orders)) dplyr::bind_rows(block, blank) else block
  }))
}

write_markdown_table <- function(tbl, output_path, table_title, table_caption, table_note = NULL) {
  ensure_parent_dir(output_path)
  clean <- tbl %>% mutate(across(everything(), ~ manuscript_table_clean_text(.x)))
  header <- paste(names(clean), collapse = " | ")
  rule <- paste(rep("---", ncol(clean)), collapse = " | ")
  body <- apply(as.data.frame(clean, stringsAsFactors = FALSE), 1, function(row) {
    paste(gsub("\\|", "\\\\|", row), collapse = " | ")
  })
  lines <- c(
    paste0("# ", table_title),
    "",
    table_caption,
    "",
    paste0("| ", header, " |"),
    paste0("| ", rule, " |"),
    paste0("| ", body, " |")
  )
  if (!is.null(table_note) && length(table_note) > 0L) {
    lines <- c(lines, "", paste0("Note: ", manuscript_restore_delta_text(table_note)))
  }
  writeLines(lines, output_path)
  normalizePath(output_path, winslash = "/", mustWork = FALSE)
}

build_supporting_analysis_nature_model_support_table <- function(all_comparisons) {
  required_cols <- c(
    "dataset_short", "response_pretty", "heterogeneity_metric", "covariate_scenario", "covariates_pretty",
    "candidate_model_pretty", "n", "selection_criterion", "delta_IC",
    "Akaike_weight", "adj_r_squared", "best_by_IC",
    "heterogeneity_term_1", "heterogeneity_term_1_estimate", "heterogeneity_term_1_p_value",
    "heterogeneity_term_1_conf_low", "heterogeneity_term_1_conf_high",
    "heterogeneity_term_2", "heterogeneity_term_2_estimate", "heterogeneity_term_2_p_value",
    "heterogeneity_term_2_conf_low", "heterogeneity_term_2_conf_high",
    "shape_support_with_uncertainty"
  )
  missing_cols <- setdiff(required_cols, names(all_comparisons))
  if (length(missing_cols) > 0L) {
    stop("Missing column(s) in delta model-comparison table: ", paste(missing_cols, collapse = ", "))
  }

  all_comparisons %>%
    filter(.data$heterogeneity_metric == "delta") %>%
    mutate(
      Dataset = manuscript_table_clean_text(.data$dataset_short),
      `Analysis / adjustment` = paste0(
        manuscript_table_clean_text(.data$response_pretty),
        "; ", manuscript_table_clean_text(.data$covariates_pretty)
      ),
      `Candidate curve` = format_model_family_table_label(.data$candidate_model_pretty),
      `Best?` = ifelse(.data$best_by_IC, "yes", ""),
      Criterion = manuscript_table_clean_text(.data$selection_criterion),
      !!paste0(manuscript_Delta_symbol(), "IC") := manuscript_fmt_num_fixed(.data$delta_IC, 2),
      Weight = manuscript_fmt_num_fixed(.data$Akaike_weight, 2),
      "Adjusted R\u00B2" = manuscript_fmt_num_fixed(.data$adj_r_squared, 2),
      `Focal heterogeneity term(s)` = manuscript_format_term_summary(
        .data$heterogeneity_term_1,
        .data$heterogeneity_term_1_estimate,
        .data$heterogeneity_term_1_conf_low,
        .data$heterogeneity_term_1_conf_high,
        .data$heterogeneity_term_1_p_value,
        .data$heterogeneity_term_2,
        .data$heterogeneity_term_2_estimate,
        .data$heterogeneity_term_2_conf_low,
        .data$heterogeneity_term_2_conf_high,
        .data$heterogeneity_term_2_p_value
      ),
      `Coefficient / shape check` = manuscript_short_shape_check(.data$shape_support_with_uncertainty),
      order_dataset = factor(.data$Dataset, levels = c(manuscript_dataset_short("topography"), manuscript_dataset_short("foliage"))),
      order_sensitivity = sensitivity_panel_order(.data$dataset_short, .data$covariate_scenario)
    ) %>%
    arrange(.data$order_dataset, .data$order_sensitivity, .data$delta_IC) %>%
    select(
      Dataset,
      `Analysis / adjustment`,
      `Candidate curve`,
      n,
      Criterion,
      dplyr::all_of(paste0(manuscript_Delta_symbol(), "IC")),
      Weight,
      dplyr::all_of("Adjusted R\u00B2"),
      `Best?`,
      `Focal heterogeneity term(s)`,
      `Coefficient / shape check`
    )
}

nature_style_flextable <- function(tbl, font_size = 8, header_size = 8, autofit = TRUE,
                                   col_widths = NULL, table_layout = "autofit") {
  if (!requireNamespace("flextable", quietly = TRUE) || !requireNamespace("officer", quietly = TRUE)) {
    stop("Writing Word tables requires the officer and flextable packages.")
  }
  header_border <- officer::fp_border(color = "#000000", width = 0.75)
  body_border <- officer::fp_border(color = "#BFBFBF", width = 0.25)
  bottom_border <- officer::fp_border(color = "#000000", width = 0.75)

  ft <- flextable::flextable(tbl)
  ft <- flextable::border_remove(ft)
  ft <- flextable::font(ft, fontname = "Arial", part = "all")
  ft <- flextable::fontsize(ft, size = font_size, part = "body")
  ft <- flextable::fontsize(ft, size = header_size, part = "header")
  ft <- flextable::bold(ft, bold = TRUE, part = "header")
  ft <- flextable::align(ft, align = "left", part = "all")
  numeric_cols <- intersect(
    c("n", "AIC", "AICc", paste0(manuscript_Delta_symbol(), "IC"), paste0(manuscript_Delta_symbol(), "AIC"), paste0(manuscript_Delta_symbol(), "AICc"), "Weight", "AIC weight", "AICc weight"),
    names(tbl)
  )
  ft <- flextable::align(ft, j = numeric_cols, align = "right", part = "all")
  ft <- flextable::valign(ft, valign = "top", part = "all")
  ft <- flextable::padding(ft, padding = 2, part = "all")
  ft <- flextable::hline_top(ft, border = header_border, part = "header")
  ft <- flextable::hline_bottom(ft, border = header_border, part = "header")
  ft <- flextable::hline(ft, border = body_border, part = "body")
  ft <- flextable::hline_bottom(ft, border = bottom_border, part = "body")
  if ("word_wrap" %in% getNamespaceExports("flextable")) {
    ft <- flextable::word_wrap(ft, j = names(tbl), value = TRUE)
  }
  if (!is.null(col_widths)) {
    if (is.null(names(col_widths))) {
      if (length(col_widths) != ncol(tbl)) {
        stop("Unnamed col_widths must have one value per table column.", call. = FALSE)
      }
      names(col_widths) <- names(tbl)
    }
    col_widths <- col_widths[names(col_widths) %in% names(tbl)]
    if (length(col_widths) > 0L) {
      ft <- flextable::width(ft, j = names(col_widths), width = as.numeric(col_widths), unit = "in")
    }
  }
  if (identical(table_layout, "fixed") && "set_table_properties" %in% getNamespaceExports("flextable")) {
    ft <- flextable::set_table_properties(ft, layout = "fixed", width = 1, align = "left")
  }
  if (isTRUE(autofit)) ft <- flextable::autofit(ft)
  ft
}

write_nature_docx_table <- function(tbl, output_path, table_title, table_caption,
                                    table_note = NULL, landscape = FALSE,
                                    font_size = 8, header_size = 8,
                                    autofit = TRUE, col_widths = NULL,
                                    table_layout = "autofit") {
  if (!requireNamespace("officer", quietly = TRUE) || !requireNamespace("flextable", quietly = TRUE)) {
    stop("Writing Word tables requires the officer and flextable packages.")
  }
  ensure_parent_dir(output_path)
  ft <- nature_style_flextable(
    tbl,
    font_size = font_size,
    header_size = header_size,
    autofit = autofit,
    col_widths = col_widths,
    table_layout = table_layout
  )

  doc <- officer::read_docx()
  margins <- officer::page_mar(top = 0.6, bottom = 0.6, left = 0.55, right = 0.55)
  section <- if (isTRUE(landscape)) {
    officer::prop_section(page_size = officer::page_size(orient = "landscape"), page_margins = margins)
  } else {
    officer::prop_section(page_size = officer::page_size(orient = "portrait"), page_margins = margins)
  }
  doc <- officer::body_set_default_section(doc, section)
  title_text <- officer::fp_text(font.family = "Arial", font.size = 13, bold = TRUE)
  doc <- officer::body_add_fpar(doc, officer::fpar(officer::ftext(table_title, prop = title_text)))
  doc <- officer::body_add_par(doc, value = table_caption, style = "Normal")
  doc <- flextable::body_add_flextable(doc, ft)
  if (!is.null(table_note) && length(table_note) > 0L) {
    for (note in table_note) {
      doc <- officer::body_add_par(doc, value = manuscript_restore_delta_text(note), style = "Normal")
    }
  }
  print(doc, target = output_path)
  normalizePath(output_path, winslash = "/", mustWork = FALSE)
}

write_legacy_nature_style_model_support_docx_tables <- function(
    main_text_table = NULL,
    all_comparisons = NULL,
    primary_support_csv_path = file.path(manuscript_deliverable_dirs$supporting_analysis, "support_primary_analysis_hierarchy_and_model_support.csv"),
    all_comparisons_csv_path = file.path(manuscript_deliverable_dirs$supporting_analysis, "support_complete_delta_model_comparisons.csv"),
    compact_support_docx_path = file.path(manuscript_deliverable_dirs$supporting_analysis, "support_primary_analysis_hierarchy_and_model_support.docx"),
    detailed_support_docx_path = file.path(manuscript_deliverable_dirs$supporting_analysis, "support_legacy_delta_model_support_table.docx")) {

  if (is.null(main_text_table)) {
    if (!file.exists(primary_support_csv_path)) {
      stop("Cannot find primary-analysis support table CSV: ", primary_support_csv_path)
    }
    main_text_table <- readr::read_csv(primary_support_csv_path, show_col_types = FALSE)
  }

  if (is.null(all_comparisons) && file.exists(all_comparisons_csv_path)) {
    all_comparisons <- readr::read_csv(all_comparisons_csv_path, show_col_types = FALSE)
  }

  main_tbl <- build_main_text_nature_model_support_table(main_text_table)
  main_caption <- "Legacy support table | Superseded by Extended Data Table 1 and Extended Data Table 2."
  main_note <- c(
    "This legacy compact support table is retained only for backward-compatible local support use. It is not a current main-text manuscript deliverable.",
    "Use Extended Data Table 1 for MacArthur & MacArthur model-comparison evidence and Extended Data Table 2 for Allouche et al. model-comparison evidence."
  )

  main_docx <- write_nature_docx_table(
    main_tbl,
    output_path = compact_support_docx_path,
    table_title = "Legacy model-support table",
    table_caption = main_caption,
    table_note = main_note,
    landscape = FALSE,
    font_size = 8,
    header_size = 8
  )

  detailed_support_docx <- NA_character_
  detailed_support_tbl <- NULL
  if (!is.null(all_comparisons)) {
    detailed_support_tbl <- build_supporting_analysis_nature_model_support_table(all_comparisons)
    detailed_support_caption <- paste0("Legacy support table | Candidate-model support for ", manuscript_delta_symbol(), "-only HDR analyses and sensitivities.")
    detailed_support_note <- c(
      paste0("Rows are candidate models within each ", manuscript_delta_symbol(), "-only analysis. ", manuscript_Delta_symbol(), "IC, Akaike weights, and adjusted R\u00B2 are computed within each analysis. The focal heterogeneity term column gives coefficient estimates with 95% confidence intervals and P values for ", manuscript_delta_symbol(), ", log(", manuscript_delta_symbol(), "), or ", manuscript_delta_symbol(), " and ", manuscript_delta_symbol(), manuscript_squared_symbol(), ", depending on the candidate model."),
      "This longer legacy support table is superseded for manuscript assembly by Extended Data Table 1 and Extended Data Table 2."
    )
    detailed_support_docx <- write_nature_docx_table(
      detailed_support_tbl,
      output_path = detailed_support_docx_path,
      table_title = "Legacy model-support table",
      table_caption = detailed_support_caption,
      table_note = detailed_support_note,
      landscape = TRUE,
      font_size = 7,
      header_size = 7
    )
  }

  list(
    main_text_table = main_tbl,
    detailed_support_table = detailed_support_tbl,
    main_docx = main_docx,
    detailed_support_docx = detailed_support_docx
  )
}

write_extended_data_table_drafts <- function(table_sets = NULL) {
  if (is.null(table_sets)) {
    table_sets <- write_extended_data_and_supporting_analysis_tables(
      assemble_delta_model_comparison_table(),
      assemble_delta_model_term_table()
    )
  }

  table_1_display <- if (!is.null(table_sets$extended_data_table_1_display)) table_sets$extended_data_table_1_display else table_sets$extended_data_table_1
  table_2_display <- if (!is.null(table_sets$extended_data_table_2_display)) table_sets$extended_data_table_2_display else table_sets$extended_data_table_2

  table_1_caption <- "Extended Data Table 1 | Model comparisons and coefficient checks for MacArthur & MacArthur (1961) dataset reanalyses (n = 13). Candidate linear, quadratic, and logarithmic HDR models were compared using AICc. The primary unadjusted, no-covariate model set is presented first, followed by the mean-foliage-height sensitivity analysis."
  table_2_caption <- "Extended Data Table 2 | Model comparisons and coefficient checks for Allouche et al. (2012) dataset reanalyses (n = 285). Candidate linear, quadratic, and logarithmic HDR models were compared using AIC. The primary unadjusted, no-covariate model set is presented first, followed by sensitivity analyses adding mean elevation alone and mean elevation plus mean elevation squared."

  markdown_outputs <- c(
    table_1 = write_markdown_table(
      table_1_display,
      file.path(manuscript_deliverable_dirs$extended_data, "extended_data_table_1_macarthur_model_comparisons_and_coefficient_checks.md"),
      "Extended Data Table 1",
      table_1_caption
    ),
    table_2 = write_markdown_table(
      table_2_display,
      file.path(manuscript_deliverable_dirs$extended_data, "extended_data_table_2_allouche_model_comparisons_and_coefficient_checks.md"),
      "Extended Data Table 2",
      table_2_caption
    )
  )

  docx_outputs <- c(table_1 = NA_character_, table_2 = NA_character_)
  if (requireNamespace("officer", quietly = TRUE) && requireNamespace("flextable", quietly = TRUE)) {
    docx_outputs["table_1"] <- write_nature_docx_table(
      table_1_display,
      output_path = file.path(manuscript_deliverable_dirs$extended_data, "extended_data_table_1_macarthur_model_comparisons_and_coefficient_checks.docx"),
      table_title = "Extended Data Table 1",
      table_caption = table_1_caption,
      landscape = FALSE,
      font_size = 8,
      header_size = 8
    )
    docx_outputs["table_2"] <- write_nature_docx_table(
      table_2_display,
      output_path = file.path(manuscript_deliverable_dirs$extended_data, "extended_data_table_2_allouche_model_comparisons_and_coefficient_checks.docx"),
      table_title = "Extended Data Table 2",
      table_caption = table_2_caption,
      landscape = FALSE,
      font_size = 8,
      header_size = 8
    )
  } else {
    warning("officer/flextable are unavailable; wrote Markdown Extended Data table drafts but skipped DOCX table drafts.")
  }

  list(markdown_outputs = markdown_outputs, docx_outputs = docx_outputs)
}


figure1_primary_analysis_ids <- function() {
  c(
    foliage = "Foliage_delta_unadjusted_gaussian",
    topography = "Topography_delta_no_mean_covariates_gaussian"
  )
}

figure1_hypothesis_rows <- function() {
  delta <- manuscript_delta_symbol()
  delta2 <- paste0(delta, manuscript_squared_symbol())

  tibble::tibble(
    display_order = 1:3,
    model_family = c("linear", "quadratic", "logarithmic"),
    formula = c(
      paste0("y ~ ", delta),
      paste0("y ~ ", delta, " + ", delta2),
      paste0("y ~ log(", delta, ")")
    ),
    hypothesis = c(
      "Habitat heterogeneity hypothesis",
      "Area\u2013heterogeneity trade-off hypothesis",
      "Diminishing returns hypothesis"
    ),
    prediction = c(
      "Biodiversity rises at a constant rate",
      "Biodiversity rises, peaks, then declines",
      "Biodiversity rises steeply, then gains taper"
    ),
    application = c(
      "Starting heterogeneity value alone does not guide targeted intervention",
      "Manage heterogeneity toward the biodiversity peak",
      "Increase heterogeneity when starting in homogeneous landscapes; when reducing heterogeneity, start in heterogeneous landscapes"
    )
  )
}

figure1_illustrative_curve_data <- function() {
  x_limits <- c(1, 300)
  x_span <- diff(x_limits)
  scaled_x <- function(x) (x - x_limits[1]) / x_span
  x <- seq(x_limits[1], x_limits[2], length.out = 201)
  functions <- list(
    linear = function(value) {
      t <- scaled_x(value)
      0.33 + 0.48 * t
    },
    logarithmic = function(value) 0.08 + 0.09642723972098036 * log(value),
    quadratic = function(value) {
      t <- scaled_x(value)
      0.28 + 1.42 * t - 1.30 * t^2
    }
  )

  dplyr::bind_rows(lapply(names(functions), function(model_family) {
    tibble::tibble(
      illustrative_family = model_family,
      point_index = seq_along(x),
      plotting_x = x,
      plotting_y = functions[[model_family]](x)
    )
  }))
}

figure1_prediction_confidence_band <- function(candidate_set, plotted_family, level = 0.95) {
  model_levels <- c("linear", "logarithmic", "quadratic")
  plotted_family <- tolower(as.character(plotted_family)[1])
  if (is.na(plotted_family) || !nzchar(plotted_family) || !plotted_family %in% model_levels) {
    stop("plotted_family must be one of: ", paste(model_levels, collapse = ", "), ".")
  }

  model <- candidate_set$models[[plotted_family]]
  if (is.null(model) || !inherits(model, "lm")) {
    stop("No fitted lm object found for plotted family '", plotted_family, "'.")
  }

  pred_grid <- candidate_set$predictions %>%
    dplyr::filter(as.character(model_family) == plotted_family) %>%
    dplyr::arrange(x)
  if (nrow(pred_grid) == 0L) {
    stop("No prediction grid rows found for plotted family '", plotted_family, "'.")
  }

  newdata <- tibble::tibble(
    x = as.numeric(pred_grid$x),
    logx = log(as.numeric(pred_grid$x))
  )
  model_terms <- setdiff(
    all.vars(stats::delete.response(stats::terms(model))),
    names(newdata)
  )
  model_frame <- tryCatch(stats::model.frame(model), error = function(e) NULL)

  for (term in model_terms) {
    if (term %in% names(pred_grid)) {
      newdata[[term]] <- pred_grid[[term]]
    } else if (!is.null(model_frame) && term %in% names(model_frame)) {
      term_values <- model_frame[[term]]
      if (is.factor(term_values)) {
        replacement <- stats::na.omit(term_values)[1]
        newdata[[term]] <- factor(
          rep(as.character(replacement), nrow(newdata)),
          levels = levels(term_values)
        )
      } else if (is.character(term_values)) {
        replacement <- stats::na.omit(term_values)[1]
        newdata[[term]] <- rep(as.character(replacement), nrow(newdata))
      } else {
        newdata[[term]] <- rep(mean(as.numeric(term_values), na.rm = TRUE), nrow(newdata))
      }
    } else {
      stop(
        "Cannot construct prediction newdata for model term '",
        term,
        "' in plotted family '",
        plotted_family,
        "'."
      )
    }
  }

  ci <- as.data.frame(stats::predict(
    model,
    newdata = newdata,
    interval = "confidence",
    level = level
  ))
  stored_fit_check <- all.equal(
    as.numeric(ci$fit),
    as.numeric(pred_grid$yhat),
    tolerance = 0,
    check.attributes = FALSE
  )
  if (!isTRUE(stored_fit_check)) {
    stop(
      "Figure 1 fitted-mean confidence calculation does not reproduce the stored prediction grid: ",
      paste(stored_fit_check, collapse = "; ")
    )
  }

  pred_grid %>%
    dplyr::mutate(
      yhat = as.numeric(ci$fit),
      conf_low = as.numeric(ci$lwr),
      conf_high = as.numeric(ci$upr)
    ) %>%
    dplyr::filter(is.finite(x), is.finite(yhat), is.finite(conf_low), is.finite(conf_high))
}


make_main_text_hdr_figure <- function(
    results = results_manuscript_analysis_hierarchy,
    save = TRUE,
    output_dir = NULL,
    file_stub = "main_text_mean_independent_heterogeneity_reveals_diminishing_return_HDRs",
    width_mm = 180,
	    height_mm = 100,
    dpi = 600,
    save_tiff = TRUE,
    save_pdf = isTRUE(capabilities("cairo"))
) {
  # Figure 1 follows this structure:
  #   a, table of hypotheses, curve shapes, predictions, and applications
	  #   b, linear, logarithmic, and unimodal alternatives for curve shape
	  #   c, primary MacArthur & MacArthur (1961) δ analysis without a mean covariate,
	  #      plotting only the candidate curve ranked first
	  #   d, primary Allouche et al. (2012) δ analysis without a mean covariate,
	  #      plotting only the candidate curve ranked first
  
  required_packages <- c("ggplot2", "dplyr", "tibble", "patchwork", "scales", "grid")
  missing_packages <- required_packages[
    !vapply(
      required_packages,
      function(pkg) suppressWarnings(requireNamespace(pkg, quietly = TRUE)),
      logical(1)
    )
  ]
  if (length(missing_packages) > 0) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
  }
  
  if (is.null(output_dir)) {
    if (exists("manuscript_deliverable_dirs", inherits = TRUE)) {
      output_dir <- manuscript_deliverable_dirs$main_text
    } else {
      output_dir <- getwd()
    }
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  primary_analysis_ids <- figure1_primary_analysis_ids()
  figure1_topography_analysis_id <- unname(primary_analysis_ids[["topography"]])
  figure1_foliage_analysis_id <- unname(primary_analysis_ids[["foliage"]])
  
  topography_set <- results$topography_delta_models$sets[[figure1_topography_analysis_id]]
  foliage_set <- results$foliage_delta_models$sets[[figure1_foliage_analysis_id]]
  
  if (is.null(topography_set)) {
    stop("Could not find ", figure1_topography_analysis_id, " in results$topography_delta_models$sets.")
  }
  if (is.null(foliage_set)) {
    stop("Could not find ", figure1_foliage_analysis_id, " in results$foliage_delta_models$sets.")
  }
  
  model_levels <- c("linear", "logarithmic", "quadratic")
  model_labels <- c(
    linear = "Linear",
    logarithmic = "Logarithmic",
    quadratic = "Quadratic"
  )
  
  # Figure 1 uses a local accessible palette and redundant linetypes. Keeping
  # these mappings local preserves the shared HDR defaults for plotting scales
  # used by other consumers.
  model_cols <- c(
    linear = "#0072B2",
    logarithmic = "#7B3294",
    quadratic = "#D55E00"
  )
  model_ltys <- c(
    linear = "longdash",
    logarithmic = "solid",
    quadratic = "dotdash"
  )
  
  model_factor <- function(x) {
    factor(as.character(x), levels = model_levels)
  }
  
  range_finite <- function(x, fallback = c(0, 1)) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(fallback)
    r <- range(x, na.rm = TRUE)
    if (!all(is.finite(r)) || diff(r) == 0) {
      mid <- ifelse(length(x) > 0L && is.finite(x[1]), x[1], mean(fallback))
      return(c(mid - 0.5, mid + 0.5))
    }
    r
  }
  
  get_covariates <- function(candidate_set) {
    covariates <- candidate_set$summary$covariates[1]
    if (is.na(covariates) || covariates %in% c("", "none")) {
      return(character(0))
    }
    covariates <- strsplit(covariates, ";", fixed = TRUE)[[1]]
    trimws(covariates[nzchar(covariates)])
  }
  
  get_sample_size <- function(candidate_set, pts = NULL) {
    w <- candidate_set$summary
    
    if (!is.null(w) && "n" %in% names(w)) {
      n_values <- unique(as.numeric(w$n[is.finite(as.numeric(w$n))]))
      if (length(n_values) > 0L) {
        return(as.integer(n_values[1]))
      }
    }
    
    if (!is.null(pts)) {
      return(as.integer(nrow(pts)))
    }
    
    NA_integer_
  }
  
  get_adjusted_points <- function(candidate_set) {
    dat <- candidate_set$data
    covariates <- get_covariates(candidate_set)
    
    plot_df <- data.frame(
      y = as.numeric(dat$y),
      x = as.numeric(dat$x)
    )
    
    if (length(covariates) > 0) {
      for (cv in covariates) {
        if (!cv %in% names(dat)) {
          stop("Covariate '", cv, "' is listed in candidate_set$summary but absent from candidate_set$data.")
        }
        plot_df[[cv]] <- as.numeric(dat[[cv]])
      }
    }
    
    finite_cols <- c("y", "x", covariates)
    keep <- apply(plot_df[, finite_cols, drop = FALSE], 1, function(z) all(is.finite(z)))
    plot_df <- plot_df[keep, , drop = FALSE]
    
    adjusted_response <- plot_df$y
    adjustment_label <- "raw response"
    
    if (length(covariates) > 0) {
      cov_formula <- stats::as.formula(
        paste("y ~", paste(covariates, collapse = " + "))
      )
      cov_fit <- stats::lm(cov_formula, data = plot_df)
      
      adjusted_response <- stats::residuals(cov_fit) +
        mean(stats::predict(cov_fit), na.rm = TRUE)
      
      adjustment_label <- paste0(
        "adjusted for ",
        format_covariate_label(paste(covariates, collapse = ";"), candidate_set$label)
      )
    }
    
    tibble::tibble(
      x = plot_df$x,
      y_raw = plot_df$y,
      adjusted_response = as.numeric(adjusted_response),
      adjustment_label = adjustment_label,
      n = nrow(plot_df)
    )
  }
  
  get_predictions <- function(candidate_set) {
    pred <- candidate_set$predictions
    
    if (!all(c("x", "yhat", "model_family") %in% names(pred))) {
      stop("candidate_set$predictions must contain x, yhat, and model_family.")
    }
    
    pred |>
      dplyr::filter(is.finite(x), is.finite(yhat)) |>
      dplyr::mutate(model_family = model_factor(model_family))
  }
  
  get_weight_inset_data <- function(candidate_set) {
    w <- candidate_set$summary
    
    if (!"Akaike_weight" %in% names(w)) {
      stop("candidate_set$summary must contain Akaike_weight.")
    }
    if (!"selection_criterion" %in% names(w)) {
      w$selection_criterion <- "IC"
    }
    if (!"delta_IC" %in% names(w)) {
      if ("IC" %in% names(w)) {
        w <- w |>
          dplyr::mutate(delta_IC = IC - min(IC, na.rm = TRUE))
      } else {
        stop("candidate_set$summary must contain delta_IC or IC.")
      }
    }
    
	    w <- w |>
	      dplyr::mutate(model_family = model_factor(model_family)) |>
	      dplyr::filter(!is.na(model_family), is.finite(Akaike_weight), is.finite(delta_IC)) |>
	      dplyr::arrange(delta_IC, model_family)
    
    if (nrow(w) == 0L) {
      stop("candidate_set$summary contains no finite IC weights for the three candidate models.")
    }
    
    criterion <- unique(as.character(w$selection_criterion))
    criterion <- criterion[!is.na(criterion) & nzchar(criterion)]
    criterion <- if (length(criterion) > 0L) criterion[1] else "IC"
    
    delta_header_label <- dplyr::case_when(
      criterion == "AIC" ~ paste0(manuscript_Delta_symbol(), "AIC"),
      criterion == "AICc" ~ paste0(manuscript_Delta_symbol(), "AICc"),
      TRUE ~ paste0(manuscript_Delta_symbol(), criterion)
    )
    
    weight_header_label <- dplyr::case_when(
      criterion == "AIC" ~ "AIC\nweights",
      criterion == "AICc" ~ "AICc\nweights",
      TRUE ~ paste0(criterion, "\nweights")
    )
    
    full_labels <- c(
      linear = "Linear",
      logarithmic = "Logarithmic",
      quadratic = "Quadratic"
    )
    
	    best_delta <- min(w$delta_IC, na.rm = TRUE)
    selected_family <- figure1_select_plotted_family(w)
    
	    rows <- w |>
	      dplyr::mutate(
	        model_label = full_labels[as.character(model_family)],
	        delta_ic_label = sprintf("%.2f", delta_IC),
	        weight_label = sprintf("%.2f", Akaike_weight),
	        is_best = abs(delta_IC - best_delta) < sqrt(.Machine$double.eps),
	        fontface = "plain"
	      )
    
	    list(
	      delta_header_label = delta_header_label,
	      weight_header_label = weight_header_label,
	      best_model_family = selected_family,
	      rows = rows
	    )
  }
  
  main_text_curve_legend_guides <- function() {
    # The main text figure uses compact inset tables in panels b–d, so no
    # separate colour/linetype legend is needed. Keeping guides disabled avoids
    # mismatches in ggplot guide key lengths after plotting only a subset of
    # candidate curves in each panel.
    ggplot2::guides(colour = "none", fill = "none", linetype = "none")
  }
  
  get_prediction_confidence_band <- function(candidate_set, plotted_family, level = 0.95) {
    figure1_prediction_confidence_band(candidate_set, plotted_family, level = level) %>%
      dplyr::mutate(model_family = model_factor(model_family))
  }
  
  nature_theme <- nature_manuscript_theme
  make_figure1_hypothesis_table_panel <- function() {
    header <- tibble::tibble(
      x = c(-0.080, 0.160, 0.321, 0.534),
	      y = 0.925,
      label = c("Hypothesis", "Model", "Prediction", "Application")
    )
    hypothesis_rows <- figure1_hypothesis_rows()
    wrap_table_cells <- function(values, width) {
      vapply(
        values,
        function(value) paste(strwrap(value, width = width), collapse = "\n"),
        character(1)
      )
    }
    body <- hypothesis_rows %>%
      dplyr::transmute(
        y = c(0.700, 0.415, 0.142),
        hypothesis = sub(" hypothesis$", "\nhypothesis", hypothesis),
        curve_shape = paste0(
          c("Linear", "Quadratic", "Logarithmic"),
          " model\n",
          formula
        ),
        prediction = wrap_table_cells(prediction, width = 28),
        application = dplyr::case_when(display_order == 1L ~ sub(" targeted intervention$", "\ntargeted intervention", application), display_order == 2L ~ application, display_order == 3L ~ sub("; ", ";\n", application, fixed = TRUE))
      )
    row_background <- tibble::tibble(
      xmin = -0.090,
      xmax = 1.060,
	      ymin = c(0.840, 0.557, 0.273, -0.010),
	      ymax = c(1.010, 0.840, 0.557, 0.273),
	      fill = c("white", "#EAF3F8", "#FBEDE7", "#F2EAF5")
    )
    row_rules <- tibble::tibble(
      x = -0.090,
      xend = 1.060,
	      y = c(0.840, 0.557, 0.273, -0.010),
	      yend = c(0.840, 0.557, 0.273, -0.010)
    )
    col_rules <- tibble::tibble(
      x = c(0.151, 0.311, 0.524),
      xend = c(0.151, 0.311, 0.524),
	      y = -0.010,
	      yend = 1.010
    )
    
    ggplot2::ggplot() +
      ggplot2::geom_rect(
        data = row_background,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
        inherit.aes = FALSE,
        colour = NA
      ) +
	      ggplot2::geom_rect(
	        data = tibble::tibble(xmin = -0.090, xmax = 1.060, ymin = -0.010, ymax = 1.010),
	        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
	        inherit.aes = FALSE,
	        fill = NA,
	        colour = "grey70",
	        linewidth = 0.18
      ) +
      ggplot2::geom_segment(
        data = row_rules,
	        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
	        inherit.aes = FALSE,
	        colour = "grey82",
	        linewidth = 0.16
      ) +
      ggplot2::geom_segment(
        data = col_rules,
	        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
	        inherit.aes = FALSE,
	        colour = "grey86",
	        linewidth = 0.14
      ) +
      ggplot2::geom_text(
        data = header,
        ggplot2::aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 0.5,
        fontface = "bold",
	        size = 2.88,
	        lineheight = 0.95,
	        colour = "grey10"
      ) +
      ggplot2::geom_text(
        data = body,
        ggplot2::aes(x = -0.080, y = y, label = hypothesis),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 0.5,
        fontface = "plain",
        size = 2.53,
	        lineheight = 0.95,
	        colour = "grey10"
      ) +
      ggplot2::geom_text(
        data = body,
        ggplot2::aes(x = 0.160, y = y, label = curve_shape),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 0.5,
        fontface = "plain",
        size = 2.53,
	        lineheight = 0.95,
	        colour = "grey10"
      ) +
      ggplot2::geom_text(
        data = body,
        ggplot2::aes(x = 0.321, y = y, label = prediction),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 0.5,
        fontface = "plain",
        size = 2.53,
	        lineheight = 0.95,
	        colour = "grey10"
      ) +
      ggplot2::geom_text(
        data = body,
        ggplot2::aes(x = 0.534, y = y, label = application),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 0.5,
        fontface = "plain",
        size = 2.53,
	        lineheight = 0.90,
	        colour = "grey10"
      ) +
      ggplot2::scale_fill_identity() +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 5),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA)
	      )
	  }
	  
	  figure1_empirical_title <- function(dataset_key) {
	    switch(
	      dataset_key,
		      foliage = "MacArthur & MacArthur dataset",
		      topography = "Allouche et al. dataset",
	      manuscript_dataset_phrase(dataset_key)
	    )
	  }
  
  figure1_panel_title_size <- 8.8
  figure1_panel_title_theme <- ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = figure1_panel_title_size,
      lineheight = 0.95,
      margin = ggplot2::margin(b = 3)
    )
  )
	  
	  make_concept_panel <- function() {
    # Schematic functions used only for plotting the literal candidate families on
    # one shared positive predictor. These deterministic coefficients are
    # illustrative: they are not fitted to or derived from any empirical model,
    # prediction, or dataset.
    schematic_x_limits <- c(1, 300)
    schematic_x_span <- diff(schematic_x_limits)
    schematic_x_from_fraction <- function(fraction) {
      schematic_x_limits[1] + fraction * schematic_x_span
    }
    concept <- figure1_illustrative_curve_data() %>%
      dplyr::transmute(
        x = plotting_x,
        y = plotting_y,
        model_family = model_factor(illustrative_family)
      )
    
    # Map the inset's normalized horizontal coordinates onto the literal
    # x = 1...300 domain.
    inset_x0 <- schematic_x_from_fraction(0.420)
    inset_x1 <- schematic_x_from_fraction(0.985)
    inset_header_x <- schematic_x_from_fraction(0.458)
    inset_line_x0 <- schematic_x_from_fraction(0.458)
    inset_line_x1 <- schematic_x_from_fraction(0.625)
    inset_label_x <- schematic_x_from_fraction(0.650)
    inset_y0 <- 0.070
    inset_y1 <- 0.350
    inset_header_y <- 0.315
    inset_rule_y <- 0.278
    concept_inset_rows <- tibble::tibble(
      model_family = model_factor(c("linear", "logarithmic", "quadratic")),
      model_label = c("Linear", "Logarithmic", "Quadratic"),
      y = c(0.230, 0.165, 0.100)
    )
    
    ggplot2::ggplot() +
	      ggplot2::geom_line(
	        data = dplyr::filter(concept, as.character(model_family) == "linear"),
	        ggplot2::aes(x, y, colour = model_family, linetype = model_family),
	        linewidth = 1.05,
	        show.legend = FALSE
	      ) +
	      ggplot2::geom_line(
	        data = dplyr::filter(concept, as.character(model_family) == "logarithmic"),
	        ggplot2::aes(x, y, colour = model_family, linetype = model_family),
	        linewidth = 1.05,
	        show.legend = FALSE
	      ) +
	      ggplot2::geom_line(
	        data = dplyr::filter(concept, as.character(model_family) == "quadratic"),
	        ggplot2::aes(x, y, colour = model_family, linetype = model_family),
        linewidth = 1.05,
        show.legend = FALSE
      ) +
      ggplot2::geom_rect(
        data = tibble::tibble(xmin = inset_x0, xmax = inset_x1, ymin = inset_y0, ymax = inset_y1),
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = "white",
        colour = "grey65",
        alpha = 0.94,
        linewidth = 0.18
      ) +
      ggplot2::geom_segment(
        data = tibble::tibble(x = inset_x0, xend = inset_x1, y = inset_rule_y, yend = inset_rule_y),
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        inherit.aes = FALSE,
        colour = "grey78",
        linewidth = 0.16
      ) +
      ggplot2::geom_text(
        data = tibble::tibble(x = inset_header_x, y = inset_header_y, label = "Model"),
        ggplot2::aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 0.5,
        size = 2.05,
        fontface = "bold",
        colour = "grey10"
      ) +
      ggplot2::geom_segment(
        data = concept_inset_rows,
        ggplot2::aes(x = inset_line_x0, xend = inset_line_x1, y = y, yend = y, colour = model_family, linetype = model_family),
        inherit.aes = FALSE,
        linewidth = 0.82,
        lineend = "butt",
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = concept_inset_rows,
        ggplot2::aes(x = inset_label_x, y = y, label = model_label),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 0.5,
        size = 1.95,
        colour = "grey10"
      ) +
      ggplot2::scale_colour_manual(
        name = "Candidate curve",
        values = model_cols,
        breaks = model_levels,
        labels = model_labels
      ) +
      ggplot2::scale_linetype_manual(
        name = "Candidate curve",
        values = model_ltys,
        breaks = model_levels,
        labels = model_labels
      ) +
      main_text_curve_legend_guides() +
      ggplot2::scale_x_continuous(
        limits = schematic_x_limits,
        breaks = c(1, 150.5, 300),
        labels = c("low", "intermediate", "high"),
        expand = ggplot2::expansion(mult = c(0.02, 0.02))
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        breaks = c(0, 1),
        labels = c("less", "more"),
        minor_breaks = NULL,
        expand = ggplot2::expansion(mult = c(0.00, 0.03))
      ) +
      ggplot2::labs(
        title = "Competing hypotheses",
        x = paste0("Heterogeneity gradient (", manuscript_delta_symbol(), ")"),
        y = "Biodiversity"
      ) +
      nature_theme() +
      figure1_panel_title_theme
  }
  
	  make_empirical_panel <- function(candidate_set, title, subtitle = NULL, x_lab, y_lab, plotted_family = NULL) {
    pts <- get_adjusted_points(candidate_set)
    pred <- get_predictions(candidate_set)
    weight_inset <- get_weight_inset_data(candidate_set)
    is_foliage <- manuscript_dataset_key(title) == "foliage" || grepl("Foliage", title, ignore.case = TRUE)
    
    # Sample sizes are reported in the figure caption rather than as panel
    # subtitles, keeping the main text figure visually compact.
    if (is.null(subtitle) || length(subtitle) == 0L) {
      subtitle <- NULL
    } else {
      subtitle <- as.character(subtitle[1])
      if (is.na(subtitle) || !nzchar(trimws(subtitle))) {
        subtitle <- NULL
      }
    }
    
    best_plotted_family <- weight_inset$best_model_family
    if (
      is.null(best_plotted_family) ||
        length(best_plotted_family) == 0L ||
        is.na(best_plotted_family[1]) ||
        !nzchar(as.character(best_plotted_family[1]))
    ) {
      best_plotted_family <- as.character(weight_inset$rows$model_family[which.min(weight_inset$rows$delta_IC)])
    }
    best_plotted_family <- as.character(best_plotted_family[1])
    
    if (is.null(plotted_family) || length(plotted_family) == 0L) {
      plotted_family <- best_plotted_family
    } else {
      plotted_family <- as.character(plotted_family)[1]
      if (is.na(plotted_family) || !nzchar(plotted_family)) {
        plotted_family <- best_plotted_family
      }
    }
    plotted_family <- tolower(plotted_family)
    if (!plotted_family %in% model_levels) {
      stop("plotted_family must be one of: ", paste(model_levels, collapse = ", "), ".")
    }
    
    pred_plotted <- pred |>
      dplyr::filter(as.character(model_family) == plotted_family)
    ci_plotted <- get_prediction_confidence_band(candidate_set, plotted_family)
    
    # The empirical y scale is trained exclusively from finite displayed data:
    # observations, the plotted fitted curve, and its confidence limits. Inset
    # geometry is calculated only after these final limits have been fixed.
    raw_xr <- range_finite(c(pts$x, pred$x))
    empirical_y_values <- c(
      pts$adjusted_response,
      pred_plotted$yhat,
      ci_plotted$conf_low,
      ci_plotted$conf_high
    )
    empirical_y_values <- empirical_y_values[is.finite(empirical_y_values)]
    if (length(empirical_y_values) == 0L) {
      stop("No finite empirical y values are available for Figure 1.")
    }
    raw_yr <- range(empirical_y_values)
    raw_xd <- diff(raw_xr)
    raw_yd <- diff(raw_yr)
    if (!is.finite(raw_yd) || raw_yd <= 0) {
      stop("Empirical Figure 1 y range must have a finite positive span.")
    }
    
    x_pad_left <- 0.050
    x_pad_right <- 0.075
    y_padding_fraction <- 0.040
    
    xr <- raw_xr + c(-x_pad_left, x_pad_right) * raw_xd
    yr <- raw_yr + c(-y_padding_fraction, y_padding_fraction) * raw_yd
    # If symmetric padding would create a negative axis for an otherwise
    # nonnegative response, clamp only that padded extension at zero.
    if (raw_yr[1] >= 0 && yr[1] < 0) {
      yr[1] <- 0
    }
    xd <- diff(xr)
    yd <- diff(yr)
    axis_y_breaks <- scales::breaks_extended()(yr)
    axis_y_breaks <- axis_y_breaks[
      is.finite(axis_y_breaks) & axis_y_breaks >= yr[1] & axis_y_breaks <= yr[2]
    ]
    
    # Be robust to either the current edited inset helper or the older helper.
    # The current helper already supplies delta_ic_label and two column headers;
    # the fallback below computes them from candidate_set$summary if needed.
    if (!"delta_ic_label" %in% names(weight_inset$rows)) {
      w_delta <- candidate_set$summary
      if (!"delta_IC" %in% names(w_delta)) {
        if ("IC" %in% names(w_delta)) {
          w_delta <- w_delta |>
            dplyr::mutate(delta_IC = IC - min(IC, na.rm = TRUE))
        } else {
          w_delta$delta_IC <- NA_real_
        }
      }
      w_delta <- w_delta |>
        dplyr::mutate(
          model_family = model_factor(model_family),
          delta_IC = as.numeric(delta_IC),
          delta_ic_label = ifelse(is.finite(delta_IC), sprintf("%.2f", delta_IC), "NA")
        ) |>
        dplyr::select(model_family, delta_IC, delta_ic_label)
      
      weight_inset$rows <- weight_inset$rows |>
        dplyr::left_join(w_delta, by = "model_family")
    }
    
    criterion <- "IC"
    if ("selection_criterion" %in% names(candidate_set$summary)) {
      criterion_values <- unique(as.character(candidate_set$summary$selection_criterion))
      criterion_values <- criterion_values[!is.na(criterion_values) & nzchar(criterion_values)]
      if (length(criterion_values) > 0L) criterion <- criterion_values[1]
    }
    
    if (!"delta_header_label" %in% names(weight_inset)) {
      weight_inset$delta_header_label <- dplyr::case_when(
        criterion == "AIC" ~ paste0(manuscript_Delta_symbol(), "AIC"),
        criterion == "AICc" ~ paste0(manuscript_Delta_symbol(), "AICc"),
        TRUE ~ paste0(manuscript_Delta_symbol(), criterion)
      )
    }
    if (!"weight_header_label" %in% names(weight_inset)) {
      weight_inset$weight_header_label <- dplyr::case_when(
        criterion == "AIC" ~ "AIC\nweights",
        criterion == "AICc" ~ "AICc\nweights",
        TRUE ~ paste0(criterion, "\nweights")
      )
    }
    weight_inset$weight_header_label <- gsub(
      "weight$",
      "weights",
      as.character(weight_inset$weight_header_label)
    )
    
    # Compact information criterion inset. Its normalized panel coordinates are
    # converted only after the limits determined by the data above are final, and the
    # explicit scale limits below prevent any inset layer from training y.
    inset_left_fraction <- if (is_foliage) 0.120 else 0.180
    inset_x0 <- xr[1] + inset_left_fraction * xd
    inset_x1 <- xr[2] - 0.012 * xd
    inset_bottom_fraction <- 0.015
    inset_height_fraction <- 0.300
    inset_y0 <- yr[1] + inset_bottom_fraction * yd
    inset_y1 <- inset_y0 + inset_height_fraction * yd
    inset_width <- inset_x1 - inset_x0
    
    # Give the header its own taller row. This keeps both lines of the IC and weight
    # header clear of the top border of the inset box.
    inset_body_rows <- nrow(weight_inset$rows)
    inset_header_height_factor <- 1.75
    inset_unit_height <- (inset_y1 - inset_y0) /
      (inset_header_height_factor + inset_body_rows)
    inset_header_height <- inset_header_height_factor * inset_unit_height
    inset_body_row_height <- inset_unit_height
    
    # Explicit column geometry. Header labels are centered over their own columns,
    # while numeric values are centred within their columns.
    line_col_left_fraction <- 0.025
    line_col_right_fraction <- 0.165
    model_col_left_fraction <- 0.190
    delta_col_left_fraction <- 0.520
    delta_col_right_fraction <- 0.720
    weight_col_left_fraction <- 0.740
    weight_col_right_fraction <- 0.985
    
    line_col_x0 <- inset_x0 + line_col_left_fraction * inset_width
    line_col_x1 <- inset_x0 + line_col_right_fraction * inset_width
    model_col_x0 <- inset_x0 + model_col_left_fraction * inset_width
    model_col_center <- inset_x0 + mean(c(model_col_left_fraction, delta_col_left_fraction - 0.035)) * inset_width
    delta_col_center <- inset_x0 + mean(c(delta_col_left_fraction, delta_col_right_fraction)) * inset_width
    weight_col_center <- inset_x0 + mean(c(weight_col_left_fraction, weight_col_right_fraction)) * inset_width
    weight_col_x1 <- inset_x0 + weight_col_right_fraction * inset_width
    
    inset_background <- tibble::tibble(
      xmin = inset_x0,
      xmax = inset_x1,
      ymin = inset_y0,
      ymax = inset_y1
    )
    
    inset_header <- tibble::tibble(
      x = c(model_col_center, delta_col_center, weight_col_center),
      y = inset_y1 - 0.5 * inset_header_height,
      label = c("Model", weight_inset$delta_header_label, weight_inset$weight_header_label)
    )
    
    inset_header_rule <- tibble::tibble(
      x = inset_x0,
      xend = inset_x1,
      y = inset_y1 - inset_header_height,
      yend = inset_y1 - inset_header_height
    )
    
    inset_rows <- weight_inset$rows |>
      dplyr::mutate(
        inset_row = dplyr::row_number(),
        y = inset_y1 - inset_header_height -
          (inset_row - 0.5) * inset_body_row_height,
        ymin = y - 0.5 * inset_body_row_height,
        ymax = y + 0.5 * inset_body_row_height,
        xmin = inset_x0,
        xmax = inset_x1,
        line_x0 = .env$line_col_x0,
        line_x1 = .env$line_col_x1,
        model_x = .env$model_col_x0,
        delta_x = .env$delta_col_center,
        weight_x = .env$weight_col_center,
        fontface = "plain"
      )
    
    point_size <- ifelse(nrow(pts) <= 20, 1.65, 0.72)
    point_alpha <- ifelse(nrow(pts) <= 20, 0.72, 0.30)
    
    ggplot2::ggplot() +
      ggplot2::geom_point(
        data = pts,
        ggplot2::aes(x, adjusted_response),
        size = point_size,
        alpha = point_alpha,
        colour = "grey20"
      ) +
      ggplot2::geom_ribbon(
        data = ci_plotted,
        ggplot2::aes(x = x, ymin = conf_low, ymax = conf_high, fill = model_family),
        inherit.aes = FALSE,
        alpha = 0.18,
        colour = NA,
        show.legend = FALSE
      ) +
      ggplot2::geom_line(
        data = pred_plotted,
        ggplot2::aes(x, yhat, colour = model_family, linetype = model_family),
        linewidth = 1.05,
        show.legend = FALSE
      ) +
      ggplot2::geom_rect(
        data = inset_background,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = "white",
        colour = NA,
        show.legend = FALSE
      ) +
      ggplot2::geom_rect(
        data = inset_background,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = NA,
        colour = "black",
        linewidth = 0.30,
        show.legend = FALSE
      ) +
      ggplot2::geom_segment(
        data = inset_header_rule,
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        inherit.aes = FALSE,
        colour = "black",
        linewidth = 0.24,
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = inset_header,
        ggplot2::aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        hjust = 0.5,
        vjust = 0.5,
        lineheight = 0.86,
        size = 2.53,
        fontface = "bold",
        show.legend = FALSE
      ) +
      ggplot2::geom_segment(
        data = inset_rows,
        ggplot2::aes(
          x = line_x0,
          xend = line_x1,
          y = y,
          yend = y,
          colour = model_family,
          linetype = model_family
        ),
        inherit.aes = FALSE,
        linewidth = 0.70,
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = inset_rows,
        ggplot2::aes(x = model_x, y = y, label = model_label),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 0.5,
        size = 2.29,
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = inset_rows,
        ggplot2::aes(x = delta_x, y = y, label = delta_ic_label),
        inherit.aes = FALSE,
        hjust = 0.5,
        vjust = 0.5,
        size = 2.29,
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = inset_rows,
        ggplot2::aes(x = weight_x, y = y, label = weight_label),
        inherit.aes = FALSE,
        hjust = 0.5,
        vjust = 0.5,
        size = 2.29,
        show.legend = FALSE
      ) +
      ggplot2::scale_colour_manual(
        name = "Candidate curve",
        values = model_cols,
        breaks = model_levels,
        labels = model_labels
      ) +
      ggplot2::scale_fill_manual(
        name = "Candidate curve",
        values = model_cols,
        breaks = model_levels,
        labels = model_labels
      ) +
      ggplot2::scale_linetype_manual(
        name = "Candidate curve",
        values = model_ltys,
        breaks = model_levels,
        labels = model_labels
      ) +
      main_text_curve_legend_guides() +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0))) +
      ggplot2::scale_y_continuous(
        limits = yr,
        breaks = axis_y_breaks,
        expand = ggplot2::expansion(mult = c(0, 0))
      ) +
      ggplot2::coord_cartesian(xlim = xr, ylim = yr, clip = "on") +
	      ggplot2::labs(
	        title = title,
	        subtitle = subtitle,
	        x = x_lab,
	        y = y_lab
	      ) +
	      nature_theme() +
	      figure1_panel_title_theme
	  }
  
  p_a <- make_figure1_hypothesis_table_panel() +
    ggplot2::labs(tag = "a")
  
  p_b <- make_concept_panel() +
    ggplot2::labs(tag = "b")
  
	  p_c <- make_empirical_panel(
	    foliage_set,
	    title = figure1_empirical_title("foliage"),
	    subtitle = NULL,
	    x_lab = paste0("Foliage height heterogeneity, ", manuscript_delta_symbol(), " (ft)"),
	    y_lab = manuscript_dataset_response("foliage")
  ) +
    ggplot2::labs(tag = "c")
  
	  p_d <- make_empirical_panel(
	    topography_set,
	    title = figure1_empirical_title("topography"),
	    subtitle = NULL,
	    x_lab = paste0("Topographic heterogeneity, ", manuscript_delta_symbol(), " (m)"),
	    y_lab = "Breeding bird species richness"
  ) +
    ggplot2::labs(tag = "d")
  
  # The empirical inset tables retain the comparison of all three models, while
  # the plotted empirical curves show only the candidate model ranked first.
  # Panel b uses a compact inset table containing only model identities, avoiding
  # a redundant bottom legend and giving the axes more room.
  plot_row <- patchwork::wrap_plots(
    p_b,
    p_c,
    p_d,
    nrow = 1,
    widths = c(0.90, 1.20, 1.40)
  )
  
	  fig <- (p_a / plot_row) +
	    patchwork::plot_layout(heights = c(0.50, 1.00))
  
  fig <- fig &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 11),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 4)),
      legend.position = "none"
    )
  
  saved_files <- character(0)
  
  # Render the patchwork figure explicitly. Passing the patchwork object directly
  # to ggsave() can leave nearly empty graphics files if patchwork or guide
  # rendering fails after opening a graphics device. Writing first to a temporary
  # file prevents a failed export from replacing an existing usable figure with a
  # blank file.
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  
  draw_patchwork_figure <- function() {
    grid::grid.newpage()
    if ("patchworkGrob" %in% getNamespaceExports("patchwork")) {
      patchwork_grob <- tryCatch(
        getExportedValue("patchwork", "patchworkGrob")(fig),
        error = function(e) NULL
      )
      if (!is.null(patchwork_grob)) {
        grid::grid.draw(patchwork_grob)
        return(invisible(TRUE))
      }
    }
    print(fig)
    invisible(TRUE)
  }
  
  open_bitmap_device <- function(path, format) {
    args <- list(
      filename = path,
      width = width_in,
      height = height_in,
      units = "in",
      res = dpi,
      bg = "white"
    )
    if (format == "tiff") {
      args$compression <- "lzw"
    }
    device_fun <- switch(
      format,
      png = grDevices::png,
      tiff = grDevices::tiff,
      stop("Unsupported bitmap graphics format: ", format, call. = FALSE)
    )
    if (isTRUE(capabilities("cairo"))) {
      tryCatch(
        do.call(device_fun, c(args, list(type = "cairo"))),
        error = function(e) do.call(device_fun, args)
      )
    } else {
      do.call(device_fun, args)
    }
  }
  
  open_pdf_device <- function(path) {
    if (isTRUE(capabilities("cairo"))) {
      grDevices::cairo_pdf(
        filename = path,
        width = width_in,
        height = height_in,
        bg = "white"
      )
    } else {
      grDevices::pdf(
        file = path,
        width = width_in,
        height = height_in,
        bg = "white"
      )
    }
  }
  
  render_figure_file <- function(path, format, required = FALSE) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    format <- tolower(as.character(format)[1])
    extension <- tools::file_ext(path)
    if (!nzchar(extension)) extension <- format
    tmp_path <- tempfile(
      pattern = paste0(tools::file_path_sans_ext(basename(path)), "_"),
      tmpdir = dirname(path),
      fileext = paste0(".", extension)
    )
    
    device_opened <- FALSE
    export_error <- NULL
    ok <- tryCatch(
      {
        if (format %in% c("png", "tiff")) {
          open_bitmap_device(tmp_path, format)
        } else if (format == "pdf") {
          open_pdf_device(tmp_path)
        } else {
          stop("Unsupported graphics format: ", format, call. = FALSE)
        }
        device_opened <- TRUE
        draw_patchwork_figure()
        grDevices::dev.off()
        device_opened <- FALSE
        
        if (!file.exists(tmp_path)) {
          stop("graphics device closed without creating the output file")
        }
        size_bytes <- file.info(tmp_path)$size
        minimum_size <- switch(format, png = 50000, tiff = 50000, pdf = 5000, 1)
        if (!is.finite(size_bytes) || size_bytes < minimum_size) {
          stop(
            "graphics export produced an unexpectedly small file (",
            size_bytes,
            " bytes), suggesting a blank or incomplete render"
          )
        }
        
        if (file.exists(path)) unlink(path)
        if (!file.rename(tmp_path, path)) {
          if (!file.copy(tmp_path, path, overwrite = TRUE)) {
            stop("could not move temporary graphics file to final output path")
          }
          unlink(tmp_path)
        }
        TRUE
      },
      error = function(e) {
        export_error <<- conditionMessage(e)
        FALSE
      },
      finally = {
        if (isTRUE(device_opened) && grDevices::dev.cur() > 1) {
          try(grDevices::dev.off(), silent = TRUE)
        }
        if (file.exists(tmp_path)) {
          try(unlink(tmp_path), silent = TRUE)
        }
      }
    )
    
    if (ok) return(path)
    
    message_text <- paste0("Graphics export failed for ", basename(path), ": ", export_error)
    if (isTRUE(required)) {
      stop(message_text, call. = FALSE)
    }
    message("Skipping optional graphics export for ", basename(path), ": ", export_error)
    NA_character_
  }
  
  if (save) {
    png_path <- file.path(output_dir, paste0(file_stub, ".png"))
    png_written <- render_figure_file(png_path, "png", required = TRUE)
    if (!is.na(png_written)) saved_files <- c(saved_files, png = png_written)
    
    if (isTRUE(save_tiff)) {
      tif_path <- file.path(output_dir, paste0(file_stub, ".tif"))
      tif_written <- render_figure_file(tif_path, "tiff", required = FALSE)
      if (!is.na(tif_written)) saved_files <- c(saved_files, tiff = tif_written)
    }
    
    if (isTRUE(save_pdf)) {
      pdf_path <- file.path(output_dir, paste0(file_stub, ".pdf"))
      pdf_written <- render_figure_file(pdf_path, "pdf", required = FALSE)
      if (!is.na(pdf_written)) saved_files <- c(saved_files, pdf = pdf_written)
    }
  }
  
  attr(fig, "saved_files") <- saved_files
  fig
}

join_analysis_registry_metadata <- function(tbl, registry) {
  joined <- tbl %>%
    left_join(registry, by = "analysis", suffix = c("", "_registry"))

  joined <- manuscript_ensure_column(joined, "dataset", NA_character_)
  joined <- manuscript_ensure_column(joined, "dataset_registry", NA_character_)
  joined <- manuscript_ensure_column(joined, "dataset_short", NA_character_)
  joined <- manuscript_ensure_column(joined, "dataset_short_registry", NA_character_)
  joined <- manuscript_ensure_column(joined, "dataset_phrase", NA_character_)
  joined <- manuscript_ensure_column(joined, "dataset_phrase_registry", NA_character_)
  joined <- manuscript_ensure_column(joined, "source_study", NA_character_)
  joined <- manuscript_ensure_column(joined, "source_study_registry", NA_character_)
  joined <- manuscript_ensure_column(joined, "heterogeneity_metric", NA_character_)
  joined <- manuscript_ensure_column(joined, "predictor", NA_character_)
  joined <- manuscript_ensure_column(joined, "covariate_scenario", NA_character_)
  joined <- manuscript_ensure_column(joined, "covariates", NA_character_)
  joined <- manuscript_ensure_column(joined, "analysis_tier", NA_character_)

  joined %>%
    mutate(
      dataset = dplyr::coalesce(dataset, dataset_registry, "unknown dataset"),
      dataset_short = dplyr::coalesce(dataset_short, dataset_short_registry, manuscript_dataset_short(dataset)),
      dataset_phrase = dplyr::coalesce(dataset_phrase, dataset_phrase_registry, manuscript_dataset_phrase(dataset)),
      source_study = dplyr::coalesce(source_study, source_study_registry, manuscript_dataset_source(dataset)),
      heterogeneity_metric = dplyr::coalesce(heterogeneity_metric, predictor),
      covariate_scenario = dplyr::coalesce(covariate_scenario, covariates, "none"),
      analysis_label = paste(
        dataset_short,
        format_metric_label(heterogeneity_metric),
        format_covariate_label(covariate_scenario),
        sep = " | "
      ),
      analysis_label = gsub("( \\| )+$", "", analysis_label)
    )
}

make_delta_model_observed_points <- function(candidate_sets, registry_tbl) {
  observed <- purrr::imap_dfr(candidate_sets, function(candidate_set, analysis_id) {
    if (is.null(candidate_set) || is.null(candidate_set$data)) return(tibble::tibble())
    if (!all(c("x", "y") %in% names(candidate_set$data))) return(tibble::tibble())
    candidate_set$data %>%
      transmute(
        analysis = .env$analysis_id,
        x = as.numeric(.data$x),
        y = as.numeric(.data$y)
      ) %>%
      filter(is.finite(.data$x), is.finite(.data$y)) %>%
      distinct()
  })
  if (nrow(observed) == 0L) return(observed)
  join_analysis_registry_metadata(observed, registry_tbl)
}

make_delta_model_sensitivity_figure <- function(summary_tbl, prediction_tbl, registry_tbl,
                                                observed_tbl = NULL,
                                                title, subtitle = NULL, ncol = 4,
                                                x_lab = paste0(manuscript_delta_symbol(), " heterogeneity"),
                                                y_lab = "Predicted response",
                                                ic_weight_label = "IC weight",
                                                require_observed = FALSE,
                                                point_alpha = 0.42,
                                                presentation = NULL) {
  custom_presentation <- !is.null(presentation)
  presentation <- utils::modifyList(
    list(
      model_colors = NULL,
      model_linetypes = NULL,
      curve_linewidth = 0.62,
      panel_tag_size = 3.05,
      axis_title_size = 7.6,
      axis_text_size = 5.8,
      show_strips = TRUE,
      show_legend = TRUE,
      table_header_size = 2.0,
      table_body_size = 1.95,
      table_swatch_linewidth = 0.78,
      table_swatch_x = 0.125,
      table_swatch_xend = 0.190,
      table_model_x = 0.220,
      plot_table_heights = c(2.35, 0.82)
    ),
    if (is.null(presentation)) list() else presentation
  )

  candidate_color_scale <- if (is.null(presentation$model_colors)) {
    hdr_scale_color("Candidate curve")
  } else {
    ggplot2::scale_color_manual(
      name = "Candidate curve",
      values = presentation$model_colors,
      limits = names(hdr_model_colors),
      drop = FALSE
    )
  }
  candidate_linetype_scale <- if (is.null(presentation$model_linetypes)) {
    hdr_scale_linetype("Candidate curve")
  } else {
    ggplot2::scale_linetype_manual(
      name = "Candidate curve",
      values = presentation$model_linetypes,
      limits = names(hdr_model_linetypes),
      drop = FALSE
    )
  }

  pred <- join_analysis_registry_metadata(prediction_tbl, registry_tbl) %>%
    filter(is.finite(x), is.finite(yhat)) %>%
    mutate(model_family = hdr_factor_model_family(model_family))
  weights <- join_analysis_registry_metadata(summary_tbl, registry_tbl) %>%
    mutate(model_family = hdr_factor_model_family(model_family))

  include_dataset <- length(unique(pred$dataset_short[!is.na(pred$dataset_short)])) > 1

  pred <- pred %>%
    mutate(
      panel_label = format_sensitivity_panel_label(dataset_short, covariate_scenario, include_dataset = include_dataset),
      panel_order = sensitivity_panel_order(dataset_short, covariate_scenario),
      dataset_order = dplyr::case_when(manuscript_dataset_key(dataset_short) == "topography" ~ 1, manuscript_dataset_key(dataset_short) == "foliage" ~ 2, TRUE ~ 99)
    )
  panel_levels <- pred %>%
    distinct(panel_label, dataset_order, panel_order) %>%
    arrange(dataset_order, panel_order) %>%
    pull(panel_label)
  panel_tags <- tibble::tibble(
    panel_label = factor(panel_levels, levels = panel_levels),
    curve_panel_tag = letters[seq_along(panel_levels)],
    weight_panel_tag = letters[length(panel_levels) + seq_along(panel_levels)]
  )

  pred <- pred %>% mutate(panel_label = factor(panel_label, levels = panel_levels))
  observed <- tibble::tibble()
  if (!is.null(observed_tbl) && nrow(observed_tbl) > 0L) {
    observed <- observed_tbl %>%
      mutate(
        panel_label = format_sensitivity_panel_label(dataset_short, covariate_scenario, include_dataset = include_dataset),
        panel_label = factor(panel_label, levels = panel_levels)
      ) %>%
      filter(!is.na(panel_label), is.finite(.data$x), is.finite(.data$y))
  }
  if (isTRUE(require_observed) && nrow(observed) == 0L) {
    stop("Observed data points are required for this figure but no finite x/y rows were available.")
  }

  weights <- weights %>%
    mutate(
      panel_label = format_sensitivity_panel_label(dataset_short, covariate_scenario, include_dataset = include_dataset),
      panel_label = factor(panel_label, levels = panel_levels),
      model_family = hdr_factor_model_family(model_family)
    ) %>%
    filter(!is.na(panel_label), is.finite(.data$delta_IC), is.finite(.data$Akaike_weight)) %>%
    arrange(.data$panel_label, .data$delta_IC, .data$model_family) %>%
    group_by(.data$panel_label) %>%
    mutate(
      support_row = dplyr::row_number(),
      support_y = 0.58 - (.data$support_row - 1) * 0.20,
      model_family_label = format_model_family_table_label(.data$model_family),
      delta_ic_label = sprintf("%.2f", .data$delta_IC),
      weight_label = sprintf("%.2f", .data$Akaike_weight)
    ) %>%
    ungroup()

  criterion_values <- unique(as.character(weights$selection_criterion))
  criterion_values <- criterion_values[!is.na(criterion_values) & nzchar(criterion_values)]
  criterion <- if (length(criterion_values) > 0L) criterion_values[1] else "IC"
  delta_header_label <- dplyr::case_when(
    criterion == "AIC" ~ paste0(manuscript_Delta_symbol(), "AIC"),
    criterion == "AICc" ~ paste0(manuscript_Delta_symbol(), "AICc"),
    TRUE ~ paste0(manuscript_Delta_symbol(), criterion)
  )
  weight_header_label <- dplyr::case_when(
    criterion == "AIC" ~ "AIC\nweight",
    criterion == "AICc" ~ "AICc\nweight",
    TRUE ~ paste0(criterion, "\nweight")
  )

  support_header <- dplyr::bind_rows(lapply(panel_levels, function(panel_label_value) {
    tibble::tibble(
      panel_label = factor(panel_label_value, levels = panel_levels),
      x = c(0.315, 0.625, 0.845),
      y = 0.835,
      label = c("Model", delta_header_label, weight_header_label)
    )
  }))
  support_background <- tibble::tibble(
    panel_label = factor(panel_levels, levels = panel_levels),
    xmin = 0.085,
    xmax = 0.980,
    ymin = 0.055,
    ymax = 0.925
  )
  support_header_rule <- tibble::tibble(
    panel_label = factor(panel_levels, levels = panel_levels),
    x = 0.085,
    xend = 0.980,
    y = 0.720,
    yend = 0.720
  )
  support_row_rules <- dplyr::bind_rows(lapply(panel_levels, function(panel_label_value) {
    tibble::tibble(
      panel_label = factor(panel_label_value, levels = panel_levels),
      x = 0.085,
      xend = 0.980,
      y = c(0.480, 0.280),
      yend = c(0.480, 0.280)
    )
  }))

  p_pred <- ggplot(pred, aes(x, yhat, color = model_family, linetype = model_family)) +
    geom_point(
      data = observed,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      color = "grey35",
      alpha = point_alpha,
      size = 0.74,
      stroke = 0
    ) +
    geom_line(linewidth = presentation$curve_linewidth) +
    geom_text(
      data = panel_tags,
      aes(label = curve_panel_tag),
      x = -Inf,
      y = Inf,
      inherit.aes = FALSE,
      hjust = -0.35,
      vjust = 1.2,
      fontface = "bold",
      size = presentation$panel_tag_size
    ) +
    facet_wrap(~panel_label, scales = "free", ncol = ncol) +
    candidate_color_scale +
    candidate_linetype_scale +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = 1.05)),
      linetype = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = 1.05))
    ) +
    labs(x = x_lab, y = y_lab) +
    nature_manuscript_theme(base_size = 7.6) +
    theme(
      strip.text = if (isTRUE(presentation$show_strips)) element_text(size = 6.8, face = "bold") else element_blank(),
      strip.background = if (isTRUE(presentation$show_strips)) NULL else element_blank(),
      axis.title = element_text(size = presentation$axis_title_size),
      axis.text = element_text(size = presentation$axis_text_size),
      axis.title.y = element_text(margin = margin(r = 2)),
      legend.position = if (isTRUE(presentation$show_legend)) "top" else "none",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.box = "horizontal",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.key.width = grid::unit(1.08, "lines"),
      legend.key.height = grid::unit(0.48, "lines"),
      plot.margin = if (custom_presentation) margin(2, 2, 0, 2, unit = "pt") else margin(1, 2, 0, 0)
    )

  p_support <- ggplot() +
    geom_rect(
      data = support_background,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "white",
      colour = "black",
      linewidth = 0.24
    ) +
    geom_segment(
      data = support_header_rule,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.22
    ) +
    geom_segment(
      data = support_row_rules,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = "grey88",
      linewidth = 0.16
    ) +
    geom_text(
      data = panel_tags,
      aes(x = 0.018, y = 0.955, label = weight_panel_tag),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      fontface = "bold",
      size = presentation$panel_tag_size
    ) +
    geom_text(
      data = support_header,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0.5,
      vjust = 0.5,
      fontface = "bold",
      lineheight = 0.86,
      size = presentation$table_header_size
    ) +
    geom_segment(
      data = weights,
      aes(
        x = presentation$table_swatch_x,
        xend = presentation$table_swatch_xend,
        y = support_y,
        yend = support_y,
        colour = model_family,
        linetype = model_family
      ),
      inherit.aes = FALSE,
      linewidth = presentation$table_swatch_linewidth
    ) +
    geom_text(
      data = weights,
      aes(x = presentation$table_model_x, y = support_y, label = model_family_label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0.5,
      size = presentation$table_body_size
    ) +
    geom_text(
      data = weights,
      aes(x = 0.625, y = support_y, label = delta_ic_label),
      inherit.aes = FALSE,
      hjust = 0.5,
      vjust = 0.5,
      size = presentation$table_body_size
    ) +
    geom_text(
      data = weights,
      aes(x = 0.845, y = support_y, label = weight_label),
      inherit.aes = FALSE,
      hjust = 0.5,
      vjust = 0.5,
      size = presentation$table_body_size
    ) +
    facet_wrap(~panel_label, ncol = ncol) +
    candidate_color_scale +
    candidate_linetype_scale +
    guides(colour = "none", linetype = "none") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void(base_size = 7.6) +
    theme(
      strip.text = element_blank(),
      legend.position = "none",
      plot.margin = if (custom_presentation) margin(0, 2, 2, 2, unit = "pt") else margin(0, 1, 1, 0),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  if (custom_presentation) {
    return(
      (patchwork::wrap_elements(p_pred) / p_support) +
        patchwork::plot_layout(heights = presentation$plot_table_heights)
    )
  }

  (patchwork::wrap_elements(p_pred) / p_support) +
    patchwork::plot_layout(heights = c(2.35, 0.82)) +
    patchwork::plot_annotation(title = title, subtitle = subtitle) &
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 8.5, margin = ggplot2::margin(0, 0, 2, 0)),
      plot.subtitle = ggplot2::element_text(size = 7, margin = ggplot2::margin(0, 0, 4, 0)),
      plot.margin = ggplot2::margin(3, 5, 3, 3)
    )
}


make_delta_model_sensitivity_set_figure <- function() {
  registry <- bind_rows(
    results_manuscript_analysis_hierarchy$topography_delta_models$registry,
    results_manuscript_analysis_hierarchy$foliage_delta_models$registry
  )
  summary_tbl <- bind_rows(
    results_manuscript_analysis_hierarchy$topography_delta_models$summary,
    results_manuscript_analysis_hierarchy$foliage_delta_models$summary
  )
  prediction_tbl <- bind_rows(
    results_manuscript_analysis_hierarchy$topography_delta_models$predictions,
    results_manuscript_analysis_hierarchy$foliage_delta_models$predictions
  )
  fig <- make_delta_model_sensitivity_figure(
    summary_tbl,
    prediction_tbl,
    registry,
    title = paste0(manuscript_delta_symbol(), " covariate sensitivity summary"),
    subtitle = NULL,
    ncol = 3
  )
  manuscript_save_png(fig, "supporting_analysis", "support_delta_candidate_model_sensitivity_set.png", width = 9.8, height = 6.8)
  fig
}



make_extended_data_delta_model_figures <- function() {
  extended_data_width_mm <- 180
  extended_data_margin_mm <- 2 * 25.4 / 72.27
  extended_data_panel_spacing_mm <- 5.5 * 25.4 / 72.27
  extended_data_swatch_x <- 0.100
  extended_data_presentation <- function(panel_count) {
    panel_width_mm <- (
      extended_data_width_mm -
        2 * extended_data_margin_mm -
        (panel_count - 1) * extended_data_panel_spacing_mm
    ) / panel_count
    swatch_span <- 10 * 1.21 / panel_width_mm

    list(
      model_colors = c(
        linear = "#0072B2",
        logarithmic = "#7B3294",
        quadratic = "#D55E00"
      ),
      model_linetypes = c(
        linear = "longdash",
        logarithmic = "solid",
        quadratic = "dotdash"
      ),
      curve_linewidth = 0.30,
      panel_tag_size = 8 / ggplot2::.pt,
      axis_title_size = 7,
      axis_text_size = 6,
      show_strips = FALSE,
      show_legend = FALSE,
      table_header_size = 6.5 / ggplot2::.pt,
      table_body_size = 6 / ggplot2::.pt,
      table_swatch_linewidth = 0.30,
      table_swatch_x = extended_data_swatch_x,
      table_swatch_xend = extended_data_swatch_x + swatch_span,
      table_model_x = extended_data_swatch_x + swatch_span + 0.02,
      plot_table_heights = c(2.55, 0.90)
    )
  }

  figure_1 <- make_delta_model_sensitivity_figure(
    results_manuscript_analysis_hierarchy$foliage_delta_models$summary,
    results_manuscript_analysis_hierarchy$foliage_delta_models$predictions,
    results_manuscript_analysis_hierarchy$foliage_delta_models$registry,
    observed_tbl = make_delta_model_observed_points(
      results_manuscript_analysis_hierarchy$foliage_delta_models$sets,
      results_manuscript_analysis_hierarchy$foliage_delta_models$registry
    ),
    title = NULL,
    subtitle = NULL,
    ncol = 2,
    x_lab = paste0("Foliage height heterogeneity, ", manuscript_delta_symbol(), " (ft)"),
    y_lab = "Bird species diversity",
    ic_weight_label = "AICc weight",
    require_observed = TRUE,
    point_alpha = 0.56,
    presentation = extended_data_presentation(2)
  )
  figure_1_path <- manuscript_save_png(
    figure_1,
    "extended_data",
    "extended_data_figure_1_macarthur_all_analyses.png",
    width = 180 / 25.4,
    height = 112 / 25.4,
    dpi = 300
  )

  figure_2 <- make_delta_model_sensitivity_figure(
    results_manuscript_analysis_hierarchy$topography_delta_models$summary,
    results_manuscript_analysis_hierarchy$topography_delta_models$predictions,
    results_manuscript_analysis_hierarchy$topography_delta_models$registry,
    observed_tbl = make_delta_model_observed_points(
      results_manuscript_analysis_hierarchy$topography_delta_models$sets,
      results_manuscript_analysis_hierarchy$topography_delta_models$registry
    ),
    title = NULL,
    subtitle = NULL,
    ncol = 3,
    x_lab = paste0("Topographic heterogeneity, ", manuscript_delta_symbol(), " (m)"),
    y_lab = "Breeding bird species richness",
    ic_weight_label = "AIC weight",
    require_observed = TRUE,
    point_alpha = 0.40,
    presentation = extended_data_presentation(3)
  )
  figure_2_path <- manuscript_save_png(
    figure_2,
    "extended_data",
    "extended_data_figure_2_allouche_all_analyses.png",
    width = 180 / 25.4,
    height = 120 / 25.4,
    dpi = 300
  )

  list(
    figure_1 = figure_1,
    figure_2 = figure_2,
    output_paths = c(
      figure_1 = figure_1_path,
      figure_2 = figure_2_path
    )
  )
}

write_manuscript_caption_files <- function() {
  main_caption <- c(
    "Logarithmic candidates ranked first in both primary three-candidate comparisons.",
    "a, Hypothesis table linking the three candidate response shapes to their predictions and applications. The applications assume that the focal heterogeneity facet is manipulable and causally affects biodiversity within the intervention range.",
    "b, Plotting-only illustrative functions show three competing response-shape hypotheses along a conceptual heterogeneity gradient: a linear increase (blue long-dashed), a logarithmic increase that flattens across the displayed range (purple solid), and a quadratic rise followed by decline (vermillion dot-dashed). The parameters were chosen for illustration and are independent of the empirical fits.",
    paste0("c, MacArthur & MacArthur (1961) dataset (n = 13): bird species diversity plotted against foliage height ", manuscript_delta_symbol(), ", where ", manuscript_delta_symbol(), " is within-site sample variance of foliage height divided by mean foliage height and is measured in feet (ft). Grey points show observed bird-diversity values against calculated foliage height ", manuscript_delta_symbol(), ". The logarithmic model ranked first by AICc within the primary three-candidate comparison without a mean covariate. The purple solid curve shows the fitted mean response; the surrounding shading gives pointwise 95% confidence intervals for that mean. The inset gives ", manuscript_Delta_symbol(), "AICc values and AICc weights for all three candidates."),
    paste0("d, Allouche et al. (2012) dataset (n = 285): breeding bird species richness plotted against topographic ", manuscript_delta_symbol(), ", where ", manuscript_delta_symbol(), " is within-cell sample variance of elevation divided by mean elevation and is measured in meters (m). Grey points show observed richness values against calculated topographic ", manuscript_delta_symbol(), ". All three candidate models were fitted with nonspatial Gaussian OLS. The logarithmic model ranked first by AIC within the primary comparison without a mean covariate. The purple solid curve shows the fitted mean response; the surrounding shading gives pointwise 95% confidence intervals for that mean. The inset gives ", manuscript_Delta_symbol(), "AIC values and AIC weights for all three candidates.")
  )
  table_caption <- c(
    "No separate main-text Table 1 is generated by the current manuscript-facing workflow.",
    "The hypothesis table is Figure 1 panel a. Model-comparison and coefficient-check evidence is reported in Extended Data Table 1 and Extended Data Table 2."
  )
  extended_data_table_captions <- c(
    "Extended Data Table 1 | MacArthur & MacArthur (1961) model comparisons and coefficient checks. Candidate linear, quadratic, and logarithmic HDR models were compared using AICc. The primary no-covariate model set is presented first, followed by the mean-foliage-height sensitivity analysis. Shape interpretation required both model support and coefficient/geometry checks.",
    "Extended Data Table 2 | Allouche et al. (2012) model comparisons and coefficient checks. Candidate linear, quadratic, and logarithmic HDR models were compared using AIC. The primary no-covariate model set is presented first, followed by sensitivity analyses adding mean elevation alone and mean elevation plus mean elevation squared. Shape interpretation required both model support and coefficient/geometry checks."
  )
  extended_data_figure_captions <- c(
    paste0(
      "Extended Data Fig. 1 | Primary and mean-foliage-height sensitivity analyses for the MacArthur & MacArthur (1961) dataset. Panels a and c show the primary comparison, and panels b and d show the mean-foliage-height sensitivity. In panels a and b, grey points are observed bird-diversity values plotted against foliage-height ",
      manuscript_delta_symbol(),
      ". Blue long-dashed, purple solid, and vermillion dot-dashed lines show the fitted mean responses from the linear, logarithmic, and quadratic candidate models, respectively. In panel b, the grey points remain the unadjusted observations, whereas the fitted lines show conditional predictions with mean foliage height held at its empirical mean; point-to-curve distances are therefore not residuals. Panels c and d report the corresponding ",
      manuscript_Delta_symbol(),
      "AICc values and AICc weights for all three candidates."
    ),
    paste0(
      "Extended Data Fig. 2 | Primary and mean-elevation sensitivity analyses for the Allouche et al. (2012) dataset. Panels a and d show the primary comparison; panels b and e add mean elevation; panels c and f add mean elevation and squared mean elevation. All three comparisons use nonspatial Gaussian ordinary least squares. In panels a–c, grey points are observed richness values plotted against topographic ",
      manuscript_delta_symbol(),
      ". Blue long-dashed, purple solid, and vermillion dot-dashed lines show the fitted mean responses from the linear, logarithmic, and quadratic candidate models, respectively. In panels b and c, the grey points remain the unadjusted observations, whereas the fitted lines show conditional predictions at the fixed elevation values described below; point-to-curve distances are therefore not residuals. In panel b, predictions hold mean elevation at its empirical mean. In panel c, predictions hold mean elevation and squared mean elevation separately at their empirical means; together, these held values reproduce the sample-average additive contribution of the two elevation terms. Panels d–f report the corresponding ",
      manuscript_Delta_symbol(),
      "AIC values and AIC weights for all three candidates."
    )
  )
  extended_data_captions <- c(
    extended_data_table_captions,
    extended_data_figure_captions
  )
  writeLines(main_caption, file.path(manuscript_deliverable_dirs$captions, "main_text_mean_independent_HDR_caption.txt"))
  writeLines(table_caption, file.path(manuscript_deliverable_dirs$captions, "main_text_table_status.txt"))
  writeLines(extended_data_table_captions, file.path(manuscript_deliverable_dirs$captions, "Extended_Data_Table_captions.txt"))
  writeLines(extended_data_figure_captions, file.path(manuscript_deliverable_dirs$captions, "Extended_Data_Figure_captions.txt"))
  writeLines(extended_data_captions, file.path(manuscript_deliverable_dirs$captions, "Extended_Data_captions.txt"))
  invisible(TRUE)
}


make_manuscript_named_deliverables_manifest <- function() {
  tibble::tribble(
    ~output_type, ~dir_key, ~filename, ~description,
    "main_text_figure", "main_text", "main_text_mean_independent_heterogeneity_reveals_diminishing_return_HDRs.png", paste0("Main-text Figure 1: panel a summarizes the three hypotheses; panel b shows plotting-only illustrative response-shape functions; panel c shows MacArthur & MacArthur (1961) observations and the fitted mean response from the AICc-first model; panel d shows Allouche et al. (2012) observations and the fitted mean response from the AIC-first nonspatial Gaussian OLS model. Panels c and d include pointwise 95% confidence intervals for fitted means and model-ranking information for the primary ", manuscript_delta_symbol(), " comparisons."),
    "extended_data_table", "extended_data", "extended_data_table_1_macarthur_model_comparisons_and_coefficient_checks.csv", "Extended Data Table 1 CSV: MacArthur & MacArthur (1961) model comparisons and coefficient checks",
    "extended_data_table", "extended_data", "extended_data_table_1_macarthur_model_comparisons_and_coefficient_checks.md", "Extended Data Table 1 Markdown draft",
    "extended_data_table", "extended_data", "extended_data_table_1_macarthur_model_comparisons_and_coefficient_checks.docx", "Extended Data Table 1 Word draft",
    "extended_data_table", "extended_data", "extended_data_table_2_allouche_model_comparisons_and_coefficient_checks.csv", "Extended Data Table 2 CSV: Allouche et al. (2012) model comparisons and coefficient checks",
    "extended_data_table", "extended_data", "extended_data_table_2_allouche_model_comparisons_and_coefficient_checks.md", "Extended Data Table 2 Markdown draft",
    "extended_data_table", "extended_data", "extended_data_table_2_allouche_model_comparisons_and_coefficient_checks.docx", "Extended Data Table 2 Word draft",
    "support_table", "supporting_analysis", "support_complete_delta_model_comparisons.csv", paste0("Support CSV: complete ", manuscript_delta_symbol(), " model comparisons"),
    "support_table", "supporting_analysis", "support_complete_delta_model_terms.csv", paste0("Support CSV: ", manuscript_delta_symbol(), " model terms"),
    "support_table", "supporting_analysis", "support_analysis_hierarchy_registry.csv", "Support CSV: analysis hierarchy registry",
    "extended_data_figure", "extended_data", "extended_data_figure_1_macarthur_all_analyses.png", paste0("Extended Data Fig. 1: MacArthur & MacArthur (1961) all-analysis ", manuscript_delta_symbol(), " model-support figure"),
    "extended_data_figure", "extended_data", "extended_data_figure_2_allouche_all_analyses.png", paste0("Extended Data Fig. 2: Allouche et al. (2012) all-analysis ", manuscript_delta_symbol(), " model-support figure"),
    "caption", "captions", "main_text_mean_independent_HDR_caption.txt", "Main-text descriptive figure caption",
    "caption", "captions", "main_text_table_status.txt", "Main-text table status note",
    "caption", "captions", "Extended_Data_Table_captions.txt", "Extended Data table captions",
    "caption", "captions", "Extended_Data_Figure_captions.txt", "Extended Data figure captions",
    "caption", "captions", "Extended_Data_captions.txt", "Extended Data captions"
  ) %>%
    mutate(file = purrr::map2_chr(dir_key, filename, ~ file.path(manuscript_deliverable_dirs[[.x]], .y)), exists = file.exists(file)) %>%
    select(output_type, description, file, exists)
}


write_manuscript_named_deliverables_manifest <- function() {
  expected <- make_manuscript_named_deliverables_manifest()
  readr::write_csv(manuscript_csv_safe(expected), file.path(manuscript_deliverable_dirs$manifests, "manuscript_named_deliverables_manifest.csv"))

  readme <- c(
    "# Manuscript named deliverables",
    "",
    paste0("This folder contains the ", manuscript_delta_symbol(), "-only HDR model-comparison figures, Extended Data Table 1 and Extended Data Table 2 drafts, support CSVs, and captions explicitly named in the Nature Ecology & Evolution Brief Communication workflow."),
    "",
    "Main-text files are in `main_text/`.",
    "Figure 1 panels c and d show fitted mean responses from the AICc-first and AIC-first models, respectively, with pointwise 95% confidence intervals for fitted means.",
    "The main manuscript workflow does not generate a separate main-text Table 1; the hypothesis table is Figure 1 panel a.",
    "Extended Data Table 1, Extended Data Table 2, Extended Data Fig. 1, and Extended Data Fig. 2 files are in `extended_data/`.",
    "Machine-readable support CSVs and analysis-hierarchy products are in `supporting_analysis/`.",
    "Captions are in `captions/`.",
    "A machine-readable manifest is in `manifests/manuscript_named_deliverables_manifest.csv`.",
    "",
    "Metric-mean diagnostic tables and figures are generated only by the optional R reimplementation of selected Pellett & Valbuena analyses."
  )
  writeLines(readme, file.path(manuscript_deliverable_dirs$manifests, "README_manuscript_named_deliverables.md"))
  expected
}



run_manuscript_deliverable_workflow <- function() {
  message("Writing manuscript-analysis figures, tables, captions, and manifests...")
  all_comparisons <- assemble_delta_model_comparison_table()
  all_terms <- assemble_delta_model_term_table()
  table_sets <- write_extended_data_and_supporting_analysis_tables(all_comparisons, all_terms)
  word_table_sets <- write_extended_data_table_drafts(table_sets)
  fig_1 <- make_main_text_hdr_figure()
  extended_data_figures <- make_extended_data_delta_model_figures()
  write_manuscript_caption_files()
  expected <- write_manuscript_named_deliverables_manifest()
  message("Manuscript-analysis deliverables complete. Output folder: ", manuscript_outputs_root)
  list(extended_data_table_drafts = word_table_sets, extended_data_tables = table_sets, figure_1 = fig_1, extended_data_figures = extended_data_figures, manifest = expected)
}
