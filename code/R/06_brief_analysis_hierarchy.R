# Defines and runs the five prespecified δ-HDR manuscript analysis sets.

# 03B MANUSCRIPT ANALYSIS HIERARCHY ####
# ============================================================ #
# The manuscript comparisons use only δ as the heterogeneity metric.
# Diagnostics of dependence on the mean across δ and other metrics remain in the
# optional Pellett & Valbuena workflow.
#
# Manuscript analysis sets:
#   1. MacArthur & MacArthur (1961) dataset:
#      Primary: bird diversity ~ f(foliage-height δ), no mean covariate.
#      Sensitivity: add mean foliage height.
#   2. Allouche et al. (2012) dataset:
#      Primary: breeding-bird richness ~ f(topographic δ), no mean covariate.
#      Sensitivities: add mean elevation; add mean elevation + mean elevation².
# Each set compares linear, logarithmic log(δ), and quadratic candidates.
# The logarithmic candidate uses log(δ) directly; the workflow does not repeat
# it after adding a small constant to δ.

analysis_hierarchy_root <- file.path(manuscript_outputs_root, "supporting_analysis", "analysis_hierarchy")
analysis_hierarchy_dirs <- list(
  primary = file.path(analysis_hierarchy_root, "primary_analyses"),
  sensitivity = file.path(analysis_hierarchy_root, "sensitivity_analyses"),
  optional_provenance = file.path(analysis_hierarchy_root, "optional_provenance"),
  registry = file.path(analysis_hierarchy_root, "analysis_registry_and_manifests")
)
workflow_setup_dirs <- c(workflow_setup_dirs, unlist(analysis_hierarchy_dirs, use.names = FALSE))
brief_output_dirs <- c(brief_output_dirs, unlist(analysis_hierarchy_dirs, use.names = FALSE))


run_topography_delta_covariate_sensitivity <- function(write_outputs = TRUE) {
  message("Running the Allouche et al. (2012) primary and sensitivity analysis sets...")
  dat <- results_full_R$allouche$allouche_model_df %>% mutate(mu2 = mu^2)

  metric_role <- "focal_mean_independent_lower_bounded_metric"
  covariate_specs <- list(
    no_mean_covariates = list(
      covariates = character(0),
      tier = "main_text_primary",
      destination = "main text / Figure 1 / Extended Data Table 2",
      manuscript_role = "primary Allouche et al. analysis",
      primary_or_sensitivity = "primary",
      dataset_display_order = 2L,
      analysis_display_order = 3L,
      rationale = "Primary landmark Allouche et al. (2012) \u03b4 response-shape analysis; no mean covariate."
    ),
    mean_elevation_linear = list(
      covariates = c("mu"),
      tier = "supporting_analysis_delta_covariate_sensitivity",
      destination = "Extended Data Table 2 / Extended Data Fig. 2",
      manuscript_role = "Allouche et al. mean-elevation sensitivity",
      primary_or_sensitivity = "sensitivity",
      dataset_display_order = 2L,
      analysis_display_order = 4L,
      rationale = "Tests whether the Allouche et al. (2012) \u03b4 response shape depends on mean elevation."
    ),
    mean_elevation_quadratic = list(
      covariates = c("mu", "mu2"),
      tier = "supporting_analysis_delta_covariate_sensitivity",
      destination = "Extended Data Table 2 / Extended Data Fig. 2",
      manuscript_role = "Allouche et al. linear-plus-quadratic mean-elevation sensitivity",
      primary_or_sensitivity = "sensitivity",
      dataset_display_order = 2L,
      analysis_display_order = 5L,
      rationale = "Tests whether the Allouche et al. (2012) \u03b4 response shape depends on linear plus quadratic mean-elevation structure."
    )
  )

  sets <- list()
  registry <- list()
  for (cov_id in names(covariate_specs)) {
    covset <- covariate_specs[[cov_id]]
    analysis_id <- paste("Topography", "delta", cov_id, "gaussian", sep = "_")
    registry[[analysis_id]] <- tibble::tibble(
      dataset = manuscript_dataset_short("topography"),
      analysis = analysis_id,
      analysis_tier = covset$tier,
      response = "breeding-bird richness",
      heterogeneity_metric = "delta",
      metric_role = metric_role,
      covariate_scenario = cov_id,
      covariates = ifelse(length(covset$covariates) == 0, "none", paste(covset$covariates, collapse = ";")),
      candidate_curve_forms = "linear; logarithmic log(\u03b4); quadratic",
      status = "executed_in_manuscript_hierarchy",
      destination = covset$destination,
      manuscript_role = covset$manuscript_role,
      primary_or_sensitivity = covset$primary_or_sensitivity,
      dataset_display_order = covset$dataset_display_order,
      analysis_display_order = covset$analysis_display_order,
      rationale = covset$rationale
    )
    sets[[analysis_id]] <- fit_delta_candidate_curve_models(
      dat,
      response = "richness",
      predictor = "delta",
      covariates = covset$covariates,
      label = analysis_id,
      predictor_role = metric_role,
      information_criterion = "AIC"
    )
  }

  out <- combine_candidate_curve_model_sets(
    sets,
    "sensitivity",
    "Topography_delta_covariate_sensitivity",
    write_outputs = write_outputs
  )
  registry_tbl <- manuscript_add_dataset_labels(bind_rows(registry))
  best <- summarize_best_curve_by_analysis(out$summary) %>%
    left_join(registry_tbl, by = "analysis", suffix = c("", "_registry")) %>%
    relocate(dataset, dataset_short, dataset_phrase, source_study, analysis_tier, destination, .before = analysis)
  write_analysis_hierarchy_csv(
    registry_tbl, "registry", "Topography_delta_covariate_sensitivity_registry.csv",
    write_outputs = write_outputs
  )
  write_analysis_hierarchy_csv(
    best, "sensitivity", "Topography_delta_covariate_sensitivity_best_models.csv",
    write_outputs = write_outputs
  )
  list(sets = out$sets, summary = out$summary, terms = out$terms,
       predictions = out$predictions, registry = registry_tbl, best = best)
}


run_foliage_delta_covariate_sensitivity <- function(write_outputs = TRUE) {
  message("Running the MacArthur & MacArthur (1961) primary and sensitivity analysis sets...")
  dat <- results_full_R$macarthur$metrics

  metric_role <- "focal_mean_independent_lower_bounded_metric"
  covariate_specs <- list(
    unadjusted = list(
      covariates = character(0),
      tier = "main_text_primary",
      destination = "main text / Figure 1 / Extended Data Table 1",
      manuscript_role = "primary MacArthur & MacArthur analysis",
      primary_or_sensitivity = "primary",
      dataset_display_order = 1L,
      analysis_display_order = 1L,
      rationale = "Primary landmark MacArthur & MacArthur (1961) \u03b4 response-shape analysis; no mean covariate."
    ),
    mean_foliage_height = list(
      covariates = c("mu"),
      tier = "supporting_analysis_delta_covariate_sensitivity",
      destination = "Extended Data Table 1 / Extended Data Fig. 1",
      manuscript_role = "MacArthur & MacArthur mean-foliage-height sensitivity",
      primary_or_sensitivity = "sensitivity",
      dataset_display_order = 1L,
      analysis_display_order = 2L,
      rationale = "Tests whether the MacArthur & MacArthur (1961) \u03b4 response shape depends on mean foliage height."
    )
  )

  sets <- list()
  registry <- list()
  for (cov_id in names(covariate_specs)) {
    covset <- covariate_specs[[cov_id]]
    analysis_id <- paste("Foliage", "delta", cov_id, "gaussian", sep = "_")
    registry[[analysis_id]] <- tibble::tibble(
      dataset = manuscript_dataset_short("foliage"),
      analysis = analysis_id,
      analysis_tier = covset$tier,
      response = "bird species diversity",
      heterogeneity_metric = "delta",
      metric_role = metric_role,
      covariate_scenario = cov_id,
      covariates = ifelse(length(covset$covariates) == 0, "none", paste(covset$covariates, collapse = ";")),
      candidate_curve_forms = "linear; logarithmic log(\u03b4); quadratic",
      status = "executed_in_manuscript_hierarchy",
      destination = covset$destination,
      manuscript_role = covset$manuscript_role,
      primary_or_sensitivity = covset$primary_or_sensitivity,
      dataset_display_order = covset$dataset_display_order,
      analysis_display_order = covset$analysis_display_order,
      rationale = covset$rationale
    )
    sets[[analysis_id]] <- fit_delta_candidate_curve_models(
      dat,
      response = "bird_species_diversity",
      predictor = "delta",
      covariates = covset$covariates,
      label = analysis_id,
      predictor_role = metric_role,
      information_criterion = "AICc"
    )
  }

  out <- combine_candidate_curve_model_sets(
    sets,
    "sensitivity",
    "Foliage_delta_covariate_sensitivity",
    write_outputs = write_outputs
  )
  registry_tbl <- manuscript_add_dataset_labels(bind_rows(registry))
  best <- summarize_best_curve_by_analysis(out$summary) %>%
    left_join(registry_tbl, by = "analysis", suffix = c("", "_registry")) %>%
    relocate(dataset, dataset_short, dataset_phrase, source_study, analysis_tier, destination, .before = analysis)
  write_analysis_hierarchy_csv(
    registry_tbl, "registry", "Foliage_delta_covariate_sensitivity_registry.csv",
    write_outputs = write_outputs
  )
  write_analysis_hierarchy_csv(
    best, "sensitivity", "Foliage_delta_covariate_sensitivity_best_models.csv",
    write_outputs = write_outputs
  )
  list(sets = out$sets, summary = out$summary, terms = out$terms,
       predictions = out$predictions, registry = registry_tbl, best = best)
}

make_manuscript_analysis_registry <- function(topography_delta_models, foliage_delta_models,
                                              write_outputs = TRUE) {
  core_registry <- tibble::tribble(
    ~dataset, ~analysis, ~analysis_tier, ~response, ~heterogeneity_metric, ~metric_role, ~covariate_scenario, ~covariates, ~candidate_curve_forms, ~status, ~destination, ~manuscript_role, ~primary_or_sensitivity, ~dataset_display_order, ~analysis_display_order, ~rationale,
    manuscript_dataset_short("foliage"), "Foliage_delta_unadjusted_gaussian", "main_text_primary", "bird species diversity", "delta", "focal_mean_independent_lower_bounded_metric", "unadjusted", "none", "linear; logarithmic log(\u03b4); quadratic", "executed_in_manuscript_hierarchy", "main text / Figure 1 / Extended Data Table 1", "primary MacArthur & MacArthur analysis", "primary", 1L, 1L, "Primary landmark MacArthur & MacArthur (1961) \u03b4 response-shape analysis; no mean covariate.",
    manuscript_dataset_short("foliage"), "Foliage_delta_mean_foliage_height_gaussian", "supporting_analysis_delta_covariate_sensitivity", "bird species diversity", "delta", "focal_mean_independent_lower_bounded_metric", "mean_foliage_height", "mu", "linear; logarithmic log(\u03b4); quadratic", "executed_in_manuscript_hierarchy", "Extended Data Table 1 / Extended Data Fig. 1", "MacArthur & MacArthur mean-foliage-height sensitivity", "sensitivity", 1L, 2L, "Tests whether the MacArthur & MacArthur (1961) \u03b4 response shape depends on mean foliage height.",
    manuscript_dataset_short("topography"), "Topography_delta_no_mean_covariates_gaussian", "main_text_primary", "breeding-bird richness", "delta", "focal_mean_independent_lower_bounded_metric", "no_mean_covariates", "none", "linear; logarithmic log(\u03b4); quadratic", "executed_in_manuscript_hierarchy", "main text / Figure 1 / Extended Data Table 2", "primary Allouche et al. analysis", "primary", 2L, 3L, "Primary landmark Allouche et al. (2012) \u03b4 response-shape analysis; no mean covariate.",
    manuscript_dataset_short("topography"), "Topography_delta_mean_elevation_linear_gaussian", "supporting_analysis_delta_covariate_sensitivity", "breeding-bird richness", "delta", "focal_mean_independent_lower_bounded_metric", "mean_elevation_linear", "mu", "linear; logarithmic log(\u03b4); quadratic", "executed_in_manuscript_hierarchy", "Extended Data Table 2 / Extended Data Fig. 2", "Allouche et al. mean-elevation sensitivity", "sensitivity", 2L, 4L, "Tests whether the Allouche et al. (2012) \u03b4 response shape depends on mean elevation.",
    manuscript_dataset_short("topography"), "Topography_delta_mean_elevation_quadratic_gaussian", "supporting_analysis_delta_covariate_sensitivity", "breeding-bird richness", "delta", "focal_mean_independent_lower_bounded_metric", "mean_elevation_quadratic", "mu;mu2", "linear; logarithmic log(\u03b4); quadratic", "executed_in_manuscript_hierarchy", "Extended Data Table 2 / Extended Data Fig. 2", "Allouche et al. linear-plus-quadratic mean-elevation sensitivity", "sensitivity", 2L, 5L, "Tests whether the Allouche et al. (2012) \u03b4 response shape depends on linear plus quadratic mean-elevation structure.",
    manuscript_dataset_short("both"), "exact_Pellett_Valbuena_replication", "optional_provenance", "outputs from the optional R reimplementation", "original Pellett and Valbuena metrics", "provenance", "original workflow", "various", "original workflow models", "executed_in_exact_replication_layer", "optional-workflow provenance or code archive", "provenance for the optional R reimplementation of selected Pellett & Valbuena analyses", "provenance", 99L, 99L, "Documents how the manuscript analyses build from the optional R reimplementation of selected Pellett & Valbuena analyses."
  )

  registry <- bind_rows(core_registry, topography_delta_models$registry, foliage_delta_models$registry) %>%
    manuscript_add_dataset_labels() %>%
    mutate(
      response_pretty = manuscript_dataset_response(dataset_short),
      heterogeneity_metric_pretty = format_metric_label(heterogeneity_metric, dataset_short),
      covariates_pretty = format_covariate_label(covariates, dataset_short),
      covariate_scenario_pretty = format_covariate_label(covariate_scenario, dataset_short)
    ) %>%
    distinct(dataset, analysis, .keep_all = TRUE) %>%
    arrange(dataset_display_order, analysis_display_order, dataset, analysis)

  write_analysis_hierarchy_csv(
    registry, "registry", "analysis_hierarchy_registry.csv",
    write_outputs = write_outputs
  )

  primary_table <- registry %>%
    filter(primary_or_sensitivity == "primary") %>%
    select(dataset, dataset_short, dataset_phrase, source_study, analysis, manuscript_role, primary_or_sensitivity, dataset_display_order, analysis_display_order, response_pretty, response, heterogeneity_metric_pretty, heterogeneity_metric, covariates_pretty, covariates, destination, rationale)
  write_analysis_hierarchy_csv(
    primary_table, "primary", "primary_analysis_hierarchy.csv",
    write_outputs = write_outputs
  )

  optional_provenance_table <- registry %>%
    filter(.data$analysis_tier == "optional_provenance") %>%
    select(dataset, dataset_short, dataset_phrase, source_study, analysis, analysis_tier, manuscript_role, primary_or_sensitivity, dataset_display_order, analysis_display_order, response_pretty, response, heterogeneity_metric_pretty, heterogeneity_metric, covariates_pretty, covariates, destination, rationale)
  write_analysis_hierarchy_csv(
    optional_provenance_table, "optional_provenance", "optional_provenance_analysis_hierarchy.csv",
    write_outputs = write_outputs
  )

  sensitivity_table <- registry %>%
    filter(grepl("^supporting_analysis", analysis_tier)) %>%
    select(dataset, dataset_short, dataset_phrase, source_study, analysis, analysis_tier, manuscript_role, primary_or_sensitivity, dataset_display_order, analysis_display_order, response_pretty, response, heterogeneity_metric_pretty, heterogeneity_metric, covariates_pretty, covariates, destination, rationale)
  write_analysis_hierarchy_csv(
    sensitivity_table, "sensitivity", "sensitivity_analysis_hierarchy.csv",
    write_outputs = write_outputs
  )

  list(registry = registry, primary_table = primary_table, optional_provenance_table = optional_provenance_table, sensitivity_table = sensitivity_table)
}


make_compact_best_model_summary <- function(topography_delta_models, foliage_delta_models,
                                            write_outputs = TRUE) {
  compact <- bind_rows(
    topography_delta_models$best %>%
      mutate(result_family = paste0(manuscript_dataset_short("topography"), " \u03b4 primary/sensitivity hierarchy")),
    foliage_delta_models$best %>%
      mutate(result_family = paste0(manuscript_dataset_short("foliage"), " \u03b4 primary/sensitivity hierarchy"))
  ) %>%
    select(
      result_family, dataset, dataset_short, dataset_phrase, source_study, analysis_tier, manuscript_role, primary_or_sensitivity, dataset_display_order, analysis_display_order, destination, analysis,
      response, heterogeneity_metric, covariates, n,
      best_model_family, best_shape_interpretation,
      best_Akaike_weight, second_best_delta_IC,
      plausible_models_delta_IC_lt_2, evidence_strength,
      best_heterogeneity_term_1,
      best_heterogeneity_term_1_estimate,
      best_heterogeneity_term_1_std_error,
      best_heterogeneity_term_1_statistic,
      best_heterogeneity_term_1_p_value,
      best_heterogeneity_term_1_conf_low,
      best_heterogeneity_term_1_conf_high,
      best_heterogeneity_term_1_clear_positive,
      best_heterogeneity_term_1_clear_negative,
      best_heterogeneity_term_2,
      best_heterogeneity_term_2_estimate,
      best_heterogeneity_term_2_std_error,
      best_heterogeneity_term_2_statistic,
      best_heterogeneity_term_2_p_value,
      best_heterogeneity_term_2_conf_low,
      best_heterogeneity_term_2_conf_high,
      best_heterogeneity_term_2_clear_positive,
      best_heterogeneity_term_2_clear_negative,
      best_quadratic_vertex_x,
      best_quadratic_vertex_inside_observed_range,
      best_quadratic_high_side_declines,
      best_heterogeneity_effect_evidence,
      best_shape_support_with_uncertainty,
      best_shape_supported_by_coefficients,
      best_shape_support_notes,
      recommended_interpretation, coefficient_qualified_interpretation, rationale
    ) %>%
    arrange(dataset_display_order, analysis_display_order, dataset, analysis)

  write_analysis_hierarchy_csv(
    compact, "registry", "compact_best_model_summary.csv",
    write_outputs = write_outputs
  )
  compact
}

make_primary_and_core_sensitivity_results <- function(topography_delta_models, foliage_delta_models,
                                                      write_outputs = TRUE) {
  primary <- bind_rows(
    foliage_delta_models$best %>%
      filter(analysis == "Foliage_delta_unadjusted_gaussian"),
    foliage_delta_models$best %>%
      filter(analysis == "Foliage_delta_mean_foliage_height_gaussian"),
    topography_delta_models$best %>%
      filter(analysis == "Topography_delta_no_mean_covariates_gaussian"),
    topography_delta_models$best %>%
      filter(analysis == "Topography_delta_mean_elevation_linear_gaussian"),
    topography_delta_models$best %>%
      filter(analysis == "Topography_delta_mean_elevation_quadratic_gaussian")
  ) %>%
    arrange(dataset_display_order, analysis_display_order, dataset, analysis) %>%
    mutate(
      interpretation_for_manuscript = dplyr::case_when(
        recommended_interpretation == "positive diminishing-return HDR" ~ "Positive diminishing-return relationship between the variance-to-mean heterogeneity metric (δ) and bird diversity or richness.",
        recommended_interpretation == "positive monotonic HDR" ~ "Positive monotonic relationship between the variance-to-mean heterogeneity metric (δ) and bird diversity or richness.",
        recommended_interpretation == "unimodal HDR candidate" ~ "Quadratic candidate estimates a within-range peak; report it as a result of the implemented candidate-set comparison and apply the stated coefficient and geometry checks.",
        TRUE ~ recommended_interpretation
      )
    ) %>%
    select(
      manuscript_role, primary_or_sensitivity, dataset_display_order, analysis_display_order, dataset, dataset_short, dataset_phrase, source_study, analysis, heterogeneity_metric, covariates,
      best_model_family, best_shape_interpretation, best_Akaike_weight,
      second_best_delta_IC, evidence_strength,
      plausible_models_delta_IC_lt_2,
      best_heterogeneity_term_1,
      best_heterogeneity_term_1_estimate,
      best_heterogeneity_term_1_std_error,
      best_heterogeneity_term_1_statistic,
      best_heterogeneity_term_1_p_value,
      best_heterogeneity_term_1_conf_low,
      best_heterogeneity_term_1_conf_high,
      best_heterogeneity_term_1_clear_positive,
      best_heterogeneity_term_1_clear_negative,
      best_heterogeneity_term_2,
      best_heterogeneity_term_2_estimate,
      best_heterogeneity_term_2_std_error,
      best_heterogeneity_term_2_statistic,
      best_heterogeneity_term_2_p_value,
      best_heterogeneity_term_2_conf_low,
      best_heterogeneity_term_2_conf_high,
      best_heterogeneity_term_2_clear_positive,
      best_heterogeneity_term_2_clear_negative,
      best_quadratic_vertex_x,
      best_quadratic_vertex_inside_observed_range,
      best_quadratic_high_side_declines,
      best_heterogeneity_effect_evidence,
      best_shape_support_with_uncertainty,
      best_shape_supported_by_coefficients,
      best_shape_support_notes,
      coefficient_qualified_interpretation,
      interpretation_for_manuscript
    )
  write_analysis_hierarchy_csv(
    primary, "primary", "primary_and_core_sensitivity_results.csv",
    write_outputs = write_outputs
  )
  primary
}


write_analysis_hierarchy_readme <- function(write_outputs = TRUE) {
  readme <- c(
    "# Manuscript-facing analysis hierarchy",
    "",
    "This folder mirrors the manuscript analysis hierarchy. Manuscript-facing HDR model comparisons use δ only. Metric-mean diagnostic tables and figures are generated only when the optional Pellett & Valbuena provenance workflow is run.",
    "",
    "## primary_analyses",
    "Contains the two primary analyses in manuscript order: MacArthur & MacArthur (1961) \u03b4 with no mean covariate, followed by Allouche et al. (2012) \u03b4 with no mean covariate.",
    "",
    "## sensitivity_analyses",
    "Contains \u03b4-only sensitivity analyses and candidate-model summaries: MacArthur & MacArthur (1961) with mean foliage height, and Allouche et al. (2012) with mean elevation or mean elevation plus mean elevation squared. It does not contain metric-mean diagnostic outputs or non-\u03b4 HDR model-comparison grids.",
    "",
    "## optional_provenance",
    "Contains the provenance record for the optional R reimplementation of selected Pellett & Valbuena analyses. This record is supporting analysis, not manuscript Extended Data.",
    "",
    "## analysis_registry_and_manifests",
    "Contains the registry and compact best-model summaries used to track which \u03b4 analyses are manuscript primaries, sensitivity analyses, or provenance-only records.",
    "",
    "## Interpretation rule",
    "Information criteria choose among linear, logarithmic, and quadratic \u03b4-HDR candidates; coefficient summaries and geometry checks then qualify the direction and certainty of the selected shape."
  )
  if (isTRUE(write_outputs)) {
    writeLines(readme, file.path(analysis_hierarchy_root, "README_analysis_hierarchy.md"))
  }
  invisible(readme)
}


run_manuscript_analysis_hierarchy <- function(write_outputs = TRUE) {
  message("Running the five prespecified δ-HDR analysis sets...")
  topography_delta_models <- run_topography_delta_covariate_sensitivity(write_outputs = write_outputs)
  foliage_delta_models <- run_foliage_delta_covariate_sensitivity(write_outputs = write_outputs)
  hierarchy <- make_manuscript_analysis_registry(
    topography_delta_models, foliage_delta_models,
    write_outputs = write_outputs
  )
  compact_best <- make_compact_best_model_summary(
    topography_delta_models, foliage_delta_models,
    write_outputs = write_outputs
  )
  primary_table <- make_primary_and_core_sensitivity_results(
    topography_delta_models, foliage_delta_models,
    write_outputs = write_outputs
  )
  write_analysis_hierarchy_readme(write_outputs = write_outputs)

  manifest <- tibble::tibble(
    output_type = c(
      "analysis_registry",
      "compact_best_model_summary",
      "primary_analysis_hierarchy",
      "primary_core_results",
      "optional_provenance_hierarchy",
      "sensitivity_analysis_hierarchy",
      "topography_delta_covariate_sensitivity",
      "foliage_delta_covariate_sensitivity",
      "readme"
    ),
    file = c(
      file.path(analysis_hierarchy_dirs$registry, "analysis_hierarchy_registry.csv"),
      file.path(analysis_hierarchy_dirs$registry, "compact_best_model_summary.csv"),
      file.path(analysis_hierarchy_dirs$primary, "primary_analysis_hierarchy.csv"),
      file.path(analysis_hierarchy_dirs$primary, "primary_and_core_sensitivity_results.csv"),
      file.path(analysis_hierarchy_dirs$optional_provenance, "optional_provenance_analysis_hierarchy.csv"),
      file.path(analysis_hierarchy_dirs$sensitivity, "sensitivity_analysis_hierarchy.csv"),
      file.path(analysis_hierarchy_dirs$sensitivity, "Topography_delta_covariate_sensitivity_model_comparison.csv"),
      file.path(analysis_hierarchy_dirs$sensitivity, "Foliage_delta_covariate_sensitivity_model_comparison.csv"),
      file.path(analysis_hierarchy_root, "README_analysis_hierarchy.md")
    )
  ) %>%
    mutate(exists = if (isTRUE(write_outputs)) file.exists(file) else NA)

  write_analysis_hierarchy_csv(
    manifest, "registry", "manuscript_analysis_hierarchy_manifest.csv",
    write_outputs = write_outputs
  )
  message("Manuscript analysis sets complete. Output location: ", analysis_hierarchy_root)
  list(topography_delta_models = topography_delta_models, foliage_delta_models = foliage_delta_models, hierarchy = hierarchy, compact_best = compact_best, primary_table = primary_table, manifest = manifest)
}
