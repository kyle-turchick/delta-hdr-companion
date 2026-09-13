# Defines file, workflow lock, and run metadata helpers without writing when sourced.

# Create each requested directory separately so vector inputs are handled safely.
ensure_dir <- function(path) {
  if (length(path) == 0L) return(invisible(character(0)))
  path <- unique(as.character(path[!is.na(path) & nzchar(path)]))
  for (p in path) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(normalizePath(path, winslash = "/", mustWork = FALSE))
}

ensure_parent_dir <- function(path) {
  ensure_dir(dirname(as.character(path)))
}

project_setup_dirs <- c(
  manifest_dir,
  script_dir,
  exact_output_root,
  hdr_output_root,
  archive_root,
  scratch_root
)

workflow_machine_name <- function() {
  machine <- Sys.info()[["nodename"]]
  if (is.na(machine) || !nzchar(machine)) {
    machine <- Sys.getenv("COMPUTERNAME", unset = Sys.getenv("HOSTNAME", unset = "unknown"))
  }
  machine
}

workflow_user_name <- function() {
  user <- Sys.info()[["user"]]
  if (is.na(user) || !nzchar(user)) {
    user <- Sys.getenv("USER", unset = Sys.getenv("USERNAME", unset = "unknown"))
  }
  user
}

read_workflow_lock <- function(lock_path = workflow_lock_path) {
  if (!file.exists(lock_path)) {
    return(character(0))
  }
  readLines(lock_path, warn = FALSE)
}

acquire_workflow_lock <- function(lock_path = workflow_lock_path) {
  ensure_dir(dirname(lock_path))
  if (file.exists(lock_path)) {
    lock_lines <- read_workflow_lock(lock_path)
    stop(
      "Workflow lock already exists at: ", lock_path,
      "\nAnother workflow session may still be running or may have stopped before releasing the lock.",
      "\nLock contents:\n", paste(lock_lines, collapse = "\n"),
      "\nConfirm that no workflow is running before deleting this lock file manually.",
      call. = FALSE
    )
  }

  lock_lines <- c(
    paste0("machine: ", workflow_machine_name()),
    paste0("user: ", workflow_user_name()),
    paste0("timestamp: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"))
  )
  writeLines(lock_lines, lock_path, useBytes = TRUE)
  invisible(lock_path)
}

release_workflow_lock <- function(lock_path = workflow_lock_path) {
  if (file.exists(lock_path)) {
    file.remove(lock_path)
  }
  invisible(TRUE)
}

project_setup_dirs <- c(
  manifest_dir,
  script_dir,
  archive_root,
  scratch_root
)

# Output folders for the optional R reimplementation of selected Pellett & Valbuena analyses.
out_dir <- exact_output_root
fig_dir <- file.path(exact_output_root, "figures")
tab_dir <- file.path(exact_output_root, "tables")
intermediate_dir <- file.path(exact_output_root, "intermediate")
replication_output_dirs <- c(exact_output_root, fig_dir, tab_dir, intermediate_dir)

# Advanced helpers for optional comparison runs between Julia and R. Defining them creates
# no files or folders; only an explicit run in audit mode calls them.
make_pv_audit_run_id <- function(time = Sys.time()) {
  format(time, "%Y%m%d_%H%M%S")
}

safe_audit_filename <- function(filename) {
  filename <- basename(as.character(filename))
  filename <- gsub("[^A-Za-z0-9._-]+", "_", filename)
  filename <- gsub("_+", "_", filename)
  filename
}

pv_audit_enabled <- function(
  audit_mode = RUN_PV_AUDIT_MODE,
  write_intermediates = WRITE_PV_AUDIT_INTERMEDIATES
) {
  isTRUE(audit_mode) && isTRUE(write_intermediates)
}

path_has_prefix <- function(path, prefix) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- normalizePath(prefix, winslash = "/", mustWork = FALSE)
  identical(path, prefix) || startsWith(path, paste0(prefix, "/"))
}

validate_pv_audit_output_root <- function(
  audit_output_root,
  project_root_path = project_root,
  exact_output_root_path = exact_output_root
) {
  if (is.null(audit_output_root) || length(audit_output_root) != 1L || !nzchar(audit_output_root)) {
    stop("PV audit output root must be a single non-empty path when audit mode is enabled.", call. = FALSE)
  }

  audit_output_root <- normalizePath(audit_output_root, winslash = "/", mustWork = FALSE)
  allowed_root <- normalizePath(file.path(exact_output_root_path, "audit_runs"), winslash = "/", mustWork = FALSE)
  project_root_path <- normalizePath(project_root_path, winslash = "/", mustWork = FALSE)
  exact_output_root_path <- normalizePath(exact_output_root_path, winslash = "/", mustWork = FALSE)

  protected_roots <- c(
    file.path(project_root_path, "data", "external", "pellett_valbuena_2025"),
    file.path(project_root_path, "results", "delta_hdr_reanalysis"),
    file.path(project_root_path, "metadata", "baselines")
  )
  protected_roots <- normalizePath(protected_roots, winslash = "/", mustWork = FALSE)

  if (basename(audit_output_root) != "r_replication_outputs") {
    stop("PV audit output root must end in r_replication_outputs: ", audit_output_root, call. = FALSE)
  }
  if (!path_has_prefix(audit_output_root, allowed_root)) {
    stop("PV audit output root must live under: ", allowed_root, call. = FALSE)
  }
  if (identical(audit_output_root, exact_output_root_path)) {
    stop("PV audit output root must not be the normal provenance output root.", call. = FALSE)
  }
  if (any(vapply(protected_roots, function(root) path_has_prefix(audit_output_root, root), logical(1)))) {
    stop("PV audit output root overlaps a protected project path: ", audit_output_root, call. = FALSE)
  }
  if (grepl("MI-HDR_NatComm", audit_output_root, fixed = TRUE)) {
    stop("PV audit output root must not point to old MI-HDR_NatComm paths.", call. = FALSE)
  }

  invisible(audit_output_root)
}

build_pv_audit_paths <- function(
  audit_output_root = PV_AUDIT_OUTPUT_ROOT,
  audit_run_id = PV_AUDIT_RUN_ID,
  create = FALSE
) {
  if (is.null(audit_run_id) || length(audit_run_id) != 1L || !nzchar(audit_run_id)) {
    audit_run_id <- make_pv_audit_run_id()
  }
  audit_run_id <- safe_audit_filename(audit_run_id)

  if (is.null(audit_output_root) || length(audit_output_root) == 0L || !nzchar(audit_output_root)) {
    audit_output_root <- file.path(
      exact_output_root,
      "audit_runs",
      audit_run_id,
      "r_replication_outputs"
    )
  }

  audit_output_root <- normalizePath(audit_output_root, winslash = "/", mustWork = FALSE)
  validate_pv_audit_output_root(audit_output_root)

  paths <- list(
    run_id = audit_run_id,
    root = audit_output_root,
    figures = file.path(audit_output_root, "figures"),
    tables = file.path(audit_output_root, "tables"),
    intermediate = file.path(audit_output_root, "intermediate"),
    manifests = file.path(audit_output_root, "manifests"),
    audit_run_manifest = file.path(audit_output_root, "audit_run_manifest.csv"),
    audit_session_info = file.path(audit_output_root, "audit_sessionInfo.txt"),
    audit_package_versions = file.path(audit_output_root, "audit_package_versions.csv"),
    audit_library_paths = file.path(audit_output_root, "audit_library_paths.csv"),
    audit_seed_records = file.path(audit_output_root, "audit_seed_records.csv"),
    audit_intermediate_file_manifest = file.path(audit_output_root, "audit_intermediate_file_manifest.csv"),
    audit_contract_export_manifest = file.path(audit_output_root, "audit_contract_export_manifest.csv")
  )

  if (isTRUE(create)) {
    ensure_dir(c(paths$root, paths$figures, paths$tables, paths$intermediate, paths$manifests))
  }

  paths
}

pv_audit_paths_if_enabled <- function(create = FALSE) {
  if (!isTRUE(RUN_PV_AUDIT_MODE)) {
    return(NULL)
  }
  build_pv_audit_paths(create = create)
}

audit_now_utc_string <- function(time = Sys.time()) {
  format(as.POSIXct(time, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

read_audit_manifest_character <- function(path) {
  if (!file.exists(path)) {
    return(tibble::tibble())
  }
  readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
}

normalize_audit_manifest_row <- function(row, schema = NULL) {
  row <- tibble::as_tibble(row)
  for (name in names(row)) {
    value <- row[[name]]
    if (inherits(value, "POSIXt")) {
      row[[name]] <- format(as.POSIXct(value, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    } else {
      row[[name]] <- as.character(value)
    }
    row[[name]][is.na(row[[name]])] <- NA_character_
  }
  if (is.null(schema)) {
    schema <- names(row)
  }
  missing_columns <- setdiff(schema, names(row))
  for (name in missing_columns) {
    row[[name]] <- NA_character_
  }
  row[, schema, drop = FALSE]
}

append_audit_csv_row <- function(row, path, schema = NULL) {
  ensure_parent_dir(path)
  existing <- read_audit_manifest_character(path)
  schema <- unique(c(schema, names(existing), names(row)))
  existing <- normalize_audit_manifest_row(existing, schema = schema)
  row <- normalize_audit_manifest_row(row, schema = schema)
  if (file.exists(path)) {
    row <- dplyr::bind_rows(existing, row)
  }
  readr::write_csv(row, path)
  invisible(path)
}

write_audit_manifest_row <- function(
  file_path,
  object_name = NA_character_,
  artifact_type = "intermediate",
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_intermediates = WRITE_PV_AUDIT_INTERMEDIATES,
  notes = NA_character_
) {
  if (!pv_audit_enabled(audit_mode, write_intermediates)) {
    return(invisible(NULL))
  }
  if (is.null(audit_paths)) {
    audit_paths <- build_pv_audit_paths(create = TRUE)
  }

  file_path <- normalizePath(file_path, winslash = "/", mustWork = FALSE)
  if (path_has_prefix(file_path, audit_paths$root)) {
    relative_path <- substring(file_path, nchar(audit_paths$root) + 2L)
  } else {
    relative_path <- file_path
  }
  info <- file.info(file_path)
  row <- tibble::tibble(
    audit_run_id = audit_paths$run_id,
    object_name = object_name,
    artifact_type = artifact_type,
    file_path = file_path,
    relative_path = relative_path,
    size_bytes = if (file.exists(file_path)) as.numeric(info$size) else NA_real_,
    recorded_time = audit_now_utc_string(),
    notes = notes
  )
  append_audit_csv_row(row, audit_paths$audit_intermediate_file_manifest)
  invisible(file_path)
}

write_audit_csv <- function(
  data,
  filename,
  subdir = "intermediate",
  object_name = tools::file_path_sans_ext(filename),
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_intermediates = WRITE_PV_AUDIT_INTERMEDIATES,
  notes = NA_character_
) {
  if (!pv_audit_enabled(audit_mode, write_intermediates)) {
    return(invisible(NULL))
  }
  if (is.null(audit_paths)) {
    audit_paths <- build_pv_audit_paths(create = TRUE)
  }

  target_dir <- switch(
    subdir,
    tables = audit_paths$tables,
    figures = audit_paths$figures,
    intermediate = audit_paths$intermediate,
    manifests = audit_paths$manifests,
    root = audit_paths$root,
    stop("Unsupported audit subdir: ", subdir, call. = FALSE)
  )
  path <- file.path(target_dir, safe_audit_filename(filename))
  ensure_parent_dir(path)
  readr::write_csv(data, path)
  write_audit_manifest_row(
    path,
    object_name = object_name,
    artifact_type = subdir,
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_intermediates,
    notes = notes
  )
  invisible(path)
}

pv_audit_export_registry <- function() {
  tibble::tribble(
    ~export_contract_id, ~target_id, ~figure_or_analysis, ~export_category, ~required_r_object_or_code_location, ~future_export_filename, ~required_columns_or_fields, ~required_sort_keys, ~comparison_method, ~priority, ~implementing_r_file, ~implementing_r_function_or_helper, ~required_from_r,
    "EC-001", "PV-F1", "Figure 1 panel a gamma densities", "plot_layer_dataframe", "run_figure1_concept_R pdf_df", "r_replication_outputs/intermediate/stage7b_F1_panel_a_density_layers_r.csv", "panel_id,mu,delta,x,density,layer_type", "panel_id,mu,delta,x", "sorted_dataframe_tolerance", "high", "code/R/04_pellett_valbuena_replication.R", "run_figure1_concept_R", "yes",
    "EC-002", "PV-F1", "Figure 1 MBH lines", "plot_layer_dataframe", "run_figure1_concept_R mbh_df", "r_replication_outputs/intermediate/stage7b_F1_panel_b_mbh_lines_r.csv", "measure,mu,heterogeneity,H,scale_note", "measure,mu", "sorted_dataframe_tolerance", "high", "code/R/04_pellett_valbuena_replication.R", "run_figure1_concept_R", "yes",
    "EC-003", "PV-F1", "Figure 1 MDR/HDR and observed HDR layers", "plot_layer_dataframe", "run_figure1_concept_R p_c_df and p_d_df", "r_replication_outputs/intermediate/stage7b_F1_panels_cd_relationship_lines_r.csv", "panel_id,relationship,measure,x,y,formula_label", "panel_id,relationship,measure,x", "sorted_dataframe_tolerance", "medium", "code/R/04_pellett_valbuena_replication.R", "run_figure1_concept_R", "yes",
    "EC-004", "PV-F2", "Figure 2 gamma measure curves", "plot_layer_dataframe", "make_gamma_measure_data gamma_df", "r_replication_outputs/intermediate/stage7b_F2_gamma_measure_curves_r.csv", "distribution,measure,mu,delta,value", "distribution,measure,delta,mu", "sorted_dataframe_tolerance", "critical", "code/R/04_pellett_valbuena_replication.R", "run_figure2_MBH_comb_R", "yes",
    "EC-005", "PV-F2", "Figure 2 beta measure curves", "plot_layer_dataframe", "make_beta_measure_data beta_df", "r_replication_outputs/intermediate/stage7b_F2_beta_measure_curves_r.csv", "distribution,measure,mu,delta,value", "distribution,measure,delta,mu", "sorted_dataframe_tolerance", "critical", "code/R/04_pellett_valbuena_replication.R", "run_figure2_MBH_comb_R", "yes",
    "EC-006", "PV-F2", "Figure 2 gamma and beta density panels", "plot_layer_dataframe", "Figure 2 gamma_pdf and beta_pdf panel data", "r_replication_outputs/intermediate/stage7b_F2_density_panel_layers_r.csv", "distribution,panel_id,mu,delta_or_alpha,x,density", "distribution,panel_id,mu,delta_or_alpha,x", "sorted_dataframe_tolerance", "high", "code/R/04_pellett_valbuena_replication.R", "run_figure2_MBH_comb_R", "yes",
    "EC-007", "PV-F2", "Figure 2 elevation empirical layer", "plot_layer_dataframe", "run_figure2_MBH_comb_R elev_long", "r_replication_outputs/intermediate/stage7b_F2_elevation_empirical_layers_r.csv", "source_row_id,mean,measure,value,raw_var,sd,cv,iod,scale_note", "source_row_id,measure,mean", "sorted_dataframe_tolerance", "critical", "code/R/04_pellett_valbuena_replication.R", "run_figure2_MBH_comb_R", "yes",
    "EC-008", "PV-F2", "Figure 2 crop empirical layer", "plot_layer_dataframe", "run_figure2_MBH_comb_R crop_long", "r_replication_outputs/intermediate/stage7b_F2_crop_empirical_layers_r.csv", "source_row_id,mean,measure,value,raw_var,sd,cv,delta,scale_note", "source_row_id,measure,mean", "sorted_dataframe_tolerance", "critical", "code/R/04_pellett_valbuena_replication.R", "run_figure2_MBH_comb_R", "yes",
    "EC-009", "PV-MBH", "Mean-biased residual-squared vectors", "statistical_test_table", "run_mean_biased_hypothesis_tests_R residual vectors", "r_replication_outputs/tables/stage7b_MBH_residual_vectors_r.csv", "test,row_id,h0_res2,h1_res2,input_mean,input_measure,source", "test,row_id", "sorted_dataframe_tolerance", "critical", "code/R/04_pellett_valbuena_replication.R", "run_mean_biased_hypothesis_tests_R", "yes",
    "EC-010", "PV-MBH", "Mean-biased test table", "statistical_test_table", "mean_biased_hypothesis_tests_R ans", "r_replication_outputs/tables/stage7b_MBH_test_statistics_r.csv", "test,t_statistic,df,p_value,mean_h0_res2,mean_h1_res2,method", "test", "p_value_table", "critical", "code/R/04_pellett_valbuena_replication.R", "run_mean_biased_hypothesis_tests_R", "yes",
    "EC-011", "PV-CAT-prep", "Catalonia raw input schema summary", "raw_input_summary", "prepare_allouche_data_R source files", "r_replication_outputs/manifests/stage7b_CAT_raw_input_schema_summary_r.csv", "input_name,path,row_count,column_count,column_names_or_raster_count,checksum_or_bundle_hash", "input_name,path", "schema", "high", "code/R/03_source_data_prep.R", "prepare_allouche_data_R", "yes",
    "EC-012", "PV-CAT-prep", "Prepared Catalonia adf dataframe", "cleaned_input_dataframe", "prepare_allouche_data_R adf", "r_replication_outputs/intermediate/stage7b_CAT_prepared_adf_r.csv", "UTM10,richness,Utm10,X_coord,Y_coord,fid,mu,range,sigma_or_sigma2,CV,delta,n,rich_per_area", "UTM10,Utm10,fid", "sorted_dataframe_tolerance", "critical", "code/R/03_source_data_prep.R", "prepare_allouche_data_R", "yes",
    "EC-013", "PV-CAT-models", "Allouche model input dataframe", "model_input_dataframe", "run_allouche_analysis_R allouche_model_df", "r_replication_outputs/intermediate/stage7b_CAT_model_input_dataframe_r.csv", "row_id,richness,mu,delta,range,mu_rs,delta_rs,range_rs,x,y", "row_id", "sorted_dataframe_tolerance", "critical", "code/R/03_source_data_prep.R", "run_allouche_analysis_R", "yes",
    "EC-014", "PV-CAT-models", "Allouche model coefficient table", "coefficient_table", "run_allouche_analysis_R models", "r_replication_outputs/tables/stage7b_CAT_model_terms_r.csv", "model,term,estimate,std_error,statistic,p_value,predictor_scale", "model,term,predictor_scale", "coefficient_table", "critical", "code/R/03_source_data_prep.R", "run_allouche_analysis_R", "yes",
    "EC-015", "PV-CAT-p-values", "Allouche labeled F-test p-values", "p_value_table", "run_allouche_analysis_R p_values", "r_replication_outputs/tables/stage7b_CAT_labeled_p_values_r.csv", "model,test_type,p_value,source_order,notes", "model,test_type", "p_value_table", "critical", "code/R/03_source_data_prep.R", "run_allouche_analysis_R", "yes",
    "EC-016", "PV-F3", "Figure 3 theoretical full grid", "prediction_grid", "theory_grid_R", "r_replication_outputs/intermediate/stage7b_F3_theory_grid_r.csv", "mu,H,richness,range_mbh,cv_mbh,var_mbh,delta_mbh,grid_n", "mu,H", "prediction_grid", "critical", "code/R/04_pellett_valbuena_replication.R", "run_figure3_corrected_HDR_R", "yes",
    "EC-017", "PV-F3", "Figure 3 theoretical model terms and fitted lines", "coefficient_table", "run_figure3_corrected_HDR_R theory_models and fitted lines", "r_replication_outputs/tables/stage7b_F3_theoretical_fits_r.csv", "measure,term,estimate,std_error,pred_x,pred_y,grid_n", "measure,term,pred_x", "coefficient_table", "critical", "code/R/04_pellett_valbuena_replication.R", "run_figure3_corrected_HDR_R", "yes",
    "EC-018", "PV-F3", "Figure 3 empirical plot-layer dataframe", "plot_layer_dataframe", "run_figure3_corrected_HDR_R emp_long", "r_replication_outputs/intermediate/stage7b_F3_empirical_layers_r.csv", "row_id,measure,x,richness,raw_mu,raw_range,CV,sigma_or_sigma2,delta,fit_y", "row_id,measure,x", "plot_layer_dataframe", "critical", "code/R/04_pellett_valbuena_replication.R", "run_figure3_corrected_HDR_R", "yes",
    "EC-019", "PV-F3", "Figure 3 Catalonia map raster layer summary", "raster_summary", "no R counterpart", "r_replication_outputs/manifests/stage7b_F3_catalonia_map_layer_status_r.csv", "path,checksum,width,height,band_count,color_range,layer_status", "path", "raster_summary", "medium", NA_character_, NA_character_, "no",
    "EC-020", "PV-S1", "SOM S1 gamma density layers", "plot_layer_dataframe", "run_som_s1_fixedCV_R df", "r_replication_outputs/intermediate/stage7b_SOM_s1_density_layers_r.csv", "mu,alpha_label,x,density,row_label", "alpha_label,mu,x", "sorted_dataframe_tolerance", "high", "code/R/04_pellett_valbuena_replication.R", "run_som_s1_fixedCV_R", "yes",
    "EC-021", "PV-S2-S5", "SOM scenario grids", "prediction_grid", "make_som_theory_scenario grids", "r_replication_outputs/intermediate/stage7b_SOM_s2_s5_scenario_grids_r.csv", "scenario,mu,H,MDR,HDR,richness,range_mbh,cv_mbh,var_mbh,grid_n", "scenario,mu,H", "prediction_grid", "high", "code/R/04_pellett_valbuena_replication.R", "run_som_s2_s5_R", "yes",
    "EC-022", "PV-S2-S5", "SOM scenario model terms", "coefficient_table", "make_som_theory_scenario terms", "r_replication_outputs/tables/stage7b_SOM_s2_s5_polyfit_terms_r.csv", "scenario,model,term,estimate,std_error,statistic,p_value,grid_n", "scenario,model,term", "coefficient_table", "high", "code/R/04_pellett_valbuena_replication.R", "run_som_s2_s5_R", "yes",
    "EC-023", "PV-S6", "SOM S6 plotted dataframe", "plot_layer_dataframe", "run_som_s6_delta_location_R adf", "r_replication_outputs/intermediate/stage7b_SOM_s6_delta_location_layers_r.csv", "row_id,delta,X_coord,Y_coord,X_scaled,Y_scaled", "row_id", "sorted_dataframe_tolerance", "medium", "code/R/04_pellett_valbuena_replication.R", "run_som_s6_delta_location_R", "yes",
    "EC-024", "PV-MAC", "MacArthur raw digitized point bundle manifest", "manifest_or_file_inventory", "load_macarthur_digitized_data source files", "r_replication_outputs/manifests/stage7b_MAC_raw_digitized_manifest_r.csv", "site,path,row_count,column_names,checksum", "site,path", "checksum", "high", "code/R/03_source_data_prep.R", "load_macarthur_digitized_data", "yes",
    "EC-025", "PV-MAC", "MacArthur cleaned digitized points", "cleaned_input_dataframe", "load_macarthur_digitized_data cleaned points", "r_replication_outputs/intermediate/stage7b_MAC_cleaned_digitized_points_r.csv", "site,row_id,x,y,edit_source", "site,row_id,x,y", "sorted_dataframe_tolerance", "critical", "code/R/03_source_data_prep.R", "load_macarthur_digitized_data", "yes",
    "EC-026", "PV-MAC", "MacArthur foliage metrics", "metric_table", "run_macarthur_analysis_R metrics", "r_replication_outputs/tables/stage7b_MAC_metric_table_r.csv", "site,mu,sigma2,delta,FHD3,FH_entropy_51,gini,range,sigma,CV,diff_entropy,BSD", "site", "metric_table", "critical", "code/R/03_source_data_prep.R", "run_macarthur_analysis_R", "yes",
    "EC-027", "PV-MAC", "MacArthur correlations", "statistical_test_table", "run_macarthur_analysis_R correlations", "r_replication_outputs/tables/stage7b_MAC_correlation_table_r.csv", "test,estimate_r,statistic,df,p_value,conf_low,conf_high,n", "test", "numeric_vector", "critical", "code/R/03_source_data_prep.R", "run_macarthur_analysis_R", "yes",
    "EC-028", "PV-MAC", "MacArthur foliage density loess layer data", "plot_layer_dataframe", "run_macarthur_analysis_R density loess layers", "r_replication_outputs/intermediate/stage7b_MAC_foliage_density_plot_layers_r.csv", "site,panel,x_density,y_height,predicted_density,loess_span,degree", "site,panel,y_height", "plot_layer_dataframe", "high", "code/R/03_source_data_prep.R", "run_macarthur_analysis_R", "yes",
    "EC-029", "PV-DTREND", "Delta trend model inputs", "model_input_dataframe", "run_delta_trend_test_R mutated elev and crop", "r_replication_outputs/intermediate/stage7b_DTREND_model_inputs_r.csv", "dataset,row_id,mean,rsvar,rssd,rsmean,rsiod,rsdelta,rscv,sqrtmean,sqrtmeandiv,x,sqrtx,sqrtxdiv", "dataset,row_id", "sorted_dataframe_tolerance", "high", "code/R/04_pellett_valbuena_replication.R", "run_delta_trend_test_R", "yes",
    "EC-030", "PV-DTREND", "Delta trend test statistics", "statistical_test_table", "run_delta_trend_test_R rows", "r_replication_outputs/tables/stage7b_DTREND_test_statistics_r.csv", "test,delta_abs_coef,delta_se,comparison_abs_coef,null_half_comparison,t_statistic,df,p_value,ratio", "test", "statistical_test_table", "critical", "code/R/04_pellett_valbuena_replication.R", "run_delta_trend_test_R", "yes",
    "EC-031", "PV-R2", "Review2 deterministic beta layer data", "plot_layer_dataframe", "run_review2_figures_R beta density layers", "r_replication_outputs/intermediate/stage7b_R2_beta_density_layers_r.csv", "figure,panel,p,q,x,density,mu,delta,phi,curve_x,curve_y", "figure,panel,p,q,x,curve_x", "plot_layer_dataframe", "medium", "code/R/04_pellett_valbuena_replication.R", "run_review2_figures_R", "yes",
    "EC-032", "PV-R2", "Review2 p/q simulation summary", "random_seed_record", "run_review2_figures_R pqdf summary", "r_replication_outputs/tables/stage7b_R2_pq_simulation_summary_r.csv", "simulation,seed_status,n,distribution,mean_mu,mean_delta,mean_var,quantiles,notes", "simulation", "numeric_vector", "high", "code/R/04_pellett_valbuena_replication.R", "run_review2_figures_R", "yes",
    "EC-033", "PV-R2", "Review2 mu/variance simulation summary", "random_seed_record", "run_review2_figures_R mudf summary", "r_replication_outputs/tables/stage7b_R2_muvar_simulation_summary_r.csv", "simulation,seed_status,n,mean_mu,mean_var,mean_delta,quantiles,notes", "simulation", "numeric_vector", "high", "code/R/04_pellett_valbuena_replication.R", "run_review2_figures_R", "yes",
    "EC-034", "PV-R2", "Review2 hexbin and loess layer data", "plot_layer_dataframe", "run_review2_figures_R hexbin and loess layers", "r_replication_outputs/intermediate/stage7b_R2_hexbin_loess_layers_r.csv", "simulation,panel,bin_x,bin_y,count,loess_x,loess_y,bins,span", "simulation,panel,bin_x,bin_y,loess_x", "plot_layer_dataframe", "medium", "code/R/04_pellett_valbuena_replication.R", "run_review2_figures_R", "yes",
    "EC-035", "PV-ENV-REPO", "Environment/runtime evidence", "environment_record", "write_audit_session_info runtime evidence", "r_replication_outputs/manifests/stage7b_environment_runtime_r.csv", "runtime,version,package,package_version,path,checksum,os,notes", "runtime,package", "environment_record", "medium", "code/R/02_io_lock_metadata.R", "write_audit_session_info", "yes"
  )
}

pv_audit_contract_row <- function(export_contract_id) {
  registry <- pv_audit_export_registry()
  row <- registry[registry$export_contract_id == export_contract_id, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("Unknown audit export contract id: ", export_contract_id, call. = FALSE)
  }
  row
}

pv_audit_split_contract_fields <- function(value) {
  fields <- strsplit(as.character(value), ",", fixed = TRUE)[[1]]
  trimws(fields[nzchar(trimws(fields))])
}

pv_audit_required_schema <- function(export_contract_id) {
  pv_audit_split_contract_fields(pv_audit_contract_row(export_contract_id)$required_columns_or_fields)
}

pv_audit_sort_keys <- function(export_contract_id) {
  pv_audit_split_contract_fields(pv_audit_contract_row(export_contract_id)$required_sort_keys)
}

pv_audit_export_subdir <- function(export_contract_id) {
  filename <- pv_audit_contract_row(export_contract_id)$future_export_filename
  if (grepl("/tables/", filename, fixed = TRUE)) return("tables")
  if (grepl("/intermediate/", filename, fixed = TRUE)) return("intermediate")
  if (grepl("/manifests/", filename, fixed = TRUE)) return("manifests")
  "root"
}

pv_audit_validate_export_schema <- function(data, export_contract_id) {
  data <- tibble::as_tibble(data)
  required <- pv_audit_required_schema(export_contract_id)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(
      "Audit export ", export_contract_id, " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  data[, required, drop = FALSE]
}

pv_audit_sort_export <- function(data, export_contract_id) {
  keys <- pv_audit_sort_keys(export_contract_id)
  keys <- keys[keys %in% names(data)]
  if (length(keys) == 0L) {
    return(data)
  }
  data <- as.data.frame(data)
  order_args <- c(data[keys], list(na.last = TRUE))
  data[do.call(order, order_args), , drop = FALSE]
}

pv_audit_record_contract_export <- function(
  export_contract_id,
  file,
  n_rows,
  n_cols,
  status,
  notes = NA_character_,
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_intermediates = WRITE_PV_AUDIT_INTERMEDIATES
) {
  if (!pv_audit_enabled(audit_mode, write_intermediates)) {
    return(invisible(NULL))
  }
  if (is.null(audit_paths)) {
    audit_paths <- build_pv_audit_paths(create = TRUE)
  }
  row <- pv_audit_contract_row(export_contract_id)
  manifest_row <- tibble::tibble(
    audit_run_id = audit_paths$run_id,
    export_contract_id = export_contract_id,
    target_id = row$target_id,
    file = normalizePath(file, winslash = "/", mustWork = FALSE),
    n_rows = n_rows,
    n_cols = n_cols,
    status = status,
    recorded_time = audit_now_utc_string(),
    notes = notes
  )
  append_audit_csv_row(manifest_row, audit_paths$audit_contract_export_manifest)
  invisible(file)
}

pv_audit_write_contract_csv <- function(
  data,
  export_contract_id,
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_intermediates = WRITE_PV_AUDIT_INTERMEDIATES,
  notes = NA_character_
) {
  if (!pv_audit_enabled(audit_mode, write_intermediates)) {
    return(invisible(NULL))
  }
  if (is.null(audit_paths)) {
    audit_paths <- build_pv_audit_paths(create = TRUE)
  }
  row <- pv_audit_contract_row(export_contract_id)
  aligned <- pv_audit_validate_export_schema(data, export_contract_id)
  aligned <- pv_audit_sort_export(aligned, export_contract_id)
  filename <- basename(row$future_export_filename)
  subdir <- pv_audit_export_subdir(export_contract_id)
  path <- write_audit_csv(
    aligned,
    filename,
    subdir = subdir,
    object_name = export_contract_id,
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_intermediates,
    notes = notes
  )
  pv_audit_record_contract_export(
    export_contract_id,
    path,
    nrow(aligned),
    ncol(aligned),
    status = "written",
    notes = notes,
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_intermediates
  )
  invisible(path)
}

write_audit_session_info <- function(
  audit_paths = NULL,
  package_names = c(
    "readr", "dplyr", "tidyr", "tibble", "purrr", "terra", "broom", "ggplot2",
    "patchwork", "hexbin", "scales", "rlang", "officer", "flextable"
  ),
  audit_mode = RUN_PV_AUDIT_MODE,
  write_intermediates = WRITE_PV_AUDIT_INTERMEDIATES
) {
  if (!isTRUE(audit_mode)) {
    return(invisible(NULL))
  }
  if (is.null(audit_paths)) {
    audit_paths <- build_pv_audit_paths(create = TRUE)
  }

  package_install_path <- function(pkg) {
    path <- tryCatch(find.package(pkg, quiet = TRUE), error = function(e) character(0))
    if (length(path) != 1L || !nzchar(path)) {
      return(NA_character_)
    }
    normalizePath(path, winslash = "/", mustWork = FALSE)
  }

  writeLines(capture.output(utils::sessionInfo()), audit_paths$audit_session_info, useBytes = TRUE)
  package_names <- unique(as.character(package_names))
  package_available <- vapply(package_names, function(pkg) requireNamespace(pkg, quietly = TRUE), logical(1))
  package_versions <- vapply(package_names, function(pkg) {
    if (!package_available[[pkg]]) {
      return(NA_character_)
    }
    as.character(utils::packageVersion(pkg))
  }, character(1))
  package_paths <- vapply(package_names, package_install_path, character(1))
  versions <- tibble::tibble(
    package = package_names,
    available = package_available,
    version = package_versions,
    loaded = paste0("package:", package_names) %in% search(),
    namespace_loaded = package_names %in% loadedNamespaces(),
    library_path_or_na = package_paths
  )
  readr::write_csv(versions, audit_paths$audit_package_versions)

  library_paths <- tibble::tibble(
    library_path_index = seq_along(.libPaths()),
    library_path = normalizePath(.libPaths(), winslash = "/", mustWork = FALSE)
  )
  readr::write_csv(library_paths, audit_paths$audit_library_paths)

  runtime_rows <- tibble::tibble(
    runtime = c("R", "Rscript", "OS", paste0("package:", package_names)),
    version = c(
      paste(R.version$major, R.version$minor, sep = "."),
      paste(R.version$major, R.version$minor, sep = "."),
      Sys.info()[["sysname"]],
      rep(NA_character_, length(package_names))
    ),
    package = c(NA_character_, NA_character_, NA_character_, package_names),
    package_version = c(NA_character_, NA_character_, NA_character_, package_versions),
    path = c(R.home(), Sys.which("Rscript"), NA_character_, package_paths),
    checksum = NA_character_,
    os = c(
      paste(Sys.info()[["sysname"]], Sys.info()[["release"]]),
      paste(Sys.info()[["sysname"]], Sys.info()[["release"]]),
      paste(Sys.info()[["sysname"]], Sys.info()[["release"]]),
      rep(paste(Sys.info()[["sysname"]], Sys.info()[["release"]]), length(package_names))
    ),
    notes = c("R runtime", "Rscript executable", "Operating system", rep("workflow-required package", length(package_names)))
  )
  pv_audit_write_contract_csv(
    runtime_rows,
    "EC-035",
    audit_paths = audit_paths,
    audit_mode = audit_mode,
    write_intermediates = write_intermediates,
    notes = "Stage 7B environment/runtime evidence"
  )

  run_manifest <- tibble::tibble(
    audit_run_id = audit_paths$run_id,
    project_root = project_root,
    source_root = repo_root,
    output_root = audit_paths$root,
    normal_output_root = exact_output_root,
    created_time = audit_now_utc_string(),
    audit_mode = TRUE,
    write_audit_intermediates = TRUE
  )
  readr::write_csv(run_manifest, audit_paths$audit_run_manifest)
  invisible(audit_paths)
}

write_audit_seed_record <- function(
  seed_label,
  seed_value,
  audit_paths = NULL,
  audit_mode = RUN_PV_AUDIT_MODE,
  write_intermediates = WRITE_PV_AUDIT_INTERMEDIATES,
  notes = NA_character_
) {
  if (!isTRUE(audit_mode)) {
    return(invisible(NULL))
  }
  if (is.null(audit_paths)) {
    audit_paths <- build_pv_audit_paths(create = TRUE)
  }
  row <- tibble::tibble(
    audit_run_id = audit_paths$run_id,
    seed_label = seed_label,
    seed_value = seed_value,
    recorded_time = audit_now_utc_string(),
    notes = notes
  )
  append_audit_csv_row(row, audit_paths$audit_seed_records)
  invisible(audit_paths$audit_seed_records)
}

activate_pv_audit_output_dirs <- function(audit_paths) {
  old_dirs <- list(
    out_dir = out_dir,
    fig_dir = fig_dir,
    tab_dir = tab_dir,
    intermediate_dir = intermediate_dir
  )
  out_dir <<- audit_paths$root
  fig_dir <<- audit_paths$figures
  tab_dir <<- audit_paths$tables
  intermediate_dir <<- audit_paths$intermediate
  invisible(old_dirs)
}

restore_pv_audit_output_dirs <- function(old_dirs) {
  if (is.null(old_dirs)) {
    return(invisible(TRUE))
  }
  out_dir <<- old_dirs$out_dir
  fig_dir <<- old_dirs$fig_dir
  tab_dir <<- old_dirs$tab_dir
  intermediate_dir <<- old_dirs$intermediate_dir
  invisible(TRUE)
}

# Output root for the manuscript HDR reanalysis.
manuscript_outputs_root <- file.path(hdr_output_root, "manuscript_tables_figures")
brief_output_dirs <- c(hdr_output_root, manuscript_outputs_root)
workflow_setup_dirs <- c(project_setup_dirs, replication_output_dirs, brief_output_dirs)

path_manifest <- tibble::tibble(
  label = c(
    "project_root", "pellett_valbuena_source_root", "exact_replication_outputs",
    "hdr_reanalysis_outputs", "manuscript_tables_figures", "archive_root",
    "local_scratch"
  ),
  path = c(
    project_root, repo_root, exact_output_root,
    hdr_output_root, manuscript_outputs_root, archive_root,
    scratch_root
  )
)

write_analysis_path_manifest <- function() {
  ensure_dir(manifest_dir)
  readr::write_csv(path_manifest, file.path(manifest_dir, "analysis_path_manifest.csv"))
  invisible(file.path(manifest_dir, "analysis_path_manifest.csv"))
}

write_project_structure_readme <- function() {
  readme_path <- file.path(manifest_dir, "README_project_structure.md")
  lines <- c(
    "# Generated project-folder guide",
    "",
    "The top-level project folder contains archived inputs, R scripts, and generated results. The main manuscript analysis and the optional R reimplementation write their results to separate folders.",
    "",
    "## Input and code folders",
    "",
    "- `data/external/pellett_valbuena_2025/`: archived inputs used by the main manuscript analysis and the optional R reimplementation, stored at the paths required by the scripts.",
    "- `code/`: R scripts for the main manuscript analysis, source-data export, and optional reimplementation.",
    "",
    "## Generated result folders",
    "",
    paste0("- Optional R reimplementation outputs: `", normalizePath(exact_output_root, winslash = "/", mustWork = FALSE), "`"),
    paste0("- Manuscript figures, tables, captions, and supporting analysis files: `", normalizePath(manuscript_outputs_root, winslash = "/", mustWork = FALSE), "`"),
    "",
    "## Generated setup records",
    "",
    "- `metadata/analysis_path_manifest.csv`: an inventory of key project paths used by the workflows.",
    "- `metadata/README_project_structure.md`: this generated guide.",
    "",
    "## Temporary and compatibility folders",
    "",
    "- `_archive_do_not_release/`: created as an empty compatibility folder; the supported workflows do not write results there.",
    "- `_local_scratch/`: temporary files used while an analysis is running, including a lock file that prevents simultaneous runs.",
    "",
    "## Analysis scope",
    "",
    "The manuscript analysis uses δ, defined as sample variance divided by the mean, for the heterogeneity–diversity relationship comparisons. Additional metric–mean diagnostics are generated only by the optional R reimplementation.",
    "",
    "## File-inventory note",
    "",
    "The generated CSV inventory retains stable machine identifiers such as `exact_replication_outputs`. These legacy names allow the scripts and checks to refer to the same paths consistently; `PELLETT_VALBUENA_WORKFLOW.md` defines the supported workflow scope.",
    "",
    "Keep `data/external/pellett_valbuena_2025/` read-only. Supported workflows write results to the output folders above, setup records to `metadata/`, and temporary run files to `_local_scratch/`."
  )
  writeLines(lines, readme_path)
  invisible(readme_path)
}

write_workflow_setup_metadata <- function() {
  write_analysis_path_manifest()
  write_project_structure_readme()
  invisible(TRUE)
}

prepare_workflow_setup <- function(write_metadata = WRITE_SETUP_METADATA, dirs = workflow_setup_dirs) {
  ensure_dir(dirs)
  if (isTRUE(write_metadata)) {
    write_workflow_setup_metadata()
  }
  invisible(TRUE)
}

# Shared path validation and numeric helpers.
check_repo_root <- function(root = repo_root) {
  required <- c("01_concept", "02_MBH_comb", "03_corrected_HDR", "04_som", "05_review")
  missing <- required[!dir.exists(file.path(root, required))]
  if (length(missing) > 0) {
    stop(
      "The archived Pellett & Valbuena input folder is incomplete. Missing: ",
      paste(missing, collapse = ", "),
      "\nChecked source root: ", root,
      "\nThe repository root is: ", project_root
    )
  }
  invisible(TRUE)
}

safe_file <- function(...) file.path(...)

rescale01 <- function(x) {
  x <- as.numeric(x)
  mn <- min(x, na.rm = TRUE)
  mx <- max(x, na.rm = TRUE)
  if (!is.finite(mn) || !is.finite(mx) || mx == mn) return(rep(NA_real_, length(x)))
  (x - mn) / (mx - mn)
}
save_png_path <- function(plot, path, width = 8, height = 5, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  suppressWarnings(
    ggplot2::ggsave(path, plot, width = width, height = height, device = "png", dpi = dpi, bg = "white")
  )
}

