# Internal runner for the main analysis, called by code/run_reanalysis.R.
# Readers should use the documented code/run_reanalysis.R command.

RUN_WORKFLOW_ON_SOURCE <- get0("RUN_WORKFLOW_ON_SOURCE", ifnotfound = TRUE, inherits = FALSE)
INSTALL_MISSING <- get0("INSTALL_MISSING", ifnotfound = FALSE, inherits = FALSE)
RUN_EXACT_REPLICATION <- get0("RUN_EXACT_REPLICATION", ifnotfound = FALSE, inherits = FALSE)
RUN_MANUSCRIPT_HIERARCHY <- get0("RUN_MANUSCRIPT_HIERARCHY", ifnotfound = TRUE, inherits = FALSE)
RUN_MANUSCRIPT_DELIVERABLES <- get0("RUN_MANUSCRIPT_DELIVERABLES", ifnotfound = TRUE, inherits = FALSE)
WRITE_SETUP_METADATA <- get0("WRITE_SETUP_METADATA", ifnotfound = isTRUE(RUN_WORKFLOW_ON_SOURCE), inherits = FALSE)
WRITE_SOURCE_PREP_OUTPUTS <- get0("WRITE_SOURCE_PREP_OUTPUTS", ifnotfound = FALSE, inherits = FALSE)

find_project_root_for_modules <- function(start = getwd(), marker = ".mi_hdr_project_root") {
  env_root <- Sys.getenv("MI_HDR_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = TRUE))
  }
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, marker))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }
  stop("Could not find the repository root containing marker: ", marker)
}

module_root <- file.path(find_project_root_for_modules(), "code", "R")
source(file.path(module_root, "00_setup_paths.R"))
source(file.path(module_root, "01_packages.R"))
source(file.path(module_root, "02_io_lock_metadata.R"))
source(file.path(module_root, "03_source_data_prep.R"))
source(file.path(module_root, "04_pellett_valbuena_replication.R"))
source(file.path(module_root, "05_delta_hdr_model_helpers.R"))
source(file.path(module_root, "06_brief_analysis_hierarchy.R"))
source(file.path(module_root, "07_brief_deliverables.R"))
source(file.path(module_root, "08_validation_helpers.R"))

run_manuscript_input_preparation <- function(root = repo_root, write_source_prep_outputs = WRITE_SOURCE_PREP_OUTPUTS) {
  check_repo_root(root)
  message("Preparing the manuscript-analysis inputs from archived files under: ", root)

  adf <- prepare_allouche_data_R(root, write_outputs = write_source_prep_outputs)
  allouche <- run_allouche_analysis_R(adf, write_outputs = write_source_prep_outputs)
  macarthur <- run_macarthur_analysis_R(root, write_outputs = write_source_prep_outputs)

  list(
    allouche = allouche,
    macarthur = macarthur,
    provenance = tibble::tibble(
      component = c("Topography inputs (Allouche et al. 2012)", "Foliage inputs (MacArthur & MacArthur 1961)"),
      status = c("prepared", "prepared")
    )
  )
}

run_delta_hdr_brief_reanalysis <- function(
    root = repo_root,
    run_exact_replication = RUN_EXACT_REPLICATION,
    run_manuscript_hierarchy = RUN_MANUSCRIPT_HIERARCHY,
    run_manuscript_deliverables = RUN_MANUSCRIPT_DELIVERABLES,
    write_source_prep_outputs = WRITE_SOURCE_PREP_OUTPUTS) {

  acquire_workflow_lock()
  workflow_completed <- FALSE
  on.exit({
    if (isTRUE(workflow_completed)) {
      release_workflow_lock()
    }
  }, add = TRUE)

  setup_dirs <- c(project_setup_dirs, brief_output_dirs)
  if (isTRUE(run_exact_replication)) {
    setup_dirs <- c(setup_dirs, replication_output_dirs)
  }
  prepare_workflow_setup(dirs = setup_dirs)

  if (isTRUE(run_exact_replication)) {
    message("Running the optional R reimplementation of selected Pellett & Valbuena analyses before the main manuscript analysis...")
    results_full_R_local <- run_all_R_replication(root)
  } else {
    message("The optional R reimplementation is disabled; running the main manuscript analysis only...")
    results_full_R_local <- run_manuscript_input_preparation(
      root,
      write_source_prep_outputs = write_source_prep_outputs
    )
  }

  assign("results_full_R", results_full_R_local, envir = globalenv())

  if (isTRUE(run_manuscript_hierarchy)) {
    results_manuscript_analysis_hierarchy_local <- run_manuscript_analysis_hierarchy()
    assign("results_manuscript_analysis_hierarchy", results_manuscript_analysis_hierarchy_local, envir = globalenv())
  } else {
    results_manuscript_analysis_hierarchy_local <- NULL
  }

  if (isTRUE(run_manuscript_deliverables)) {
    if (is.null(results_manuscript_analysis_hierarchy_local)) {
      stop("RUN_MANUSCRIPT_DELIVERABLES = TRUE requires RUN_MANUSCRIPT_HIERARCHY = TRUE.")
    }
    results_manuscript_named_deliverables_local <- run_manuscript_deliverable_workflow()
    assign("results_manuscript_named_deliverables", results_manuscript_named_deliverables_local, envir = globalenv())
  } else {
    results_manuscript_named_deliverables_local <- NULL
  }

  run_manifest <- tibble::tibble(
    component = c("input_preparation", "analysis_hierarchy", "named_deliverables"),
    status = c(
      "completed",
      ifelse(is.null(results_manuscript_analysis_hierarchy_local), "skipped", "completed"),
      ifelse(is.null(results_manuscript_named_deliverables_local), "skipped", "completed")
    ),
    output_root = c(hdr_output_root, analysis_hierarchy_root, manuscript_outputs_root)
  )
  readr::write_csv(run_manifest, file.path(manifest_dir, "split_brief_run_manifest.csv"))

  message("Main manuscript workflow complete.")
  message("Manuscript-analysis output folder: ", manuscript_outputs_root)

  result <- list(
    results_full_R = results_full_R_local,
    results_manuscript_analysis_hierarchy = results_manuscript_analysis_hierarchy_local,
    results_manuscript_named_deliverables = results_manuscript_named_deliverables_local,
    run_manifest = run_manifest
  )
  workflow_completed <- TRUE
  result
}

if (isTRUE(RUN_WORKFLOW_ON_SOURCE)) {
  workflow_results_delta_hdr_brief <- run_delta_hdr_brief_reanalysis()
}
