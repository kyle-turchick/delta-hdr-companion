# Configures workflow controls, finds the repository root, and defines paths.

# 00 PROJECT CONFIGURATION, PATHS, PACKAGES ####
# ============================================================ #
# 00A user toggles
# A caller can override these defaults before sourcing this file, for example:
#   RUN_EXACT_REPLICATION <- TRUE

RUN_WORKFLOW_ON_SOURCE <- get0("RUN_WORKFLOW_ON_SOURCE", ifnotfound = TRUE, inherits = FALSE)
INSTALL_MISSING <- get0("INSTALL_MISSING", ifnotfound = FALSE, inherits = FALSE)
# The legacy name RUN_EXACT_REPLICATION is retained for compatibility. It is
# FALSE by default; TRUE runs the optional Pellett & Valbuena reimplementation
# before the main manuscript analysis.
RUN_EXACT_REPLICATION <- get0("RUN_EXACT_REPLICATION", ifnotfound = FALSE, inherits = FALSE)
RUN_MANUSCRIPT_HIERARCHY <- get0("RUN_MANUSCRIPT_HIERARCHY", ifnotfound = TRUE, inherits = FALSE)
RUN_MANUSCRIPT_DELIVERABLES <- get0("RUN_MANUSCRIPT_DELIVERABLES", ifnotfound = TRUE, inherits = FALSE)
WRITE_SETUP_METADATA <- get0("WRITE_SETUP_METADATA", ifnotfound = isTRUE(RUN_WORKFLOW_ON_SOURCE), inherits = FALSE)

# These advanced controls write comparison evidence only when audit mode is enabled.
RUN_PV_AUDIT_MODE <- get0("RUN_PV_AUDIT_MODE", ifnotfound = FALSE, inherits = FALSE)
PV_AUDIT_OUTPUT_ROOT <- get0("PV_AUDIT_OUTPUT_ROOT", ifnotfound = NULL, inherits = FALSE)
WRITE_PV_AUDIT_INTERMEDIATES <- get0("WRITE_PV_AUDIT_INTERMEDIATES", ifnotfound = FALSE, inherits = FALSE)
PV_AUDIT_RUN_ID <- get0("PV_AUDIT_RUN_ID", ifnotfound = NULL, inherits = FALSE)

# This setup does not delete existing outputs from the manuscript analysis.

RUN_EXPENSIVE_GLOBAL_RESAMPLING <- get0("RUN_EXPENSIVE_GLOBAL_RESAMPLING", ifnotfound = FALSE, inherits = FALSE)
GLOBAL_RESAMPLING_N <- get0("GLOBAL_RESAMPLING_N", ifnotfound = 10000000L, inherits = FALSE)
REVIEW2_SIM_N <- get0("REVIEW2_SIM_N", ifnotfound = 1000000L, inherits = FALSE)

# 00B project and source paths
# The repository root is the top-level folder containing .mi_hdr_project_root.
project_root_marker <- ".mi_hdr_project_root"

discover_project_root <- function(start = getwd(), marker = project_root_marker) {
  env_root <- Sys.getenv("MI_HDR_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) {
    if (!dir.exists(env_root)) {
      stop("MI_HDR_PROJECT_ROOT is set but does not exist: ", env_root)
    }
    root <- normalizePath(env_root, winslash = "/", mustWork = TRUE)
    if (!file.exists(file.path(root, marker))) {
      stop(
        "MI_HDR_PROJECT_ROOT is set to: ", root,
        "\nExpected marker file not found: ", marker
      )
    }
    return(root)
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

  stop(
    "Could not find the repository root. Set MI_HDR_PROJECT_ROOT to the ",
    "top-level repository folder, or run this command from inside it. ",
    "Expected marker file: ", marker
  )
}

project_root <- discover_project_root()

# This archived Pellett & Valbuena input folder remains read-only.
pv_source_rel_path <- file.path("data", "external", "pellett_valbuena_2025")
pv_source_root <- file.path(project_root, pv_source_rel_path)
if (!all(dir.exists(file.path(pv_source_root, c("01_concept", "02_MBH_comb", "03_corrected_HDR", "04_som", "05_review"))))) {
  stop(
    "Could not find the archived Pellett & Valbuena input folders under: ",
    file.path(project_root, pv_source_rel_path),
    "\nThe source folder should contain 01_concept, 02_MBH_comb, 03_corrected_HDR, 04_som, and 05_review."
  )
}
repo_root <- normalizePath(pv_source_root, winslash = "/", mustWork = TRUE)
# The legacy name repo_root refers to the archived Pellett & Valbuena input
# folder, not the repository's top-level folder.

# Output folder variable names remain stable for functions that use tab_dir,
# fig_dir, and intermediate_dir.
manifest_dir <- file.path(project_root, "metadata")
script_dir <- file.path(project_root, "code")
exact_output_root <- file.path(project_root, "provenance", "pellett_valbuena_replication_outputs")
hdr_output_root <- file.path(project_root, "results", "delta_hdr_reanalysis")
archive_root <- file.path(project_root, "_archive_do_not_release")
scratch_root <- file.path(project_root, "_local_scratch")
workflow_lock_path <- file.path(scratch_root, ".workflow_run_lock")
