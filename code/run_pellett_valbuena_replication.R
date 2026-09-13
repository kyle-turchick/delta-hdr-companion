# Runs the optional R reimplementation of selected Pellett & Valbuena analyses.
# From the repository root, run:
#   Rscript --vanilla code/run_pellett_valbuena_replication.R

RUN_WORKFLOW_ON_SOURCE <- get0("RUN_WORKFLOW_ON_SOURCE", ifnotfound = TRUE, inherits = FALSE)
INSTALL_MISSING <- get0("INSTALL_MISSING", ifnotfound = FALSE, inherits = FALSE)
WRITE_SETUP_METADATA <- get0("WRITE_SETUP_METADATA", ifnotfound = isTRUE(RUN_WORKFLOW_ON_SOURCE), inherits = FALSE)

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
source(file.path(module_root, "08_validation_helpers.R"))

run_pellett_valbuena_replication <- function(root = repo_root) {
  acquire_workflow_lock()
  workflow_completed <- FALSE
  on.exit({
    if (isTRUE(workflow_completed)) {
      release_workflow_lock()
    }
  }, add = TRUE)

  prepare_workflow_setup(dirs = c(project_setup_dirs, replication_output_dirs))
  result <- run_all_R_replication(root)
  workflow_completed <- TRUE
  result
}

if (isTRUE(RUN_WORKFLOW_ON_SOURCE)) {
  workflow_results_pellett_valbuena_replication <- run_pellett_valbuena_replication()
}
