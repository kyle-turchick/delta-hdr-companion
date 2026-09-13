# Default entry point for the main manuscript analysis. The optional Pellett &
# Valbuena R reimplementation remains disabled unless RUN_EXACT_REPLICATION is
# set to TRUE; it can also be run with code/run_pellett_valbuena_replication.R.

RUN_WORKFLOW_ON_SOURCE <- get0("RUN_WORKFLOW_ON_SOURCE", ifnotfound = TRUE, inherits = FALSE)
RUN_EXACT_REPLICATION <- get0("RUN_EXACT_REPLICATION", ifnotfound = FALSE, inherits = FALSE)
RUN_MANUSCRIPT_HIERARCHY <- get0("RUN_MANUSCRIPT_HIERARCHY", ifnotfound = TRUE, inherits = FALSE)
RUN_MANUSCRIPT_DELIVERABLES <- get0("RUN_MANUSCRIPT_DELIVERABLES", ifnotfound = TRUE, inherits = FALSE)
INSTALL_MISSING <- get0("INSTALL_MISSING", ifnotfound = FALSE, inherits = FALSE)

source(file.path("code", "run_delta_hdr_brief_reanalysis.R"))
