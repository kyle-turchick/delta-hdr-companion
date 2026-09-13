# Run instructions

## Requirements

The [environment guide](environment/README.md) records R 4.6.0, the add-on package versions, and the operating system used for validation. Using the recorded R and package versions provides the closest match to those runs, although other versions may also work.

Install [R from the Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/), then open its interactive console. At the R prompt, shown as `>`, run:

```r
install.packages(c("broom", "dplyr", "flextable", "ggplot2", "hexbin", "officer", "patchwork", "purrr", "readr", "rlang", "scales", "terra", "tibble", "tidyr"), repos = "https://cloud.r-project.org")
```

An R package is an add-on library used by the analysis. This command obtains packages from CRAN. It installs the versions currently available there, which may differ from the tested versions. The workflow checks for required packages but does not install them automatically.

## Run the main manuscript analysis

The main manuscript workflow uses the 402 inputs marked `canonical_and_optional` in [`data/INPUT_MANIFEST.csv`](data/INPUT_MANIFEST.csv). That stable machine value means that the inputs are shared by the main and optional workflows.

### Option 1: command line

The command line is an interface where you type commands. Open a terminal application, such as Terminal on macOS, Command Prompt, or PowerShell on Windows, and change to the repository root: the extracted main project folder containing `README.md`, `code/`, `data/`, `environment/`, and `results/`. The command `cd` means "change directory." For example, use `cd "/full/path/to/extracted-project-folder"` with the folder's actual path. In that terminal window, run `Rscript --version`. If it prints an R version, `Rscript` is available. If the command is not recognized, follow the command line setup instructions for your R installation before continuing. Then run:

```sh
Rscript --vanilla code/run_reanalysis.R
```

`Rscript` runs an R script from the command line. `--vanilla` starts R without loading user startup files or a saved workspace.

### Option 2: R or RStudio Console

If you prefer R or RStudio, use the Console, where you enter R commands. The working directory is the folder where R looks for project files. Set it to the repository root: the extracted main project folder containing `README.md`, `.mi_hdr_project_root`, `code/`, `data/`, `environment/`, and `results/`. In RStudio, choose **Session > Set Working Directory > Choose Directory...** and select that folder. In another R interface, use its control for setting the working directory or run `setwd("/full/path/to/extracted-project-folder")` with the folder's actual path. Then paste this complete block into the Console:

```r
launch_delta_hdr_main_workflow <- function() {
  if (!file.exists(".mi_hdr_project_root")) {
    stop("Set the working directory to the top-level project folder first.")
  }

  project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
  previous_root <- Sys.getenv("MI_HDR_PROJECT_ROOT", unset = NA_character_)

  on.exit({
    if (is.na(previous_root)) {
      Sys.unsetenv("MI_HDR_PROJECT_ROOT")
    } else {
      Sys.setenv(MI_HDR_PROJECT_ROOT = previous_root)
    }
  }, add = TRUE)

  Sys.setenv(MI_HDR_PROJECT_ROOT = project_root)

  rscript <- file.path(
    R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  )

  system2(
    rscript,
    c("--vanilla", file.path("code", "run_reanalysis.R"))
  )
}

status <- launch_delta_hdr_main_workflow()

if (status != 0L) {
  stop("The workflow ended with exit status ", status, ".")
}
```

The marker check detects a wrong working directory before the workflow starts. The block finds the `Rscript` program in the active R installation and starts the same `code/run_reanalysis.R` entry point in a separate `--vanilla` process. It temporarily assigns the selected repository root for the child process and restores any previous `MI_HDR_PROJECT_ROOT` value afterward. Objects and control variables in the open Console therefore do not change which workflow components run.

Do not use `source("code/run_reanalysis.R")` as an equivalent shortcut. `source()` evaluates the entry script inside the current interactive workspace, where existing variables can enable the optional workflow, skip required components, or change how packages are installed.

Either launch route prepares the archived inputs and runs five prespecified analysis sets, each comparing linear, logarithmic, and quadratic candidate models. It creates manuscript figures, tables, captions, and supporting analysis files beneath `results/delta_hdr_reanalysis/manuscript_tables_figures/`. The five included figures and tables are listed in the [results guide](results/README.md).

A successful run finishes without an execution error (exit status 0), removes `_local_scratch/.workflow_run_lock`, prints `Main manuscript workflow complete.`, and reports the output folder for the manuscript analysis. A run may finish with warnings. Because warning counts can vary across software environments, review the warning text together with the completion message and expected outputs. Treat warnings about missing inputs or packages, failed calculations, model convergence, or files that could not be written as problems to resolve.

The lock file prevents simultaneous repository workflows. Run only one workflow at a time. If `_local_scratch/.workflow_run_lock` remains after an interrupted run, inspect it and confirm that no R workflow is active before removing it manually.

If the workflow stops before the completion message, read the final error and keep the complete terminal output. Correct only the reported setup problem. For example, install a named required package, return to the repository root, or obtain a fresh project copy if an included input is missing. Do not edit archived inputs or analysis code simply to make the run continue. If the cause is unclear, include the complete terminal output when seeking help. If the workflow lock remains, follow the instructions above for checking the lock file.

## Export figure and table source data

The repository already includes five CSV files containing source data: structured tables containing values for Figure 1b–d, Extended Data Fig. 1, Extended Data Fig. 2, Extended Data Table 1, and Extended Data Table 2. Figure 1a is a conceptual hypothesis table rather than a direct display of the CSV values. The included CSV files can be inspected without rerunning an analysis.

To recreate all five files, run this separate command from the repository root:

```sh
Rscript --vanilla code/export_figure_table_source_data.R
```

The exporter prepares the five prespecified analysis sets, validates a complete temporary copy, and replaces all five CSV files together beneath [`results/delta_hdr_reanalysis/source_data/`](results/delta_hdr_reanalysis/source_data/). If replacement fails, it restores the prior complete set. It does not write the manuscript figures and tables or products from the optional workflow. The [source data guide](results/delta_hdr_reanalysis/source_data/README.md) maps the files to the empirical and illustrative figure panels and Extended Data tables and defines every field and unit. The exporter uses the same lock that prevents simultaneous workflows.

A successful export returns exit status 0, removes the workflow lock, and prints `Figure and table source-data export complete. Files:` followed by the source data folder. Review any unexpected warnings.

## Run the optional Pellett & Valbuena workflow

The main manuscript analysis does not require the optional workflow. In this separate workflow, selected Pellett & Valbuena procedures originally written in Julia were written again in R. The [optional workflow guide](PELLETT_VALBUENA_WORKFLOW.md) is the authoritative source for its scientific purpose, five additional archived inputs, ways to start it from the command line or R/RStudio Console, signs of successful completion, output locations, comparison scope, troubleshooting, attribution, and reuse.
