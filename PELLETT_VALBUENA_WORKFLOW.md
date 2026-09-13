# Optional Pellett & Valbuena workflow

This optional workflow reimplements selected procedures from Pellett & Valbuena (2025) in R. A workflow is an ordered set of steps that reads inputs, performs calculations, and creates outputs. An R reimplementation means that selected procedures originally written in another programming language were written again in R. The corresponding Pellett & Valbuena workflow was written in Julia.

The present manuscript reanalyzes the MacArthur & MacArthur and Allouche et al. datasets to compare possible shapes of the relationship between environmental heterogeneity and bird diversity. The main workflow recreates the present manuscript's computational analyses and five included products: Figure 1, Extended Data Figs. 1–2, and Extended Data Tables 1–2. Pellett & Valbuena’s work helps explain the methods and where the archived data came from.

The optional workflow creates a separate collection of R tables and figures for conceptual examples and analyses based on the MacArthur and Catalonia data. Selected results from this R workflow were also compared with corresponding Julia results.

## Why this workflow is optional

The main workflow alone recreates these computational analyses and five included products. You do not need to run the optional workflow to reproduce or understand those results.

The optional workflow provides more detail about the methods, the origins of the archived data, and selected comparisons between R and Julia. It does not replace or alter the main manuscript analysis, and it does not run Julia. Its products are written separately under `provenance/pellett_valbuena_replication_outputs/`. Main Table 1 and Extended Data Table 3 are interpretive summaries assembled for the manuscript, not outputs of either workflow.

## Requirements

Install the packages listed in the [main run instructions](RUN_INSTRUCTIONS.md). Automatic package installation is disabled. The [environment guide](environment/README.md) records the tested software environment: the R version, package versions, and operating system used for validation.

The project already contains all archived inputs required by the optional workflow. It uses the 402 inputs shared with the main workflow and these five additional files:

| Path | Role |
|---|---|
| `data/external/pellett_valbuena_2025/02_MBH_comb/data/elevation/heterogeneity_inds.csv` | Derived elevation summary used in calculations |
| `data/external/pellett_valbuena_2025/02_MBH_comb/data/crop_cover/heterogeneity_inds.csv` | Derived crop cover summary used in calculations |
| `data/external/pellett_valbuena_2025/02_MBH_comb/data/elevation/lowres_full_world.tif` | Global elevation map context |
| `data/external/pellett_valbuena_2025/02_MBH_comb/data/crop_cover/lowres_crop.tif` | Global crop cover map context |
| `data/external/pellett_valbuena_2025/03_corrected_HDR/empirical/01_carnicer/data/elevation/catalunya_lowres.tif` | Catalonia elevation map context |

The two CSV files are archived calculation summaries. A CSV file is a table stored as plain text. The three TIFF files store gridded map data that provide geographic context. The supported workflow reads the archived summaries instead of regenerating them from the original global map files.

[`data/INPUT_MANIFEST.csv`](data/INPUT_MANIFEST.csv) is a manifest, or structured file list, for all 407 archived inputs. It records each file’s path, byte size, and SHA-256 file fingerprint. A SHA-256 fingerprint is a digital identifier calculated from a file’s bytes to confirm exact file identity. The included inputs were verified against the cited deposits; readers do not need to calculate these fingerprints before running the workflow. A matching fingerprint verifies file identity, not agreement between scientific analyses.

Run only one repository workflow at a time. A workflow lock is a temporary file that prevents two repository workflows from running simultaneously. This workflow creates the lock at `_local_scratch/.workflow_run_lock`. A successful run removes it when finished.

## Run the optional workflow

### Option 1: command line

A command line is an interface in which you type commands. Examples include Terminal, Command Prompt, and PowerShell.

Install R. Then open one terminal application and run:

```sh
Rscript --version
```

If this prints an R version, `Rscript` is available. If the command is not recognized, follow the command line setup instructions for your R installation.

In the same terminal, change to the extracted main project folder. This folder contains `README.md`, `.mi_hdr_project_root`, `code/`, `data/`, `environment/`, and `results/`. The command `cd` means “change directory.” Replace the example path below with the project’s actual path:

```sh
cd "/full/path/to/project folder"
```

Quotation marks allow the path to contain spaces. In Windows Command Prompt, use `cd /d` instead of `cd` if the project folder is on a different drive.

From the main project folder, run:

```sh
Rscript --vanilla code/run_pellett_valbuena_replication.R
```

`Rscript` runs an R script from the command line. `--vanilla` starts R without loading personal startup files or a saved workspace.

### Option 2: R or RStudio Console

The working directory is the folder where R looks for the project files. Set it to the extracted main project folder containing `README.md`, `.mi_hdr_project_root`, `code/`, `data/`, `environment/`, and `results/`. The file `.mi_hdr_project_root` is an empty marker file used to confirm that R is working in the main project folder.

In RStudio, choose **Session > Set Working Directory > Choose Directory...** and select that folder. In another R interface, use its working directory control or run the following command with the folder’s actual path:

```r
setwd("/full/path/to/project folder")
```

On Windows, you can use forward slashes in a path entered in R, as shown above.

Then paste this complete block into the Console:

```r
launch_delta_hdr_optional_workflow <- function() {
  if (!file.exists(".mi_hdr_project_root")) {
    stop("Set the working directory to the main project folder first.")
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

  if (!file.exists(rscript)) {
    stop("Could not find Rscript in the active R installation: ", rscript)
  }

  workflow_script <- file.path(
    "code",
    "run_pellett_valbuena_replication.R"
  )

  if (!file.exists(workflow_script)) {
    stop("Could not find the optional workflow script: ", workflow_script)
  }

  message("Using project folder: ", project_root)

  status <- system2(
    rscript,
    c("--vanilla", workflow_script)
  )

  as.integer(status)
}

status <- launch_delta_hdr_optional_workflow()

if (status != 0L) {
  stop("The optional workflow ended with exit status ", status, ".")
}
```

The marker check detects an error caused by selecting the wrong folder before the workflow starts. The block identifies the selected project folder, finds `Rscript` or `Rscript.exe` in the active R installation, and starts `code/run_pellett_valbuena_replication.R` in a new, independent R session using `--vanilla`. It temporarily supplies the full project path to that session and restores any previous `MI_HDR_PROJECT_ROOT` value afterward. The quoted `setwd()` path and the relative workflow script path allow the project folder’s path to contain spaces.

Do not use `source("code/run_pellett_valbuena_replication.R")` as an equivalent shortcut. `source()` would run the script inside the current R session, where data or settings already loaded in that session could influence the workflow. Starting a new R session provides a cleaner and more repeatable route.

## Successful completion and outputs

An exit status is the number returned when a process ends; exit status 0 means that no execution error was reported.

A successful run:

- returns exit status 0;
- prints messages beginning `Optional tables written to:` and `Optional figures written to:`;
- creates the expected output groups; and
- removes `_local_scratch/.workflow_run_lock` when it finishes.

The run creates 59 files in the reader’s project copy. These generated files are not included among the files distributed with the companion.

| Location | Contents |
|---|---|
| `provenance/pellett_valbuena_replication_outputs/tables/` | 14 result tables in CSV format |
| `provenance/pellett_valbuena_replication_outputs/intermediate/` | Four supporting CSV files used to construct results |
| `provenance/pellett_valbuena_replication_outputs/figures/` | 38 PNG figures |
| `provenance/pellett_valbuena_replication_outputs/R_replication_manifest.csv` | One CSV file recording the status of selected parts of the R reimplementation |
| `metadata/analysis_path_manifest.csv` | Key project paths used during the run |
| `metadata/README_project_structure.md` | A generated guide to project folders |

The first four locations together contain 19 CSV files and 38 PNG files. The two files under `metadata/` are shared setup records. These products are separate from the manuscript figures, tables, and source data files.

A run may print an aggregate notice that warnings occurred without printing each warning. A warning count alone is not a success criterion. Evaluate the run using exit status 0, both completion messages, the expected output groups, and removal of the workflow lock. Preserve and review any specific warning or error text that is actually printed. Treat a nonzero exit status, either missing completion message, a missing expected output group, or a retained workflow lock as a problem.

## What validation established

### Repeatability between R runs

Using the recorded R and package versions, all 19 CSV outputs matched exactly across the validated R runs, and 32 of the 38 PNG figures matched exactly. The other six were produced by two figure routines that randomly select displayed points and do not force the same selection on every run. They retained their dimensions, structure, panels, labels, overall distribution of displayed points, and curve shapes. Review found no substantive visual difference.

### Comparisons between R and Julia

Separate checks paired selected R outputs with corresponding Julia outputs created from the Julia files preserved in Pellett & Valbuena’s Zenodo version 6 deposit. For selected paired outputs, numerical results that did not depend on random sampling agreed within the predefined tolerance for each comparison. A numerical tolerance is a small permitted difference used when software implementations may calculate or store values slightly differently.

Outputs based on simulation were assessed for similarity rather than identical random draws. Figures were assessed with criteria defined separately for each comparison.

Together, these checks support the selected R reimplementation. They do not cover every Julia procedure, intermediate file, or figure.

## Scientific and technical boundaries

The optional workflow uses the archived 10,000 × 13 MacArthur foliage table as a fixed input. The data guide explains the table’s documented source and role in the workflow. See [`data/README.md`](data/README.md).

Some generated filenames and fields use identifiers from the original workflow because the software depends on their exact names. In those identifiers, `exact` labels an export intended for comparison with Julia. It does not mean that the complete Julia workflow was reproduced with identical file bytes.

## Troubleshooting

If the workflow reports that it cannot find `.mi_hdr_project_root`, select the extracted main project folder as the working directory and try again.

If `Rscript` or a required package is missing, install or configure the named requirement using the [run instructions](RUN_INSTRUCTIONS.md) and [environment guide](environment/README.md).

If the workflow reports that an included input is missing, confirm that the project was fully extracted. Obtain a fresh project copy if necessary. Do not edit an archived input or analytical script merely to make the workflow continue.

The workflow lock prevents simultaneous repository workflows. If `_local_scratch/.workflow_run_lock` remains after an interrupted run, first confirm that no R workflow is active. Remove the retained lock only after that confirmation.

If the cause of a failure is unclear, preserve the complete Console or terminal output when seeking help.

## Data sources, attribution, and reuse

The five inputs used only by the optional workflow are redistributed unchanged from the Pellett & Valbuena Zenodo deposits and are listed in [`data/INPUT_MANIFEST.csv`](data/INPUT_MANIFEST.csv). See [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for the versioned deposit records, source attribution, and reuse terms.

The project license does not replace the separate reuse terms that apply to the 407 archived inputs. Consult the project [`LICENSE`](LICENSE) to identify the terms that apply to each file group.
