# Code guide

This folder contains the R scripts that run the main manuscript analysis, export the values underlying its figures and tables, and run the optional workflow that rewrites selected Pellett & Valbuena procedures in R.

## Run the main manuscript analysis

From the repository root (the main project folder containing `README.md` and `.mi_hdr_project_root`, an empty marker file that helps the scripts locate this folder), run:

```sh
Rscript --vanilla code/run_reanalysis.R
```

See the [run instructions](../RUN_INSTRUCTIONS.md) for R and package requirements.

## Export figure and table source data

Run the source data exporter separately from the main manuscript workflow and the optional workflow:

```sh
Rscript --vanilla code/export_figure_table_source_data.R
```

It creates the five structured CSV files containing values for the empirical and illustrative figure panels and Extended Data Tables 1–2, as described in the [source data guide](../results/delta_hdr_reanalysis/source_data/README.md). It does not write the main manuscript figures and tables or products from the optional workflow.

## Main scripts

- [`run_reanalysis.R`](run_reanalysis.R) starts the supported manuscript analysis.
- [`export_figure_table_source_data.R`](export_figure_table_source_data.R) creates the separate source data package for figures and tables.
- [`run_delta_hdr_brief_reanalysis.R`](run_delta_hdr_brief_reanalysis.R) carries out the manuscript steps used by that command. It is not a second user command.
- [`run_pellett_valbuena_replication.R`](run_pellett_valbuena_replication.R) starts the separate optional workflow for selected Pellett & Valbuena procedures rewritten in R. Use the [optional workflow guide](../PELLETT_VALBUENA_WORKFLOW.md) for its scientific purpose and complete command line and R/RStudio Console instructions.

## Numbered R files

The main scripts load these files in order. Each file groups related tasks and is not run on its own.

| Order | File | Purpose |
|---:|---|---|
| 00 | [`R/00_setup_paths.R`](R/00_setup_paths.R) | Find the repository root and define paths and run settings. |
| 01 | [`R/01_packages.R`](R/01_packages.R) | Check and load required packages. |
| 02 | [`R/02_io_lock_metadata.R`](R/02_io_lock_metadata.R) | Read and write files, prevent overlapping runs, and record setup information. |
| 03 | [`R/03_source_data_prep.R`](R/03_source_data_prep.R) | Prepare source data used by the analyses. |
| 04 | [`R/04_pellett_valbuena_replication.R`](R/04_pellett_valbuena_replication.R) | Define the optional Pellett & Valbuena analysis. |
| 05 | [`R/05_delta_hdr_model_helpers.R`](R/05_delta_hdr_model_helpers.R) | Define δ heterogeneity–diversity relationship (HDR) models and helpers for manuscript outputs. |
| 06 | [`R/06_brief_analysis_hierarchy.R`](R/06_brief_analysis_hierarchy.R) | Fit and compare the primary and sensitivity manuscript models. |
| 07 | [`R/07_brief_deliverables.R`](R/07_brief_deliverables.R) | Produce manuscript figures, tables, and prospective caption files. |
| 08 | [`R/08_validation_helpers.R`](R/08_validation_helpers.R) | Check files and compare results with reference values. |

## Separate workflows

The manuscript command loads files 00 through 08. File 04 defines tools for the optional analysis, but loading those definitions does not run it. The optional script loads files 00 through 04 and 08, then starts the separate analysis. Use the documented entry commands rather than running numbered files directly.
