# Reshaping landmark heterogeneity–diversity relationships

This companion contains the data, R code, and results for a reanalysis of two classic datasets of bird diversity. The main manuscript workflow consists of ordered R steps that prepare data, fit models, and create results. It uses these observational data to estimate heterogeneity–diversity relationships (HDRs), which are associations between bird diversity and variation in environmental conditions. For both datasets, `δ = sample variance / mean`; sample variance measures how spread out the environmental values are. Before fitting the models, we specified linear, logarithmic, and quadratic candidate model families to represent three alternative response shapes: a steady increase, an increase with diminishing returns, and a rise followed by a decline.

The logarithmic candidate ranked first in both primary analyses and all three sensitivity analyses. This support is relative to the specified candidate models in these two datasets; it does not establish a universal response or demonstrate a causal intervention effect.

## Download and run the main manuscript analysis

You do not need Git to download this project. If you can view its GitHub page, select **Code**, choose **Download ZIP**, and extract the downloaded file.

Before running code, you can open the included [figures and tables](#figures-and-tables) and inspect the [CSV files containing source data](results/delta_hdr_reanalysis/source_data/README.md) directly.

Install R and the required add-on packages as described in the [run instructions](RUN_INSTRUCTIONS.md). The [environment guide](environment/README.md) records the tested R and package versions.

Prefer to begin in R or RStudio? The [R or RStudio Console route](RUN_INSTRUCTIONS.md#option-2-r-or-rstudio-console) starts the same entry script in a separate R process without loading a saved workspace or personal startup files.

Open a terminal application, then change to the extracted main project folder. The command `cd` means "change directory." For example, use `cd "/full/path/to/extracted-project-folder"` with the folder's actual path. This folder is the repository root; it contains `README.md`, `code/`, `data/`, `environment/`, and `results/`. From that folder, run:

```sh
Rscript --vanilla code/run_reanalysis.R
```

This command starts the main manuscript workflow. The [run instructions](RUN_INSTRUCTIONS.md) explain the required setup, how to check that `Rscript` is available, what the command creates, and how to recognize successful completion.

## Inputs

The repository contains 407 archived input files. The main manuscript workflow uses 402 inputs: 386 TIFF raster files and 16 comma-separated value (CSV) tables. The optional Pellett & Valbuena R workflow uses five additional archived inputs: two CSV summaries required for calculations and three TIFF rasters that provide map context.

The [input manifest](data/INPUT_MANIFEST.csv) is the file list used to identify and verify these inputs. It records each file's workflow scope, byte size, and SHA-256 checksum (a file fingerprint used to confirm that the file has not changed). The retained machine value `canonical_and_optional` identifies inputs shared by the main and optional workflows; `optional_only` identifies the five inputs used only by the optional workflow.

A clean release copy also contains `RELEASE_MANIFEST.csv`. This file inventories the other 452 files, bringing the total to 453. For each file, it records the relative path, byte size, and SHA-256 fingerprint. It excludes itself because a file cannot stably record its own fingerprint.

## Repository contents

- [`code/`](code/) contains the R scripts and a [code guide](code/README.md).
- [`data/`](data/) contains the archived inputs, their manifest, and a [data guide](data/README.md).
- [`environment/`](environment/) records the R version and add-on package versions used for testing and includes an [environment guide](environment/README.md).
- [`results/`](results/) contains five included products from the manuscript analysis (three figures and two tables), CSV files containing source data for Figure 1b–d, Extended Data Fig. 1, Extended Data Fig. 2, Extended Data Table 1, and Extended Data Table 2, and a [results guide](results/README.md).

## Source datasets

The manuscript analysis uses MacArthur & MacArthur data on foliage height and forest bird diversity and Allouche et al. data on Catalonia elevation and breeding birds. The [data guide](data/README.md) describes their variables, processing, and provenance.

## Figures and tables

- [Figure 1](results/delta_hdr_reanalysis/manuscript_tables_figures/main_text/main_text_mean_independent_heterogeneity_reveals_diminishing_return_HDRs.png)
- [Extended Data Fig. 1](results/delta_hdr_reanalysis/manuscript_tables_figures/extended_data/extended_data_figure_1_macarthur_all_analyses.png)
- [Extended Data Fig. 2](results/delta_hdr_reanalysis/manuscript_tables_figures/extended_data/extended_data_figure_2_allouche_all_analyses.png)
- [Extended Data Table 1](results/delta_hdr_reanalysis/manuscript_tables_figures/extended_data/extended_data_table_1_macarthur_model_comparisons_and_coefficient_checks.md)
- [Extended Data Table 2](results/delta_hdr_reanalysis/manuscript_tables_figures/extended_data/extended_data_table_2_allouche_model_comparisons_and_coefficient_checks.md)

Main Table 1 and Extended Data Table 3 are interpretive summaries assembled for the manuscript, not generated computational outputs. The [results guide](results/README.md#manuscript-to-archive-map) explains their relationship to the five included products.

## Optional workflow

The [optional workflow guide](PELLETT_VALBUENA_WORKFLOW.md) explains how selected Pellett & Valbuena procedures originally written in Julia were written again in R. This optional workflow provides more detail about the methods, the origins of the archived data, and selected comparisons between R and Julia. The main workflow recreates the manuscript's computational analyses and five included products: Figure 1, Extended Data Figs. 1–2, and Extended Data Tables 1–2. The detailed guide explains how to start the workflow from the command line or R/RStudio Console, how to recognize success, where to find its outputs, and what the comparisons covered. The five additional inputs are already included. Generated optional products remain only in the reader’s project copy.

## Citation, attribution, and reuse

[`CITATION.cff`](CITATION.cff) provides the companion’s current machine-readable citation metadata. It currently records the companion title and type, author, and affiliation; final release and publication details will be added when authoritative values exist. Because reuse terms differ among file groups, the project [`LICENSE`](LICENSE) explains which terms apply.

The 407 archived inputs are redistributed unchanged from the Pellett & Valbuena data and code deposits. Project-authored material and redistributed inputs have separate terms. Before reuse, consult the [data guide](data/README.md), [third-party notices](THIRD_PARTY_NOTICES.md), project [`LICENSE`](LICENSE), and included [CC BY 4.0 legal code](LICENSES/CC-BY-4.0.txt).
