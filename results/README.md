# Results guide

This folder contains five included computational products from the manuscript analysis. Readers can open them directly before running any code. The map below distinguishes these products from the manuscript's two interpretive summary tables.

## Included products

### Manuscript-to-archive map

| Manuscript item | Relationship to the archive |
|---|---|
| [Figure 1](delta_hdr_reanalysis/manuscript_tables_figures/main_text/main_text_mean_independent_heterogeneity_reveals_diminishing_return_HDRs.png) | Included computational figure. Panel a presents conceptual hypotheses, predictions, and applications; panel b shows illustrative curves; panels c–d show the primary MacArthur and Catalonia fits. |
| Main Table 1 | Interpretive summary assembled for the manuscript; not a generated computational output or an included archive artifact. |
| [Extended Data Fig. 1](delta_hdr_reanalysis/manuscript_tables_figures/extended_data/extended_data_figure_1_macarthur_all_analyses.png) | Included computational figure: the MacArthur primary analysis and sensitivity analysis that includes mean foliage height. |
| [Extended Data Fig. 2](delta_hdr_reanalysis/manuscript_tables_figures/extended_data/extended_data_figure_2_allouche_all_analyses.png) | Included computational figure: the Catalonia primary analysis and sensitivity analyses that include mean elevation or mean elevation plus squared mean elevation. |
| [Extended Data Table 1](delta_hdr_reanalysis/manuscript_tables_figures/extended_data/extended_data_table_1_macarthur_model_comparisons_and_coefficient_checks.md) | Included computational table: MacArthur model comparisons and checks of response shape. |
| [Extended Data Table 2](delta_hdr_reanalysis/manuscript_tables_figures/extended_data/extended_data_table_2_allouche_model_comparisons_and_coefficient_checks.md) | Included computational table: Catalonia model comparisons and checks of response shape. |
| Extended Data Table 3 | Interpretive summary assembled for the manuscript; not a generated computational output or an included archive artifact. |

### Reading the figures and tables

Figure 1a summarizes conceptual hypotheses and their conditional applications. The applications assume that the focal heterogeneity facet is manipulable and causally affects biodiversity within the intervention range.

Figure 1b shows illustrative plotting coordinates, not curves fitted to observations. Figure 1c–d show the primary observations and selected logarithmic fitted means. Their pointwise 95% confidence intervals describe uncertainty in the fitted mean at each prediction coordinate, not prediction intervals for individual observations.

The MacArthur plots use foliage-height δ in feet and a dimensionless Shannon index of bird species diversity for 13 sites. The Catalonia plots use topographic δ in meters and breeding bird species richness, a count, for 285 retained 10 × 10 km grid cells. In both datasets, δ is sample variance divided by the environmental mean.

Extended Data Fig. 1 and Extended Data Fig. 2 use the same model styling as Figure 1: blue long-dashed lines for linear models, purple solid lines for logarithmic models, and vermilion dot-dashed lines for quadratic models. The tables of model rankings beneath the plots repeat these mappings with labeled swatches.

In Extended Data Fig. 1, panels a and c present the primary MacArthur analysis, and panels b and d present the sensitivity analysis that includes mean foliage height. In Extended Data Fig. 2, panels a and d present the primary Catalonia analysis, panels b and e add mean elevation, and panels c and f add mean elevation and squared mean elevation.

The curves in sensitivity analyses are conditional predictions, with adjustment variables held fixed. Mean foliage height or mean elevation is held at its observed sample mean. For the two-term Catalonia adjustment, the squared covariate is held at the mean of squared cell-level mean elevations, not the square of the mean elevation across cells. The points remain unadjusted observations, so distances from points to these conditional curves are not residuals. The [source data guide](delta_hdr_reanalysis/source_data/README.md) defines the exact fields, analysis identifiers, units, fixed covariate settings, and pointwise confidence intervals. Extended Data Figs. 1–2 show no confidence bands.

Extended Data Table 1 and Extended Data Table 2 use Akaike’s information criterion (AIC) or its small-sample correction (AICc) to rank candidate models within each analysis. Lower values rank higher; ΔAIC or ΔAICc is each candidate’s difference from the lowest value, and Akaike weights summarize relative support within that candidate set. In the tables, `n` is sample size, and bracketed values in the focal term column are 95% confidence intervals.

## Figure and table source data

The repository also includes five CSV tables containing source data for Figure 1b–d, Extended Data Fig. 1, Extended Data Fig. 2, Extended Data Table 1, and Extended Data Table 2. Figure 1a is a conceptual hypothesis table rather than a direct display of these CSV values. The source data files can be opened directly in spreadsheet or statistical software and are generated separately from the main and optional workflows:

- [`observations.csv`](delta_hdr_reanalysis/source_data/observations.csv) records the 13 MacArthur and 285 Catalonia observations and their environmental summaries.
- [`fitted_mean_curves.csv`](delta_hdr_reanalysis/source_data/fitted_mean_curves.csv) records curves of the fitted mean response for all five analyses and three candidate families, including the pointwise 95% intervals for fitted means used only in Figure 1c–d. Pointwise means that each reported interval applies separately at its prediction coordinate.
- [`model_comparisons.csv`](delta_hdr_reanalysis/source_data/model_comparisons.csv) records candidate formulas, information criteria, differences, weights, ranks, and checks of response shape.
- [`model_terms.csv`](delta_hdr_reanalysis/source_data/model_terms.csv) records coefficient estimates, uncertainty, and roles of model terms.
- [`figure1_illustrative_curves.csv`](delta_hdr_reanalysis/source_data/figure1_illustrative_curves.csv) records the plotting coordinates for the illustrative curves in Figure 1b.

The [source data guide](delta_hdr_reanalysis/source_data/README.md) maps these files to the empirical and illustrative figure panels and Extended Data tables and defines every field and unit. To recreate all five files, run this separate command from the repository root:

```sh
Rscript --vanilla code/export_figure_table_source_data.R
```

## Regenerate the products from the manuscript analysis

From the repository root, run the main manuscript command:

```sh
Rscript --vanilla code/run_reanalysis.R
```

The [run instructions](../RUN_INSTRUCTIONS.md) explain the required R packages, tested environment, output locations, and completion checks.

## Additional generated files

The main manuscript command also creates analysis inventories, captions, alternative formats, and setup records. These additional files are created only in your local copy and are not among the files provided in the repository; `.gitignore` identifies them. The five included computational products in the map above are the only products of the manuscript analysis provided in the repository.

Standalone full caption files are generated by this command under `results/delta_hdr_reanalysis/manuscript_tables_figures/captions/`. The reading notes in this guide are reader documentation, not new generated captions or a replacement for the manuscript.

The optional Pellett & Valbuena command creates a separate output tree beneath `provenance/`. Those files are created only in your local copy and are not among the files provided in the repository. The optional command does not produce or replace the five included products from the manuscript analysis. See the [optional workflow guide](../PELLETT_VALBUENA_WORKFLOW.md).
