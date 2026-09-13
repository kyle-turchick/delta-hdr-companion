# Figure and table source data

These five comma-separated value (CSV) files connect the included analyses to Figure 1b–d, Extended Data Fig. 1, Extended Data Fig. 2, Extended Data Table 1, and Extended Data Table 2. A CSV is a table stored as plain text that can be opened in spreadsheet or statistical software. Figure 1a is a conceptual hypothesis table and is not derived directly from these CSVs.

Readers can inspect the included files before running any code:

| File | Rows | What it contains | Unique key |
|---|---:|---|---|
| `observations.csv` | 298 | Observed responses and calculated heterogeneity values | `dataset_id`, `observation_id` |
| `fitted_mean_curves.csv` | 3,000 | Prediction coordinates and fitted mean responses | `analysis_id`, `model_family`, `prediction_index` |
| `model_comparisons.csv` | 15 | Candidate model rankings and checks of response shape | `analysis_id`, `model_family` |
| `model_terms.csv` | 47 | Coefficient estimates and uncertainty | `analysis_id`, `model_family`, `term_order` |
| `figure1_illustrative_curves.csv` | 603 | Plotting coordinates for Figure 1b's conceptual curves | `illustrative_family`, `point_index` |

A **unique key** is the smallest combination of fields that identifies one row. **Lineage** fields connect a row to its prepared site, raster cell, or analysis. The heterogeneity metric δ (delta) equals the sample variance divided by the environmental mean.

## Datasets and analyses

The MacArthur dataset contains 13 sites (A–M). Its environmental variable is foliage height, measured in feet (`ft`); sample variance is in square feet (`ft^2`), foliage-height δ is in feet, and the response is a dimensionless Shannon index for bird species diversity, calculated with the natural logarithm. The Catalonia dataset contains 285 retained raster cells. Its environmental variable is elevation, measured in meters (`m`); sample variance is in square meters (`m^2`), topographic δ is in meters, and the response is breeding bird species richness, recorded as a count.

Each of the five prespecified analysis sets compares linear, logarithmic, and quadratic candidate models using Gaussian ordinary least squares (OLS), which fits a mean response by minimizing squared differences between observed and fitted values. The Catalonia comparisons are nonspatial: the fitted models do not include a spatial correlation structure. A **primary** analysis addresses the main comparison; a **sensitivity** analysis repeats it with stated adjustment terms.

| Order | Analysis ID | Dataset | Role | Adjustment | n | Criterion |
|---:|---|---|---|---|---:|---|
| 1 | `Foliage_delta_unadjusted_gaussian` | `macarthur_1961` | `primary` | `unadjusted` | 13 | AICc |
| 2 | `Foliage_delta_mean_foliage_height_gaussian` | `macarthur_1961` | `sensitivity` | `mean_foliage_height` | 13 | AICc |
| 3 | `Topography_delta_no_mean_covariates_gaussian` | `allouche_2012` | `primary` | `no_mean_covariates` | 285 | AIC |
| 4 | `Topography_delta_mean_elevation_linear_gaussian` | `allouche_2012` | `sensitivity` | `mean_elevation_linear` | 285 | AIC |
| 5 | `Topography_delta_mean_elevation_quadratic_gaussian` | `allouche_2012` | `sensitivity` | `mean_elevation_quadratic` | 285 | AIC |

Here, `n` is sample size. Akaike's information criterion (AIC) and its small-sample correction (AICc) rank candidates within an analysis; lower values rank first. A ΔAIC or ΔAICc value is a candidate's difference from the minimum, and an Akaike weight is its normalized relative support within the candidate set.

## Figure 1a conceptual mapping

Figure 1a presents these hypotheses directly; it does not use rows from the five CSVs containing source data. Its display order in the panel differs from the analytical candidate order above. In the formulas, `y` is the biodiversity response, `~` means "is modeled as a function of," δ is heterogeneity, δ² is squared heterogeneity, and `log(δ)` is the natural logarithm of heterogeneity.

| Display order | Family | Formula | Hypothesis | Prediction | Application |
|---:|---|---|---|---|---|
| 1 | linear | `y ~ δ` | Habitat heterogeneity hypothesis | Biodiversity rises at a constant rate | Starting heterogeneity value alone does not guide targeted intervention |
| 2 | quadratic | `y ~ δ + δ²` | Area–heterogeneity trade-off hypothesis | Biodiversity rises, peaks, then declines | Manage heterogeneity toward the biodiversity peak |
| 3 | logarithmic | `y ~ log(δ)` | Diminishing returns hypothesis | Biodiversity rises steeply, then gains taper | Increase heterogeneity when starting in homogeneous landscapes; when reducing heterogeneity, start in heterogeneous landscapes |

The applications assume that the focal heterogeneity facet is manipulable and causally affects biodiversity within the intervention range, as defined in the Figure 1 caption.

## Figure and table map

| Product | Analysis or content | Source tables |
|---|---|---|
| Figure 1a | Conceptual hypothesis, model, prediction, and application table | Defined in the source that generates the figure, not these CSVs |
| Figure 1b | Illustrative curves used only for plotting | `figure1_illustrative_curves.csv` |
| Figure 1c | MacArthur primary observations, selected fitted mean and interval, and ranking of all candidates | `observations.csv`, `fitted_mean_curves.csv`, `model_comparisons.csv` |
| Figure 1d | Catalonia primary observations, selected fitted mean and interval, and ranking of all candidates | `observations.csv`, `fitted_mean_curves.csv`, `model_comparisons.csv` |
| Extended Data Fig. 1, panels a–b | MacArthur primary curves and sensitivity curves that include mean foliage height | `observations.csv`, `fitted_mean_curves.csv` |
| Extended Data Fig. 1, panels c–d | Corresponding MacArthur candidate rankings | `model_comparisons.csv` |
| Extended Data Fig. 2, panels a–c | Catalonia primary curves and sensitivity curves adjusted for elevation | `observations.csv`, `fitted_mean_curves.csv` |
| Extended Data Fig. 2, panels d–f | Corresponding Catalonia candidate rankings | `model_comparisons.csv` |
| Extended Data Table 1 | MacArthur comparisons and coefficient checks | `model_comparisons.csv`, `model_terms.csv` |
| Extended Data Table 2 | Catalonia comparisons and coefficient checks | `model_comparisons.csv`, `model_terms.csv` |

Figure 1c–d show pointwise 95% confidence intervals for the selected fitted means. "Pointwise" means that each prediction coordinate has its own interval; the intervals are not a simultaneous band for the entire curve. Exactly 400 rows in `fitted_mean_curves.csv` contain intervals: 200 for each primary logarithmic curve. Extended Data Fig. 1 and Extended Data Fig. 2 show no confidence bands.

The sensitivity panels plot the unadjusted observed responses at their calculated δ coordinates.

Their curves show conditional predictions with adjustment variables held at specified values. Predictions that include mean foliage height or mean elevation hold the relevant adjustment variable at its empirical mean. The sensitivity analysis with two elevation terms holds mean elevation and squared mean elevation separately at their empirical means; together, these held values reproduce the average additive contribution of the two elevation terms across the sample.

A residual is an observed response minus its fitted value under the same model conditions. Because the plotted points are unadjusted and the curves use fixed adjustment values, distances from the points to the curves are not residuals.

## Field dictionary: `observations.csv`

| Field | Meaning |
|---|---|
| `dataset_id` | Stable machine identifier for the MacArthur or Catalonia dataset. |
| `observation_id` | Identifier unique within a dataset: site A–M for MacArthur or the retained source-cell identifier for Catalonia. |
| `site_id` | MacArthur site letter A–M; `NA` for Catalonia rows. |
| `utm10` | Catalonia source-cell identifier retained from the prepared data; `NA` for MacArthur rows. |
| `source_fid` | Catalonia source-raster lineage identifier; `NA` for MacArthur rows. |
| `raster_filename` | Filename of the Catalonia raster associated with the retained cell; `NA` for MacArthur rows. |
| `environmental_variable` | Environmental quantity summarized for the observation: foliage height or elevation. |
| `environmental_mean` | Arithmetic mean of the retained environmental values. |
| `environmental_mean_unit` | Unit of `environmental_mean`: `ft` or `m`. |
| `environmental_sample_variance` | Sample variance of the retained environmental values. |
| `environmental_variance_unit` | Squared unit of the sample variance: `ft^2` or `m^2`. |
| `delta` | Calculated heterogeneity value: sample variance divided by environmental mean. |
| `delta_unit` | Unit of δ: `ft` or `m`. |
| `response_name` | Scientific response represented by `response_value`. |
| `response_value` | Observed Shannon index for bird species diversity or observed breeding bird species richness. |
| `response_unit` | `dimensionless` for the Shannon index or `count` for richness. |
| `retained_value_count` | Number of finite environmental values used to calculate the observation's mean and sample variance. |

## Field dictionary: `fitted_mean_curves.csv`

| Field | Meaning |
|---|---|
| `analysis_id` | Stable machine identifier for one of the five analysis sets. |
| `analysis_display_order` | Order of the five analysis sets in repository outputs. |
| `dataset_id` | Stable machine identifier for the dataset. |
| `analysis_role` | `primary` for a main comparison or `sensitivity` for an adjustment analysis. |
| `adjustment_scenario` | Stable machine identifier for the covariate adjustment setting. |
| `model_family` | Candidate response shape: `linear`, `logarithmic`, or `quadratic`. |
| `candidate_display_order` | Controlled candidate order: linear, logarithmic, quadratic. |
| `prediction_index` | Row number within the prediction grid for an analysis and model family, ordered by δ. |
| `delta` | δ coordinate at which the fitted mean was calculated. |
| `delta_unit` | Unit of δ: `ft` or `m`. |
| `fitted_mean_response` | Mean response estimated by the model at the prediction coordinate and any fixed adjustment values. |
| `response_name` | Scientific response represented by the fitted mean. |
| `response_unit` | `dimensionless` for the Shannon index or `count` for richness. |
| `fixed_mean_foliage_height_ft` | Empirical mean foliage height held fixed for the MacArthur sensitivity prediction; otherwise `NA`. |
| `fixed_mean_elevation_m` | Empirical mean elevation held fixed for the relevant Catalonia sensitivity predictions; otherwise `NA`. |
| `fixed_mean_elevation_squared_m2` | Empirical mean of squared mean elevation held fixed for the Catalonia sensitivity analysis that includes mean elevation and squared mean elevation; otherwise `NA`. |
| `confidence_interval_included` | `TRUE` when the row contains a confidence interval for the fitted mean and `FALSE` otherwise. |
| `confidence_interval_type` | `pointwise_fitted_mean` for included intervals; otherwise `NA`. |
| `confidence_level` | Interval confidence level, `0.95`, when an interval is included; otherwise `NA`. |
| `confidence_lower` | Lower endpoint of the included confidence interval for the fitted mean; otherwise `NA`. |
| `confidence_upper` | Upper endpoint of the included confidence interval for the fitted mean; otherwise `NA`. |

Before reading the two model tables, note their formula notation: `y` is the response, `x` is δ, `logx` is the natural logarithm of δ, and `I(x^2)` is squared δ. The adjustment term `mu` is mean foliage height or mean elevation, depending on the dataset, and `mu2` is squared mean elevation. A plus sign adds a model term.

## Field dictionary: `model_comparisons.csv`

| Field | Meaning |
|---|---|
| `analysis_id` | Stable machine identifier for one of the five analysis sets. |
| `analysis_display_order` | Order of the five analysis sets in repository outputs. |
| `dataset_id` | Stable machine identifier for the dataset. |
| `analysis_role` | `primary` for a main comparison or `sensitivity` for an adjustment analysis. |
| `adjustment_scenario` | Stable machine identifier for the covariate adjustment setting. |
| `model_family` | Candidate response shape: `linear`, `logarithmic`, or `quadratic`. |
| `candidate_display_order` | Controlled candidate order: linear, logarithmic, quadratic. |
| `formula` | Fitted model formula using the workflow's standardized predictor names. |
| `n` | Number of observations fitted in the model. |
| `k` | Number of estimated parameters used in the information criterion. |
| `information_criterion` | Ranking criterion for the analysis: `AIC` or `AICc`. |
| `information_criterion_value` | AIC or AICc value at full precision. |
| `delta_ic` | Difference between this candidate's criterion value and the minimum within the analysis. |
| `akaike_weight` | Normalized relative support for the candidate within the analysis; the three weights sum to one apart from numerical rounding. |
| `rank` | Criterion rank within the analysis, with 1 denoting the minimum and exact ties sharing a rank. |
| `best_by_ic` | `TRUE` when the candidate has the exact minimum criterion value within the analysis. |
| `predictor_min` | Minimum observed δ used to fit the analysis. |
| `predictor_max` | Maximum observed δ used to fit the analysis. |
| `quadratic_vertex` | δ coordinate of the quadratic curve's vertex; `NA` for other families. |
| `quadratic_vertex_inside_observed_range` | Whether the quadratic vertex lies within the observed δ range; `NA` for other families. |
| `quadratic_high_side_declines` | Whether the fitted quadratic declines toward the upper end of the observed δ range; `NA` for other families. |
| `shape_supported_by_coefficients` | Whether coefficient signs, P values, and confidence intervals support the interpretation of the fitted response shape. For quadratic models, this field also incorporates the implemented geometry checks. |
| `shape_support_with_uncertainty` | Machine-readable label summarizing coefficient and geometry evidence for the fitted shape. |
| `shape_interpretation` | Machine-readable interpretation of the fitted response shape. |

### Machine-coded shape values

| Value | Meaning |
|---|---|
| `positive_monotonic` | Fitted response increases across the observed δ range. |
| `positive_diminishing_return` | Fitted response increases but its gains weaken as δ increases. |
| `unimodal_in_observed_range` | Fitted quadratic rises and then declines within the observed δ range. |
| `positive_monotonic_supported_by_coefficient_CI` | Coefficient signs, P values, and confidence intervals support a positive monotonic shape. |
| `positive_diminishing_return_supported_by_coefficient_CI` | Coefficient signs, P values, and confidence intervals support a positive diminishing-return shape. |
| `unimodal_supported_by_coefficients_CI_and_geometry` | Coefficient checks plus the implemented vertex and upper-range decline checks support a quadratic rise and decline. |
| `IC_selected_quadratic_shape_but_coefficient_or_geometry_uncertain` | Legacy machine label for a quadratic shape that rises and declines but has uncertain coefficient or geometry evidence; it does not indicate that the quadratic candidate ranked first in the comparison of three candidates. |

## Field dictionary: `model_terms.csv`

| Field | Meaning |
|---|---|
| `analysis_id` | Stable machine identifier for one of the five analysis sets. |
| `analysis_display_order` | Order of the five analysis sets in repository outputs. |
| `dataset_id` | Stable machine identifier for the dataset. |
| `analysis_role` | `primary` for a main comparison or `sensitivity` for an adjustment analysis. |
| `adjustment_scenario` | Stable machine identifier for the covariate adjustment setting. |
| `model_family` | Candidate response shape: `linear`, `logarithmic`, or `quadratic`. |
| `candidate_display_order` | Controlled candidate order: linear, logarithmic, quadratic. |
| `formula` | Fitted model formula using the workflow's standardized predictor names. |
| `term_order` | Coefficient order in the fitted model matrix. |
| `term` | Exact term name in the model matrix. |
| `term_role` | Machine-readable description of the term's scientific role. |
| `estimate` | Estimated coefficient. |
| `std_error` | Standard error of the coefficient estimate. |
| `statistic` | Coefficient test statistic reported by the fitted model. |
| `p_value` | P value for the coefficient test. |
| `conf_low` | Lower endpoint of the 95% coefficient confidence interval. |
| `conf_high` | Upper endpoint of the 95% coefficient confidence interval. |

The `term_role` values translate as follows: `intercept` is the fitted baseline; `heterogeneity_linear_term` and `heterogeneity_logarithmic_term` are the focal δ effects; `heterogeneity_quadratic_linear_term` and `heterogeneity_quadratic_squared_term` are the two focal quadratic terms; and `covariate_or_other_term` is an adjustment term.

## Field dictionary: `figure1_illustrative_curves.csv`

| Field | Meaning |
|---|---|
| `illustrative_family` | Conceptual curve identity: linear, logarithmic, or quadratic. |
| `point_index` | Order of the plotting point within the conceptual curve. |
| `plotting_x` | Unitless coordinate along the conceptual heterogeneity gradient. |
| `plotting_y` | Unitless plotting height chosen to illustrate the hypothesis about response shape. |

The Figure 1b coordinates come from functions used only for plotting. They are not empirical observations, fitted values, or parameter estimates.

## Missing values and precision

The files use `NA`, not blank text, when a field does not apply. They contain no `NaN` values (undefined numeric results) or infinite values. Computers store many decimal values as close binary approximations rather than exact decimal numbers. The exporter writes decimal values with 17 significant digits in the C locale, which uses a period as the decimal mark, to make repeated exports deterministic. This makes some approximations visible: for example, intended values 0.95 and 0.639 appear as `0.94999999999999996` and `0.63900000000000001`. The exporter then rereads every column using its intended type and requires the table structure to match exactly. Text, whole numbers, `TRUE`/`FALSE` values, missingness, and zero values must match exactly. For each finite, nonzero decimal value, it requires `abs(observed - expected) / abs(expected) <= 1e-12`. This threshold verifies preservation through CSV writing and rereading; it is not a tolerance for differences among analytical results.

## Regenerate the source data CSV files

Regeneration is optional for readers who only want to inspect the included source data. To regenerate them, open a terminal application, change to the repository root (the main project folder containing `README.md`), and run:

```sh
Rscript --vanilla code/export_figure_table_source_data.R
```

`Rscript` runs an R file from the command line, and `--vanilla` starts R without loading or saving a personal workspace or startup profile. The exporter prepares and fits the five analysis sets once, validates the tables in memory, and writes and rereads a temporary complete set. It then replaces all five files together. If replacement or final validation fails, it restores the prior complete set; if no prior set existed, it removes any partial new set. The exporter runs separately from the main manuscript analysis and the optional R reimplementation of selected Pellett & Valbuena analyses. It does not write figures, manuscript tables, or products from the optional workflow.

## Attribution and reuse

The exporter is project software covered by the MIT grant described in the repository [`LICENSE`](../../../LICENSE). This README and the five derived CSVs fall within the project's CC BY 4.0 scope only to the extent described in that license. Underlying archived values retain the third-party treatment recorded in [`THIRD_PARTY_NOTICES.md`](../../../THIRD_PARTY_NOTICES.md). The [data guide](../../../data/README.md) and [`INPUT_MANIFEST.csv`](../../../data/INPUT_MANIFEST.csv) identify the archived inputs and their provenance.
