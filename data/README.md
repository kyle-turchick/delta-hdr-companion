# Data guide

The repository contains 407 archived input files: 402 used by the main manuscript workflow and five used only by the optional Pellett & Valbuena R workflow. A comma-separated value (CSV) file is a tabular text file; a TIFF file is a gridded raster image.

[`INPUT_MANIFEST.csv`](INPUT_MANIFEST.csv) is the file list used to identify and verify every input. Each row represents one file.

| Workflow scope | Files | Role |
|---|---:|---|
| Main and optional (`canonical_and_optional`) | 386 TIFF and 16 CSV | Inputs shared by the main manuscript and optional workflows |
| Optional only (`optional_only`) | 3 TIFF and 2 CSV | Archived calculation summaries and map context used only by the optional R reimplementation |

The scope strings in parentheses are stable machine values and must not be renamed.

## Input manifest columns

- `path`: the file's location relative to the repository root.
- `deposit_member_path`: the matching file location inside the archived Pellett & Valbuena deposit.
- `workflow_scope`: the stable machine value identifying which workflow uses the file.
- `attribution_key`: the source group used to connect the file with its attribution entry: `macarthur_foliage`, `catalonia_bird_atlas`, `srtm_elevation`, or `copernicus_crop_cover`.
- `file_type`: the recorded file format, CSV or TIFF.
- `size_bytes`: the exact file size in bytes.
- `sha256`: the SHA-256 file fingerprint used to confirm unchanged bytes.
- `zenodo_v5_doi`: the Digital Object Identifier (DOI) for the matching Zenodo v5 deposit.
- `zenodo_v6_doi`: the DOI for the matching Zenodo v6 deposit.
- `deposit_record_declared_license`: the license treatment declared by the deposit record.
- `redistribution_status`: whether the included file bytes are unchanged from the verified deposit member.

The repository contains several kinds of material, not a single set of untouched raw data.

| Material | Source and role |
|---|---|
| MacArthur & MacArthur (1961) | Original study |
| Foliage profiles A–M | Digitized from the original study |
| `foliage_density.csv` | Archived processed input for site heterogeneity |
| Bird species diversity values in R | Shannon index values embedded in the analysis code |
| Catalonia bird and coordinate files | Bird data and coordinates from the Catalan Breeding Bird Atlas |
| Catalonia elevation tiles | Archived elevation inputs derived from SRTM |
| Optional heterogeneity summaries | Archived derived calculation inputs for the optional workflow |
| Optional map rasters | Archived context layers for optional figures |
| Calculated statistics, models, figures, and tables | Created by the R workflows |

## Heterogeneity measure

For both manuscript datasets, `δ = sample variance / mean`. δ is designed to reduce dependence on the environmental mean and retains the source variable's units: feet for foliage height and meters for elevation. The corresponding sample variances have square feet and square meters as their units.

## MacArthur & MacArthur forest sites

The analysis covers 13 forest sites, labeled A–M.

- [`foliage_density.csv`](external/pellett_valbuena_2025/05_review/data/foliage_density.csv) has one column per site and 10,000 processed values of foliage height in feet per column. Pellett & Valbuena describe these values as samples from smoothed curves fitted to digitized profiles. The main manuscript workflow reads this deposited table directly, preserving the exact input for all sites used in the reported analyses. From it, the workflow calculates each site’s mean foliage height, sample variance, and foliage-height δ. The digitized point files described below support provenance and the optional figures of foliage profiles.
- [`digitised_points/`](external/pellett_valbuena_2025/05_review/data/digitised_points/) contains `A.csv` through `M.csv`, digitized from Figure 1 of MacArthur and MacArthur (1961). In each file, `x` is a digitized value of foliage density, expressed as square feet of leaf silhouette per cubic foot of space, and `y` is a digitized value of foliage height in feet. The optional R reimplementation applies the documented digitization corrections and boundary points before creating the optional figures of foliage profiles.
- [`code/R/03_source_data_prep.R`](../code/R/03_source_data_prep.R) embeds the 13 dimensionless Shannon index values for bird species diversity reported for sites A through M by MacArthur and MacArthur (1961). Their index uses the natural logarithm, and Pellett and Valbuena reused the published values. The script retains the published site order and pairs the values with the foliage columns in the same order. The R source stores these response values directly rather than reading them from a separate input file.

## Allouche et al. Catalonia grid cells

The analysis re-examines Allouche et al. (2012). Bird records come from the Catalan Breeding Bird Atlas 1999–2002, collected by the Institut Català d'Ornitologia and supplied to Pellett & Valbuena by Jofre Carnicer and Lluís Brotons. Each analytical unit is a 10 × 10 km Universal Transverse Mercator (UTM) grid cell.

- [`birds.csv`](external/pellett_valbuena_2025/03_corrected_HDR/empirical/01_carnicer/data/bird/birds.csv) contains the `UTM10` cell identifier and species columns. Richness is the number of species columns with a numeric value greater than zero and is a count.
- [`coordinates.csv`](external/pellett_valbuena_2025/03_corrected_HDR/empirical/01_carnicer/data/bird/coordinates.csv) has 386 rows linking `Utm10` identifiers to `X_coord` and `Y_coord`. Metadata retained with the source rasters identify the coordinate reference system as World Geodetic System 1984 / UTM zone 31N (`EPSG:32631`), with easting and northing in meters. The workflow uses the supplied coordinate values without transforming them.
- [`elevation/elevgrid/`](external/pellett_valbuena_2025/03_corrected_HDR/empirical/01_carnicer/data/elevation/elevgrid/) contains elevation tiles `0.tif` through `385.tif`, derived from the Shuttle Radar Topography Mission (SRTM), with values in meters. Pellett & Valbuena cite Jarvis et al. (2008) for the elevation source.

The workflow keeps elevation values that are not missing or infinite and are greater than zero. For each tile it calculates the number of retained pixels, mean, sample variance, and topographic δ. It pairs tile summaries with coordinate rows in numeric filename order, then matches bird richness to the tile summaries using the `UTM10` and `Utm10` cell identifiers. Cells with all required matched values and more than 14,000 retained elevation pixels are kept, yielding 285 cells. Pixels measure environmental heterogeneity within a cell; they are not observations of bird diversity.

## Inputs for the optional workflow

The optional workflow uses five additional archived inputs at their locations in the source project:

| Input | Role |
|---|---|
| `02_MBH_comb/data/elevation/heterogeneity_inds.csv` | Derived elevation summary used in calculations |
| `02_MBH_comb/data/crop_cover/heterogeneity_inds.csv` | Derived crop cover summary used in calculations |
| `02_MBH_comb/data/elevation/lowres_full_world.tif` | Global elevation context for an optional figure |
| `02_MBH_comb/data/crop_cover/lowres_crop.tif` | Global crop cover context for an optional figure |
| `03_corrected_HDR/empirical/01_carnicer/data/elevation/catalunya_lowres.tif` | Catalonia elevation context for an optional figure |

The two CSV files are archived derived calculation inputs. The three TIFF files provide map context. They are not collectively described as original raw data. See the [optional workflow guide](../PELLETT_VALBUENA_WORKFLOW.md) for the command and generated products.

## Required paths and ordering

Inputs are under [`data/external/pellett_valbuena_2025/`](external/pellett_valbuena_2025/README.md). The numbered folder names come from the source project, and the scripts expect those paths to remain in place.

Do not rename required paths, reorder `coordinates.csv`, or rename or reorder the numbered elevation rasters. The manuscript workflow pairs coordinate rows and rasters by position.

## Provenance and reuse

All 407 input files are redistributed unchanged. Their exact byte sizes and SHA-256 file fingerprints match corresponding members of both versioned Pellett & Valbuena records: [Zenodo v5](https://doi.org/10.5281/zenodo.15438378) and [Zenodo v6](https://doi.org/10.5281/zenodo.17207705). A SHA-256 fingerprint is used to confirm that a file has not changed. Both records declare CC BY 4.0. The [third-party notices](../THIRD_PARTY_NOTICES.md) record the license treatment declared by the deposits, source attribution, and project license boundary.

Cite Pellett & Valbuena (2025), the applicable versioned Zenodo record, and the underlying source appropriate to each input group:

- MacArthur, R. H. & MacArthur, J. W. “On Bird Species Diversity.” *Ecology* **42**, 594–598 (1961). https://doi.org/10.2307/1932254.
- Allouche, O. et al. “Area–heterogeneity tradeoff and the diversity of ecological communities.” *PNAS* **109**, 17495–17500 (2012).
- Estrada, J., Pedrocchi, V., Brotons, L. & Herrando, S. *Atles Dels Ocells Nidificants de Catalunya 1999–2002*. Institut Català d'Ornitologia/Lynx Editions, Barcelona (2004).
- Carnicer, J., Brotons, L., Herrando, S. & Sol, D. “Improved empirical tests of area-heterogeneity tradeoffs.” *PNAS* **110**, E2858–E2860 (2013). https://doi.org/10.1073/pnas.1222681110.
- Jarvis, A., Reuter, H. I., Nelson, A. & Guevara, E. *Hole-filled seamless SRTM data V4*. International Centre for Tropical Agriculture (2008).
- Buchhorn, M. et al. *Copernicus Global Land Service: Land Cover 100m: collection 3: epoch 2019: Globe* (2020). https://doi.org/10.5281/zenodo.3939050.
- Pellett, C. & Valbuena, R. “Disentangling dispersion from mean reveals true heterogeneity-diversity relationships.” *Nature Communications* **16**, 8532 (2025). https://doi.org/10.1038/s41467-025-64287-0.

The project license does not relicense the 407 files listed in `INPUT_MANIFEST.csv`. Project-authored documentation stored beside those inputs remains covered by the applicable license for project documentation. Source data for figures and tables that were derived by the project use the project license only for the project's selection, arrangement, annotations, and newly derived analytical results; that grant does not replace the terms attached to underlying archived values.
