# Third-party notices

## Scope

This repository redistributes 407 archived input files unchanged. [`data/INPUT_MANIFEST.csv`](data/INPUT_MANIFEST.csv) is the file list used to identify and verify them; its 407 rows define the exact file-level scope:

- 14 MacArthur foliage CSV files: `A.csv` through `M.csv` and `foliage_density.csv`;
- two Catalonia bird and coordinate CSV files;
- 386 Catalonia elevation TIFF files, `0.tif` through `385.tif`;
- two heterogeneity summary CSV files used only by the optional workflow; and
- three TIFF files used only by the optional workflow to provide map context.

The exclusion of these 407 files from the project license applies to these specific files. It does not apply to project-authored documentation stored beside the inputs, including `data/external/pellett_valbuena_2025/README.md`.

## Pellett & Valbuena deposits

The 407 files retain their verified bytes, filenames, and required folder locations within the repository. Their exact byte sizes and SHA-256 checksums (file fingerprints used to confirm unchanged files) match corresponding members of both:

- Pellett, C. & Valbuena, R. (2025). *Data and code – Disentangeling dispersion from mean reveals true heterogeneity-diversity relationships* (Version 5) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.15438378.
- Pellett, C. & Valbuena, R. (2025). *Data and code – Disentangling dispersion from mean reveals true heterogeneity-diversity relationships* (Version 6) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.17207705.

Both cited Zenodo deposit records declare the deposited material under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/). This repository relies on that deposit-declared treatment and does not claim that the project independently reconstructed every upstream permission. Cite the related article:

Pellett, C. & Valbuena, R. “Disentangling dispersion from mean reveals true heterogeneity-diversity relationships.” *Nature Communications* **16**, 8532 (2025). https://doi.org/10.1038/s41467-025-64287-0.

## Attribution for underlying sources

Retain the attribution appropriate to each input group:

- MacArthur foliage inputs: MacArthur, R. H. & MacArthur, J. W. “On Bird Species Diversity.” *Ecology* **42**, 594–598 (1961). https://doi.org/10.2307/1932254.
- Catalonia bird and coordinate inputs: Estrada, J., Pedrocchi, V., Brotons, L. & Herrando, S. *Atles Dels Ocells Nidificants de Catalunya 1999–2002*. Institut Català d'Ornitologia/Lynx Editions, Barcelona (2004); and Carnicer, J., Brotons, L., Herrando, S. & Sol, D. “Improved empirical tests of area-heterogeneity tradeoffs.” *PNAS* **110**, E2858–E2860 (2013). https://doi.org/10.1073/pnas.1222681110. The deposited source documentation acknowledges Jofre Carnicer, Lluís Brotons, and the Institut Català d'Ornitologia.
- Elevation inputs derived from SRTM: Jarvis, A., Reuter, H. I., Nelson, A. & Guevara, E. *Hole-filled seamless SRTM data V4*. International Centre for Tropical Agriculture (2008). https://srtm.csi.cgiar.org/srtmdata/.
- Crop cover inputs: Buchhorn, M. et al. *Copernicus Global Land Service: Land Cover 100m: collection 3: epoch 2019: Globe* (2020). https://doi.org/10.5281/zenodo.3939050.

## License boundary

The project's license grant does not relicense the 407 input files listed in `data/INPUT_MANIFEST.csv`. Those files remain attributed to their sources and are redistributed in reliance on the CC BY 4.0 treatment declared by the cited Pellett & Valbuena deposits.

Project-authored R code, repository configuration, documentation, the five included products from the manuscript analysis, environment records, and derived source data for figures and tables are covered only by the separately stated project license scopes. For source data derived by the project, that grant covers the project's selection, arrangement, annotations, and newly derived analytical results. It does not replace the source terms applying to underlying archived values.

This repository verified that all 407 files match their deposited counterparts exactly and preserves the supplied source citations. It did not independently reconstruct every upstream permission.
