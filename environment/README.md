# Tested R environment

This folder records the software environment used to test the R workflows: the R version, add-on package versions, and detailed session information, including the operating system.

- [`R_VERSION.txt`](R_VERSION.txt) records R 4.6.0.
- [`PACKAGE_VERSIONS.csv`](PACKAGE_VERSIONS.csv) records the 14 directly required R packages and the version tested for each one. In that table, `package` is the add-on library name and `tested_version` is the tested version.
- [`SESSION_INFO.txt`](SESSION_INFO.txt) is the detailed technical record of the tested platform, attached packages, and loaded namespaces. A namespace is the internal package environment loaded by R.

## Setup

Install R and the packages listed in the [run instructions](../RUN_INSTRUCTIONS.md). R 4.6.0 is the tested version and the best choice for matching the recorded validation; other versions may work but were not established by that validation. `PACKAGE_VERSIONS.csv` records tested package versions but does not install or lock them. The workflow checks required packages and does not install them automatically.

The `terra` package and graphics or document packages may require system libraries specific to the operating system. Differences in fonts and graphics libraries can change rendered file bytes even when analytical values are unchanged.

## Demonstrated repeatability

Validation used R 4.6.0 with the recorded package versions on Apple Silicon macOS. In that environment, the main workflow created all 42 expected outputs in an isolated copy. All 7,698 numerical entries read from 25 workflow CSV files matched each of two earlier R runs exactly, and the five included figure and table files and five CSV files containing source data remained unchanged at the file level. The [optional workflow guide](../PELLETT_VALBUENA_WORKFLOW.md) describes separate checks of repeatability within R and comparisons of corresponding R and Julia outputs.

## Troubleshooting

Run commands from the repository root: the main project folder containing `README.md` and the marker file `.mi_hdr_project_root`. The marker helps the scripts find the required folders. Run only one workflow at a time. `_local_scratch/.workflow_run_lock` prevents simultaneous runs; if a process ends unexpectedly, inspect that file and confirm that no R workflow is active before removing it manually.
