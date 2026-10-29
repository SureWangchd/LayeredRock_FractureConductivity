# Beyond Bulk Electrical Equivalence

Models, data and plotting scripts for *Beyond Bulk Electrical Equivalence: Structural Control of Fracture-Induced Conductivity Enhancement in Layered Rocks*.

Version 1.1.0: [10.5281/zenodo.23768041](https://doi.org/10.5281/zenodo.23768041). The unchanged version 1.0.0 remains available at [10.5281/zenodo.23755034](https://doi.org/10.5281/zenodo.23755034).

### Numerical models

`solver/` contains the MATLAB and Python source for layered two-dimensional models, finite three-dimensional ellipsoids, and CT-based sandstone models. The calculations use COMSOL Multiphysics 6.3.0.290 with LiveLink for MATLAB R2025b. The saved environment record is [`runtime_versions.json`](simulationCase/parameters_and_results/runs/round7b_localfit/runtime_versions.json); its executable paths describe the original machine and should be adjusted on another installation. Entry functions and their required inputs are listed in `solver/entry_points.csv`.

Representative solved COMSOL models are supplied as four ZIP files in the [v1.0.0 release](https://github.com/SureWangchd/LayeredRock_FractureConductivity/releases/tag/v1.0.0): the two-dimensional configuration and intervention pairs, four finite three-dimensional states, and two refined CT pairs. The four archives contain 13 solved MPH files: five two-dimensional, four synthetic three-dimensional, and four CT files. The original `model_inventory.csv` retains its v1.0.0 figure labels; `model_inventory_manuscript_20261029.csv` adds the current manuscript labels without changing the models.

### Data and figures

`simulationCase/` contains the numerical tables, sampled fields and MATLAB plotting scripts. `simulationCase/figure_index.csv` maps these directories to the twelve figures in the 2026-10-29 manuscript layout. `figure_index_v1.0.0.csv` preserves the released index, and `figure_number_map_20261029.csv` documents the old-to-current numbering. For example, run `simulationCase/Fig02/run_replot.m` in MATLAB to reproduce the connection-intervention figure (current Figure 3). Each figure directory contains its input data and plotting entry point. The three figures in `Shared_Fig01_04_07/` share a plotting script and are current Figures 1, 6 and 8. Directory and script names are historical identifiers, not current manuscript figure numbers; current Figure 2 uses `Fig05/run_replot.m`. Source paths and numerical data are retained. The label in current Figure 6e is corrected to "Peripheral power (NN)".

`simulationCase/parameters_and_results/` contains the case definitions and numerical response tables used in the study. The S1 and S5 source images are available from the Imperial College Consortium on Pore-scale Imaging and Modelling at [10.6084/m9.figshare.1189274.v1](https://doi.org/10.6084/m9.figshare.1189274.v1) and [10.6084/m9.figshare.1189281.v1](https://doi.org/10.6084/m9.figshare.1189281.v1).

The original 80 prospective prediction JSON records, batch-ready markers and their source mapping are described in [`PREDICTION_RECORDS.md`](simulationCase/parameters_and_results/PREDICTION_RECORDS.md). Selected core-mesh comparisons, finite-element trace checks and small-effect diagnostics are available in [`numerical_clarifications_20261029`](simulationCase/parameters_and_results/numerical_clarifications_20261029/README.md). The original `file_inventory.csv` and `content_summary.json` describe the v1.0.0 assembly; the new version's package inventory identifies the complete updated contents. See [`CHANGES_v1.1.0.md`](CHANGES_v1.1.0.md) for this revision.

The twelve current manuscript caption sources are in [`manuscript_captions_20261029`](simulationCase/manuscript_captions_20261029/README.md). The bounded [CT fixed-mesh review](simulationCase/parameters_and_results/round7b_review_20261029/README.md) preserves the actual status and logs: S1 produced **0/4** accepted new stricter states, so its four-state solver-change and full combined scales remain unknown; S5 has **4/4** previously accepted stricter states and a combined observed scale of 0.004527284519 microS. Completing the review or packaging the records does not mean that S1 passed a stricter numerical check. No new multi-GB MPH files are added.

### License and support

Project code is distributed under GPL-3.0-only; project data, prediction records and MPH model files are distributed under CC BY 4.0. Original CT images retain their source attribution and license conditions. COMSOL and MATLAB are required to run the numerical models and are not included.

This work was supported by the National Natural Science Foundation of China (Grant 42374154).
