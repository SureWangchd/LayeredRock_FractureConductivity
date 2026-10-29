# Beyond Bulk Electrical Equivalence

Models, data and plotting scripts for *Beyond Bulk Electrical Equivalence: Structural Control of Fracture-Induced Conductivity Enhancement in Layered Rocks*.

Version 1.0.0: [10.5281/zenodo.23755034](https://doi.org/10.5281/zenodo.23755034).

### Numerical models

`solver/` contains the MATLAB and Python source for layered two-dimensional models, finite three-dimensional ellipsoids, and CT-based sandstone models. The calculations use COMSOL Multiphysics 6.3 with LiveLink for MATLAB. Entry functions and their required inputs are listed in `solver/entry_points.csv`.

Representative solved COMSOL models are supplied as four ZIP files in the [v1.0.0 release](https://github.com/SureWangchd/LayeredRock_FractureConductivity/releases/tag/v1.0.0): the two-dimensional configuration and intervention pairs, four finite three-dimensional states, and two refined CT pairs. The archive inventory identifies each model and its parameters.

### Data and figures

`simulationCase/` contains the numerical tables, sampled fields and MATLAB plotting scripts. `simulationCase/figure_index.csv` maps these directories to the twelve manuscript figures. For example, run `simulationCase/Fig02/run_replot.m` in MATLAB to reproduce the connection-intervention figure. Each figure directory contains its input data and plotting entry point. The three figures in `Shared_Fig01_04_07/` share a plotting script.

`simulationCase/parameters_and_results/` contains the case definitions and numerical response tables used in the study. The S1 and S5 source images are available from the Imperial College Consortium on Pore-scale Imaging and Modelling at [10.6084/m9.figshare.1189274.v1](https://doi.org/10.6084/m9.figshare.1189274.v1) and [10.6084/m9.figshare.1189281.v1](https://doi.org/10.6084/m9.figshare.1189281.v1).

### License and support

Project code is distributed under GPL-3.0; project data and model files are distributed under CC BY 4.0. Original CT images retain their source attribution. COMSOL and MATLAB are required to run the numerical models and are not included.

This work was supported by the National Natural Science Foundation of China (Grant 42374154).
