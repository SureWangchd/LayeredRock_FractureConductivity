# Numerical source entry points

`scripts/` contains project-owned model-building, solving and numerical-postprocessing source. COMSOL is not distributed. Steady electric-current models used COMSOL 6.3 with MATLAB LiveLink. The exact source entry points and expected run directories are listed in `entry_points.csv`.

The compact stage manifests, parameters and original response CSVs are in `../simulationCase/parameters_and_results/`. They retain the original relative `runs/<stage>` names. Figure reproduction uses the saved input groups and requires no COMSOL solve. Continuing a numerical stage also requires its prerequisite MPH, FE operators or CT mesh as stated in the entry table; these compact records alone do not constitute a fresh full-study solver workspace.

For a solver function accepting `root`, supply a working project root containing its listed `runs/` directory and add `solver/scripts` to the MATLAB or Python path. Legacy runners that derive their root from their own location expect `runs/` alongside `scripts/`. Configure the local COMSOL installation and an available LiveLink connection before using those runners. Source paths in JSON/CSV records identify the original computation.

Final CT results use `r7b_fix_*` local conforming cutting and fixed-interface refinement. `r7_ct_native_solver.m` and its phase-validation helper are retained because the final CT solver calls these generic imported-mesh electrical routines. The retired all-surface Union/FreeTet and earlier centroid-labelled fracture launch routes are not included.

The small `round5_taskpack/` and `round6_taskpack/` directories contain the analytical coefficient/interval modules and numerical design JSONs imported by the project scripts, plus their required algebraic self-tests. They contain no task prose or literature materials.
