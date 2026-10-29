# Numerical clarifications for the 2026-10-29 revision

These tables make selected evidence behind the numerical qualifications directly available. They do not add new target solutions or refit predictions. The original computation files remain unchanged; `source_index.csv` records the project-relative extraction source for each table.

The core-mesh and quadratic-trace tables distinguish the R5 and R6 cases: common physical core geometry does not by itself imply the same finite-element discretization. The observed R5 predictor and target core meshes differ, while the checked R6 core quadratic elements and interface traces coincide. The tables, rather than a blanket common-mesh assertion, define the comparisons.

The other tables retain all 64 R5 M8 single-side checks, the eight R6 joint contrasts, the small B cancellation at 0.18 S/m, the synthetic 3D marker scale and the specified base/refined mesh sizes. Observed changes on mesh refinement are numerical resolution diagnostics, not rigorous continuum error bounds. Prediction timestamps are saved workflow records, not independent preregistration.

`modal_detail_extract.py` uses the Python standard library. `modal_evidence_audit.py` additionally requires NumPy, SciPy and h5py. Both take `--project-root PATH --output-dir PATH`, read saved outputs and do not launch COMSOL. The full core connectivity check requires original `reference/mesh.mat`, `reference.mat`, target `source_background/mesh.mat` and `background_responses.mat` files. The compact public release contains the resulting comparison tables and representative solved MPH files, not every auxiliary matrix from the full project. Missing source files therefore raise an error; they are not treated as matching meshes.

Example against a complete local computation tree:

```text
python modal_evidence_audit.py --project-root /path/to/full/project --output-dir /path/to/check-output
python modal_detail_extract.py --project-root /path/to/full/project --output-dir /path/to/check-output
```

Source paths inside the copied tables are historical local provenance. The data are CC BY 4.0 and the scripts are GPL-3.0-only, consistent with the repository licenses.
