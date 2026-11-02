# Public derivative scope

This v1.2.0 directory is derived from the completed 2026-11-02 calculation package. The original private delivery is retained locally. The public derivative includes the FE scripts, original numerical inputs and computed outputs, meshes, fixed sensitivity scenarios, measurement-tolerance propagation and Figure S1 reproduction material. It does not contain user task plans, correspondence drafts, manuscript-insertion drafts, third-party full-text manuals/articles, long publisher excerpts or downloaded website code.

The input CSV was relocated from an internal task-bundle directory to `data/`; its numerical content is unchanged and matches the independently extracted public Table 1. `run_validation.py` and `validate_outputs.py` use that public path. The validator no longer requires an internal manuscript-insertion draft. No FE equations, solver settings, target values, predictions or fitted parameters were changed for publication. Source documents remain available through the DOI and official links in `DATA_SOURCES.md`.

The original `INPUT_PROTOCOL.json`, saved freeze record and run amendments remain unchanged. These document the actual retrospective execution, not independent historical blinding. Existing early-stage check records are retained with their original scope; the final comparison contains 39 rows. `PUBLIC_FILE_MANIFEST.csv` records the included file paths and sizes. The package digest is supplied in the release receipt, outside the ZIP, to avoid a self-referential digest.

The full v1.2.0 repository package also retains the previous layered-rock calculation archive. This standalone directory can reproduce the external comparison with its stated Python environment and does not require COMSOL or MATLAB.
