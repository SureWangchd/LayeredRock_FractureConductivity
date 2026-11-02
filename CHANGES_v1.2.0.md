# Version 1.2.0 — external electrical comparison

Version DOI: `10.5281/zenodo.23879386`. The v1.1.0 and v1.0.0 releases, their sources and 13 representative MPH files remain accessible and are not overwritten.

This version adds the completed `external_electrical_validation/` calculation: published Han2022 measurement inputs; fitted-to-boundary tetrahedral FE code; saved cell and cylinder meshes/results; nominal, fixed-scenario and measurement-tolerance comparisons; four-level mesh trends; and Figure S1 reproduction files. The external calculation uses Python/Gmsh/SciPy, separate from the original COMSOL studies. No new FE solve or target adjustment was performed merely to package these records.

This is a public derivative of the completed local calculation package. The numerical input CSV is relocated to `data/`; two scripts use that path, and the output checker no longer requires an internal manuscript-insertion draft. Numerical algorithms, calibration/holdout definitions, input values and results are unchanged. User task plans, correspondence drafts, manuscript work notes, publisher excerpts, full third-party manuals and downloaded website code are omitted. Authoritative DOI/manual links and the independently checked numerical transcription are retained. The original private delivery remains local.

The comparison is retrospective and has the limits stated in its README and assessment: no same-core pre-fracture measurement is supplied, and the strict A/B peripheral-configuration claims are not directly tested experimentally. Earlier S1 CT check limitations and the 80 original prospective numerical prediction states are retained unchanged.

The standalone external ZIP and the complete v1.2.0 source/data ZIP contain the same public external directory. Their published download digests are recorded in the release verification receipt. `package_inventory_v1.2.0.csv` and `package_summary_v1.2.0.json` describe this version; earlier inventories remain historical records.
