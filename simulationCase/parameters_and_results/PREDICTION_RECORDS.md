# Original prospective prediction records

This directory supplements the aggregate response tables with the original 48 round-5 and 32 round-6 per-state `prediction.json` files. The records are under `runs/round5_modal_closure/predictions/` and `runs/round6_applicability/predictions/`. Their method files, complete-batch `predictions_ready.json` markers and completion markers are included beside them. Round 6 also includes its saved field-readiness marker and pre-truth side/rank intervals.

`prediction_record_sources.csv` maps every included record back to its project-relative source path, states its original timestamp where present, and gives its byte size. JSON contents are preserved as saved; original local paths inside a record are historical provenance, not portable entry points. Use the archive-relative paths in the manifest to locate files. Local timestamps without an explicit UTC offset use Asia/Shanghai (UTC+08:00).

The round-5 complete batch was marked ready at `2026-10-26 01:07:38`; the round-6 batch at `2026-10-26T16:38:39.080619+08:00`. These saved workflow records document the intended prediction-before-target sequence. They are not an independently timestamped preregistration or a proof immune to later editing. A target completion marker records completion, not the starting time of every individual target solve.

Only the initial prospective batches are counted here: 48 plus 32 states. Development cases and adaptive crossover follow-ups are not relabeled as independent first-batch predictions. Original R5/R6 method versions, finite-model limits and failures remain unchanged. This archival supplement does not rerun predictions, refit methods or change numerical results.

The complete batches and batch markers are added together for the versioned archive update. Four R6 prediction files were already present among the figure-specific inputs in v1.0.0; those historical plotting copies remain in place and are not four additional independent states. The canonical batch paths above contain 80 distinct states. Code is GPL-3.0-only; prediction records and other project data are CC BY 4.0 as described in the repository license files.
