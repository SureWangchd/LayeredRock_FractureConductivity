"""Summarize existing round-3 results; never solve, fit, or modify raw data.

Run from any directory:
    python scripts/summarize_round3.py
    python scripts/summarize_round3.py --root F:/岩层裂缝/runs/round3_predict

Only successful truth states enter the output tables. Predictions without a
matching truth stay in their original raw/frozen files, not placeholder rows.
Relative increment errors are blank near zero; absolute errors are retained.
Shared background and basis costs are repeated metadata, not additive costs.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
from typing import Any


HOMOGENIZED_DEVELOPMENT_CASE_IDS = (
    "core_h20_L60_A_m3_full",
    "core_h20_L60_B_m3_full",
    "core_h20_L60_C_m3_full",
    "strict_h10_L60_A_m3_full",
    "strict_h10_L60_B_m3_full",
)


def number(value: Any) -> float | None:
    """Keep missing/nonfinite values missing; never turn them into zero."""
    if value is None or value == "" or isinstance(value, bool):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def key(case_id: Any, sigma: Any) -> tuple[str, str] | None:
    sf = number(sigma)
    if not case_id or sf is None:
        return None
    return str(case_id), format(sf, ".12g")


def ratio(numerator: Any, denominator: Any) -> float | None:
    a, b = number(numerator), number(denominator)
    return a / b if a is not None and b is not None and b != 0 else None


def flatten_scalars(data: dict[str, Any], prefix: str) -> dict[str, Any]:
    """Expose scalar local diagnostics; arrays remain in their source JSON."""
    result: dict[str, Any] = {}
    for name, value in data.items():
        field = prefix + name
        if isinstance(value, dict):
            result.update(flatten_scalars(value, field + "_"))
        elif isinstance(value, (str, int, float, bool)) or value is None:
            result[field] = value
    return result


def read_json(path: Path, warnings: list[str]) -> dict[str, Any] | None:
    try:
        data = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError) as exc:
        warnings.append(f"Cannot read {path}: {exc}")
        return None
    if not isinstance(data, dict):
        warnings.append(f"Expected a JSON object: {path}")
        return None
    return data


def read_csv(path: Path, warnings: list[str]) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            return [row for row in csv.DictReader(handle) if any(row.values())]
    except (OSError, csv.Error) as exc:
        warnings.append(f"Cannot read {path}: {exc}")
        return []


def shared_background_metadata(
    folder: Path, warnings: list[str]
) -> dict[str, Any]:
    path = folder / "background_summary.json"
    if not path.is_file():
        return {}
    data = read_json(path, warnings)
    if data is None:
        return {}
    out = {
        "background_endcap_dV_V": data.get("cap_dV_V"),
        "background_dofs": data.get("dofs"),
        "background_free_dofs": data.get("free_dofs"),
        "shared_background_cost_repeated_across_targets": True,
    }
    out.update(flatten_scalars(data.get("cost", {}), "shared_background_cost_"))
    return out


def load_truths(
    root: Path, scope: str, warnings: list[str]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    shared: dict[Path, dict[str, Any]] = {}
    seen: dict[tuple[str, str], Path] = {}
    for path in sorted((root / scope).glob("*/sf_*/SUCCESS.json")):
        data = read_json(path, warnings)
        if data is None:
            continue
        if data.get("status", "success").lower() != "success":
            warnings.append(f"Skipping non-success record: {path}")
            continue
        required = ("G0_S", "G_S", "deltaG_S", "U_S", "R_S", "cross_S")
        ident = key(data.get("case_id"), data.get("sigma_f_S_m"))
        if ident is None or any(number(data.get(name)) is None for name in required):
            warnings.append(f"Skipping incomplete success record: {path}")
            continue
        if ident in seen:
            raise ValueError(
                f"Ambiguous truth for {ident}: {seen[ident]} and {path}; "
                "use distinct case IDs instead of silently choosing a result."
            )
        seen[ident] = path
        folder = path.parent.parent
        if folder not in shared:
            shared[folder] = shared_background_metadata(folder, warnings)
        row = {name: value for name, value in data.items() if not isinstance(value, (dict, list))}
        for name in required + ("sigma_f_S_m",):
            row[name] = number(data[name])
        row.update(shared[folder])
        local = data.get("local", {})
        if isinstance(local, dict):
            row.update(flatten_scalars(local, "local_"))
        row.update(scope=scope, truth_success_file=path.relative_to(root).as_posix())
        rows.append(row)
    return sorted(rows, key=lambda row: (row["case_id"], row["sigma_f_S_m"]))


def energy_rows(truths: list[dict[str, Any]], near_zero: float) -> list[dict[str, Any]]:
    rows = []
    for truth in truths:
        row = {k: v for k, v in truth.items() if not k.startswith("local_")}
        actual, upper, correction = truth["deltaG_S"], truth["U_S"], truth["R_S"]
        residual = upper - correction - actual
        cross_residual = truth["cross_S"] - actual
        row.update(
            U_minus_R_recomputed_S=upper - correction,
            identity_residual_S=residual,
            identity_absolute_error_S=abs(residual),
            identity_relative_error_unfloored=abs(residual) / abs(actual) if abs(actual) > near_zero else None,
            cross_residual_S=cross_residual,
            cross_absolute_error_S=abs(cross_residual),
            cross_relative_error_unfloored=abs(cross_residual) / abs(actual) if abs(actual) > near_zero else None,
            deltaG_near_zero=abs(actual) <= near_zero,
            relative_error_near_zero_threshold_S=near_zero,
            cancellation_factor_U_R_over_delta=(abs(upper) + abs(correction)) / abs(actual) if abs(actual) > near_zero else None,
        )
        rows.append(row)
    return rows


def response_rows(truths: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Attach backward adjacent-point d ln(deltaG) / d ln(sigma_f)."""
    rows = []
    previous: dict[str, dict[str, Any]] = {}
    for truth in truths:
        row = dict(truth)
        row["relative_gain"] = ratio(truth["deltaG_S"], truth["G0_S"])
        row["relative_gain_percent"] = 100 * row["relative_gain"] if row["relative_gain"] is not None else None
        prev = previous.get(truth["case_id"])
        if prev is None:
            row["adjacent_log_sensitivity_status"] = "first_available_state"
        else:
            s0, s1 = prev["sigma_f_S_m"], truth["sigma_f_S_m"]
            d0, d1 = prev["deltaG_S"], truth["deltaG_S"]
            row.update(
                previous_sigma_f_S_m=s0,
                previous_deltaG_S=d0,
                adjacent_deltaG_change_S=d1 - d0,
            )
            if s1 > s0 > 0 and d0 > 0 and d1 > 0:
                row["adjacent_dln_deltaG_dln_sigma_f"] = math.log(d1 / d0) / math.log(s1 / s0)
                row["adjacent_log_sensitivity_status"] = "computed_between_available_positive_states"
            else:
                row["adjacent_log_sensitivity_status"] = "undefined_for_nonpositive_gain_or_sigma"
        previous[truth["case_id"]] = truth
        rows.append(row)
    return rows


def comparison_rows(
    predictions: list[tuple[dict[str, Any], Path]],
    truths: list[dict[str, Any]],
    root: Path,
    near_zero: float,
    warnings: list[str],
) -> tuple[list[dict[str, Any]], int]:
    lookup = {key(row["case_id"], row["sigma_f_S_m"]): row for row in truths}
    rows, unmatched = [], 0
    for prediction, path in predictions:
        ident = key(prediction.get("case_id"), prediction.get("sigma_f_S_m"))
        estimate = number(prediction.get("predicted_deltaG_S"))
        if ident is None or estimate is None:
            warnings.append(f"Skipping incomplete prediction row in {path}")
            continue
        truth = lookup.get(ident)
        if truth is None:
            unmatched += 1
            continue
        actual, g0 = truth["deltaG_S"], truth["G0_S"]
        error = estimate - actual
        row = dict(prediction)  # Preserve method/version and additional cost fields.
        row.update(
            case_id=truth["case_id"],
            sigma_f_S_m=truth["sigma_f_S_m"],
            predicted_deltaG_S=estimate,
            truth_deltaG_S=actual,
            truth_G0_S=g0,
            truth_G_S=truth["G_S"],
            truth_U_S=truth["U_S"],
            truth_R_S=truth["R_S"],
            increment_error_S=error,
            increment_absolute_error_S=abs(error),
            increment_signed_relative_error=error / abs(actual) if abs(actual) > near_zero else None,
            increment_relative_error=abs(error) / abs(actual) if abs(actual) > near_zero else None,
            increment_relative_error_percent=100 * abs(error) / abs(actual) if abs(actual) > near_zero else None,
            increment_error_over_background=ratio(error, abs(g0)),
            increment_absolute_error_over_background=ratio(abs(error), abs(g0)),
            truth_relative_gain=ratio(actual, g0),
            predicted_relative_gain=ratio(estimate, g0),
            deltaG_near_zero=abs(actual) <= near_zero,
            relative_error_near_zero_threshold_S=near_zero,
            prediction_minus_truth_S=error,
            upper_U_minus_prediction_S=truth["U_S"] - estimate,
            truth_solve_seconds=truth.get("solve_seconds"),
            truth_postprocess_seconds=truth.get("postprocess_seconds"),
            truth_reused_seen_sigma10=truth.get("reused_seen_sigma10"),
            truth_source=truth.get("source"),
            truth_success_file=truth["truth_success_file"],
            prediction_source_file=path.relative_to(root).as_posix(),
            comparison_status="matched_existing_success",
        )
        if number(prediction.get("G0_S")) is not None:
            row["prediction_G0_minus_truth_G0_S"] = number(prediction["G0_S"]) - g0
        for name, value in truth.items():
            if name.startswith("shared_background_cost_"):
                row[name] = value
        rows.append(row)
    return rows, unmatched


def homogenized_comparison_rows(
    root: Path,
    development: list[dict[str, Any]],
    near_zero: float,
    warnings: list[str],
) -> tuple[list[dict[str, Any]], int]:
    """Reuse one solved anisotropic-reference prediction for the five cases.

    This is a deliberately common baseline, not target-specific calibration.
    The fixed list prevents accidentally applying this matching assumption to
    later wavy backgrounds or unrelated development geometries.
    """
    path = root / "homogenized_reference" / "homogenized_predictions.csv"
    if not path.is_file():
        return [], 0
    references: dict[str, dict[str, Any]] = {}
    for row in read_csv(path, warnings):
        ident = key("common_reference", row.get("sigma_f_S_m"))
        if ident is None or number(row.get("predicted_deltaG_S")) is None:
            warnings.append(f"Skipping incomplete homogenized reference row in {path}")
            continue
        if ident[1] in references:
            raise ValueError(f"More than one common homogenized prediction at sigma_f={ident[1]} in {path}")
        references[ident[1]] = row
    predictions = []
    eligible = [row for row in development if row["case_id"] in HOMOGENIZED_DEVELOPMENT_CASE_IDS]
    for truth in eligible:
        sf_key = key(truth["case_id"], truth["sigma_f_S_m"])[1]
        reference = references.get(sf_key)
        if reference is None:
            continue
        prediction = dict(reference)
        prediction.update(
            case_id=truth["case_id"],
            method="homogenized_directional_reference",
            common_reference_case_id="common_homogenized_reference",
            identical_reference_prediction_for_all_five_backgrounds=True,
            reference_solve_seconds_shared=reference.get("solve_seconds"),
            reference_background_solve_seconds_shared=reference.get("background_solve_seconds_shared"),
            directional_matching_scope="Common x/z responses of the five prescribed flat development backgrounds; not a full tensor match and no target-specific calibration.",
            prediction_cost_scope="Full COMSOL homogeneous-anisotropic reference solve, reused across explicit backgrounds; not a background reduced-order solve. Shared reference cost must not be counted five times.",
        )
        predictions.append((prediction, path))
    rows, _ = comparison_rows(predictions, eligible, root, near_zero, warnings)
    return rows, len(references)


def write_nonempty_csv(path: Path, rows: list[dict[str, Any]], preferred: list[str]) -> int:
    if not rows:
        return 0
    available = set().union(*(row.keys() for row in rows))
    fields = [name for name in preferred if name in available]
    fields.extend(sorted(available - set(fields)))
    temporary = path.with_suffix(path.suffix + ".tmp")
    try:
        with temporary.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return len(rows)


def summarize(root: Path, near_zero: float = 1e-12) -> dict[str, Any]:
    warnings: list[str] = []
    development = load_truths(root, "development", warnings)
    prediction_files = sorted((root / "development").glob("*/predictor_development_raw.csv"))
    prediction_files += sorted((root / "development").glob("*/predictor_enriched_raw.csv"))
    prediction_files += sorted((root / "development").glob("*/predictor_enriched_rankguard.csv"))
    raw = [(row, path) for path in prediction_files for row in read_csv(path, warnings)]
    comparisons, unmatched = comparison_rows(raw, development, root, near_zero, warnings)
    truth_preferred = ["case_id", "sigma_f_S_m", "G0_S", "G_S", "deltaG_S", "relative_gain", "U_S", "R_S", "U_minus_R_S", "cross_S", "identity_relative_error"]
    prediction_preferred = ["case_id", "method", "segments", "candidate_columns", "rank", "sigma_f_S_m", "predicted_deltaG_S", "truth_deltaG_S", "increment_error_S", "increment_relative_error", "increment_absolute_error_over_background"]
    outputs = {
        "energy_decomposition.csv": energy_rows(development, near_zero),
        "conductivity_response.csv": response_rows(development),
        "predictor_development.csv": comparisons,
    }
    counts = {}
    for filename, rows in outputs.items():
        preferred = prediction_preferred if filename.startswith("predictor_") else truth_preferred
        counts[filename] = write_nonempty_csv(root / filename, rows, preferred)
        if not rows and (root / filename).exists():
            warnings.append(f"No current rows for {filename}; an older output exists and was not changed.")

    homogenized_rows, reference_count = homogenized_comparison_rows(root, development, near_zero, warnings)
    counts["homogenized_comparison.csv"] = write_nonempty_csv(
        root / "homogenized_comparison.csv", homogenized_rows, prediction_preferred
    )
    if not homogenized_rows and (root / "homogenized_comparison.csv").exists():
        warnings.append("No current matched homogeneous-reference rows; existing homogenized_comparison.csv was not changed.")

    # Test truths are deliberately joined only to frozen predictions, not raw
    # development predictions. This summary does not certify freeze timing.
    test = load_truths(root, "test", warnings)
    frozen_files = [root / "predictions_frozen.csv", root / "test" / "predictions_frozen.csv"]
    frozen_files += sorted((root / "test").glob("*/predictions_frozen.csv"))
    frozen = [(row, path) for path in frozen_files if path.is_file() for row in read_csv(path, warnings)]
    test_rows, test_unmatched = comparison_rows(frozen, test, root, near_zero, warnings)
    for row in test_rows:
        row["prediction_provenance_note"] = "Joined from predictions_frozen.csv; this summarizer does not independently certify freeze-before-truth timing."
    counts["prediction_test.csv"] = write_nonempty_csv(root / "prediction_test.csv", test_rows, prediction_preferred)
    if not test_rows and (root / "prediction_test.csv").exists():
        warnings.append("No current matched test rows; existing prediction_test.csv was not changed.")
    return {
        "root": str(root),
        "written_rows": counts,
        "development_success_states": len(development),
        "development_prediction_rows_without_truth": unmatched,
        "homogenized_reference_states_available": reference_count,
        "test_success_states": len(test),
        "frozen_prediction_rows_without_truth": test_unmatched,
        "empty_outputs": "Not created; no placeholder result rows.",
        "sensitivity_definition": "Backward adjacent available states: ln(deltaG_i/deltaG_previous)/ln(sigma_i/sigma_previous); only for positive gains and conductivities.",
        "cost_note": "Shared background costs and per-basis preparation costs repeat across target rows; do not sum them across sigma values.",
        "warnings": warnings,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1] / "runs" / "round3_predict")
    parser.add_argument("--near-zero-deltaG-S", type=float, default=1e-12, dest="near_zero")
    args = parser.parse_args()
    if not args.root.is_dir():
        parser.error(f"Result directory does not exist: {args.root}")
    if args.near_zero < 0 or not math.isfinite(args.near_zero):
        parser.error("The near-zero increment threshold must be finite and nonnegative.")
    report = summarize(args.root.resolve(), args.near_zero)
    # ASCII escaping keeps console output portable on Windows; CSV is UTF-8.
    print(json.dumps(report, ensure_ascii=True, indent=2))


if __name__ == "__main__":
    main()
