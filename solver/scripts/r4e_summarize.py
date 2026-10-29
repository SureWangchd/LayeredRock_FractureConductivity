"""Summarize the bounded round-4 factorial extension; no solver or old writes.

Missing records remain explicit. Observed mesh changes are diagnostics, not
rigorous FEM error estimates. Uses only Python standard-library modules.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / "runs" / "round4_extension_20261024"
LAYOUTS = ("NN", "NF", "FN", "FF")
QUANTITIES = ("gA", "gB", "D", "SfA", "SfB", "M")
EFFECTS = {
    "lower_near_at_upper_near": {"NN": 1, "FN": -1},
    "lower_near_at_upper_far": {"NF": 1, "FF": -1},
    "upper_near_at_lower_near": {"NN": 1, "NF": -1},
    "upper_near_at_lower_far": {"FN": 1, "FF": -1},
    "lower_average": {"NN": .5, "FN": -.5, "NF": .5, "FF": -.5},
    "upper_average": {"NN": .5, "NF": -.5, "FN": .5, "FF": -.5},
    "interaction": {"NN": 1, "NF": -1, "FN": -1, "FF": 1},
    "joint_NN_minus_FF": {"NN": 1, "FF": -1},
}
MESH_NOTE = (
    "Absolute observed fine-minus-coarse changes summed over participating "
    "underlying A/B states with absolute contrast weights; cancellation is not "
    "used as an error bound. One mesh trend only, not a rigorous FEM error bound."
)


def read(path: Path, notes: list[str]) -> dict:
    if not path.is_file():
        return {}
    try:
        obj = json.loads(path.read_text(encoding="utf-8-sig"))
        if not isinstance(obj, dict):
            raise ValueError("expected JSON object")
        return obj
    except (OSError, ValueError) as exc:
        notes.append(f"Record temporarily unavailable: {path}: {exc}")
        return {}


def finite(value):
    try:
        val = float(value)
    except (TypeError, ValueError):
        return None
    return val if math.isfinite(val) else None


def difference(left, right):
    return left - right if left is not None and right is not None else None


def ratio(numerator, denominator):
    return numerator / denominator if numerator is not None and denominator not in (None, 0) else None


def token(sf):
    return "sf_" + format(sf, ".10g").replace(".", "p")


def write_csv(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_response(layout, role, sf, folder, reused, case_id, out, notes):
    target = folder / token(sf)
    success = read(target / "SUCCESS.json", notes)
    marginal = read(target / "marginal_scalar.json", notes)
    prediction = read(target / "prediction.json", notes)
    valid = success.get("status") == "success" and all(finite(success.get(k)) is not None for k in ("sigma_f_S_m", "G0_S", "G_S", "deltaG_S"))
    if valid and not math.isclose(float(success["sigma_f_S_m"]), sf, rel_tol=1e-13, abs_tol=1e-14):
        raise ValueError(f"sigma mismatch: {target}")
    if success and success.get("case_id") != case_id:
        raise ValueError(f"case mismatch: {target}: {success.get('case_id')} != {case_id}")
    row = {
        "layout": layout, "role": role, "case_id": case_id, "sigma_f_S_m": sf,
        "status": "success" if valid else "pending_or_incomplete",
        "reused_round4_state": reused,
        "source_folder": folder.as_posix(), "truth_file": (target / "SUCCESS.json").as_posix(),
        "G0_S": finite(success.get("G0_S")) if valid else None,
        "G_S": finite(success.get("G_S")) if valid else None,
        "deltaG_S": finite(success.get("deltaG_S")) if valid else None,
        "relative_gain": finite(success.get("relative_gain")) if valid else None,
        "Sf_m": finite(marginal.get("Sf_m")) if valid else None,
        "U_S": finite(success.get("U_S")) if valid else None,
        "R_S": finite(success.get("R_S")) if valid else None,
        "reaction_flux_relative_error": finite(success.get("reaction_flux_relative_error")) if valid else None,
        "energy_relative_error": finite(success.get("energy_relative_error")) if valid else None,
        "identity_relative_error": finite(success.get("identity_relative_error")) if valid else None,
        "solve_seconds": finite(success.get("solve_seconds")) if valid else None,
        "scope": "Known round-4 results informed this mechanism experiment; not an independent blind geometry test.",
        "Sf_scope": "Actual target-field posthoc derivative; not background prediction.",
    }
    fine_folder = out / "precision" / f"{role}_{layout}_{token(sf)}"
    fine = read(fine_folder / "SUCCESS.json", notes)
    fmarginal = read(fine_folder / "marginal_scalar.json", notes)
    fine_valid = fine.get("status") == "success" and finite(fine.get("fine_deltaG_S")) is not None
    if fine_valid:
        if not math.isclose(float(fine["sigma_f_S_m"]), sf, rel_tol=1e-13, abs_tol=1e-14):
            raise ValueError(f"fine sigma mismatch: {fine_folder}")
        if fine.get("case_id") != case_id:
            raise ValueError(f"fine case mismatch: {fine_folder}")
    row.update({
        "fine_status": "success" if fine_valid else "not_available",
        "fine_deltaG_S": finite(fine.get("fine_deltaG_S")) if fine_valid else None,
        "fine_G0_S": finite(fine.get("fine_G0_S")) if fine_valid else None,
        "fine_G_S": finite(fine.get("fine_G_S")) if fine_valid else None,
        "fine_Sf_m": finite(fmarginal.get("Sf_m", fine.get("fine_Sf_m"))) if fine_valid else None,
        "fine_triangles": finite(fine.get("fine_triangles")) if fine_valid else None,
        "fine_file": (fine_folder / "SUCCESS.json").as_posix() if fine_valid else None,
    })
    return row, prediction


def get_quantity(index, layout, sf, name, fine=False):
    field = ("fine_" if fine else "") + ("Sf_m" if name in ("SfA", "SfB", "M") else "deltaG_S")
    a, b = index[(layout, "A", sf)].get(field), index[(layout, "B", sf)].get(field)
    if name in ("gA", "SfA"):
        return a
    if name in ("gB", "SfB"):
        return b
    return difference(b, a)


def observed_corner_change(index, layout, sf, name):
    field = "Sf_m" if name in ("SfA", "SfB", "M") else "deltaG_S"
    roles = ("A",) if name in ("gA", "SfA") else ("B",) if name in ("gB", "SfB") else ("A", "B")
    changes = [difference(index[(layout, role, sf)].get("fine_" + field), index[(layout, role, sf)].get(field)) for role in roles]
    return sum(abs(x) for x in changes) if all(x is not None for x in changes) else None


def contrast(values, weights):
    return sum(weights[k] * values[k] for k in weights) if all(values[k] is not None for k in weights) else None


def rank_sign(value):
    if value is None:
        return "unavailable"
    return "B_greater" if value > 0 else "A_greater" if value < 0 else "equal_point_values"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=DEFAULT, help="Extension output directory")
    args = parser.parse_args()
    out = args.root.resolve()
    notes = []
    manifest = read(out / "manifest.json", notes)
    if not manifest:
        raise SystemExit("manifest.json is required")
    sigmas = [float(v) for v in manifest["sigma_f_S_m"]]
    lookup = {}
    for old in manifest["old_references"]:
        lookup[(old["factorial_layout"], old["configuration_role"])] = (Path(old["source_folder"]), True, old["case_id"])
    for c in manifest["cases"]:
        lookup[(c["layout"], c["configuration_role"])] = (out / "cases" / c["id"], False, c["id"])
    responses, predictions = [], []
    for sf in sigmas:
        for layout in LAYOUTS:
            for role in ("A", "B"):
                folder, reused, ident = lookup[(layout, role)]
                row, prediction = load_response(layout, role, sf, folder, reused, ident, out, notes)
                responses.append(row)
                if not reused:
                    estimated = finite(prediction.get("predicted_deltaG_S"))
                    actual = row["deltaG_S"]
                    err = difference(estimated, actual)
                    lower = finite(prediction.get("certified_lower_deltaG_S"))
                    predictions.append({
                        "layout": layout, "role": role, "case_id": ident, "sigma_f_S_m": sf,
                        "status": "complete" if actual is not None and estimated is not None else "prediction_saved_truth_pending" if estimated is not None else "pending",
                        "G0_S": row["G0_S"], "actual_deltaG_S": actual,
                        "predicted_deltaG_S": estimated, "signed_error_S": err,
                        "absolute_error_S": abs(err) if err is not None else None,
                        "relative_increment_error": ratio(abs(err), abs(actual)) if err is not None and actual is not None else None,
                        "certified_lower_deltaG_S": lower, "certified_upper_deltaG_S": estimated,
                        "actual_inside_discrete_bounds": lower <= actual <= estimated if all(v is not None for v in (lower, actual, estimated)) else None,
                        "prediction_saved_local": prediction.get("prediction_saved_local"),
                        "phase": prediction.get("phase"), "method_version": prediction.get("method_version"),
                        "independent_blind_geometry_test": False,
                        "prediction_file": (folder / token(sf) / "prediction.json").as_posix(),
                        "truth_file": row["truth_file"],
                        "certificate_scope": "Current discrete FE system only; does not include mesh error or natural-rock uncertainty.",
                    })
    index = {(r["layout"], r["role"], r["sigma_f_S_m"]): r for r in responses}
    for row in responses:
        row["D_B_minus_A_S"] = get_quantity(index, row["layout"], row["sigma_f_S_m"], "D")
        row["M_B_minus_A_m"] = get_quantity(index, row["layout"], row["sigma_f_S_m"], "M")
        row["D_point_rank"] = rank_sign(row["D_B_minus_A_S"])
    pindex = {(r["layout"], r["role"], r["sigma_f_S_m"]): r for r in predictions}
    for row in predictions:
        a, b = [pindex[(row["layout"], role, row["sigma_f_S_m"])] for role in ("A", "B")]
        dp = difference(b["predicted_deltaG_S"], a["predicted_deltaG_S"])
        dt = difference(b["actual_deltaG_S"], a["actual_deltaG_S"])
        low = difference(b["certified_lower_deltaG_S"], a["certified_upper_deltaG_S"])
        high = difference(b["certified_upper_deltaG_S"], a["certified_lower_deltaG_S"])
        row.update({
            "pair_predicted_D_S": dp, "pair_actual_D_S": dt,
            "pair_point_ranking_correct": rank_sign(dp) == rank_sign(dt) if dp is not None and dt is not None else None,
            "pair_certified_D_lower_S": low, "pair_certified_D_upper_S": high,
            "pair_certificate_rank": "B_greater" if low is not None and low > 0 else "A_greater" if high is not None and high < 0 else "overlap_zero_cannot_rank" if low is not None and high is not None else "unavailable",
            "pair_fields_repeated_for_roles": True,
        })
    effects, refined = [], []
    for sf in sigmas:
        for name in QUANTITIES:
            unit = "m" if name in ("SfA", "SfB", "M") else "S"
            vals = {layout: get_quantity(index, layout, sf, name) for layout in LAYOUTS}
            fvals = {layout: get_quantity(index, layout, sf, name, fine=True) for layout in LAYOUTS}
            scales = {layout: observed_corner_change(index, layout, sf, name) for layout in LAYOUTS}
            for effect, weights in EFFECTS.items():
                value = contrast(vals, weights)
                row = {
                    "sigma_f_S_m": sf, "quantity": name, "unit": unit, "effect": effect,
                    "status": "complete" if value is not None else "partial",
                    **{f"{layout}_value": vals[layout] for layout in LAYOUTS},
                    **{f"weight_{layout}": weights.get(layout, 0) for layout in LAYOUTS},
                    "effect_value": value,
                    "interpretation": "Controlled finite position contrast; averages are a symmetric algebraic allocation, not unique independent physical shares.",
                }
                effects.append(row)
                if not math.isclose(sf, .159375, rel_tol=1e-13, abs_tol=1e-14):
                    continue
                fine_value = contrast(fvals, weights)
                weighted_scale = sum(abs(weight) * scales[layout] for layout, weight in weights.items()) if all(scales[k] is not None for k in weights) else None
                all_corner_scale = sum(scales.values()) if all(v is not None for v in scales.values()) else None
                available = all(v is not None for v in (value, fine_value, weighted_scale))
                stable = value * fine_value > 0 if available else None
                resolved = min(abs(value), abs(fine_value)) > 3 * weighted_scale if available else None
                status = "not_available" if not available else "sign_changed" if not stable and value != 0 and fine_value != 0 else "same_sign_exceeds_observed_scale" if stable and resolved else "near_observed_mesh_change_scale"
                refined.append({
                    "sigma_f_S_m": sf, "quantity": name, "unit": unit, "effect": effect,
                    "status": status, "coarse_value": value, "fine_value": fine_value,
                    "fine_minus_coarse": difference(fine_value, value),
                    "observed_weighted_sum_underlying_mesh_changes": weighted_scale,
                    "observed_all_four_corner_underlying_mesh_changes": all_corner_scale,
                    **{f"{layout}_observed_underlying_mesh_change": scales[layout] for layout in LAYOUTS},
                    **{f"{layout}_fine_value": fvals[layout] for layout in LAYOUTS},
                    "diagnostic_factor": 3, "coarse_fine_same_sign": stable,
                    "both_magnitudes_exceed_three_times_observed_scale": resolved,
                    "mesh_change_interpretation": MESH_NOTE,
                    "missing_fine_Sf_note": "Fine Sf is read only from a saved target integral; never replaced by a coarse value or inferred from fine G." if name in ("SfA", "SfB", "M") else None,
                })
    write_csv(out / "responses.csv", responses)
    write_csv(out / "factorial_effects.csv", effects)
    write_csv(out / "prediction_comparison.csv", predictions)
    write_csv(out / "fine_factorial_comparison.csv", refined)
    summary = {
        "generated_local": datetime.now().isoformat(timespec="seconds"),
        "expected_response_rows": 24, "response_rows": len(responses),
        "completed_states": sum(r["status"] == "success" for r in responses),
        "completed_new_states": sum(r["status"] == "success" and not r["reused_round4_state"] for r in responses),
        "reused_states": sum(r["status"] == "success" and r["reused_round4_state"] for r in responses),
        "complete_new_prediction_comparisons": sum(r["status"] == "complete" for r in predictions),
        "fine_g_states": sum(r["fine_deltaG_S"] is not None for r in responses),
        "fine_Sf_states": sum(r["fine_Sf_m"] is not None for r in responses),
        "factorial_rows": len(effects), "fine_factorial_rows": len(refined),
        "notes": notes,
        "scope": "Read current extension states plus exact 12 old NN/FF references; no solver, no old writes, no predictive-method adjustment.",
    }
    (out / "summary_status.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
