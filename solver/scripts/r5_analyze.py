# -*- coding: utf-8 -*-
"""Round-5 development order selection and prospective validation tables.

Reads saved predictions and appropriate COMSOL responses; never solves a
model or changes an algorithm. Development selection uses only the 24 old
states. Prospective target reads require an explicit --phase request and the
saved all-batch prediction gate. No method is frozen by this script.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/round5_modal_closure"
ORDERS = [0, 1, 2, 4, 8, 16, 32, 64]
DEV_EFFECTS = {
    "lower_near_at_upper_near": {"NN": 1, "FN": -1},
    "lower_near_at_upper_far": {"NF": 1, "FF": -1},
    "upper_near_at_lower_near": {"NN": 1, "NF": -1},
    "upper_near_at_lower_far": {"FN": 1, "FF": -1},
    "lower_average": {"NN": .5, "NF": .5, "FN": -.5, "FF": -.5},
    "upper_average": {"NN": .5, "NF": -.5, "FN": .5, "FF": -.5},
    "interaction": {"NN": 1, "NF": -1, "FN": -1, "FF": 1},
    "joint_NN_minus_FF": {"NN": 1, "FF": -1},
}
GATE_EFFECTS = set(list(DEV_EFFECTS)[:4])
NEW_EFFECTS = {
    "lower_7_to_1_at_upper_1": {"P11": 1, "P71": -1},
    "lower_7_to_1_at_upper_7": {"P17": 1, "P77": -1},
    "upper_7_to_1_at_lower_1": {"P11": 1, "P17": -1},
    "upper_7_to_1_at_lower_7": {"P71": 1, "P77": -1},
    "lower_4_to_1_at_upper_1": {"P11": 1, "P41": -1},
    "lower_7_to_4_at_upper_1": {"P41": 1, "P71": -1},
    "upper_4_to_1_at_lower_1": {"P11": 1, "P14": -1},
    "upper_7_to_4_at_lower_1": {"P14": 1, "P17": -1},
    "joint_P11_minus_P77": {"P11": 1, "P77": -1},
    "interaction_1_7": {"P11": 1, "P17": -1, "P71": -1, "P77": 1},
}


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def tag(value):
    return "sf_" + format(float(value), ".10g").replace(".", "p")


def sf_key(value):
    return format(float(value), ".10g")


def sign(value):
    return 1 if value > 0 else -1 if value < 0 else 0


def relative(error, reference):
    return None if reference == 0 else abs(error / reference)


def methods(prediction):
    value = prediction["methods"]
    if isinstance(value, dict):
        value = [value] if "K" in value else list(value.values())
    result = {}
    for method in value:
        k = int(method["K"])
        if k in result:
            raise ValueError(f"Duplicate prediction order K={k}")
        for field in ["G0_S", "G_S", "deltaG_S"]:
            if not math.isfinite(float(method[field])):
                raise ValueError("Nonfinite saved prediction: " + field)
        result[k] = method
    return result


def load_states(design, phase):
    saved, missing = {}, []
    for state in design["development_states" if phase == "development" else "new_states"]:
        role, layout, sf = state["role"], state["layout"], float(state["sigma_f_S_m"])
        if phase == "development":
            prediction_file = OUT / "development" / f"{role}_{layout}" / tag(sf) / "prediction.json"
            truth_file = Path(state["existing_target_directory"]) / "SUCCESS.json"
        else:
            prediction_file = OUT / "predictions" / state["case_id"] / tag(sf) / "prediction.json"
            truth_file = OUT / "truth" / state["case_id"] / tag(sf) / "SUCCESS.json"
        if not prediction_file.is_file() or not truth_file.is_file():
            missing.append({"state_id": state["state_id"], "prediction_exists": prediction_file.is_file(),
                            "truth_exists": truth_file.is_file()})
            continue
        p, truth = read_json(prediction_file), read_json(truth_file)
        if str(p["role"]) != role or str(p["layout"]) != layout:
            raise ValueError("Prediction role/layout mismatch: " + str(prediction_file))
        if not math.isclose(float(p["sigma_f_S_m"]), sf, rel_tol=1e-12, abs_tol=1e-14):
            raise ValueError("Prediction conductivity mismatch: " + str(prediction_file))
        if not math.isclose(float(truth["sigma_f_S_m"]), sf, rel_tol=1e-12, abs_tol=1e-14):
            raise ValueError("Target conductivity mismatch: " + str(truth_file))
        key = (role, layout, sf_key(sf))
        if key in saved:
            raise ValueError("Duplicate state key: " + repr(key))
        saved[key] = {"state": state, "prediction": p, "methods": methods(p), "truth": truth,
                      "prediction_file": prediction_file.as_posix(), "truth_file": truth_file.as_posix()}
    return saved, missing


def comparison_rows(saved, phase):
    rows = []
    for (role, layout, sf), entry in saved.items():
        truth = entry["truth"]
        for k, method in sorted(entry["methods"].items()):
            error = float(method["deltaG_S"]) - float(truth["deltaG_S"])
            rel = relative(error, float(truth["deltaG_S"]))
            rows.append({
                "phase": phase, "case_id": entry["state"]["case_id"], "role": role, "layout": layout,
                "sigma_f_S_m": float(sf), "method": method.get("method", "M0" if k == 0 else "MK"), "K": k,
                "G0_prediction_S": method["G0_S"], "G0_COMSOL_S": truth["G0_S"],
                "G_prediction_S": method["G_S"], "G_COMSOL_S": truth["G_S"],
                "g_prediction_S": method["deltaG_S"], "g_COMSOL_S": truth["deltaG_S"],
                "g_signed_error_S": error, "g_absolute_error_S": abs(error),
                "g_relative_error": rel, "g_relative_error_percent": None if rel is None else 100*rel,
                "truth_source": "existing_round4_development" if phase == "development" else "new_full_COMSOL_post_prediction",
                "prediction_file": entry["prediction_file"], "truth_file": entry["truth_file"],
            })
    return rows


def fine_evidence():
    path = ROOT / "runs/round4_extension_20261024/fine_factorial_comparison.csv"
    rows = read_csv(path)
    return {(r["quantity"], r["effect"], sf_key(r["sigma_f_S_m"])): r for r in rows}


def effect_value(saved, quantity, sf, weights, k=None):
    roles = [quantity[1:]] if quantity in ["gA", "gB"] else ["A", "B"]
    total = 0.0
    for role in roles:
        role_weight = -1 if quantity == "D" and role == "A" else 1
        for layout, weight in weights.items():
            entry = saved.get((role, layout, sf))
            if entry is None or (k is not None and k not in entry["methods"]):
                return None
            value = entry["truth"]["deltaG_S"] if k is None else entry["methods"][k]["deltaG_S"]
            total += role_weight * weight * float(value)
    return total


def effect_rows(saved, phase):
    effects = DEV_EFFECTS if phase == "development" else NEW_EFFECTS
    evidence = fine_evidence() if phase == "development" else {}
    orders = sorted({k for entry in saved.values() for k in entry["methods"]})
    points = sorted({sf for _, _, sf in saved}, key=float)
    rows = []
    for sf in points:
        for quantity in ["gA", "gB", "D"]:
            for name, weights in effects.items():
                truth = effect_value(saved, quantity, sf, weights)
                if truth is None:
                    continue
                fine = evidence.get((quantity, name, sf))
                resolved = bool(fine and fine["status"] == "same_sign_exceeds_observed_scale")
                gate = (phase == "development" and quantity in ["gA", "gB"] and name in GATE_EFFECTS
                        and sf == sf_key(.159375) and resolved)
                fvalue = float(fine["fine_value"]) if fine else None
                for k in orders:
                    prediction = effect_value(saved, quantity, sf, weights, k)
                    if prediction is None:
                        continue
                    error = prediction - truth
                    ferror = None if fvalue is None else prediction - fvalue
                    # Tiny cancellation/interactions retain absolute errors without
                    # being silently promoted to strong relative-error evidence.
                    unresolved_cancellation = name.startswith("joint_") or name.startswith("interaction")
                    relative_meaningful = resolved or not unresolved_cancellation
                    rel = relative(error, truth) if relative_meaningful else None
                    frel = relative(ferror, fvalue) if resolved else None
                    rows.append({
                        "phase": phase, "sigma_f_S_m": float(sf), "quantity": quantity, "effect": name, "K": k,
                        "weights": json.dumps(weights, sort_keys=True), "prediction_S": prediction,
                        "COMSOL_S": truth, "signed_error_S": error, "absolute_error_S": abs(error),
                        "relative_error": rel, "relative_error_percent": None if rel is None else 100*rel,
                        "point_sign_correct": sign(prediction) == sign(truth),
                        "fine_COMSOL_S": fvalue, "fine_signed_error_S": ferror,
                        "fine_relative_error": frel, "fine_relative_error_percent": None if frel is None else 100*frel,
                        "fine_point_sign_correct": None if fvalue is None else sign(prediction) == sign(fvalue),
                        "fine_mesh_evidence_status": fine["status"] if fine else "no_selected_refinement",
                        "fine_observed_mesh_scale_S": float(fine["observed_weighted_sum_underlying_mesh_changes"]) if fine else None,
                        "required_for_development_K_selection": gate,
                        "relative_error_interpretation": "resolved_by_saved_mesh_trend" if resolved else "descriptive_no_mesh_resolution_claim" if relative_meaningful else "omitted_for_unresolved_cancellation_use_absolute_error",
                    })
    return rows


def ranking_rows(saved, phase):
    rows = []
    layouts = sorted({layout for _, layout, _ in saved})
    points = sorted({sf for _, _, sf in saved}, key=float)
    orders = sorted({k for entry in saved.values() for k in entry["methods"]})
    for layout in layouts:
        for sf in points:
            a, b = saved.get(("A", layout, sf)), saved.get(("B", layout, sf))
            if a is None or b is None:
                continue
            td = b["truth"]["deltaG_S"] - a["truth"]["deltaG_S"]
            for k in orders:
                if k not in a["methods"] or k not in b["methods"]:
                    continue
                pd = b["methods"][k]["deltaG_S"] - a["methods"][k]["deltaG_S"]
                rows.append({"phase": phase, "layout": layout, "sigma_f_S_m": float(sf), "K": k,
                             "D_prediction_S": pd, "D_COMSOL_S": td, "D_signed_error_S": pd-td,
                             "D_absolute_error_S": abs(pd-td), "point_sign_correct": sign(pd) == sign(td),
                             "interpretation": "Point-value rank only; near-tie significance requires separate mesh evidence."})
    return rows


def choose_order(saved, comparisons, sides, missing):
    candidates = []
    for k in ORDERS:
        gains = [r for r in comparisons if r["K"] == k]
        relevant = [r for r in sides if r["K"] == k and r["required_for_development_K_selection"]]
        complete = len(gains) == 24 and len(relevant) == 8
        max_gain = max((r["g_relative_error"] for r in gains), default=None)
        max_side = max((r["fine_relative_error"] for r in relevant), default=None)
        signs = all(r["fine_point_sign_correct"] for r in relevant) if relevant else False
        gains_ok = complete and max_gain <= .01
        sides_ok = complete and max_side <= .10 and signs
        candidates.append({
            "K": k, "available_gain_states": len(gains), "required_side_effects_available": len(relevant),
            "complete": complete, "gain_max_relative_error": max_gain,
            "gain_median_relative_error": statistics.median(r["g_relative_error"] for r in gains) if gains else None,
            "required_fine_side_max_relative_error": max_side, "required_fine_side_all_signs_correct": signs,
            "gain_requirement_passed": gains_ok, "side_requirement_passed": sides_ok,
            "both_requirements_passed": gains_ok and sides_ok,
        })
    usable = [r["K"] for r in candidates if 1 <= r["K"] <= 32 and r["both_requirements_passed"]]
    gain_only = [r["K"] for r in candidates if r["gain_requirement_passed"]]
    selected = min(usable) if usable else None
    stability = {"status": "incomplete", "Mref_suggested": None}
    if all(32 in e["methods"] and 64 in e["methods"] for e in saved.values()) and len(saved) == 24:
        gain_changes = [abs((e["methods"][64]["deltaG_S"] - e["methods"][32]["deltaG_S"])/e["methods"][64]["deltaG_S"]) for e in saved.values()]
        side_changes = []
        for r in sides:
            if r["K"] != 32 or not r["required_for_development_K_selection"]:
                continue
            matching = [s for s in sides if s["K"] == 64 and s["quantity"] == r["quantity"] and s["effect"] == r["effect"] and s["sigma_f_S_m"] == r["sigma_f_S_m"]]
            assert len(matching) == 1
            side_changes.append(abs((matching[0]["prediction_S"]-r["prediction_S"])/r["fine_COMSOL_S"]))
        assert len(side_changes) == 8
        r32 = next(r for r in candidates if r["K"] == 32)
        stable = max(gain_changes) <= .0001 and max(side_changes) <= .01
        stability = {"status": "complete", "g_32_to_64_max_relative_change": max(gain_changes),
                     "resolved_side_32_to_64_max_change_over_fine_effect": max(side_changes),
                     "gain_change_threshold": .0001, "side_change_threshold": .01,
                     "changes_pass_thresholds": stable, "order32_passes_truth_accuracy_requirements": r32["both_requirements_passed"],
                     "Mref_suggested": 32 if stable and r32["both_requirements_passed"] else 64}
    return {
        "status": "development_complete_pending_root_review" if not missing and all(c["complete"] for c in candidates) else "development_incomplete",
        "scope": "Only 24 previously solved development states; no prospective target used for order selection.",
        "method_frozen_written": False, "candidate_order_results": candidates,
        "gain_requirement": {"max_relative_error": .01, "state_count": 24},
        "side_requirement": {"max_relative_error": .10, "require_sign": True, "count": 8,
                             "sigma_f_S_m": .159375, "reference": "Saved fine_factorial_comparison.csv fine values; only gA/gB four single-side effects verified by prior refinement."},
        "minimum_K_for_gain_only": min(gain_only) if gain_only else None,
        "recommended_K_star": selected, "low_order_meets_both": selected is not None,
        "fallback": None if selected is not None else "No K in 1..32 satisfies both; inspect 64 and discretization before claiming few modes suffice.",
        "Mref_32_64_check": stability, "missing_states": missing,
        "precision_limits": "Gain accuracy and tiny side-difference accuracy are distinct. Prior refinement diagnostics are observed mesh trends, not rigorous error bounds. Joint cancellation and interactions are not forced through relative-error gates.",
    }


def merge_delivery_tables():
    comparisons, sides, ranks = [], [], []
    for phase in ["development", "prospective"]:
        for suffix, destination in [("comparison", comparisons), ("side_effects", sides), ("ranking", ranks)]:
            source = OUT / f"{phase}_{suffix}.csv"
            if source.is_file():
                destination.extend(read_csv(source))
    crossover = OUT / "crossover_gain_comparison.csv"
    if crossover.is_file():
        comparisons.extend(read_csv(crossover))
    write_csv(OUT / "closure_comparison.csv", comparisons)
    write_csv(OUT / "side_effects.csv", sides)
    write_csv(OUT / "ranking_comparison.csv", ranks)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["development", "prospective", "all"], default="development")
    args = parser.parse_args()
    design = read_json(OUT / "design.json")
    summaries = {}
    for phase in (["development", "prospective"] if args.phase == "all" else [args.phase]):
        if phase == "prospective":
            if not (OUT / "predictions_ready.json").is_file() or not (OUT / "method_frozen.json").is_file():
                raise RuntimeError("Do not read new targets before the all-predictions gate and method freeze are saved")
        saved, missing = load_states(design, phase)
        comparisons = comparison_rows(saved, phase)
        sides = effect_rows(saved, phase)
        ranks = ranking_rows(saved, phase)
        write_csv(OUT / f"{phase}_comparison.csv", comparisons)
        write_csv(OUT / f"{phase}_side_effects.csv", sides)
        write_csv(OUT / f"{phase}_ranking.csv", ranks)
        if phase == "development":
            summary = choose_order(saved, comparisons, sides, missing)
            write_json(OUT / "development_order_selection.json", summary)
        else:
            summaries_k = []
            for k in sorted({r["K"] for r in comparisons}):
                rows = [r for r in comparisons if r["K"] == k]
                rankrows = [r for r in ranks if r["K"] == k]
                summaries_k.append({"K": k, "states": len(rows),
                                    "gain_median_relative_error": statistics.median(r["g_relative_error"] for r in rows),
                                    "gain_max_relative_error": max(r["g_relative_error"] for r in rows),
                                    "gain_max_absolute_error_S": max(r["g_absolute_error_S"] for r in rows),
                                    "point_rank_correct": sum(r["point_sign_correct"] for r in rankrows),
                                    "point_rank_count": len(rankrows)})
            summary = {"status": "complete" if not missing else "incomplete", "states": len(saved), "missing": missing,
                       "by_order": summaries_k, "selection_used_new_truth": False,
                       "scope": "Prospective validation within the fixed strict-core flat-periphery family, not arbitrary-geometry generalization."}
            write_json(OUT / "prospective_summary.json", summary)
        summaries[phase] = {"saved_states": len(saved), "missing_states": len(missing), "comparison_rows": len(comparisons),
                            "side_rows": len(sides), "ranking_rows": len(ranks),
                            "recommended_K_star": summary.get("recommended_K_star")}
    merge_delivery_tables()
    print(json.dumps(summaries, ensure_ascii=False))


if __name__ == "__main__":
    main()
