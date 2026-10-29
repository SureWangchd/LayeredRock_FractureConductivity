"""Fifth-round modal diagnostics from saved PREDICTIONS only.

Reads prediction.json and analytic port_coefficients.csv. No target truth,
target fields, fitted coefficients, MATLAB, or COMSOL calls are involved.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "runs" / "round5_modal_closure"
SIDES = ("lower", "upper")
PAIR_DESIGNS = {
    "development": [
        ("FF", "NF", "lower_near_upper_far"),
        ("FF", "FN", "upper_near_lower_far"),
        ("FF", "NN", "both_near"),
        ("FN", "NN", "lower_near_upper_near"),
        ("NF", "NN", "upper_near_lower_near"),
    ],
    "new": [
        ("P77", "P17", "lower_7_to_1_upper_7"),
        ("P77", "P71", "upper_7_to_1_lower_7"),
        ("P77", "P11", "both_7_to_1"),
        ("P71", "P11", "lower_7_to_1_upper_1"),
        ("P17", "P11", "upper_7_to_1_lower_1"),
        ("P71", "P41", "lower_7_to_4_upper_1"),
        ("P41", "P11", "lower_4_to_1_upper_1"),
        ("P17", "P14", "upper_7_to_4_lower_1"),
        ("P14", "P11", "upper_4_to_1_lower_1"),
    ],
}
PAIR_DESIGNS["precision"] = PAIR_DESIGNS["development"] + [
    ("P77", "P11", "both_7_to_1")
]
for phase in ("limited_crossover_prediction_only", "limited_crossover_prediction_reuse"):
    PAIR_DESIGNS[phase] = [
        ("P77", "P17", "lower_7_to_1_upper_7"),
        ("P77", "P11", "both_7_to_1"),
        ("P17", "P11", "upper_7_to_1_lower_1"),
    ]


def load_coefficients(path: Path) -> dict:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    coefficients = {}
    for row in rows:
        key = (row["layout"], row["side"], int(row["n"]))
        if key in coefficients:
            raise ValueError(f"Duplicate analytic coefficient: {key}")
        numeric = {
            key: float(row[key])
            for key in (
                "k_per_m", "Y_actual_S_m2", "Y_reference_S_m2",
                "norm_N_m", "width_m", "out_of_plane_thickness_m",
            )
        }
        numeric["layer_source_case_id"] = row.get("layer_source_case_id", "")
        coefficients[(row["layout"], row["side"], int(row["n"]))] = numeric
    return coefficients


def model_y(coeff: dict, n: int, order: int) -> float:
    # The engine changes only nonzero modes 1..K. Its zero mode and tail
    # remain the shared M0 FE realization, including the external drive.
    return coeff["Y_actual_S_m2"] if 1 <= n <= order else coeff["Y_reference_S_m2"]


def load_predictions(out: Path) -> list[dict]:
    states = []
    for folder in ("development", "predictions", "prediction_precision", "crossover_predictions"):
        for path in sorted((out / folder).glob("*/*/prediction.json")):
            value = json.loads(path.read_text(encoding="utf-8-sig"))
            if value.get("target_truth_used") is not False:
                raise ValueError(f"Prediction provenance missing/unsafe: {path}")
            methods = value["methods"]
            if isinstance(methods, dict):
                methods = [methods]
            for method in methods:
                if method.get("target_truth_used", False):
                    raise ValueError(f"Method uses target truth: {path}")
                if len(method["port_v"]) != 2 or any(len(v) != 65 for v in method["port_v"]):
                    raise ValueError(f"Expected 2 by 65 predicted modes: {path}")
                if not all(math.isfinite(v) for side in method["port_v"] for v in side):
                    raise ValueError(f"Nonfinite predicted mode: {path}")
                states.append({
                    "phase": value["phase"], "mesh": value["mesh"],
                    "role": value["role"], "layout": value["layout"],
                    "sigma_f_S_m": float(value["sigma_f_S_m"]),
                    "K": int(method["K"]), "method": method,
                    "method_version": value["method_version"],
                    "prediction_saved_local": value["prediction_saved_local"],
                    "prediction_source": str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
                })
    return states


def metadata(state: dict) -> dict:
    return {key: state[key] for key in (
        "phase", "mesh", "role", "layout", "sigma_f_S_m", "K",
        "method_version", "prediction_saved_local", "prediction_source",
    )}


def port_rows(states: list[dict], coefficients: dict) -> list[dict]:
    rows = []
    for state in states:
        for si, side in enumerate(SIDES):
            for n, v in enumerate(state["method"]["port_v"][si]):
                coeff = coefficients[state["layout"], side, n]
                k = coeff["k_per_m"]
                y = model_y(coeff, n, state["K"])
                row = metadata(state)
                row.update({
                    "side": side, "n": n, "k_per_m": k,
                    "transverse_scale_1_over_k_mm": 1000 / k if n else "",
                    "v_predicted_V": v, "Y_actual_S_m2": coeff["Y_actual_S_m2"],
                    "Y_reference_S_m2": coeff["Y_reference_S_m2"],
                    "Y_model_S_m2": y,
                    "actual_position_information_retained": 1 <= n <= state["K"],
                    "norm_N_m": coeff["norm_N_m"],
                    "out_of_plane_thickness_m": coeff["out_of_plane_thickness_m"],
                    "layer_source_case_id": coeff["layer_source_case_id"],
                    "field_provenance": "selfconsistent_prediction_no_target_field",
                    "mode_units": "V", "admittance_units": "S/m^2",
                })
                rows.append(row)
    return rows


def pair_rows(states: list[dict], coefficients: dict) -> tuple[list[dict], list[dict]]:
    groups = {}
    for state in states:
        key = (state["phase"], state["mesh"], state["role"],
               state["sigma_f_S_m"], state["K"], state["method_version"])
        layouts = groups.setdefault(key, {})
        if state["layout"] in layouts:
            raise ValueError(f"Duplicate prediction state {key} / {state['layout']}")
        layouts[state["layout"]] = state
    terms, summaries = [], []
    for key, layouts in sorted(groups.items()):
        phase, mesh, role, sf, order, version = key
        for before, after, comparison in PAIR_DESIGNS.get(phase, []):
            if before not in layouts or after not in layouts:
                continue  # A partially completed batch remains rerunnable.
            p, q = layouts[before], layouts[after]
            mp, mq = p["method"], q["method"]
            dv = float(mp.get("DV_V", 1.0))
            if not math.isclose(dv, float(mq.get("DV_V", 1.0)), rel_tol=1e-14):
                raise ValueError("Compared models have different external voltages")
            common = {
                "phase": phase, "mesh": mesh, "role": role, "sigma_f_S_m": sf,
                "K": order, "layout_p": before, "layout_q": after,
                "comparison_q_minus_p": comparison, "method_version": version,
                "prediction_p_saved_local": p["prediction_saved_local"],
                "prediction_q_saved_local": q["prediction_saved_local"],
                "prediction_p_source": p["prediction_source"],
                "prediction_q_source": q["prediction_source"],
            }
            modal = {side: [] for side in SIDES}
            detail = []
            for si, side in enumerate(SIDES):
                for n in range(order + 1):
                    cp, cq = coefficients[before, side, n], coefficients[after, side, n]
                    for check in ("norm_N_m", "out_of_plane_thickness_m", "Y_reference_S_m2"):
                        if not math.isclose(cp[check], cq[check], rel_tol=1e-12, abs_tol=1e-15):
                            raise ValueError(f"Compared models have different shared {check}")
                    vp, vq = mp["port_v"][si][n], mq["port_v"][si][n]
                    delta_y = model_y(cq, n, order) - model_y(cp, n, order)
                    term = cp["out_of_plane_thickness_m"] * cp["norm_N_m"] * delta_y * vp * vq / dv**2
                    modal[side].append(term)
                    k = cp["k_per_m"]
                    detail.append({**common,
                        "side": side, "n": n, "k_per_m": k,
                        "transverse_scale_1_over_k_mm": 1000 / k if n else "",
                        "Y_p_model_S_m2": model_y(cp, n, order),
                        "Y_q_model_S_m2": model_y(cq, n, order),
                        "delta_Y_model_S_m2": delta_y,
                        "v_p_predicted_V": vp, "v_q_predicted_V": vq,
                        "v_p_times_v_q_V2": vp * vq,
                        "bilinear_effect_S": term,
                        "cumulative_side_through_n_S": math.fsum(modal[side]),
                        "cumulative_both_sides_through_n_S": 0.0,
                        "interpretation": "two_predicted_fields_finite_change_identity_not_unique_causal_share",
                    })
            for row in detail:
                through = row["n"] + 1
                row["cumulative_both_sides_through_n_S"] = math.fsum(
                    math.fsum(modal[side][:through]) for side in SIDES)
            terms.extend(detail)
            lower, upper = (math.fsum(modal[side]) for side in SIDES)
            total = math.fsum((lower, upper))
            delta_g = float(mq["deltaG_S"]) - float(mp["deltaG_S"])
            delta_G = float(mq["G_S"]) - float(mp["G_S"])
            delta_G0 = float(mq["G0_S"]) - float(mp["G0_S"])
            residual_G = total - delta_G
            residual_g = total - delta_G0 - delta_g
            scale = math.fsum(abs(v) for side in SIDES for v in modal[side])
            row = {**common,
                "predicted_G_q_minus_p_S": delta_G,
                "predicted_g_q_minus_p_S": delta_g,
                "predicted_G0_q_minus_p_S": delta_G0,
                "lower_mode_effect_S": lower, "upper_mode_effect_S": upper,
                "mode_sum_effect_S": total,
                "sum_abs_individual_mode_effects_S": scale,
                "identity_residual_vs_G_difference_S": residual_G,
                "identity_residual_vs_g_difference_after_background_S": residual_g,
                "relative_identity_residual_on_G_difference": abs(residual_G / delta_G) if abs(delta_G) > 1e-14 else "",
                "relative_identity_residual_on_abs_mode_sum": abs(residual_G / scale) if scale > 1e-14 else "",
                "zero_mode_effect_S": modal["lower"][0] + modal["upper"][0],
                "field_provenance": "two_selfconsistent_predictions_no_target_fields",
                "tail_treatment": "common_reference_above_K_cancels_in_same_K_layout_difference",
            }
            for cutoff in (1, 2, 4, 8, 16, 32, 64):
                row[f"cumulative_through_{cutoff}_S"] = math.fsum(
                    math.fsum(modal[side][: min(order, cutoff) + 1]) for side in SIDES)
            summaries.append(row)
    return terms, summaries


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def mode_scale_rows(terms: list[dict]) -> list[dict]:
    """Selected prediction-only 1..4/5..8/9..16/17..32/33..64 bands.

    t*N*deltaY has units S, and the voltage-product mean has units V^2.
    Their product / DV^2 supplies the band effect if all coefficients share
    one sign; the signed individual terms always supply the exact effect.
    """
    groups = {}
    for term in terms:
        if term["mesh"] != "base":
            continue
        selected = (
            term["phase"] == "development"
            and math.isclose(term["sigma_f_S_m"], .159375, rel_tol=1e-12)
            and term["K"] in (32, 64)
            and term["layout_p"] == "FF" and term["layout_q"] == "NN"
        ) or (
            term["phase"] == "new"
            and any(math.isclose(term["sigma_f_S_m"], x, rel_tol=1e-12) for x in (.18, 10.0))
            and term["K"] == 32 and term["layout_p"] == "P77"
            and term["layout_q"] in ("P11", "P17", "P71")
        )
        if not selected:
            continue
        key = tuple(term[k] for k in (
            "phase", "mesh", "role", "sigma_f_S_m", "K", "layout_p", "layout_q", "side"
        ))
        groups.setdefault(key, []).append(term)
    rows = []
    bands = ((1, 4), (5, 8), (9, 16), (17, 32), (33, 64))
    for key, available in sorted(groups.items()):
        phase, mesh, role, sf, order, before, after, side = key
        for first, last in bands:
            block = [term for term in available if first <= term["n"] <= last]
            # All requested physical modes use the common W/2=.06 m,
            # t=1 m and DV=1 V of the saved prediction engine.
            weights = [.06 * term["delta_Y_model_S_m2"] for term in block]
            absolute_weight = math.fsum(abs(v) for v in weights)
            product_mean = (
                math.fsum(abs(w) * term["v_p_times_v_q_V2"] for w, term in zip(weights, block)) / absolute_weight
                if absolute_weight else ""
            )
            signs = {1 if w > 0 else -1 for w in weights if w}
            rows.append({
                "phase": phase, "mesh": mesh, "role": role,
                "sigma_f_S_m": sf, "K": order, "layout_p": before,
                "layout_q": after, "side": side, "n_first": first, "n_last": last,
                "scale_1_over_k_max_mm": 120 / (math.pi * first),
                "scale_1_over_k_min_mm": 120 / (math.pi * last),
                "retained_modes_in_band": len(block),
                "band_has_actual_position_information": first <= order,
                "signed_response_weight_sum_S": math.fsum(weights),
                "absolute_response_weight_sum_S": absolute_weight,
                "response_weighted_predicted_vp_vq_V2": product_mean,
                "response_coefficient_sign_consistent": len(signs) <= 1,
                "bilinear_band_effect_S": math.fsum(term["bilinear_effect_S"] for term in block),
                "sum_abs_bilinear_band_terms_S": math.fsum(abs(term["bilinear_effect_S"]) for term in block),
                "field_provenance": "two_selfconsistent_predictions_no_target_fields",
                "unretained_band_meaning": "common_reference_tail_cancels_not_insulating_not_zero_field",
                "interpretation": "signed_response_times_activated_voltage_product_not_unique_causal_share",
            })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    out = args.out.resolve()
    coefficients = load_coefficients(out / "port_coefficients.csv")
    states = load_predictions(out)
    modes = port_rows(states, coefficients)
    terms, summaries = pair_rows(states, coefficients)
    write_csv(out / "predicted_port_modes.csv", modes)
    write_csv(out / "predicted_mode_effects.csv", terms)
    write_csv(out / "mode_effects_summary.csv", summaries)
    scales = mode_scale_rows(terms)
    write_csv(out / "predicted_mode_scale_summary.csv", scales)
    print(json.dumps({
        "prediction_method_states": len(states), "port_mode_rows": len(modes),
        "finite_change_mode_rows": len(terms), "pair_summaries": len(summaries),
        "selected_mode_scale_rows": len(scales),
        "phases_present": sorted({state["phase"] for state in states}),
        "maximum_absolute_identity_residual_S": max(
            (abs(row["identity_residual_vs_G_difference_S"]) for row in summaries), default=None),
        "target_files_read": 0, "output": str(out),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
