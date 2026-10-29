"""Summarize only Round 5 new geometry and numerical precision (stdlib).

No COMSOL calls, no predictor changes, and no old-case directory sweep.
Pending states are reported as pending rather than treated as zero results.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def optional(path: Path) -> dict:
    return read(path) if path.is_file() else {}


def sf_token(value: float) -> str:
    return "sf_" + format(value, ".10g").replace(".", "p")


def write_csv(path: Path, rows: list, first: list) -> None:
    fields = list(first)
    for row in rows:
        fields.extend(k for k in row if k not in fields)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def count_elements(stats: dict, kind: str):
    types, counts = stats.get("types", []), stats.get("numelem", [])
    if isinstance(types, str):
        types = [types]
    if isinstance(counts, (int, float)):
        counts = [counts]
    return counts[types.index(kind)] if kind in types else None


def meta(case: dict) -> dict:
    return {"case_id": case["id"], "role": case.get("configuration_role"),
            "layout": case.get("layout", case.get("layout_id"))}


def geometry_rows(out: Path, design: dict) -> list:
    cases = {c["id"]: c for c in design["cases"]}
    ids = list(dict.fromkeys(s["case_id"] for s in design["new_states"]))
    rows = []
    for case_id in ids:
        case = cases[case_id]
        folder = out / "truth" / case_id
        source = folder / "source_background"
        path = source / "SUCCESS.json"
        row = dict(meta(case), status="pending", actual_geometry_source=str(path))
        if not path.is_file():
            rows.append(row)
            continue
        result = read(path)
        bg = optional(folder / "background_summary.json")
        stats = optional(source / "mesh_stats.json")
        row.update(status=result.get("status"), sigma_low_S_m=case["sigma_low_S_m"],
                   sigma_high_S_m=case["sigma_high_S_m"], L_mm=case["L_mm"],
                   a_mm=case["a_mm"], b_mm=case["b_mm"],
                   theta_deg=case["theta_deg"], cx_mm=case["cx_mm"], cz_mm=case["cz_mm"],
                   out_of_plane_thickness_m=case["out_of_plane_thickness_m"],
                   geometry_partition_interfaces_mm=json.dumps(case["z_interfaces_mm"], separators=(",", ":")),
                   background_layer_phases=json.dumps(case["layer_materials"], separators=(",", ":")))
        for key in ["background_high_area_mm2", "background_low_area_mm2", "original_ellipse_area_mm2",
                    "original_ellipse_replaced_high_mm2", "original_ellipse_replaced_low_mm2",
                    "capminus_area_mm2", "capplus_area_mm2", "capminus_high_fraction", "capplus_high_fraction",
                    "triangles", "min_mesh_quality", "mean_mesh_quality", "solve_seconds", "elapsed_seconds",
                    "comsol_version", "analytic_expected_S_m", "analytic_relative_error",
                    "reaction_flux_relative_error", "reaction_power_relative_error"]:
            row[key] = result.get(key)
        row["G0_S"] = bg.get("G0_S", result.get("sigma_eff_S_m"))
        row["background_scalar_voltage_dofs"] = bg.get("dofs")
        row["mesh_edge_elements"] = count_elements(stats, "edg")
        row["mesh_vertex_elements_not_node_count"] = count_elements(stats, "vtx")
        row["mesh_node_count"] = result.get("mesh_nodes", result.get("nodes"))
        row["mesh_node_count_note"] = "Blank when not exported; scalar FE DOFs and vertex elements are not node counts."
        row["target_high_area_mm2"] = result["background_high_area_mm2"] - result["original_ellipse_replaced_high_mm2"]
        row["target_low_area_mm2"] = result["background_low_area_mm2"] - result["original_ellipse_replaced_low_mm2"]
        row["target_fracture_area_mm2"] = result["original_ellipse_area_mm2"]
        row["target_area_basis"] = "Actual COMSOL background/ellipse domain integrals on unchanged target geometry; final matrix areas by subtraction."
        row["background_area_closure_error_mm2"] = result["background_high_area_mm2"] + result["background_low_area_mm2"] - 14400
        row["target_area_closure_error_mm2"] = row["target_high_area_mm2"] + row["target_low_area_mm2"] + row["target_fracture_area_mm2"] - 14400
        row["high_background_area_error_mm2"] = result["background_high_area_mm2"] - 7200
        row["low_background_area_error_mm2"] = result["background_low_area_mm2"] - 7200
        row["high_replacement_error_mm2"] = result["original_ellipse_replaced_high_mm2"] - math.pi * 30 * 1.5 / 2
        row["low_replacement_error_mm2"] = result["original_ellipse_replaced_low_mm2"] - math.pi * 30 * 1.5 / 2
        rows.append(row)
    return rows


def target_states(out: Path, design: dict) -> list:
    states = [(s, "initial_48") for s in design["new_states"]]
    path = out / "crossover_targets.json"
    if path.is_file():
        data = read(path)
        values = data.get("targets", []) if isinstance(data, dict) else data
        states.extend((s, "limited_crossover") for s in values)
    unique = {}
    for state, phase in states:
        key = (state["case_id"], float(state["sigma_f_S_m"]))
        unique.setdefault(key, (state, phase))
    return list(unique.values())


def numerical_rows(out: Path, design: dict, geometry: list) -> list:
    cases = {c["id"]: c for c in design["cases"]}
    geom = {r["case_id"]: r for r in geometry}
    rows = []
    for state, phase in target_states(out, design):
        case_id, sf = state["case_id"], state["sigma_f_S_m"]
        folder = out / "truth" / case_id / sf_token(sf)
        path = folder / "SUCCESS.json"
        row = dict(meta(cases[case_id]), kind="full_COMSOL_target", phase=phase,
                   sigma_f_S_m=sf, mesh="base", status="pending", source=str(path))
        if path.is_file():
            result = read(path)
            for key in ["status", "G0_S", "G_S", "deltaG_S", "relative_gain", "solve_seconds", "postprocess_seconds",
                        "reaction_flux_relative_error", "energy_relative_error", "identity_relative_error",
                        "cross_relative_error", "matrix_energy_relative_error"]:
                row[key] = result.get(key)
            row["scalar_voltage_dofs"] = geom[case_id].get("background_scalar_voltage_dofs")
            for key in ["triangles", "min_mesh_quality", "mean_mesh_quality"]:
                row[key] = geom[case_id].get(key)
            for name in ["solution.mph", "fields.mat", "parameters.json"]:
                row[name.replace(".", "_") + "_saved"] = (folder / name).is_file()
            row["power_source_W"] = result.get("G_S")  # DV=1 V and t=1 m, fixed and saved.
            row["power_note"] = "DV=1 V; G numerically equals source power. Signed volume-power difference is not exported by r3_truth; relative balance is recorded."
        rows.append(row)
    for row_geom in geometry:
        case_id = row_geom["case_id"]
        path = out / "truth" / case_id / "background_summary.json"
        if not path.is_file():
            continue
        result = read(path)
        p = result.get("provenance", {})
        rows.append(dict(meta(cases[case_id]), kind="full_COMSOL_passive_background", phase="initial_48",
                         mesh="base", status="success", G_S=result["G0_S"],
                         scalar_voltage_dofs=result.get("dofs"),
                         triangles=row_geom.get("triangles"),
                         reaction_flux_relative_error=result.get("reaction_flux_relative_error"),
                         energy_relative_error=result.get("energy_relative_error"),
                         solve_seconds=p.get("background_solve_seconds"),
                         restored_solve_seconds=p.get("restored_background_solve_seconds"), source=str(path)))
    for path in sorted((out / "precision").glob("*/SUCCESS.json")):
        result = read(path)
        case_id = result["case_id"]
        case = cases[case_id]
        for state_name in ["target", "background"]:
            metrics = result[state_name]
            row = dict(meta(case), kind="full_COMSOL_" + state_name, phase="target_mesh_refinement",
                       mesh="fine", status=result.get("status"),
                       sigma_f_S_m=result["sigma_f_S_m"] if state_name == "target" else None,
                       source=str(path))
            for key in ["G_S", "reaction_positive_A", "reaction_ground_A", "power_source_W", "power_volume_W",
                        "reaction_flux_relative_error", "energy_relative_error", "triangles", "solve_seconds"]:
                row[key] = metrics.get(key)
            row["min_mesh_quality"] = result.get("mesh", {}).get("fine_min_mesh_quality")
            row["mean_mesh_quality"] = result.get("mesh", {}).get("fine_mean_mesh_quality")
            rows.append(row)
    for group in ["predictions", "prediction_precision", "prediction_precision_base", "crossover_predictions"]:
        for path in sorted((out / group).glob("*/sf_*/prediction.json")):
            prediction = read(path)
            for method in prediction["methods"]:
                row = {"case_id": prediction["case_id"], "role": prediction["role"], "layout": prediction["layout"],
                       "kind": "self_consistent_combined_prediction", "phase": group,
                       "sigma_f_S_m": prediction["sigma_f_S_m"], "mesh": prediction.get("mesh"),
                       "K": method["K"], "status": "saved_prediction", "source": str(path)}
                for key in ["G_S", "G0_S", "deltaG_S", "flux_balance_relative", "reaction_power_balance_relative",
                            "affine_power_balance_relative", "mean_current_balance_relative", "free_residual_relative",
                            "minimum_passivity_eigenvalue", "passivity_cholesky_flag", "ndofs", "nfree",
                            "total_power_W", "modal_correction_power_W"]:
                    row[key] = method.get(key)
                rows.append(row)
    return rows


def coarse_prediction_path(out: Path, prediction: dict) -> Path:
    folder = sf_token(prediction["sigma_f_S_m"])
    short = prediction["role"] + "_" + prediction["layout"]
    case_id = prediction["case_id"]
    candidates = [out / "prediction_precision_base" / short / folder / "prediction.json",
                  out / "prediction_precision_base" / case_id / folder / "prediction.json",
                  out / "development" / short / folder / "prediction.json",
                  out / "predictions" / case_id / folder / "prediction.json",
                  out / "crossover_predictions" / case_id / folder / "prediction.json"]
    return next((p for p in candidates if p.is_file()), candidates[0])


def mesh_rows(out: Path) -> list:
    rows = []
    for path in sorted((out / "precision").glob("*/SUCCESS.json")):
        result = read(path)
        row = {"kind": "full_target_mesh", "case_id": result["case_id"], "sigma_f_S_m": result["sigma_f_S_m"],
               "status": result["status"], "source_fine": str(path),
               "source_coarse": str(Path(result["source_target_folder"]) / "SUCCESS.json"),
               "scope": "Full physical COMSOL target and its own same-fine-mesh passive background; no predictor comparison."}
        for key in ["coarse_G_S", "fine_G_S", "coarse_G0_S", "fine_G0_S", "coarse_deltaG_S", "fine_deltaG_S",
                    "fine_minus_coarse_G_S", "fine_minus_coarse_G0_S", "fine_minus_coarse_deltaG_S",
                    "G_relative_change", "G0_relative_change", "deltaG_relative_change", "source_triangles", "fine_triangles"]:
            row[key] = result.get(key)
        row["gain_mesh_change_uS"] = result["fine_minus_coarse_deltaG_S"] * 1e6
        row["gain_mesh_change_percent"] = result["deltaG_relative_change"] * 100
        rows.append(row)
    for path in sorted((out / "prediction_precision").glob("*/sf_*/prediction.json")):
        fine = read(path)
        coarse_path = coarse_prediction_path(out, fine)
        coarse = optional(coarse_path)
        if coarse:
            assert coarse["case_id"] == fine["case_id"]
            assert math.isclose(coarse["sigma_f_S_m"], fine["sigma_f_S_m"], rel_tol=1e-13, abs_tol=1e-14)
        coarse_methods = {m["K"]: m for m in coarse.get("methods", [])}
        for method in fine["methods"]:
            row = {"kind": "combined_core_reference_mesh", "case_id": fine["case_id"], "role": fine["role"], "layout": fine["layout"],
                   "sigma_f_S_m": fine["sigma_f_S_m"], "K": method["K"], "status": "missing_corresponding_base_method",
                   "fine_G_S": method["G_S"], "fine_G0_S": method["G0_S"], "fine_deltaG_S": method["deltaG_S"],
                   "fine_scalar_voltage_dofs": method.get("ndofs"), "source_fine": str(path), "source_coarse": str(coarse_path),
                   "scope": "Same retained modes and analytic exterior; finite-element core/common-reference mesh changed. Separate from full-target mesh error."}
            base = coarse_methods.get(method["K"])
            if base:
                row.update(status="matched", coarse_G_S=base["G_S"], coarse_G0_S=base["G0_S"], coarse_deltaG_S=base["deltaG_S"],
                           source_scalar_voltage_dofs=base.get("ndofs"))
                for quantity in ["G", "G0", "deltaG"]:
                    difference = method[quantity + "_S"] - base[quantity + "_S"]
                    row["fine_minus_coarse_" + quantity + "_S"] = difference
                    row[quantity + "_relative_change"] = abs(difference) / max(abs(base[quantity + "_S"]), 1e-12)
                row["gain_mesh_change_uS"] = row["fine_minus_coarse_deltaG_S"] * 1e6
                row["gain_mesh_change_percent"] = row["deltaG_relative_change"] * 100
            elif method["K"] == 64 and fine["layout"] in {"P11", "P77"} and math.isclose(fine["sigma_f_S_m"], 10):
                row["status"] = "fine64_truncation_check_no_base_required"
                row["scope"] = "This extra fine M64 is used only against fine M32 for retained-mode convergence; no base M64 was required or computed."
            rows.append(row)
    return rows


def modal_truncation_rows(out: Path) -> list:
    rows = []
    for group, label in [("development", "development_base_mesh"),
                         ("prediction_precision", "selected_fine_mesh")]:
        for path in sorted((out / group).glob("*/sf_*/prediction.json")):
            prediction = read(path)
            methods = {m["K"]: m for m in prediction["methods"]}
            if not {32, 64}.issubset(methods):
                continue
            lower, higher = methods[32], methods[64]
            difference = higher["deltaG_S"] - lower["deltaG_S"]
            rows.append({"case_id": prediction["case_id"], "role": prediction["role"],
                         "layout": prediction["layout"], "sigma_f_S_m": prediction["sigma_f_S_m"],
                         "comparison_scope": label, "mesh": prediction.get("mesh"),
                         "K_lower": 32, "K_higher": 64,
                         "deltaG_M32_S": lower["deltaG_S"], "deltaG_M64_S": higher["deltaG_S"],
                         "M64_minus_M32_deltaG_S": difference,
                         "M64_minus_M32_deltaG_uS": difference * 1e6,
                         "gain_relative_change_denominator_M32": abs(difference) / max(abs(lower["deltaG_S"]), 1e-12),
                         "gain_change_percent": abs(difference) / max(abs(lower["deltaG_S"]), 1e-12) * 100,
                         "M64_minus_M32_G_S": higher["G_S"] - lower["G_S"],
                         "scalar_voltage_dofs": higher.get("ndofs"),
                         "meaning": "Same FE mesh and independently solved interface voltage; modal-order trend only, not a target-mesh error bound.",
                         "source": str(path)})
    return rows


def selected_effect_rows(out: Path) -> list:
    refined = {}
    for path in sorted((out / "precision").glob("*/SUCCESS.json")):
        result = read(path)
        if result.get("status") != "success":
            continue
        parts = result["case_id"].split("_")
        role, layout = parts[-2:]
        refined[(role, layout, float(result["sigma_f_S_m"]))] = (result, path)
    # First position is distance from lower core interface, second from upper.
    # Every side effect below holds the other side fixed. Both-side net and
    # factorial interaction explicitly use distinct comparison types.
    formulas = [
        ("upper_near_at_lower_far", "single_upper", [("P71", 1), ("P77", -1)]),
        ("lower_near_at_upper_far", "single_lower", [("P17", 1), ("P77", -1)]),
        ("upper_near_at_lower_near", "single_upper", [("P11", 1), ("P17", -1)]),
        ("lower_near_at_upper_near", "single_lower", [("P11", 1), ("P71", -1)]),
        ("both_near_minus_both_far", "net_both_sides", [("P11", 1), ("P77", -1)]),
        ("factorial_interaction", "interaction", [("P11", 1), ("P17", -1), ("P71", -1), ("P77", 1)]),
    ]
    candidates = []
    sigmas = sorted({key[2] for key in refined})
    for role in ["A", "B"]:
        for sf in sigmas:
            for effect, kind, terms in formulas:
                indexed = [((role, layout, sf), weight) for layout, weight in terms]
                if all(key in refined for key, weight in indexed):
                    candidates.append((effect, kind, role, "", sf, indexed))
    for sf in sigmas:
        for layout in ["P11", "P77", "P17", "P71", "P41", "P14"]:
            indexed = [(("B", layout, sf), 1), (("A", layout, sf), -1)]
            if all(key in refined for key, weight in indexed):
                candidates.append(("B_minus_A", "AB_ranking", "B_minus_A", layout, sf, indexed))
    rows = []
    for effect, kind, role, layout, sf, terms in candidates:
        coarse = sum(weight * refined[key][0]["coarse_deltaG_S"] for key, weight in terms)
        fine = sum(weight * refined[key][0]["fine_deltaG_S"] for key, weight in terms)
        scale = sum(abs(weight) * abs(refined[key][0]["fine_minus_coarse_deltaG_S"]) for key, weight in terms)
        same_sign = fine * coarse > 0
        resolved = same_sign and abs(fine) > 3 * scale
        rows.append({"effect": effect, "comparison_kind": kind, "role": role, "layout": layout,
                     "sigma_f_S_m": sf, "coarse_effect_S": coarse, "fine_effect_S": fine,
                     "coarse_effect_uS": coarse * 1e6, "fine_effect_uS": fine * 1e6,
                     "fine_minus_coarse_effect_uS": (fine - coarse) * 1e6,
                     "sum_abs_underlying_gain_mesh_changes_S": scale,
                     "sum_abs_underlying_gain_mesh_changes_uS": scale * 1e6,
                     "fine_signal_to_observed_mesh_scale": abs(fine) / scale if scale > 0 else None,
                     "same_sign_on_two_meshes": same_sign,
                     "resolved_on_3x_observed_mesh_scale": resolved,
                     "diagnostic_classification": "resolved_on_observed_scale" if resolved else "near_tie_or_not_resolved_on_observed_scale",
                     "diagnostic_rule": "Same sign on both meshes and |fine effect| > 3*sum(|underlying fine-minus-coarse gain|). Observed diagnostic scale, NOT a rigorous FEM error bound.",
                     "linear_combination": "; ".join(f"{weight:+d}*{key[0]}_{key[1]}" for key, weight in terms),
                     "sources": "; ".join(str(refined[key][1]) for key, weight in terms)})
    return rows


def max_abs(rows: list, key: str):
    values = [abs(row[key]) for row in rows if isinstance(row.get(key), (int, float)) and math.isfinite(row[key])]
    return max(values) if values else None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    out = args.root / "runs" / "round5_modal_closure"
    design = read(out / "design.json")
    geometry = geometry_rows(out, design)
    numerical = numerical_rows(out, design, geometry)
    mesh = mesh_rows(out)
    truncation = modal_truncation_rows(out)
    selected_effects = selected_effect_rows(out)
    write_csv(out / "actual_geometry.csv", geometry, ["case_id", "role", "layout", "status"])
    write_csv(out / "numerical_validation.csv", numerical, ["case_id", "role", "layout", "kind", "phase", "sigma_f_S_m", "mesh", "K", "status"])
    write_csv(out / "mesh_comparison.csv", mesh, ["kind", "case_id", "role", "layout", "sigma_f_S_m", "K", "status"])
    write_csv(out / "modal_truncation_comparison.csv", truncation, ["case_id", "role", "layout", "sigma_f_S_m", "comparison_scope", "mesh"])
    write_csv(out / "selected_mesh_effects.csv", selected_effects, ["effect", "comparison_kind", "role", "layout", "sigma_f_S_m"])
    complete = [r for r in numerical if r["kind"] == "full_COMSOL_target" and r.get("status") == "success"]
    initial = [r for r in complete if r["phase"] == "initial_48"]
    combined = [r for r in numerical if r["kind"] == "self_consistent_combined_prediction"]
    full_mesh = [r for r in mesh if r["kind"] == "full_target_mesh"]
    core_mesh = [r for r in mesh if r["kind"] == "combined_core_reference_mesh" and r["status"] == "matched"]
    summary = {"status": "initial_48_complete" if len(initial) == 48 else "partial_not_all_initial_targets_complete",
               "initial_targets_completed": len(initial), "initial_targets_expected": 48,
               "actual_geometry_cases_completed": sum(r["status"] == "success" for r in geometry),
               "full_target_refinements_completed": len(full_mesh), "combined_core_mesh_comparisons": len(core_mesh),
               "combined_core_missing_base_comparisons": sum(r["status"] == "missing_corresponding_base_method" for r in mesh),
               "fine64_extra_truncation_states_no_base_required": sum(r["status"] == "fine64_truncation_check_no_base_required" for r in mesh),
               "modal_truncation_development_comparisons": sum(r["comparison_scope"] == "development_base_mesh" for r in truncation),
               "modal_truncation_fine_comparisons": sum(r["comparison_scope"] == "selected_fine_mesh" for r in truncation),
               "max_development_M32_to_M64_gain_change_uS": max_abs([r for r in truncation if r["comparison_scope"] == "development_base_mesh"], "M64_minus_M32_deltaG_uS"),
               "max_fine_M32_to_M64_gain_change_uS": max_abs([r for r in truncation if r["comparison_scope"] == "selected_fine_mesh"], "M64_minus_M32_deltaG_uS"),
               "selected_fullmesh_effect_comparisons": len(selected_effects),
               "selected_effects_not_resolved_on_observed_mesh_scale": sum(not r["resolved_on_3x_observed_mesh_scale"] for r in selected_effects),
               "initial_target_max_reaction_flux_relative_error": max_abs(initial, "reaction_flux_relative_error"),
               "initial_target_max_energy_relative_error": max_abs(initial, "energy_relative_error"),
               "all_completed_target_max_reaction_flux_relative_error": max_abs(complete, "reaction_flux_relative_error"),
               "all_completed_target_max_energy_relative_error": max_abs(complete, "energy_relative_error"),
               "combined_prediction_max_flux_balance_relative": max_abs(combined, "flux_balance_relative"),
               "combined_prediction_max_power_balance_relative": max_abs(combined, "reaction_power_balance_relative"),
               "combined_prediction_max_free_residual_relative": max_abs(combined, "free_residual_relative"),
               "actual_background_max_area_closure_error_mm2": max_abs(geometry, "background_area_closure_error_mm2"),
               "actual_target_max_area_closure_error_mm2": max_abs(geometry, "target_area_closure_error_mm2"),
               "full_target_max_gain_mesh_change_uS": max_abs(full_mesh, "gain_mesh_change_uS"),
               "full_target_max_gain_mesh_change_percent": max_abs(full_mesh, "gain_mesh_change_percent"),
               "combined_core_max_gain_mesh_change_uS": max_abs(core_mesh, "gain_mesh_change_uS"),
               "combined_core_max_gain_mesh_change_percent": max_abs(core_mesh, "gain_mesh_change_percent"),
               "limits": ["Mesh trends are not rigorous error bounds.",
                          "Full-target mesh changes and combined-model core/reference mesh changes are separate numerical effects.",
                          "Good whole-domain current/power balance does not establish local flux accuracy or resolve tiny differences.",
                          "Final target phase areas use actual unchanged-geometry background/ellipse integrals by subtraction.",
                          "No unavailable node counts or unsolved target results are inferred."]}
    (out / "numerical_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
