# -*- coding: utf-8 -*-
"""Summarize actual Round-7 COMSOL outputs without solving or changing data.

Default completion requires the prescribed synthetic target/background and four
refinement states.  --partial writes explicitly incomplete diagnostic tables.
CT input availability and CT target completion are tracked separately; neither
existence of image files nor a frozen CT plan counts as a solved CT target.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "runs/round7_3d_ct"


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def write_csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", encoding="utf-8-sig", newline="") as handle:
        if fields:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)


def finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def safe_ratio(a, b):
    return a/b if finite(a) and finite(b) and b != 0 else None


def relative_error(a, b):
    return abs(a-b)/abs(b) if finite(a) and finite(b) and b != 0 else None


def token(sf):
    return "sf_"+f"{sf:g}".replace(".", "p")


def path_string(path):
    return Path(path).resolve().as_posix()


def successful(path):
    if not Path(path).is_file():
        return None
    value = read_json(path)
    return value if value.get("status") == "success" else None


def clean_json(obj):
    if isinstance(obj, dict):
        return {k: clean_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [clean_json(v) for v in obj]
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--partial", action="store_true", help="Write incomplete diagnostic tables; never marks research complete")
    args = parser.parse_args()
    out = args.out.resolve()
    design = read_json(out/"design.json")
    cases = {c["id"]: c for c in [*design["cases"], *design["m0_cases"]]}
    missing, invalid = [], []
    state_rows, background_rows, geometry_rows, mesh_rows, model_rows = [], [], [], [], []
    states = {}
    geometry = {}
    backgrounds = {}
    refinements = {}

    def require(path, note):
        value = successful(path)
        if value is None:
            missing.append({"required": note, "path": path_string(path)})
        return value

    def index_model(folder, r, level, background):
        mph = folder/"solution.mph"
        if not mph.is_file():
            invalid.append({"issue": "success_without_solved_mph", "path": path_string(folder)})
        model_rows.append({
            "case_id": r.get("case_id"), "role": r.get("role"), "layout": r.get("layout"), "kind": r.get("kind"),
            "dimension": 3, "load": r.get("load"), "sigma_f_S_m": None if background else r.get("sigma_f_S_m"),
            "is_background": background, "mesh_level": level, "mph_path": path_string(mph),
            "success_json_path": path_string(folder/"SUCCESS.json"), "parameters_json_path": path_string(folder/"parameters.json"),
            "slices_mat_path": path_string(folder/"slices.mat"),
            "slice_y45mm_csv_path": path_string(folder/"slice_y45mm.csv"), "slice_y60mm_csv_path": path_string(folder/"slice_y60mm.csv"), "slice_y75mm_csv_path": path_string(folder/"slice_y75mm.csv"),
            "entry_script": path_string(ROOT/"scripts/r7_synthetic.m"), "reused_previous_round_solution": False,
            "comsol_version": r.get("comsol_version"), "solve_seconds": r.get("solve_seconds"), "completed_local": r.get("completed_local"),
        })

    for c in cases.values():
        branch = "m0" if c["kind"] == "synthetic3d_M0" else "synthetic3d"
        folder = out/branch/c["id"]
        bg = require(folder/"background_z/SUCCESS.json", c["id"]+" true z background")
        if bg:
            backgrounds[c["id"]] = bg
        axes = ["z", "x"] if branch == "synthetic3d" else ["z"]
        if c["id"] == design["cases"][0]["id"]:
            axes.append("y")
        # M0 x is an optional extra background written by the existing runner.
        if branch == "m0" and successful(folder/"background_x/SUCCESS.json"):
            axes.append("x")
        for axis in axes:
            r = bg if axis == "z" else require(folder/f"background_{axis}/SUCCESS.json", c["id"]+f" {axis} background")
            if not r:
                continue
            expected_sigma = c["geometry_checks"][f"sigma_{axis}_background_analytic_S_m"]
            expected_G = expected_sigma*0.12
            err = relative_error(r["G_S"], expected_G)
            if err is None or err > 1e-5:
                invalid.append({"issue": "background_analytic_mismatch", "case_id": c["id"], "axis": axis, "relative_error": err})
            background_rows.append({
                "case_id": c["id"], "role": c["role"], "layout": c["layout"], "kind": c["kind"], "load": axis,
                "G0_S": r["G_S"], "G0_analytic_S": expected_G, "G0_analytic_relative_error": err,
                "sigma_app_S_m": r.get("sigma_app_S_m"), "sigma_analytic_S_m": expected_sigma,
                "reaction_flux_relative_error": r.get("reaction_flux_relative_error"), "energy_relative_error": r.get("energy_relative_error"),
                "gradient_flux_relative_error": r.get("gradient_flux_relative_error"), "tetrahedra": r.get("tetrahedra"),
                "solve_seconds": r.get("solve_seconds"), "value_kind": "COMSOL actual passive background",
            })
            index_model(folder/f"background_{axis}", r, "base", True)

        geo_path = folder/"actual_geometry.json"
        if geo_path.is_file():
            geo = read_json(geo_path)
            geometry[c["id"]] = geo
            analytic = c["geometry_checks"]
            m0 = branch == "m0"
            row = {
                "case_id": c["id"], "role": c["role"], "layout": c["layout"], "kind": c["kind"],
                "Lx_mm": 120, "Ly_mm": 120, "Lz_mm": 120, "a_mm": 30, "b_mm": 30, "c_mm": 1.5,
                "sample_volume_mm3_actual": geo.get("sample_volume_mm3"),
                "ellipsoid_volume_mm3_actual": geo.get("ellipsoid_volume_mm3"), "ellipsoid_volume_mm3_analytic": analytic["ellipsoid_volume_mm3"],
                "ellipsoid_volume_relative_error": geo.get("ellipsoid_relative_volume_error"),
                "background_low_volume_mm3_actual": None if m0 else geo.get("original_low_volume_mm3"),
                "background_high_volume_mm3_actual": None if m0 else geo.get("original_high_volume_mm3"),
                "replaced_low_volume_mm3_actual": geo.get("replaced_low_volume_mm3"), "replaced_high_volume_mm3_actual": geo.get("replaced_high_volume_mm3"),
                "replaced_low_volume_mm3_analytic": analytic["replaced_low_volume_mm3"], "replaced_high_volume_mm3_analytic": analytic["replaced_high_volume_mm3"],
                "replaced_low_relative_error": relative_error(geo.get("replaced_low_volume_mm3"), analytic["replaced_low_volume_mm3"]),
                "replaced_high_relative_error": relative_error(geo.get("replaced_high_volume_mm3"), analytic["replaced_high_volume_mm3"]),
                "remaining_low_volume_mm3_actual": None if m0 else geo.get("remaining_low_volume_mm3"),
                "remaining_high_volume_mm3_actual": None if m0 else geo.get("remaining_high_volume_mm3"),
                "uniform_lower_volume_mm3_actual": geo.get("uniform_lower_volume_mm3"), "uniform_upper_volume_mm3_actual": geo.get("uniform_upper_volume_mm3"),
                "cap_minus_high_fraction_analytic": analytic["endcap_minus"]["high_fraction_analytic"],
                "cap_plus_high_fraction_analytic": analytic["endcap_plus"]["high_fraction_analytic"],
                "cap_volume_each_mm3_analytic": analytic["endcap_minus"]["volume_mm3_analytic"],
                "cap_volume_minus_mm3_actual_quadrature": bg.get("caps", {}).get("minus_volume_mm3") if bg else None,
                "cap_volume_plus_mm3_actual_quadrature": bg.get("caps", {}).get("plus_volume_mm3") if bg else None,
                "minimum_core_clearance_mm_analytic": analytic["minimum_core_interface_clearance_mm"],
                "layer_source_json": path_string(ROOT/c["source_design"]),
                "actual_geometry_json": path_string(geo_path), "domain_volumes_csv": path_string(folder/"domain_volumes.csv"),
                "volume_note": "M0 homogeneous exterior phase volumes are undefined; only explicit core replacement volumes and homogeneous-zone volumes are physical." if m0 else "Actual COMSOL order-8 phase-volume integration. Finite-cap phase labels are analytically proven, cap quadrature volume recorded separately.",
            }
            for axis, bounds in zip("xyz", geo.get("ellipse_mesh_bbox_mm", [[], [], []])):
                if len(bounds) == 2:
                    row[f"ellipse_{axis}_min_mm_mesh"] = bounds[0]
                    row[f"ellipse_{axis}_max_mm_mesh"] = bounds[1]
            geometry_rows.append(row)
        else:
            missing.append({"required": c["id"]+" actual geometry integrations", "path": path_string(geo_path)})

        for sf in [0.18, 10.0]:
            target_folder = folder/token(sf)
            r = require(target_folder/"SUCCESS.json", c["id"]+" "+token(sf)+" target")
            if not r:
                continue
            if not bg:
                invalid.append({"issue": "target_without_own_background", "case_id": c["id"], "sigma_f_S_m": sf})
                continue
            g = r["G_S"]-bg["G_S"]
            if finite(r.get("g_S")) and abs(g-r["g_S"]) > 1e-14:
                invalid.append({"issue": "target_saved_gain_mismatch", "case_id": c["id"], "sigma_f_S_m": sf})
            states[(c["id"], sf)] = {**r, "g_S": g, "G0_S": bg["G_S"], "relative_gain": g/bg["G_S"]}
            cap = r.get("caps", {})
            state_rows.append({
                "state_id": c["id"]+"_"+token(sf), "case_id": c["id"], "role": c["role"], "layout": c["layout"], "kind": c["kind"],
                "sigma_f_S_m": sf, "G0_S": bg["G_S"], "Gf_S": r["G_S"], "g_S": g, "eta": g/bg["G_S"],
                "sigma0_app_S_m": bg["G_S"]/.12, "sigmaf_app_S_m": r["G_S"]/.12, "delta_sigma_app_S_m": g/.12,
                "G_to_sigma_factor_per_m": 1/.12,
                "reaction_flux_relative_error": r.get("reaction_flux_relative_error"), "energy_relative_error": r.get("energy_relative_error"),
                "gradient_flux_relative_error": r.get("gradient_flux_relative_error"), "source_power_W": r.get("source_power_W"), "volume_dissipation_W": r.get("volume_dissipation_W"),
                "fracture_dissipation_W": r.get("fracture_dissipation_W"),
                "Jy_net_volume_integral_A_m": r.get("Jy_net_volume_integral_A_m"), "Jy_abs_volume_integral_A_m": r.get("Jy_abs_volume_integral_A_m"),
                "Jy_squared_volume_integral_A2_per_m": r.get("Jy_squared_volume_integral_A2_per_m"), "all_J_squared_volume_integral_A2_per_m": r.get("all_J_squared_volume_integral_A2_per_m"),
                "Jy_squared_fraction": r.get("Jy_squared_fraction"),
                "cap_minus_mean_V": cap.get("minus_mean_V"), "cap_plus_mean_V": cap.get("plus_mean_V"),
                "cap_mean_potential_difference_V": cap["plus_mean_V"]-cap["minus_mean_V"] if finite(cap.get("plus_mean_V")) and finite(cap.get("minus_mean_V")) else None,
                "tetrahedra": r.get("tetrahedra"), "min_mesh_quality": r.get("min_mesh_quality"), "solve_seconds": r.get("solve_seconds"),
                "success_json": path_string(target_folder/"SUCCESS.json"),
                "comparison_scope": "M0 reference; original exterior/electrode-zone material has been homogenized" if branch == "m0" else "Same-role perimeter shifts hold core, finite ellipsoid and explicit electrode buffers fixed",
            })
            index_model(target_folder, r, "base", False)

    for request in design["refinement_targets"]:
        cid, sf = request["case_id"], request["sigma_f_S_m"]
        folder = out/"precision"/(cid+"_"+token(sf))
        r = require(folder/"target/SUCCESS.json", cid+" "+token(sf)+" refined target")
        b = require(folder/"background/SUCCESS.json", cid+" "+token(sf)+" refined background")
        base = states.get((cid, sf))
        if not r or not b or not base:
            continue
        gf = r["G_S"]-b["G_S"]
        refinements[(cid, sf)] = {"target": r, "background": b, "g_S": gf, "delta_g_S": gf-base["g_S"]}
        geo_path = folder/"actual_geometry.json"
        fine_geo = read_json(geo_path) if geo_path.is_file() else {}
        base_geo = geometry.get(cid, {})
        mesh_rows.append({
            "case_id": cid, "role": base["role"], "layout": base["layout"], "sigma_f_S_m": sf,
            "base_G0_S": base["G0_S"], "fine_G0_S": b["G_S"], "base_Gf_S": base["G_S"], "fine_Gf_S": r["G_S"],
            "base_g_S": base["g_S"], "fine_g_S": gf, "gain_change_S_fine_minus_base": gf-base["g_S"], "gain_relative_change": relative_error(base["g_S"], gf),
            "base_tetrahedra": base.get("tetrahedra"), "fine_tetrahedra": r.get("tetrahedra"),
            "base_min_mesh_quality": base.get("min_mesh_quality"), "fine_min_mesh_quality": r.get("min_mesh_quality"),
            "base_ellipsoid_volume_mm3": base_geo.get("ellipsoid_volume_mm3"), "fine_ellipsoid_volume_mm3": fine_geo.get("ellipsoid_volume_mm3"),
            "base_ellipsoid_volume_relative_error": base_geo.get("ellipsoid_relative_volume_error"), "fine_ellipsoid_volume_relative_error": fine_geo.get("ellipsoid_relative_volume_error"),
            "base_replaced_low_volume_mm3": base_geo.get("replaced_low_volume_mm3"), "fine_replaced_low_volume_mm3": fine_geo.get("replaced_low_volume_mm3"),
            "base_replaced_high_volume_mm3": base_geo.get("replaced_high_volume_mm3"), "fine_replaced_high_volume_mm3": fine_geo.get("replaced_high_volume_mm3"),
            "fine_reaction_flux_relative_error": r.get("reaction_flux_relative_error"), "fine_energy_relative_error": r.get("energy_relative_error"),
            "selection_reason": request["reason"], "scope": "Discretization comparison at fixed ideal CAD geometry; finite-element representation also changes, no image-resolution issue in synthetic models.",
        })
        index_model(folder/"target", r, "fine", False)
        index_model(folder/"background", b, "fine", True)

    # Pairwise comparison uses the own-background gain, never raw conductivity
    # as a substitute and never the two-dimensional one-metre thickness rule.
    g_lookup = {(r["role"], r["layout"], sf): r for (_, sf), r in states.items() if r["kind"] == "synthetic3d_explicit"}
    effects = []
    max_abs_mesh_delta = max([abs(v["delta_g_S"]) for v in refinements.values()], default=None)

    def effect(name, role, sf, terms, scope, reference):
        selected = [(weight, g_lookup.get((r, layout, sf))) for weight, r, layout in terms]
        if any(state is None for _, state in selected):
            return
        value = sum(weight*state["g_S"] for weight, state in selected)
        fine = [(weight, refinements.get((state["case_id"], sf))) for weight, state in selected]
        have_fine = all(r is not None for _, r in fine)
        fine_value = sum(weight*r["g_S"] for weight, r in fine) if have_fine else None
        scale = sum(abs(weight) for weight, _ in selected)*max_abs_mesh_delta if max_abs_mesh_delta is not None else None
        source = g_lookup.get(reference)
        effects.append({
            "effect_id": name, "role": role, "sigma_f_S_m": sf, "definition": "; ".join(f"{weight:+g}*g({r},{layout})" for weight, r, layout in terms),
            "controlled_gain_difference_S": value, "controlled_delta_sigma_difference_S_m": value/.12,
            "relative_to_reference_gain": safe_ratio(value, source["g_S"] if source else None), "reference_role_layout": f"{reference[0]}_{reference[1]}",
            "comparison_scope": scope, "directly_refined_contrast": have_fine,
            "fine_controlled_gain_difference_S": fine_value,
            "controlled_difference_change_S_fine_minus_base": fine_value-value if have_fine else None,
            "absolute_sum_of_paired_gain_changes_S": sum(abs(weight*r["delta_g_S"]) for weight, r in fine) if have_fine else None,
            "empirical_global_mesh_scale_S": scale,
            "above_empirical_mesh_scale": abs(value) > scale if scale is not None else None,
            "mesh_scale_scope": "Empirical screen only: sum of absolute weights times largest change among four preselected refinements. Not a certified error bound or a direct refinement of every side contrast.",
        })

    contrasts = [
        ("lower_near_minus_far_upper_near", [(1, "P11"), (-1, "P71")]),
        ("lower_near_minus_far_upper_far", [(1, "P17"), (-1, "P77")]),
        ("upper_near_minus_far_lower_near", [(1, "P11"), (-1, "P17")]),
        ("upper_near_minus_far_lower_far", [(1, "P71"), (-1, "P77")]),
        ("joint_near_minus_far", [(1, "P11"), (-1, "P77")]),
        ("interaction", [(1, "P11"), (-1, "P17"), (-1, "P71"), (1, "P77")]),
    ]
    for role in ("A", "B"):
        for sf in (.18, 10.0):
            for name, terms in contrasts:
                effect(name, role, sf, [(w, role, layout) for w, layout in terms], "Pure prescribed peripheral relocation within one fixed core; actual CAD/mesh volume changes must still be inspected.", (role, "P77", sf))
    for layout in ("P11", "P17", "P71", "P77"):
        for sf in (.18, 10.0):
            effect("B_minus_A_"+layout, "B-A", sf, [(1, "B", layout), (-1, "A", layout)], "Core layer arrangement differs, finite-cap labels and analytical replacement volumes match; not a pure peripheral-only comparison.", ("A", layout, sf))

    # Actual A/B and same-role replacement-volume differences are not hidden
    # behind the analytically matched design.
    volume_pairs = []
    for layout in ("P11", "P17", "P71", "P77"):
        ids = [f"r7_strict_{role}_{layout}" for role in ("A", "B")]
        if all(cid in geometry for cid in ids):
            a, b = [geometry[cid] for cid in ids]
            ve = cases[ids[0]]["geometry_checks"]["ellipsoid_volume_mm3"]
            volume_pairs.append({
                "pair": "B_minus_A_"+layout, "case_minus": ids[0], "case_plus": ids[1],
                "low_replacement_difference_mm3": b["replaced_low_volume_mm3"]-a["replaced_low_volume_mm3"],
                "high_replacement_difference_mm3": b["replaced_high_volume_mm3"]-a["replaced_high_volume_mm3"],
                "low_replacement_difference_over_analytic_ellipsoid_volume": (b["replaced_low_volume_mm3"]-a["replaced_low_volume_mm3"])/ve,
                "high_replacement_difference_over_analytic_ellipsoid_volume": (b["replaced_high_volume_mm3"]-a["replaced_high_volume_mm3"])/ve,
                "analytical_replacement_difference_mm3": 0.0, "scope": "Actual integration residual reported separately; no forced numerical phase matching.",
            })
    for role in ("A", "B"):
        reference_id = f"r7_strict_{role}_P77"
        if reference_id not in geometry:
            continue
        a = geometry[reference_id]
        ve = cases[reference_id]["geometry_checks"]["ellipsoid_volume_mm3"]
        for layout in ("P11", "P17", "P71"):
            cid = f"r7_strict_{role}_{layout}"
            if cid not in geometry:
                continue
            b = geometry[cid]
            volume_pairs.append({
                "pair": role+"_"+layout+"_minus_P77", "case_minus": reference_id, "case_plus": cid,
                "low_replacement_difference_mm3": b["replaced_low_volume_mm3"]-a["replaced_low_volume_mm3"],
                "high_replacement_difference_mm3": b["replaced_high_volume_mm3"]-a["replaced_high_volume_mm3"],
                "low_replacement_difference_over_analytic_ellipsoid_volume": (b["replaced_low_volume_mm3"]-a["replaced_low_volume_mm3"])/ve,
                "high_replacement_difference_over_analytic_ellipsoid_volume": (b["replaced_high_volume_mm3"]-a["replaced_high_volume_mm3"])/ve,
                "analytical_replacement_difference_mm3": 0.0, "scope": "Same core and inclusion CAD; actual mesh-integration residual is retained.",
            })

    m0_comparisons = []
    for (role, layout, sf), full in g_lookup.items():
        m0 = states.get((f"r7_strict_{role}_M0", sf))
        if not m0:
            continue
        m0_comparisons.append({
            "role": role, "layout": layout, "sigma_f_S_m": sf,
            "explicit_G0_S": full["G0_S"], "M0_G0_S": m0["G0_S"], "background_relative_mismatch": relative_error(m0["G0_S"], full["G0_S"]),
            "explicit_Gf_S": full["G_S"], "M0_Gf_S": m0["G_S"], "explicit_g_S": full["g_S"], "M0_g_S": m0["g_S"],
            "M0_gain_error_S": m0["g_S"]-full["g_S"], "M0_gain_relative_error": relative_error(m0["g_S"], full["g_S"]),
            "M0_total_G_relative_error": relative_error(m0["G_S"], full["G_S"]),
            "same_M0_for_all_four_peripheral_layouts": True,
            "interpretation": "M0 has no peripheral-position label; its within-role position effects are exactly zero. This checks total-gain simplification, not a new 3D modal predictor.",
        })

    frozen_rows = []
    for (role, layout, sf), state in g_lookup.items():
        geo = geometry.get(state["case_id"])
        if not geo:
            continue
        # A passive, flat, z-loaded layer stack has spatially constant J0z.
        # Its background solution was independently checked against this exact
        # one-dimensional result. Actual replacement volumes come from COMSOL.
        current_density = state["G0_S"]/.0144  # DeltaV=1 V, electrode area=.0144 m^2
        E_low, E_high = current_density/.01, current_density/.1
        potential_low = geo["replaced_low_volume_mm3"]*1e-9*(sf-.01)*E_low**2
        potential_high = geo["replaced_high_volume_mm3"]*1e-9*(sf-.1)*E_high**2
        potential = potential_low+potential_high
        frozen_rows.append({
            "case_id": state["case_id"], "role": role, "layout": layout, "sigma_f_S_m": sf,
            "true_background_G0_S": state["G0_S"], "electrode_area_m2": .0144, "DeltaV_V": 1.0,
            "background_J0z_magnitude_A_m2": current_density, "background_E_low_magnitude_V_m": E_low, "background_E_high_magnitude_V_m": E_high,
            "replaced_low_volume_mm3_actual": geo["replaced_low_volume_mm3"], "replaced_high_volume_mm3_actual": geo["replaced_high_volume_mm3"],
            "frozen_potential_low_contribution_S": potential_low, "frozen_potential_high_contribution_S": potential_high,
            "frozen_potential_gain_S": potential, "actual_gain_S": state["g_S"],
            "redistribution_correction_inferred_S": potential-state["g_S"],
            "actual_gain_over_frozen_potential": safe_ratio(state["g_S"], potential),
            "correction_independently_field_integrated_this_round": False,
            "is_prospective_prediction": False,
            "basis": "Posthoc fixed-voltage energy identity: g_pot=int_DeltaSigma |E0|^2 dV/DeltaV^2; g_pot-g is inferred redistribution correction. E0 follows validated exact flat-layer background; actual g is read from full target. This is not a new target-field energy-integral verification or predictive solver.",
        })

    dimension_rows = []
    old_reference = out/"two_dimensional_reference.csv"
    if old_reference.is_file():
        with old_reference.open(encoding="utf-8-sig", newline="") as handle:
            old_rows = list(csv.DictReader(handle))
        for old in old_rows:
            role, layout, sf = old["role"], old["layout"], float(old["sigma_f_S_m"])
            current = g_lookup.get((role, layout, sf))
            if not current:
                continue
            eta2 = float(old["eta"])
            ds2 = float(old["delta_sigma_app_S_m"])
            geo = geometry.get(current["case_id"], {})
            fraction3 = safe_ratio(geo.get("ellipsoid_volume_mm3"), geo.get("sample_volume_mm3"))
            item = {
                "role": role, "layout": layout, "sigma_f_S_m": sf,
                "two_dimensional_source_case": old["case_id"], "finite_three_dimensional_case": current["case_id"],
                "delta_sigma_2d_app_S_m": ds2, "delta_sigma_3d_app_S_m": current["g_S"]/.12,
                "eta_2d": eta2, "eta_3d": current["relative_gain"],
                "eta_3d_over_2d": safe_ratio(current["relative_gain"], eta2),
                "inclusion_fraction_2d_actual": float(old["ellipse_area_fraction_actual"]),
                "inclusion_fraction_3d_actual": fraction3,
                "inclusion_fraction_3d_over_2d": safe_ratio(fraction3, float(old["ellipse_area_fraction_actual"])),
                "inclusion_fraction_ratio_analytic": float(old["finite3d_to_extruded2d_inclusion_fraction_ratio"]),
                "new_2d_COMSOL_solves": 0,
                "two_dimensional_truth_source": old["truth_file"],
                "interpretation_limit": "Finite y footprint and inclusion fraction change together with 3D bypass; amplitude ratios do not isolate a pure dimensional effect. Raw G is not compared across dimensions.",
            }
            lower_pair = ("P11", "P71") if layout in ("P11", "P71") else ("P17", "P77")
            upper_pair = ("P11", "P17") if layout in ("P11", "P17") else ("P71", "P77")
            for name, pair in [("lower_near", lower_pair), ("upper_near", upper_pair), ("joint_near", ("P11", "P77"))]:
                near, far = (g_lookup.get((role, value, sf)) for value in pair)
                item[name+"_pair"] = pair[0]+" minus "+pair[1]
                item[name+"_2d_delta_gain_sigma_app_S_m"] = float(old[name+"_delta_gain_sigma_app_S_m"])
                item[name+"_2d_gain_change_relative_to_far_gain"] = float(old[name+"_gain_change_relative_to_far_gain"])
                item[name+"_3d_delta_gain_sigma_app_S_m"] = (near["g_S"]-far["g_S"])/.12 if near and far else None
                item[name+"_3d_gain_change_relative_to_far_gain"] = (near["g_S"]-far["g_S"])/far["g_S"] if near and far else None
            dimension_rows.append(item)

    # Preserve CT evidence boundaries and do not overwrite the CT worker's CSV.
    ct_manifest = out/"ct_data_manifest.json"
    ct_response = out/"ct_response.csv"
    ct_summary = {
        "manifest_exists": ct_manifest.is_file(), "manifest_path": path_string(ct_manifest) if ct_manifest.is_file() else None,
        "response_csv_exists": ct_response.is_file(), "response_csv_path": path_string(ct_response) if ct_response.is_file() else None,
        "status": "CT must be reported from its own actual target/background records; this synthetic summary does not certify CT completion.",
    }
    if ct_response.is_file():
        with ct_response.open(encoding="utf-8-sig", newline="") as f:
            ct_summary["csv_rows"] = len(list(csv.DictReader(f)))
    complete = not missing and not invalid
    summary = {
        "schema": "round7_synthetic_actual_summary_v1", "synthetic_prescribed_work_complete": complete,
        "partial_diagnostic_mode": args.partial, "target_count_actual": len(states),
        "explicit_targets_actual": sum(s["kind"] == "synthetic3d_explicit" for s in states.values()),
        "M0_targets_actual": sum(s["kind"] == "synthetic3d_M0" for s in states.values()),
        "background_direction_states_actual": len(background_rows), "refined_target_background_pairs_actual": len(refinements),
        "explicit_background_direction_states_actual": sum(r["kind"] == "synthetic3d_explicit" for r in background_rows),
        "M0_background_direction_states_actual": sum(r["kind"] == "synthetic3d_M0" for r in background_rows),
        "solved_model_count_indexed": len(model_rows), "missing_required": missing, "invalid_checks": invalid,
        "background_max_relative_analytic_error": max([r["G0_analytic_relative_error"] for r in background_rows if finite(r["G0_analytic_relative_error"])], default=None),
        "max_reaction_flux_relative_error": max([r["reaction_flux_relative_error"] for r in state_rows if finite(r["reaction_flux_relative_error"])], default=None),
        "max_energy_relative_error": max([r["energy_relative_error"] for r in state_rows if finite(r["energy_relative_error"])], default=None),
        "max_gradient_flux_relative_error": max([r["gradient_flux_relative_error"] for r in state_rows if finite(r["gradient_flux_relative_error"])], default=None),
        "max_fine_gain_relative_change": max([r["gain_relative_change"] for r in mesh_rows if finite(r["gain_relative_change"])], default=None),
        "max_fine_absolute_gain_change_S": max_abs_mesh_delta,
        "M0_gain_relative_error_min": min([r["M0_gain_relative_error"] for r in m0_comparisons if finite(r["M0_gain_relative_error"])], default=None),
        "M0_gain_relative_error_max": max([r["M0_gain_relative_error"] for r in m0_comparisons if finite(r["M0_gain_relative_error"])], default=None),
        "explicit_Jy_squared_fraction_min": min([r["Jy_squared_fraction"] for r in state_rows if r["kind"] == "synthetic3d_explicit" and finite(r["Jy_squared_fraction"])], default=None),
        "explicit_Jy_squared_fraction_max": max([r["Jy_squared_fraction"] for r in state_rows if r["kind"] == "synthetic3d_explicit" and finite(r["Jy_squared_fraction"])], default=None),
        "finite_3D_diagnostic": "Nonzero Jy_squared_fraction or integral(abs(Jy)) demonstrates out-of-plane currents even when net Jy cancels by symmetry; this is not a conductivity tensor claim.",
        "mesh_scope": "Four preselected target/background refinements, two directly refined joint peripheral contrasts. Global mesh scale is empirical, not a guaranteed uncertainty bound.",
        "CT_evidence_status": ct_summary,
        "units": {"G_G0_g": "S", "sigma": "S/m", "G_to_sigma_factor": "1/0.12 m", "Jy_volume_integral": "A*m", "Jy_squared_volume_integral": "A^2/m", "volume": "mm^3"},
    }
    # A non-partial premature run leaves a machine-readable status but no new
    # success-labelled report. Partial data are explicitly tagged throughout.
    write_json(out/"synthetic_summary_status.json", clean_json(summary))
    if not complete and not args.partial:
        raise SystemExit(f"Round 7 synthetic results incomplete: {len(missing)} missing, {len(invalid)} invalid; use --partial only for diagnostics. Status saved.")
    completeness = "complete_prescribed_synthetic_work" if complete else "INCOMPLETE_DIAGNOSTIC_ONLY"
    outputs = {
        "synthetic3d_response.csv": state_rows,
        "background_direction_checks.csv": background_rows,
        "controlled_effects.csv": effects,
        "fine_joint_comparison.csv": [r for r in effects if r["effect_id"] == "joint_near_minus_far" and r["directly_refined_contrast"]],
        "geometry_manifest.csv": geometry_rows,
        "volume_matching_actual.csv": volume_pairs,
        "mesh_comparison.csv": mesh_rows,
        "M0_comparison.csv": m0_comparisons,
        "retrospective_frozen_field.csv": frozen_rows,
        "two_vs_three_dimensional_comparison.csv": dimension_rows,
        "model_index.csv": model_rows,
    }
    for name, rows in outputs.items():
        for row in rows:
            row["synthetic_completion_status"] = completeness
        write_csv(out/name, rows)
    if complete:
        write_json(out/"synthetic_summary.json", clean_json(summary))
    print(json.dumps({"synthetic_complete": complete, "actual_targets": len(states), "fine_pairs": len(refinements), "csv_outputs": list(outputs), "missing_required_count": len(missing), "invalid_count": len(invalid)}, ensure_ascii=True))


if __name__ == "__main__":
    main()
