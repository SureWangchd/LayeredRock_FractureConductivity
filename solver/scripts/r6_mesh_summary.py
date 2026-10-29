# -*- coding: utf-8 -*-
"""Summarize only existing Round-6 mesh evidence and solved model files.

This script never invokes a solver, changes predictions, or creates missing
measurements. It separates target-mesh change, composite-reference mesh change,
finite-modal truncation, interface quadrature and the numerical domination test.
"""
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "runs/round6_applicability"
ORDERS = (0, 8, 32, 64)


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def write_csv(path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def stamp():
    return datetime.now(timezone(timedelta(hours=8))).isoformat()


def sf_tag(sf):
    return "sf_"+format(float(sf), ".10g").replace(".", "p")


def state_id(state):
    # This is the criterion-table identifier, deliberately different from the
    # filesystem-safe design state_id which prefixes r6 and substitutes p.
    return f"L{int(state['full_length_mm'])}_{state['role']}_{state['layout']}_sf{float(state['sigma_f_S_m']):.10g}"


def meta(state):
    return {"state_id":state_id(state), "case_id":state["case_id"], "full_length_mm":state["full_length_mm"],
            "role":state["role"], "layout":state["layout"], "sigma_f_S_m":state["sigma_f_S_m"]}


def methods(prediction):
    values = prediction["methods"]
    if isinstance(values, dict):
        values = [values] if "K" in values else list(values.values())
    result = {int(v["K"]): v for v in values}
    assert set(result) == set(ORDERS), "A completed prediction requires all four frozen orders"
    return result


def optional(path, missing):
    if Path(path).is_file():
        return read_json(path)
    missing.append(str(path))
    return None


def existing_path(value, out):
    path = Path(value)
    return path if path.is_absolute() else out/path


def geometry_rows(out, design, missing):
    result = []
    fields = ["background_high_area_mm2", "background_low_area_mm2", "original_ellipse_area_mm2",
              "original_ellipse_replaced_high_mm2", "original_ellipse_replaced_low_mm2", "capminus_area_mm2",
              "capplus_area_mm2", "capminus_high_area_mm2", "capminus_low_area_mm2", "capplus_high_area_mm2",
              "capplus_low_area_mm2", "capminus_high_fraction", "capplus_high_fraction"]
    for case in design["cases"]:
        path = out/"truth"/case["id"]/"actual_geometry_and_background.json"
        data = optional(path, missing)
        if data is None:
            continue
        row = {"case_id":case["id"], "role":case["configuration_role"], "full_length_mm":case["full_length_mm"],
               "layout":case["layout"], "ellipse_axial_endcap_length_mm":case["endcap_length_mm"],
               "core_clearance_analytic_mm":case["geometry_only_checks"]["clearance_to_fixed_core_interface_mm"],
               "source":path.as_posix(), "area_units":"mm^2", "area_source":"COMSOL integrated domains on actual new geometry"}
        row.update({k:data.get(k) for k in fields})
        row["target_final_high_area_from_integrated_replacement_mm2"] = data["background_high_area_mm2"]-data["original_ellipse_replaced_high_mm2"]
        row["target_final_low_area_from_integrated_replacement_mm2"] = data["background_low_area_mm2"]-data["original_ellipse_replaced_low_mm2"]
        row["target_final_fracture_area_from_integrated_geometry_mm2"] = data["original_ellipse_area_mm2"]
        row["area_closure_error_mm2"] = sum(row[k] for k in ["target_final_high_area_from_integrated_replacement_mm2", "target_final_low_area_from_integrated_replacement_mm2", "target_final_fracture_area_from_integrated_geometry_mm2"])-14400
        row["replacement_high_minus_analytic_mm2"] = data["original_ellipse_replaced_high_mm2"]-case["expected_replaced_high_area_mm2"]
        row["replacement_low_minus_analytic_mm2"] = data["original_ellipse_replaced_low_mm2"]-case["expected_replaced_low_area_mm2"]
        row["endcap_matching_claim"] = "No matching claim across A/B; current finite cap fractions recorded separately"
        result.append(row)
    return result


def index_models(out):
    result = []
    for path in sorted(out.rglob("*.mph")):
        rel = path.relative_to(out)
        if rel.parts[0] not in {"truth", "precision", "reference"}:
            continue
        category = "other_new_round6_model"
        state = "saved_model"
        if rel.parts[0] == "reference":
            if path.name == "background_m0.mph":
                category = "COMSOL_M0_reference_background"
            else:
                category = "physical_passive_ellipse_source_background_for_reference"
        elif rel.parts[0] == "truth":
            if path.name == "solution.mph":
                category = "full_COMSOL_target"
            elif "source_background" in rel.parts:
                category = "full_physical_background_geometry_source"
            else:
                category = "full_physical_no_fracture_background"
        elif path.name == "fine_target.mph":
            category = "refined_full_COMSOL_target"
        elif path.name == "fine_background.mph":
            category = "refined_same_mesh_physical_background"
        companion = None
        for name in ["SUCCESS.json", "reference_summary.json", "fine_target_metrics.json", "fine_background_metrics.json", "background_summary.json"]:
            candidate = path.parent/name
            if candidate.is_file():
                companion = candidate
                state = "solved_with_saved_metrics"
                break
        params = path.parent/"parameters.json"
        p = read_json(params) if params.is_file() else {}
        result.append({"model_category":category,"case_id":p.get("id", ""),"full_length_mm":p.get("full_length_mm", ""),
                       "sigma_f_S_m":p.get("sigma_frac_S_m", "") if "target" in category else "",
                       "model_file":path.resolve().as_posix(),"relative_model_file":rel.as_posix(),
                       "model_status":state,"companion_metrics_file":companion.as_posix() if companion else "",
                       "size_MB":path.stat().st_size/1e6,
                       "interpretation":"Composite modal predictions are MATLAB fields plus reference MPH, not this MPH alone" if "reference" in category else "Saved actual COMSOL stationary model"})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--allow-partial", action="store_true", help="Write only observed data and mark incomplete; never fill absent results")
    args = parser.parse_args()
    out = args.out.resolve()
    design = read_json(out/"design.json")
    request = read_json(out/"refinement_targets.json")
    targets = request["targets"] if isinstance(request,dict) else request
    assert len(design["new_states"]) == 32 and len(design["cases"]) == 16 and len(targets) <= 6
    states = {(s["case_id"],float(s["sigma_f_S_m"])):s for s in design["new_states"]}
    missing, rows, per_state, base_predictions, target_results = [], [], {}, {}, {}
    trunc_base, trunc_fine = [], []
    for state in design["new_states"]:
        sid, cid, token = state_id(state),state["case_id"],sf_tag(state["sigma_f_S_m"])
        predicted_path = out/"predictions"/cid/token/"prediction.json"
        truth_path = out/"truth"/cid/token/"SUCCESS.json"
        p = optional(predicted_path,missing)
        truth = optional(truth_path,missing)
        if truth is not None:
            assert truth["status"] == "success" and truth["case_id"] == cid
            target_results[sid] = truth
        per_state[sid] = {"target_gain_change_S":None,"reference_gain_change_S":None,
                          "reference_gain_change_by_order_S":{},"target_mesh_available":False,"reference_mesh_available":False}
        if p is None:
            continue
        mm = methods(p)
        base_predictions[sid] = mm
        change = mm[64]["deltaG_S"]-mm[32]["deltaG_S"]
        trunc_base.append(abs(change))
        per_state[sid]["base_M64_minus_M32_gain_S"] = change
        rows.append({"evidence_type":"finite_modal_order_change",**meta(state),"mesh":"base","method_from":"M32","method_to":"M64",
                     "from_gain_S":mm[32]["deltaG_S"],"to_gain_S":mm[64]["deltaG_S"],"signed_change_S":change,
                     "relative_change":abs(change/mm[64]["deltaG_S"]),"source":predicted_path.as_posix(),
                     "scope":"Finite modal truncation change on one mesh; not a target mesh error"})
    target_changes, prediction_changes = [], {k:[] for k in ORDERS}
    completed_refinements = 0
    for entry in targets:
        state = states[(entry["case_id"],float(entry["sigma_f_S_m"]))]
        sid = state_id(state)
        target_folder = existing_path(entry.get("output",f"precision/{state['case_id']}_{sf_tag(state['sigma_f_S_m'])}"),out)
        fine_path = target_folder/"SUCCESS.json"
        fine = optional(fine_path,missing)
        if fine is not None:
            assert fine["status"] == "success" and fine["case_id"] == state["case_id"]
            change = fine["fine_minus_coarse_deltaG_S"]
            assert abs(change-(fine["fine_deltaG_S"]-fine["coarse_deltaG_S"])) < 1e-14
            completed_refinements += 1
            target_changes.append(abs(change))
            per_state[sid].update({"target_gain_change_S":change,"target_mesh_available":True,
                "coarse_target_gain_S":fine["coarse_deltaG_S"],"fine_target_gain_S":fine["fine_deltaG_S"],
                "target_refinement_file":fine_path.as_posix()})
            rows.append({"evidence_type":"full_COMSOL_target_mesh_change",**meta(state),"mesh":"base_to_fine","method_from":"full_COMSOL_base","method_to":"full_COMSOL_fine",
                         "from_gain_S":fine["coarse_deltaG_S"],"to_gain_S":fine["fine_deltaG_S"],"signed_change_S":change,
                         "relative_change":fine["deltaG_relative_change"],"base_G0_S":fine["coarse_G0_S"],"fine_G0_S":fine["fine_G0_S"],
                         "base_G_S":fine["coarse_G_S"],"fine_G_S":fine["fine_G_S"],"base_triangles":fine["source_triangles"],"fine_triangles":fine["fine_triangles"],
                         "target_flux_relative_error":fine["target"]["reaction_flux_relative_error"],"target_energy_relative_error":fine["target"]["energy_relative_error"],
                         "background_flux_relative_error":fine["background"]["reaction_flux_relative_error"],"background_energy_relative_error":fine["background"]["energy_relative_error"],
                         "source":fine_path.as_posix(),"scope":"Full-target gain uses its own same-fine-mesh no-fracture background"})
        pred_path = out/"prediction_precision"/state["case_id"]/sf_tag(state["sigma_f_S_m"])/"prediction.json"
        p = optional(pred_path,missing)
        if p is None or sid not in base_predictions:
            continue
        mm, base = methods(p),base_predictions[sid]
        per_state[sid]["reference_mesh_available"] = True
        for order in ORDERS:
            change = mm[order]["deltaG_S"]-base[order]["deltaG_S"]
            prediction_changes[order].append(abs(change))
            per_state[sid]["reference_gain_change_by_order_S"][f"M{order}"] = change
            if order == 64:
                per_state[sid]["reference_gain_change_S"] = change
            rows.append({"evidence_type":"composite_reference_core_mesh_change",**meta(state),"mesh":"base_to_fine","method_from":f"M{order}_base","method_to":f"M{order}_fine",
                         "from_gain_S":base[order]["deltaG_S"],"to_gain_S":mm[order]["deltaG_S"],"signed_change_S":change,
                         "relative_change":abs(change/base[order]["deltaG_S"]),"base_G0_S":base[order]["G0_S"],"fine_G0_S":mm[order]["G0_S"],
                         "base_G_S":base[order]["G_S"],"fine_G_S":mm[order]["G_S"],"source":pred_path.as_posix(),
                         "scope":"Composite model reference/core mesh trend, separate from physical full-target refinement"})
        change = mm[64]["deltaG_S"]-mm[32]["deltaG_S"]
        trunc_fine.append(abs(change))
        per_state[sid]["fine_M64_minus_M32_gain_S"] = change
        rows.append({"evidence_type":"finite_modal_order_change",**meta(state),"mesh":"fine","method_from":"M32","method_to":"M64",
                     "from_gain_S":mm[32]["deltaG_S"],"to_gain_S":mm[64]["deltaG_S"],"signed_change_S":change,
                     "relative_change":abs(change/mm[64]["deltaG_S"]),"source":pred_path.as_posix(),
                     "scope":"Finite modal truncation on fine composite mesh, separate from full COMSOL"})
    quadrature, margins = [], []
    for group in design["source_groups"]:
        for level in ["base","fine"]:
            path = out/"reference"/(group["source_key"]+"_"+level)/"reference_summary.json"
            summary = optional(path,missing)
            if summary is not None:
                trace = summary["trace"]
                quadrature.append(trace["gauss20_vs40_relative_fro"])
                rows.append({"evidence_type":"interface_projection_quadrature", "source_key":group["source_key"],"role":group["role"],"full_length_mm":group["full_length_mm"],"mesh":level,
                             "ndofs":summary["ndofs"],"gauss20_vs40_relative_fro":trace["gauss20_vs40_relative_fro"],"gauss20_vs40_max_abs_m":trace["gauss20_vs40_max_abs_m"],
                             "nonzero_background_max_abs_V":trace["nonzero_background_max_abs_V"],"source":path.as_posix(),
                             "scope":"Same trace basis,20-point vs40-point interface Gauss integration; not field/truth accuracy"})
            for sf in [.18,10.]:
                path = out/"predicted_backgrounds"/(group["source_key"]+"_"+level)/(sf_tag(sf)+"_conditions.json")
                # Precision orchestration only computes selected source/sigma groups.
                required = level=="base" or any(e["role"]==group["role"] and e["full_length_mm"]==group["full_length_mm"] and e["sigma_f_S_m"]==sf for e in targets)
                if not required and not path.is_file():
                    continue
                info = optional(path,missing)
                if info is None:
                    continue
                margins.append(info["domination_spectral_margin"])
                rows.append({"evidence_type":"finite_matrix_domination_condition","source_key":group["source_key"],"role":group["role"],"full_length_mm":group["full_length_mm"],"mesh":level,"sigma_f_S_m":sf,
                             "ndofs":info["ndofs"],"maximum_reference_trace_energy_eigenvalue":info["maximum_reference_trace_energy_eigenvalue"],
                             "domination_spectral_margin":info["domination_spectral_margin"],"numerical_domination_pass":info["numerical_domination_pass"],
                             "physical_matrix_relative_skew":info["physical_matrix_relative_skew"],"source":path.as_posix(),
                             "scope":"Observed finite-matrix sufficient condition, not interval-arithmetic or continuum proof"})
    geometry = geometry_rows(out,design,missing)
    if missing and not args.allow_partial:
        raise RuntimeError("Round6 evidence incomplete; no summary files written. Missing:\n"+"\n".join(missing))
    max_target = max(target_changes,default=None)
    global_scale = 2*max_target if max_target is not None else None
    evidence = {"generated_local":stamp(),"status":"complete" if not missing else "partial_observed_data_only",
        "full_target_count":len(target_results),"expected_full_target_count":32,"planned_refinement_count":len(targets),"completed_target_refinement_count":completed_refinements,
        "global_conservative_gain_resolution_S":global_scale,
        "resolution_definition":"2 times maximum observed absolute full-target fine-minus-base gain change among the prospectively selected refinements; contrasts multiply this per-state scale by sum(abs(weights))",
        "resolution_scope":"Descriptive conservative observed scale only; not a proven discretization error bound, not incorporated in finite-modal Q/B intervals",
        "maximum_abs_target_mesh_gain_change_S":max_target,
        "maximum_abs_composite_reference_mesh_gain_change_by_order_S":{f"M{k}":max(v,default=None) for k,v in prediction_changes.items()},
        "maximum_abs_base_M64_minus_M32_gain_S":max(trunc_base,default=None),
        "maximum_abs_fine_M64_minus_M32_gain_S":max(trunc_fine,default=None),
        "maximum_gauss20_vs40_relative_fro":max(quadrature,default=None),"minimum_domination_spectral_margin":min(margins,default=None),
        "per_state_reference_gain_change_definition":"M64 fine-minus-base; four individual orders retained in reference_gain_change_by_order_S",
        "per_state":per_state,"missing_files":missing,
        "separation":"Targetmesh, composite/coremesh, finite-modal order, interface quadrature and spectral condition are separate observations; no combined fulltruth guarantee"}
    write_csv(out/"mesh_comparison.csv",rows)
    write_json(out/"mesh_evidence.json",evidence)
    write_csv(out/"actual_geometry.csv",geometry)
    model_rows = index_models(out)
    write_csv(out/"model_index.csv",model_rows)
    print(json.dumps({"status":evidence["status"],"targets":len(target_results),"target_refinements":completed_refinements,
                      "composite_refinement_states":sum(x["reference_mesh_available"] for x in per_state.values()),
                      "mesh_rows":len(rows),"actual_geometry_rows":len(geometry),"model_rows":len(model_rows),
                      "observed_gain_resolution_S":global_scale,"missing_files":len(missing)},ensure_ascii=False))


if __name__ == "__main__":
    main()
