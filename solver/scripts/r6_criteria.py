# -*- coding: utf-8 -*-
"""Prospective R6 output-specific applicability diagnostics and validation.

No COMSOL call, target field, fitted coefficient, or changed solver is used.
`predict` uses reference traces only to form Q/B and freezes all 32 criteria
before the full target gate. `validate` reads target scalars only after that
gate. Intervals concern a specified finite-modal model on its common FE space,
not the complete COMSOL model or a continuum/mesh error certificate.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone, timedelta
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "runs/round6_applicability"
EFFECTS = {
    "lower_7_to_1_at_upper_1": {"P11": 1, "P71": -1},
    "lower_7_to_1_at_upper_7": {"P17": 1, "P77": -1},
    "upper_7_to_1_at_lower_1": {"P11": 1, "P17": -1},
    "upper_7_to_1_at_lower_7": {"P71": 1, "P77": -1},
    "lower_average": {"P11": .5, "P17": .5, "P71": -.5, "P77": -.5},
    "upper_average": {"P11": .5, "P17": -.5, "P71": .5, "P77": -.5},
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
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def stamp():
    return datetime.now(timezone(timedelta(hours=8))).isoformat()


def sign(x):
    return 1 if x > 0 else -1 if x < 0 else 0


def rel(error, reference):
    return abs(error / reference) if reference else None


def sf_tag(value):
    return "sf_" + format(float(value), ".10g").replace(".", "p")


def key(length, role, layout, sf):
    return (int(float(length)), str(role), str(layout), float(sf))


def state_id(k):
    return f"L{k[0]}_{k[1]}_{k[2]}_sf{format(k[3], '.10g')}"


def method_dict(p):
    methods = p["methods"]
    if isinstance(methods, dict):
        methods = [methods] if "K" in methods else list(methods.values())
    result = {int(m["K"]): m for m in methods}
    if set(result) != {0, 8, 32, 64}:
        raise ValueError("Expected exactly M0, M8, M32, M64")
    for m in result.values():
        if len(m["port_v"]) != 2 or any(len(v) != 65 for v in m["port_v"]):
            raise ValueError("Expected 2 x 65 reference trace")
        vals = [m[x] for x in ("G0_S", "G_S", "deltaG_S")] + [v for side in m["port_v"] for v in side]
        if not all(math.isfinite(v) for v in vals):
            raise ValueError("Nonfinite prediction")
        if abs(m["G_S"] - m["G0_S"] - m["deltaG_S"]) > 1e-11:
            raise ValueError("Prediction gain/background inconsistent")
    return result


def load_predictions(out):
    design = read_json(ROOT / "round6_taskpack/round6_design.json")
    expected = {key(s["full_length_mm"], s["role"], s["layout"], s["sigma_f_S_m"]): s for s in design["states"]}
    result = {}
    for path in sorted((out / "predictions").glob("*/sf_*/prediction.json")):
        p = read_json(path)
        length = p.get("length_mm", p.get("full_length_mm"))
        if length is None:
            import re
            length = int(re.search(r"L(40|64)", p["case_id"])[1])
        k = key(length, p["role"], p["layout"], p["sigma_f_S_m"])
        if k not in expected:
            raise ValueError("Prediction outside frozen geometry: " + str(path))
        if k in result:
            raise ValueError("Duplicate prediction state: " + str(k))
        result[k] = {"p": p, "methods": method_dict(p), "path": path, "design": expected[k]}
    if set(result) != set(expected):
        raise ValueError(f"Need all 32 saved predictions, found {len(result)}; missing {sorted(set(expected)-set(result))}")
    return result


def coefficients(out):
    path = out / "port_coefficients.csv"
    if not path.exists():
        path = ROOT / "runs/round5_modal_closure/port_coefficients.csv"
    result = {}
    for r in read_csv(path):
        result[r["layout"], r["side"], int(r["n"])] = (float(r["Y_reference_S_m2"]), float(r["Y_actual_S_m2"]))
    return result, path


def interval(reference_gain, trace, coeff, layout, reference_order, target_order, background_difference=0.):
    """Pure-Python independent evaluation of supplied Q/B, with bg correction.

    The two backgrounds may have small saved numerical differences. The bound
    first applies to G; subtracting (G0_target-G0_reference) yields the g bound.
    Target solved traces and target gain do not enter this calculation.
    """
    Q, B, absolute, sides = 0., 0., 0., []
    for side_index, side in enumerate(("lower", "upper")):
        ym0, yt0 = coeff[layout, side, 0]
        if min(ym0, yt0) <= 0 or not math.isclose(ym0, yt0, rel_tol=1e-10, abs_tol=1e-13):
            raise ValueError("Different/nonpositive zero modes invalidate the common drive formula")
        qs, bs, terms = 0., 0., []
        for n in range(1, target_order + 1):
            ym, yt = coeff[layout, side, n]
            if min(ym, yt) <= 0:
                raise ValueError("Full peripheral admittance must be positive")
            yr = yt if n <= reference_order else ym
            dy = yt - yr
            weight = .06 * trace[side_index][n] ** 2  # t=DV=1, W/2=.06 m
            q, b = weight * dy, weight * dy * dy / yt
            qs += q
            bs += b
            terms.append(q)
        Q += qs
        B += bs
        absolute += sum(abs(v) for v in terms)
        sides.append({"side": side, "Q_S": qs, "B_S": bs, "terms_S": terms})
    low, high = Q - B - background_difference, Q - background_difference
    lower, upper = reference_gain + low, reference_gain + high
    bound = max(abs(low), abs(high))
    return {
        "Q_S": Q, "B_S": B, "background_difference_S": background_difference,
        "correction_lower_S": low, "correction_upper_S": high,
        "gain_lower_S": lower, "gain_upper_S": upper,
        "absolute_reference_error_bound_S": bound,
        "relative_reference_error_bound": bound/lower if lower > 0 else None,
        "absolute_weighted_response_sum_S": absolute,
        "lower_Q_S": sides[0]["Q_S"], "upper_Q_S": sides[1]["Q_S"],
        "lower_B_S": sides[0]["B_S"], "upper_B_S": sides[1]["B_S"],
    }


def base_meta(k, e):
    return {"state_id": state_id(k), "case_id": e["p"]["case_id"], "full_length_mm": k[0],
            "role": k[1], "layout": k[2], "sigma_f_S_m": k[3]}


def make_criteria(saved, coeff):
    rows, indexed, diagnostics = [], {}, []
    for k, e in saved.items():
        for kr in (0, 8):
            reference = e["methods"][kr]
            for n in (32, 64):
                target = e["methods"][n]
                d = interval(reference["deltaG_S"], reference["port_v"], coeff, k[2], kr, n,
                             target["G0_S"] - reference["G0_S"])
                r = {**base_meta(k, e), "reference_K": kr, "target_N": n,
                     "reference_gain_S": reference["deltaG_S"], **d,
                     "primary_tolerance": .02, "bound_scope": f"M{kr}_to_M{n}_same_FE_space_only",
                     "target_field_used": False, "target_gain_used_in_criterion": False,
                     "reference_trace_file": str(e["path"]),
                     "reference_stationarity_relative": reference.get("free_residual_relative"),
                     "reference_minimum_passivity_eigenvalue": reference.get("minimum_passivity_eigenvalue")}
                for tol in (1, 2, 5):
                    r[f"accepted_{tol}pct"] = d["relative_reference_error_bound"] is not None and d["relative_reference_error_bound"] <= tol/100
                r["label_2pct"] = "usable_for_specified_finite_model" if r["accepted_2pct"] else "cannot_confirm"
                rows.append(r)
                indexed[k, kr, n] = r
                diagnostics.append({**base_meta(k, e), "reference_K": kr, "target_N": n,
                                    "target_finite_gain_S": target["deltaG_S"],
                                    "gain_lower_S": d["gain_lower_S"], "gain_upper_S": d["gain_upper_S"],
                                    "finite_model_contained": d["gain_lower_S"]-1e-12 <= target["deltaG_S"] <= d["gain_upper_S"]+1e-12,
                                    "diagnostic_not_used_to_construct_bound": True})
    return rows, indexed, diagnostics


def contrast_row(saved, intervals, keys_weights, kr, n, name, category, meta):
    lower = upper = point = 0.
    for k, a in keys_weights:
        r = intervals[k, kr, n]
        lower += a * (r["gain_lower_S"] if a >= 0 else r["gain_upper_S"])
        upper += a * (r["gain_upper_S"] if a >= 0 else r["gain_lower_S"])
        point += a * saved[k]["methods"][kr]["deltaG_S"]
    distance = min(abs(lower), abs(upper)) if lower * upper > 0 else 0.
    absolute_bound = max(abs(point-lower), abs(point-upper))
    eta = absolute_bound/distance if distance else None
    return {**meta, "contrast": name, "category": category, "reference_K": kr, "target_N": n,
            "reference_contrast_S": point, "contrast_lower_S": lower, "contrast_upper_S": upper,
            "sign_decision": "positive" if lower > 0 else "negative" if upper < 0 else "cannot_determine",
            "absolute_reference_error_bound_S": absolute_bound,
            "relative_amplitude_error_bound": eta,
            "amplitude_10pct_confirmed": eta is not None and eta <= .1 if category == "single_side" else None,
            "amplitude_tolerance": .1 if category == "single_side" else None,
            "contrast_terms": json.dumps([{ "state_id": state_id(k), "weight": a} for k,a in keys_weights]),
            "scope": "specified_finite_model; no assumed cancellation of common error",
            "net_sign_will_not_be_chased": category in ("joint_net", "interaction")}


def make_contrasts(saved, intervals):
    rows = []
    for length in (40, 64):
        for sf in (.18, 10.):
            for kr in (0, 8):
                for n in (32, 64):
                    for role in ("A", "B"):
                        for name, weights in EFFECTS.items():
                            category = "joint_net" if name.startswith("joint") else "interaction" if name.startswith("interaction") else "single_side"
                            rows.append(contrast_row(saved, intervals,
                                [(key(length, role, layout, sf), a) for layout,a in weights.items()], kr,n,name,category,
                                {"full_length_mm": length, "role": role, "sigma_f_S_m": sf, "layout": "multiple"}))
                    for layout in ("P11", "P77", "P17", "P71"):
                        rows.append(contrast_row(saved, intervals,
                            [(key(length,"B",layout,sf),1), (key(length,"A",layout,sf),-1)], kr,n,"B_minus_A","AB_rank",
                            {"full_length_mm": length, "role": "B_minus_A", "sigma_f_S_m": sf, "layout": layout}))
    return rows


def weighted_quantile_mode(weights, fraction):
    total = sum(weights)
    if total <= 0:
        return None
    cumulative = 0.
    for n, w in enumerate(weights, 1):
        cumulative += w
        if cumulative >= fraction * total:
            return n
    return len(weights)


def source_rows(saved, coeff):
    rows = []
    for k,e in saved.items():
        for kr in (0,8):
            for i,side in enumerate(("lower","upper")):
                v = e["methods"][kr]["port_v"][i]
                energy, weighted, terms, bs = [], [], [], []
                for n in range(1,65):
                    yr0,yt = coeff[k[2],side,n]
                    yr = yt if n <= kr else yr0
                    energy.append(.06*yr*v[n]**2)
                    q = .06*(yt-yr)*v[n]**2
                    weighted.append(abs(q)); terms.append(q)
                    bs.append(.06*(yt-yr)**2/yt*v[n]**2)
                n50,n90 = weighted_quantile_mode(energy,.5), weighted_quantile_mode(energy,.9)
                nq90 = weighted_quantile_mode(weighted,.9)
                trace_weight=[.06*x*x for x in v[1:]]
                trace_n50=weighted_quantile_mode(trace_weight,.5)
                trace_n90=weighted_quantile_mode(trace_weight,.9)
                dl,du=e["design"]["lower_distance_mm"],e["design"]["upper_distance_mm"]
                lower_near=key(k[0],k[1],f"P1{du}",k[3]); lower_far=key(k[0],k[1],f"P7{du}",k[3])
                upper_near=key(k[0],k[1],f"P{dl}1",k[3]); upper_far=key(k[0],k[1],f"P{dl}7",k[3])
                model_gain=lambda state:saved[state]["methods"][kr]["deltaG_S"]
                r = {**base_meta(k,e), "reference_K": kr, "side": side,
                     "clearance_to_fixed_core_interface_mm": e["design"]["clearance_to_fixed_core_interface_mm"],
                     "lower_distance_mm": e["design"]["lower_distance_mm"], "upper_distance_mm": e["design"]["upper_distance_mm"],
                     "nonuniform_trace_rms_V": math.sqrt(.5*sum(x*x for x in v[1:])),
                     "trace_v_squared_n50":trace_n50,"trace_v_squared_n90":trace_n90,
                     "trace_half_wavelength_at_n90_mm":120/trace_n90 if trace_n90 else None,
                     "reference_nonzero_mode_energy_S": sum(energy), "energy_n50": n50, "energy_n90": n90,
                     "energy_half_wavelength_at_n90_mm": 120/n90 if n90 else None,
                     "weighted_response_n90": nq90, "Q_side_S": sum(terms), "B_side_S": sum(bs),
                     "absolute_Q_side_sum_S": sum(weighted),
                     "reference_lower_near_effect_at_this_upper_S":model_gain(lower_near)-model_gain(lower_far),
                     "reference_upper_near_effect_at_this_lower_S":model_gain(upper_near)-model_gain(upper_far),
                     "max_absolute_delta_Y_S_m2": max(abs(coeff[k[2],side,n][1]-coeff[k[2],side,n][0]) for n in range(1,65)),
                     "interpretation": "reference excitation/response diagnostic; modal terms are not unique causal contributions"}
                for start,end in ((1,2),(3,8),(9,16),(17,32),(33,64)):
                    r[f"Q_n{start}_{end}_S"] = sum(terms[start-1:end])
                    r[f"absolute_Q_n{start}_{end}_S"] = sum(weighted[start-1:end])
                    r[f"energy_n{start}_{end}_S"] = sum(energy[start-1:end])
                rows.append(r)
    return rows


def prospective(out):
    out.mkdir(parents=True,exist_ok=True)
    if (out/"criteria_ready.json").exists():
        print(json.dumps(read_json(out/"criteria_ready.json"), ensure_ascii=False,indent=2))
        return
    if list((out/"truth").glob("*/sf_*/SUCCESS.json")):
        raise RuntimeError("Refuse creating prospective criteria after target truths exist")
    saved = load_predictions(out)
    coeff,path = coefficients(out)
    rows,indexed,diagnostics = make_criteria(saved,coeff)
    contrasts = make_contrasts(saved,indexed)
    write_csv(out/"applicability_predictions.csv",rows)
    write_json(out/"applicability_predictions.json",rows)
    write_csv(out/"side_and_rank_intervals_pretruth.csv",contrasts)
    write_json(out/"side_and_rank_intervals_pretruth.json",contrasts)
    sources=source_rows(saved,coeff)
    write_csv(out/"source_scale_response.csv",sources)
    write_csv(out/"source_scale_response_pretruth.csv",sources)
    write_csv(out/"finite_model_interval_diagnostic.csv",diagnostics)
    summary = {"saved_local":stamp(),"complete":True,"states":len(saved),"interval_rows":len(rows),
               "contrast_rows":len(contrasts),"target_truth_read":False,"coefficients_file":str(path),
               "scope":"Finite-modal same-space criteria; not full-COMSOL error certification",
               "settings":{"reference_orders":[0,8],"target_orders":[32,64],"primary_tolerance":.02,
                           "secondary_tolerances":[.01,.05],"single_side_amplitude_tolerance":.1,
                           "width_m":.12,"thickness_m":1,"voltage_V":1},
               "finite_model_containment_diagnostic":{"contained":sum(r["finite_model_contained"] for r in diagnostics),"total":len(diagnostics)},
               "acceptance_counts":{f"M{k}_to_M{n}":{f"{tol}pct":sum(r[f"accepted_{tol}pct"] for r in rows if r["reference_K"]==k and r["target_N"]==n) for tol in (1,2,5)} for k in (0,8) for n in (32,64)}}
    write_json(out/"criteria_ready.json",summary)
    print(json.dumps(summary,ensure_ascii=False,indent=2))


def truth_entries(out,saved):
    if not (out/"predictions_ready.json").exists() or not (out/"criteria_ready.json").exists():
        raise RuntimeError("Prediction/criteria all-batch gate required before reading targets")
    result = {}
    for k,e in saved.items():
        path = out/"truth"/e["p"]["case_id"]/sf_tag(k[3])/"SUCCESS.json"
        if not path.exists():
            raise ValueError("Missing full target: "+str(path))
        value = read_json(path)
        if value.get("status") != "success" or abs(value["sigma_f_S_m"]-k[3])>1e-12:
            raise ValueError("Invalid target: "+str(path))
        if path.stat().st_mtime < (out/"criteria_ready.json").stat().st_mtime:
            raise ValueError("Target predates prospective criteria: "+str(path))
        result[k] = {"truth":value,"path":path}
    return result


def load_mesh_evidence(out):
    """Optional root-generated case-level evidence; no fake uncertainty bounds.

    root may write mesh_evidence.json with global_conservative_gain_resolution_S
    and per_state[state_id] {target_gain_change_S,reference_gain_change_S}.
    A global observed scale is descriptive, not a rigorous error bound.
    """
    path=out/"mesh_evidence.json"
    return read_json(path) if path.exists() else {}


def validation(out):
    saved=load_predictions(out)
    truths=truth_entries(out,saved)
    criteria=read_json(out/"applicability_predictions.json")
    contrasts=read_json(out/"side_and_rank_intervals_pretruth.json")
    mesh=load_mesh_evidence(out)
    byid={state_id(k):k for k in saved}
    comparisons=[]
    for k,e in saved.items():
        t=truths[k]["truth"]
        for order,m in sorted(e["methods"].items()):
            error=m["deltaG_S"]-t["deltaG_S"]
            comparisons.append({**base_meta(k,e),"K":order,"G0_prediction_S":m["G0_S"],"G0_COMSOL_S":t["G0_S"],
                "G_prediction_S":m["G_S"],"G_COMSOL_S":t["G_S"],"gain_prediction_S":m["deltaG_S"],"gain_COMSOL_S":t["deltaG_S"],
                "gain_signed_error_S":error,"gain_absolute_error_S":abs(error),"gain_relative_error":rel(error,t["deltaG_S"]),
                "M64_minus_M32_gain_S":e["methods"][64]["deltaG_S"]-e["methods"][32]["deltaG_S"],
                "truth_file":str(truths[k]["path"]),"prediction_file":str(e["path"])})
    write_csv(out/"closure_new_source_comparison.csv",comparisons)
    validation_rows=[]
    for r in criteria:
        k=byid[r["state_id"]]; e=saved[k]; t=truths[k]["truth"]
        kr,n=r["reference_K"],r["target_N"]
        reference=e["methods"][kr]["deltaG_S"]
        finite=e["methods"][n]["deltaG_S"]
        actual=rel(reference-t["deltaG_S"],t["deltaG_S"])
        mf=mesh.get("per_state",{}).get(r["state_id"],{})
        for tol in (1,2,5):
            accepted=r[f"accepted_{tol}pct"]; passed=actual<=tol/100
            status="accepted_and_passed" if accepted and passed else "false_accept_full_COMSOL" if accepted else "conservative_unconfirmed_but_passed" if passed else "unconfirmed_and_failed_tolerance"
            validation_rows.append({**{x:r[x] for x in ("state_id","case_id","full_length_mm","role","layout","sigma_f_S_m")},
                "reference_K":kr,"target_N":n,"tolerance":tol/100,"prospective_accepted":accepted,
                "COMSOL_relative_gain_error":actual,"COMSOL_pass":passed,"validation_class":status,
                "finite_model_gain_S":finite,"finite_model_contained":r["gain_lower_S"]-1e-12<=finite<=r["gain_upper_S"]+1e-12,
                "gain_lower_S":r["gain_lower_S"],"gain_upper_S":r["gain_upper_S"],"COMSOL_gain_S":t["deltaG_S"],
                "full_COMSOL_contained_descriptive_only":r["gain_lower_S"]<=t["deltaG_S"]<=r["gain_upper_S"],
                "M64_minus_M32_gain_S":e["methods"][64]["deltaG_S"]-e["methods"][32]["deltaG_S"],
                "finite_model_minus_COMSOL_gain_S":finite-t["deltaG_S"],
                "target_mesh_gain_change_S":mf.get("target_gain_change_S"),"reference_mesh_gain_change_S":mf.get("reference_gain_change_S"),
                "interpretation":"No full-COMSOL guarantee was asserted; truncation and mesh differences are separate diagnostics"})
    write_csv(out/"applicability_validation.csv",validation_rows)
    contrast_rows=[]
    observed_scale=mesh.get("global_conservative_gain_resolution_S")
    for r in contrasts:
        terms=json.loads(r["contrast_terms"])
        weights=[(byid[x["state_id"]],x["weight"]) for x in terms]
        truth=sum(a*truths[k]["truth"]["deltaG_S"] for k,a in weights)
        finite=sum(a*saved[k]["methods"][r["target_N"]]["deltaG_S"] for k,a in weights)
        error=r["reference_contrast_S"]-truth
        scale=sum(abs(a) for _,a in weights)*observed_scale if observed_scale is not None else None
        resolved=abs(truth)>scale if scale is not None else None
        decision=r["sign_decision"]
        direction_correct=sign(truth)==(1 if decision=="positive" else -1) if decision!="cannot_determine" else None
        contrast_rows.append({**r,"COMSOL_contrast_S":truth,"finite_model_contrast_S":finite,
            "finite_model_contained":r["contrast_lower_S"]-1e-12<=finite<=r["contrast_upper_S"]+1e-12,
            "full_COMSOL_contained_descriptive_only":r["contrast_lower_S"]<=truth<=r["contrast_upper_S"],
            "reference_absolute_error_S":abs(error),"reference_relative_amplitude_error":rel(error,truth) if r["category"]=="single_side" else None,
            "sign_decision_correct_against_base_COMSOL":direction_correct,
            "observed_mesh_resolution_scale_S":scale,"resolved_against_observed_mesh_scale":resolved,
            "amplitude_10pct_actual_base_COMSOL":rel(error,truth)<=.1 if r["category"]=="single_side" and truth else None,
            "resolution_note":"Observed conservative scale from selected meshes; not a universal certified error bound" if scale is not None else "mesh scale not supplied; no resolution claim"})
    write_csv(out/"side_and_rank_intervals.csv",contrast_rows)
    sources=read_csv(out/"source_scale_response_pretruth.csv")
    for r in sources:
        k=byid[r["state_id"]]
        dl,du=int(r["lower_distance_mm"]),int(r["upper_distance_mm"])
        gain=lambda layout:truths[key(k[0],k[1],layout,k[3])]["truth"]["deltaG_S"]
        r["COMSOL_lower_near_effect_at_this_upper_S"]=gain(f"P1{du}")-gain(f"P7{du}")
        r["COMSOL_upper_near_effect_at_this_lower_S"]=gain(f"P{dl}1")-gain(f"P{dl}7")
        r["COMSOL_joint_near_effect_S"]=gain("P11")-gain("P77")
        r["COMSOL_scalar_added_only_after_prediction_gate"]=True
    write_csv(out/"source_scale_response.csv",sources)
    summaries={}
    for kr in (0,8):
        for n in (32,64):
            for tol in (1,2,5):
                rr=[r for r in validation_rows if r["reference_K"]==kr and r["target_N"]==n and r["tolerance"]==tol/100]
                summaries[f"M{kr}_to_M{n}_{tol}pct"]={"n":len(rr),"accepted":sum(r["prospective_accepted"] for r in rr),"actual_pass":sum(r["COMSOL_pass"] for r in rr),
                    "false_accepts":sum(r["validation_class"]=="false_accept_full_COMSOL" for r in rr),"conservative":sum(r["validation_class"]=="conservative_unconfirmed_but_passed" for r in rr),
                    "finite_model_contained":sum(r["finite_model_contained"] for r in rr)}
    summary={"validated_local":stamp(),"new_target_states":len(saved),"new_geometries":16,"independent_source_lengths_mm":[40,64],
        "criteria_scope":"same-space specified finite-modal models; full target is external validation",
        "gain_errors":{f"M{o}":{"median_relative_error":median(r["gain_relative_error"] for r in comparisons if r["K"]==o),
                                 "max_relative_error":max(r["gain_relative_error"] for r in comparisons if r["K"]==o)} for o in (0,8,32,64)},
        "classification":summaries,"mesh_evidence_supplied":bool(mesh),"contrast_rows":len(contrast_rows),
        "all_finite_contrast_intervals_contain":all(r["finite_model_contained"] for r in contrast_rows)}
    write_json(out/"applicability_validation_summary.json",summary)
    print(json.dumps(summary,ensure_ascii=False,indent=2))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action",choices=["predict","validate"])
    parser.add_argument("--out",type=Path,default=DEFAULT_OUT)
    args=parser.parse_args()
    (prospective if args.action=="predict" else validation)(args.out.resolve())


if __name__=="__main__":
    main()
