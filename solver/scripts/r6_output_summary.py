# -*- coding: utf-8 -*-
"""Summarize frozen R6 criteria against a completed 32-state validation batch.

Reads completed analysis tables, never changes the frozen criterion, launches
a solver, fits a parameter, or constructs a summary from a partial batch.
The primary comparison is M0/M8 -> M64; M32 remains a separate diagnostic.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from statistics import median

import r6_criteria as io


def number(v):
    return None if v in (None, "", "None") else float(v)


def boolean(v):
    if v in (True, "True", "true", "1", 1):
        return True
    if v in (False, "False", "false", "0", 0):
        return False
    if v in (None, "", "None"):
        return None
    raise ValueError("Invalid Boolean value: "+repr(v))


def count(rows, field, expected=True):
    return sum(boolean(r[field]) is expected for r in rows)


def contrast_id(r):
    return (f"L{int(float(r['full_length_mm']))}_{r['role']}_{r['layout']}_"
            f"sf{float(r['sigma_f_S_m']):.10g}__{r['contrast']}")


def contrast_details(r):
    value=number(r["COMSOL_contrast_S"])
    point=number(r["reference_contrast_S"])
    low,high=number(r["contrast_lower_S"]),number(r["contrast_upper_S"])
    return {"contrast_id":contrast_id(r),"reference_K":int(r["reference_K"]),
        "full_length_mm":int(float(r["full_length_mm"])),"role":r["role"],"layout":r["layout"],
        "sigma_f_S_m":number(r["sigma_f_S_m"]),"contrast":r["contrast"],"category":r["category"],
        "COMSOL_contrast_S":value,"reference_contrast_S":point,"absolute_error_S":abs(point-value),
        "interval_lower_S":low,"interval_upper_S":high,"interval_width_S":high-low,
        "COMSOL_outside_interval_distance_S":max(low-value,value-high,0.),
        "reference_relative_amplitude_error":number(r["reference_relative_amplitude_error"]),
        "sign_decision":r["sign_decision"],"mesh_resolved":boolean(r["resolved_against_observed_mesh_scale"]),
        "observed_mesh_resolution_scale_S":number(r["observed_mesh_resolution_scale_S"]),
        "contrast_terms":r["contrast_terms"]}


def summarize_gains(validation, comparisons):
    groups,failures=[],[]
    indexed={(r["state_id"],int(r["K"])):r for r in comparisons}
    for kr in (0,8):
        for tolerance in (.01,.02,.05):
            rows=[r for r in validation if int(r["reference_K"])==kr and int(r["target_N"])==64 and math.isclose(float(r["tolerance"]),tolerance)]
            if len(rows)!=32:
                raise ValueError("Incomplete gain validation group")
            accepted=[r for r in rows if boolean(r["prospective_accepted"])]
            unknown=[r for r in rows if not boolean(r["prospective_accepted"])]
            false=[r for r in accepted if not boolean(r["COMSOL_pass"])]
            conservative=[r for r in unknown if boolean(r["COMSOL_pass"])]
            errors=[float(r["COMSOL_relative_gain_error"]) for r in rows]
            groups.append({"reference_K":kr,"target_N":64,"output_group":"total_gain","tolerance":tolerance,
                "n":32,"prospectively_usable":len(accepted),"prospectively_unknown":len(unknown),
                "actual_within_tolerance":count(rows,"COMSOL_pass"),"actual_exceeds_tolerance":count(rows,"COMSOL_pass",False),
                "false_accept_against_full_COMSOL":len(false),"conservative_unknown_but_actual_pass":len(conservative),
                "finite_model_interval_contains":count(rows,"finite_model_contained"),
                "full_COMSOL_interval_contains_descriptive":count(rows,"full_COMSOL_contained_descriptive_only"),
                "median_relative_error":median(errors),"maximum_relative_error":max(errors),
                "unknown_state_ids":[r["state_id"] for r in unknown],
                "false_accept_state_ids":[r["state_id"] for r in false],
                "conservative_state_ids":[r["state_id"] for r in conservative]})
            for r in rows:
                if boolean(r["COMSOL_pass"]):
                    continue
                c=indexed[r["state_id"],kr]
                failures.append({"output_group":"total_gain","reference_K":kr,"target_N":64,"tolerance":tolerance,
                    "state_id":r["state_id"],"classification":r["validation_class"],
                    "reference_gain_S":float(c["gain_prediction_S"]),"COMSOL_gain_S":float(c["gain_COMSOL_S"]),
                    "absolute_error_S":float(c["gain_absolute_error_S"]),"relative_error":float(c["gain_relative_error"]),
                    "interval_lower_S":float(r["gain_lower_S"]),"interval_upper_S":float(r["gain_upper_S"]),
                    "M64_minus_M32_gain_S":float(r["M64_minus_M32_gain_S"]),
                    "target_mesh_gain_change_S":number(r["target_mesh_gain_change_S"]),
                    "reference_mesh_gain_change_S":number(r["reference_mesh_gain_change_S"])})
    return groups,failures


def contrast_groups(contrasts):
    predicates={
        "single_side_all48":lambda r:r["category"]=="single_side",
        "single_side_controlled_edges32":lambda r:r["category"]=="single_side" and not r["contrast"].endswith("average"),
        "single_side_averages16":lambda r:r["category"]=="single_side" and r["contrast"].endswith("average"),
        "AB_ranks16":lambda r:r["category"]=="AB_rank",
        "joint_net8":lambda r:r["category"]=="joint_net",
        "interaction8":lambda r:r["category"]=="interaction",
    }
    expected={"single_side_all48":48,"single_side_controlled_edges32":32,"single_side_averages16":16,
              "AB_ranks16":16,"joint_net8":8,"interaction8":8}
    groups,details=[],{}
    for kr in (0,8):
        for name,predicate in predicates.items():
            rows=[r for r in contrasts if int(r["reference_K"])==kr and int(r["target_N"])==64 and predicate(r)]
            if len(rows)!=expected[name]:
                raise ValueError("Incomplete contrast group: "+name)
            confirmed=[r for r in rows if r["sign_decision"]!="cannot_determine"]
            unresolved=[r for r in rows if boolean(r["resolved_against_observed_mesh_scale"]) is False]
            resolved=[r for r in rows if boolean(r["resolved_against_observed_mesh_scale"]) is True]
            wrong=[r for r in confirmed if boolean(r["sign_decision_correct_against_base_COMSOL"]) is False]
            amplitude=name.startswith("single_side")
            amp_confirmed=[r for r in rows if boolean(r["amplitude_10pct_confirmed"])] if amplitude else []
            amp_false=[r for r in amp_confirmed if boolean(r["amplitude_10pct_actual_base_COMSOL"]) is False]
            amp_failed=[r for r in rows if boolean(r["amplitude_10pct_actual_base_COMSOL"]) is False] if amplitude else []
            entry={"reference_K":kr,"target_N":64,"output_group":name,"n":len(rows),
                "interval_excludes_zero":len(confirmed),"interval_includes_zero":len(rows)-len(confirmed),
                "confirmed_sign_matches_base_COMSOL":count(confirmed,"sign_decision_correct_against_base_COMSOL"),
                "confirmed_sign_wrong_base_COMSOL":len(wrong),
                "mesh_resolved_total":len(resolved),"mesh_unresolved_total":len(unresolved),
                "mesh_resolution_unavailable":len(rows)-len(resolved)-len(unresolved),
                "confirmed_sign_matches_and_mesh_resolved":sum(boolean(r["sign_decision_correct_against_base_COMSOL"]) is True and boolean(r["resolved_against_observed_mesh_scale"]) is True for r in confirmed),
                "confirmed_sign_wrong_and_mesh_resolved":sum(boolean(r["resolved_against_observed_mesh_scale"]) is True for r in wrong),
                "confirmed_sign_but_mesh_unresolved":sum(boolean(r["resolved_against_observed_mesh_scale"]) is False for r in confirmed),
                "finite_model_interval_contains":count(rows,"finite_model_contained"),
                "full_COMSOL_interval_contains_descriptive":count(rows,"full_COMSOL_contained_descriptive_only"),
                "maximum_absolute_reference_error_S":max(float(r["reference_absolute_error_S"]) for r in rows),
                "amplitude_tolerance":.1 if amplitude else None,
                "amplitude_confirmed_10pct":len(amp_confirmed) if amplitude else None,
                "amplitude_actual_within_10pct":count(rows,"amplitude_10pct_actual_base_COMSOL") if amplitude else None,
                "amplitude_actual_exceeds_10pct":len(amp_failed) if amplitude else None,
                "amplitude_actual_exceeds_10pct_and_resolved":sum(boolean(r["resolved_against_observed_mesh_scale"]) is True for r in amp_failed) if amplitude else None,
                "amplitude_false_accept_full_COMSOL":len(amp_false) if amplitude else None,
                "amplitude_false_accept_and_resolved":sum(boolean(r["resolved_against_observed_mesh_scale"]) is True for r in amp_false) if amplitude else None,
                "median_relative_amplitude_error":median(float(r["reference_relative_amplitude_error"]) for r in rows) if amplitude else None,
                "maximum_relative_amplitude_error":max(float(r["reference_relative_amplitude_error"]) for r in rows) if amplitude else None}
            groups.append(entry)
            details[f"M{kr}_{name}"]={
                "confirmed_wrong_signs":[contrast_details(r) for r in wrong],
                "interval_includes_zero":[contrast_details(r) for r in rows if r["sign_decision"]=="cannot_determine"],
                "mesh_unresolved":[contrast_details(r) for r in unresolved],
                "amplitude_over_10pct":[contrast_details(r) for r in amp_failed],
                "amplitude_false_accept":[contrast_details(r) for r in amp_false],
                "full_COMSOL_outside_finite_interval_descriptive":[contrast_details(r) for r in rows if boolean(r["full_COMSOL_contained_descriptive_only"]) is False],
            }
    return groups,details


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=io.DEFAULT_OUT)
    args=parser.parse_args();out=args.out.resolve()
    completion=io.read_json(out/"truth_initial_complete.json")
    summary=io.read_json(out/"applicability_validation_summary.json")
    gate=io.read_json(out/"predictions_ready.json")
    if completion.get("status")!="success" or completion.get("count")!=32 or summary.get("new_target_states")!=32 or gate.get("count")!=32:
        raise RuntimeError("Require completed32 full targets, complete validation and the prior prediction gate")
    validation=io.read_csv(out/"applicability_validation.csv")
    comparisons=io.read_csv(out/"closure_new_source_comparison.csv")
    contrasts=io.read_csv(out/"side_and_rank_intervals.csv")
    if len(validation)!=384 or len(comparisons)!=128 or len(contrasts)!=320:
        raise RuntimeError("Incomplete machine tables; no summary was written")
    if len({r["state_id"] for r in comparisons})!=32:
        raise RuntimeError("Comparison table does not contain32 unique states")
    mesh_path=out/"mesh_evidence.json"
    mesh=io.read_json(mesh_path) if mesh_path.exists() else {}
    if mesh and not summary.get("mesh_evidence_supplied"):
        raise RuntimeError("Mesh evidence is newer than the validation tables; rerun r6_criteria.py validate first")
    gains,gain_failures=summarize_gains(validation,comparisons)
    contrast_stats,contrast_details=contrast_groups(contrasts)
    records=gains+contrast_stats
    result={"generated_local":io.stamp(),"full_target_count":32,"geometry_count":16,"new_lengths_mm":[40,64],
        "primary_comparison_N":64,"finite_model_error_bounds_are_not_full_COMSOL_certificates":True,
        "gain_summary":gains,"contrast_summary":contrast_stats,"gain_actual_tolerance_failures":gain_failures,
        "contrast_details":contrast_details,"mesh_evidence_available":bool(mesh),
        "global_observed_gain_resolution_S":mesh.get("global_conservative_gain_resolution_S"),
        "interpretation_notes":[
            "Unknown means not certified by the specified finite-model interval; it is not a declaration of failure.",
            "Counts for controlled edges and averages are separate; averaging is not an additional independent geometry.",
            "Sign matches against base COMSOL and observed mesh resolution are reported separately.",
            "Single-side amplitude uses the frozen10% requirement; total gain uses1%,2%,5% separately.",
            "Joint and interaction results use absolute differences and uncertainty; no small net sign is pursued.",
            "Mesh resolution is an observed conservative scale from six selected cases, not a rigorous global bound.",
            "Detailed full-COMSOL interval containment is descriptive because the proven interval concerns M64."],
        "inputs":["applicability_validation.csv","closure_new_source_comparison.csv","side_and_rank_intervals.csv","mesh_evidence.json if available"]}
    io.write_json(out/"output_specific_summary.json",result)
    io.write_csv(out/"output_specific_summary.csv",[{k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in r.items()} for r in records])
    io.write_csv(out/"output_specific_gain_failures.csv",gain_failures)
    print("Completed output-specific summary:32 states; M0/M8 -> M64; total gain,32 controlled edges,16 averages,16 AB,8 joint,8 interaction per method.")
    for r in gains:
        if r["tolerance"]==.02:
            print(f"M{r['reference_K']} gain2%: accepted={r['prospectively_usable']}, unknown={r['prospectively_unknown']}, actual_fail={r['actual_exceeds_tolerance']}, false_accept={r['false_accept_against_full_COMSOL']}, conservative={r['conservative_unknown_but_actual_pass']}")
    for r in contrast_stats:
        if r["reference_K"]==8:
            print(f"M8 {r['output_group']}: sign_confirmed={r['interval_excludes_zero']}/{r['n']}, wrong={r['confirmed_sign_wrong_base_COMSOL']}, mesh_resolved={r['mesh_resolved_total']}, amp_fail={r['amplitude_actual_exceeds_10pct']}")


if __name__=="__main__":
    main()
