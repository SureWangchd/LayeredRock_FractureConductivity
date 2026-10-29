"""Summarize existing round-4 crossover results without running any solver.

Reads only runs/round4_crossover. Earlier results are used only through their
round-4 copied SUCCESS/reuse records; original files are not re-audited.

Outputs conductivity_response.csv, crossover_summary.csv, and
prediction_vs_comsol.csv. Missing states or predictions are never filled.
All paired D quantities mean deltaG_B - deltaG_A, not total-G differences.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Any


# JSON and MATLAB CSV serialize the same binary sigma with different final
# digits. Match numerically without rounding values or changing saved results.
SIGMA_MATCH_RTOL = 1e-13
SIGMA_MATCH_ATOL = 1e-14


def num(value: Any) -> float | None:
    if value is None or value == "" or isinstance(value, bool):
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def truthy(value: Any) -> bool:
    return str(value).lower() in {"true", "1", "1.0"}


def sigma_key(value: Any) -> str | None:
    value = num(value)
    return format(value, ".12g") if value is not None else None


def sign(value: Any, floor: float) -> str:
    value = num(value)
    if value is None:
        return "unavailable"
    if abs(value) <= floor:
        return "numerical_tie"
    return "B_greater" if value > 0 else "A_greater"


def read_json(path: Path, warnings: list[str], optional: bool = False) -> dict:
    if optional and not path.is_file():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError) as exc:
        warnings.append(f"Unreadable current record {path}: {exc}")
        return {}
    if not isinstance(value, dict):
        warnings.append(f"Expected JSON object: {path}")
        return {}
    return value


def read_csv(path: Path, warnings: list[str]) -> list[dict]:
    try:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))
    except (OSError, csv.Error) as exc:
        warnings.append(f"Unreadable current table {path}: {exc}")
        return []


def scalar_fields(data: dict, prefix: str = "") -> dict:
    return {prefix+k: v for k, v in data.items() if isinstance(v, (str, int, float, bool)) or v is None}


def relative(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def load_refinements(root: Path, warnings: list[str]) -> dict[tuple, dict]:
    candidates: dict[tuple, list[dict]] = defaultdict(list)
    for path in sorted((root / "precision").glob("**/refinement_comparison.csv")):
        for row in read_csv(path, warnings):
            ident = (row.get("case_id"), num(row.get("sigma_f_S_m")))
            if not all(ident) or num(row.get("fine_deltaG_S")) is None:
                warnings.append(f"Incomplete refinement row: {path}")
                continue
            candidates[ident].append({**row, "source_file": relative(path, root)})
    selected = {}
    for ident, rows in candidates.items():
        # If more than one finer mesh is saved, explicitly use the finest.
        row = max(rows, key=lambda r: num(r.get("fine_triangles")) or -1)
        selected[ident] = {**row, "available_refinement_rows": len(rows)}
    return selected


def match_refinement(refinements: dict, case_id: str, sf: float) -> dict | None:
    matches = [(saved_sf, row) for (saved_case, saved_sf), row in refinements.items()
               if saved_case == case_id and math.isclose(saved_sf, sf,
                   rel_tol=SIGMA_MATCH_RTOL, abs_tol=SIGMA_MATCH_ATOL)]
    if len(matches) > 1:
        raise ValueError(f"Ambiguous numeric refinement match for {case_id}, sigma={sf!r}; "
                         f"candidate sigmas={[item[0] for item in matches]}")
    if not matches:
        return None
    saved_sf, row = matches[0]
    return {**row, "sigma_match_absolute_difference_S_m": abs(saved_sf-sf),
            "sigma_match_relative_tolerance": SIGMA_MATCH_RTOL,
            "sigma_match_absolute_tolerance_S_m": SIGMA_MATCH_ATOL}


def load_states(root: Path, refinements: dict, warnings: list[str]) -> list[dict]:
    states = []
    seen = set()
    paths = sorted((root / "original").glob("*/sf_*/SUCCESS.json"))
    paths += sorted((root / "periphery").glob("*/sf_*/SUCCESS.json"))
    for path in paths:
        result = read_json(path, warnings)
        if not result:
            continue
        if str(result.get("status", "success")).lower() != "success":
            continue
        required = ("G0_S", "G_S", "deltaG_S", "sigma_f_S_m")
        if any(num(result.get(k)) is None for k in required):
            warnings.append(f"Incomplete SUCCESS: {path}")
            continue
        state_folder, background_folder = path.parent, path.parent.parent
        params = read_json(background_folder / "parameters.json", warnings, optional=True)
        if not params:
            params = read_json(state_folder / "parameters.json", warnings, optional=True)
        prediction = read_json(state_folder / "prediction.json", warnings, optional=True)
        original = path.relative_to(root).parts[0] == "original"
        role = prediction.get("role") or result.get("role") or params.get("role") or params.get("configuration_role")
        layout = prediction.get("layout") or result.get("layout") or params.get("layout")
        if original:
            role = role or background_folder.name
            layout = layout or "original"
        if role not in {"A", "B"} or not layout:
            warnings.append(f"Missing explicit role/layout; no filename guessing for new geometry: {path}")
            continue
        case_id = result.get("case_id") or params.get("id")
        if not case_id:
            warnings.append(f"Missing case_id: {path}")
            continue
        sf = num(result["sigma_f_S_m"])
        ident = (layout, role, sigma_key(sf))
        if ident in seen:
            raise ValueError(f"Duplicate layout/role/sigma state {ident}; distinguish the geometry explicitly.")
        seen.add(ident)
        if prediction and sigma_key(prediction.get("sigma_f_S_m")) != sigma_key(sf):
            warnings.append(f"Prediction sigma does not match SUCCESS, prediction ignored: {state_folder}")
            prediction = {}
        phase = result.get("round4_phase") or result.get("phase") or prediction.get("phase") or "unrecorded"
        marginal = read_json(state_folder / "marginal_scalar.json", warnings, optional=True)
        row = {
            **scalar_fields(result),
            "case_id": case_id, "layout": layout, "role": role,
            "sigma_f_S_m": sf, "G0_S": num(result["G0_S"]),
            "G_S": num(result["G_S"]), "deltaG_S": num(result["deltaG_S"]),
            "relative_gain": num(result["deltaG_S"])/num(result["G0_S"]) if num(result["G0_S"]) else None,
            "phase": phase, "prediction_phase": prediction.get("phase"),
            "scope": "original_seen_geometry" if original else "new_peripheral_geometry",
            "truth_source_file": relative(path, root),
            "round4_reused": truthy(result.get("round4_reused", False)),
            "round4_source_folder": result.get("round4_source_folder"),
            "Sf_integral_m": num(marginal.get("Sf_m")),
            "marginal_source_file": relative(state_folder / "marginal_scalar.json", root) if marginal else None,
            "marginal_scope": "Posthoc t/DV^2 integral ellipse |E|^2; dG/dsigma_f, unit m; not a background prediction." if marginal else None,
            "_prediction": prediction,
        }
        row["relative_gain_percent"] = 100*row["relative_gain"] if row["relative_gain"] is not None else None
        local = result.get("local", {})
        if isinstance(local, dict):
            endcap = local.get("endcap", {})
            if isinstance(endcap, dict):
                row["finite_endcap_dV_V"] = num(endcap.get("dV_V"))
            for k in ("trace_difference_relative_to_absolute_exchange", "net_to_absolute_exchange"):
                row["local_"+k] = num(local.get(k))
        fine = match_refinement(refinements, case_id, sf)
        if fine:
            row.update(
                fine_deltaG_S=num(fine["fine_deltaG_S"]),
                fine_G0_S=num(fine.get("fine_G0_S")), fine_G_S=num(fine.get("fine_G_S")),
                fine_minus_current_deltaG_S=num(fine["fine_deltaG_S"])-row["deltaG_S"],
                fine_source_file=fine["source_file"], fine_triangles=num(fine.get("fine_triangles")),
                fine_available_refinement_rows=fine["available_refinement_rows"],
                fine_sigma_match_absolute_difference_S_m=fine["sigma_match_absolute_difference_S_m"],
                fine_sigma_match_relative_tolerance=fine["sigma_match_relative_tolerance"],
                fine_sigma_match_absolute_tolerance_S_m=fine["sigma_match_absolute_tolerance_S_m"],
            )
        states.append(row)
    states.sort(key=lambda r: (r["layout"], r["role"], r["sigma_f_S_m"]))
    previous = {}
    for row in states:
        ident = (row["layout"], row["role"])
        prior = previous.get(ident)
        if prior is not None:
            ds = row["sigma_f_S_m"]-prior["sigma_f_S_m"]
            row.update(previous_sigma_f_S_m=prior["sigma_f_S_m"],
                       adjacent_secant_dDeltaG_dsigma_m=(row["deltaG_S"]-prior["deltaG_S"])/ds)
        previous[ident] = row
    return states


def pair_states(states: list[dict], tie_floor: float, mesh_factor: float) -> list[dict]:
    indexed: dict[tuple, dict] = defaultdict(dict)
    for row in states:
        indexed[(row["layout"], sigma_key(row["sigma_f_S_m"]))][row["role"]] = row
    pairs = []
    for (layout, _), roles in indexed.items():
        if set(roles) != {"A", "B"}:
            continue
        a, b = roles["A"], roles["B"]
        difference = b["deltaG_S"]-a["deltaG_S"]
        pair = {
            "layout": layout, "sigma_f_S_m": a["sigma_f_S_m"],
            "case_A": a["case_id"], "case_B": b["case_id"],
            "G0_A_S": a["G0_S"], "G0_B_S": b["G0_S"],
            "deltaG_A_S": a["deltaG_S"], "deltaG_B_S": b["deltaG_S"],
            "D_B_minus_A_S": difference, "D_point_sign": sign(difference,tie_floor),
            "Sf_A_m": a.get("Sf_integral_m"), "Sf_B_m": b.get("Sf_integral_m"),
            "phase_A": a["phase"], "phase_B": b["phase"],
            "source_A": a["truth_source_file"], "source_B": b["truth_source_file"],
            "mesh_rank_status": "not_checked",
            "mesh_near_tie_factor": mesh_factor,
            "numerical_tie_threshold_S": tie_floor,
        }
        if num(pair["Sf_A_m"]) is not None and num(pair["Sf_B_m"]) is not None:
            pair["dD_dsigma_integral_B_minus_A_m"] = pair["Sf_B_m"]-pair["Sf_A_m"]
        af, bf = num(a.get("fine_deltaG_S")), num(b.get("fine_deltaG_S"))
        if af is not None and bf is not None:
            fine_d = bf-af
            observed = abs(af-a["deltaG_S"])+abs(bf-b["deltaG_S"])
            same = difference*fine_d > 0
            near = min(abs(difference),abs(fine_d)) <= max(tie_floor,mesh_factor*observed)
            status = "near_tie_at_observed_mesh_change" if near else "same_sign_on_finer_mesh" if same else "sign_changed_on_finer_mesh"
            pair.update(
                fine_deltaG_A_S=af, fine_deltaG_B_S=bf, fine_D_B_minus_A_S=fine_d,
                fine_D_sign=sign(fine_d,tie_floor), fine_D_minus_current_D_S=fine_d-difference,
                observed_pair_mesh_change_scale_S=observed,
                mesh_rank_status=status,
                coarse_fine_sign_agree=same,
                mesh_change_scope="Sum of observed absolute A/B increment changes, multiplied by an explicit diagnostic factor for near-tie labeling; not a rigorous FEM error bound.",
                fine_source_A=a.get("fine_source_file"),fine_source_B=b.get("fine_source_file"),
            )
        elif af is not None or bf is not None:
            pair["mesh_rank_status"] = "one_role_refined_only"
        pairs.append(pair)
    pairs.sort(key=lambda r:(r["layout"],r["sigma_f_S_m"]))
    previous = {}
    for pair in pairs:
        prior = previous.get(pair["layout"])
        if prior is not None:
            pair["previous_paired_sigma_f_S_m"] = prior["sigma_f_S_m"]
            pair["adjacent_secant_dD_dsigma_m"] = (pair["D_B_minus_A_S"]-prior["D_B_minus_A_S"])/(pair["sigma_f_S_m"]-prior["sigma_f_S_m"])
        previous[pair["layout"]] = pair
    lookup = {(r["layout"],sigma_key(r["sigma_f_S_m"])):r for r in pairs}
    for row in states:
        pair = lookup.get((row["layout"],sigma_key(row["sigma_f_S_m"])))
        if pair:
            row.update({k:v for k,v in pair.items() if k not in {"layout","sigma_f_S_m"}})
            row["paired_values_repeated_for_A_and_B"] = True
        else:
            row["pair_status"] = "counterpart_not_yet_available"
    return pairs


def crossover_rows(pairs: list[dict], tie_floor: float) -> list[dict]:
    layouts: dict[str,list[dict]] = defaultdict(list)
    for row in pairs:
        layouts[row["layout"]].append(row)
    output = []
    for layout, part in sorted(layouts.items()):
        brackets = []
        for lo, hi in zip(part,part[1:]):
            if lo["D_B_minus_A_S"]*hi["D_B_minus_A_S"] >= 0:
                continue
            sl, sh = lo["sigma_f_S_m"], hi["sigma_f_S_m"]
            row = {
                "layout":layout,"entry_type":"observed_sign_change_bracket",
                "crossing_index":len(brackets)+1,
                "sigma_lo_S_m":sl,"sigma_hi_S_m":sh,
                "bracket_width_S_m":sh-sl,"bracket_width_over_midpoint":(sh-sl)/((sh+sl)/2),
                "width_within_one_percent":(sh-sl)/((sh+sl)/2)<=.01,
                "D_lo_S":lo["D_B_minus_A_S"],"D_hi_S":hi["D_B_minus_A_S"],
                "sign_lo":lo["D_point_sign"],"sign_hi":hi["D_point_sign"],
                "deltaG_A_lo_S":lo["deltaG_A_S"],"deltaG_B_lo_S":lo["deltaG_B_S"],
                "deltaG_A_hi_S":hi["deltaG_A_S"],"deltaG_B_hi_S":hi["deltaG_B_S"],
                "fine_D_lo_S":lo.get("fine_D_B_minus_A_S"),"fine_D_hi_S":hi.get("fine_D_B_minus_A_S"),
                "fine_sign_lo":lo.get("fine_D_sign"),"fine_sign_hi":hi.get("fine_D_sign"),
                "mesh_status_lo":lo["mesh_rank_status"],"mesh_status_hi":hi["mesh_rank_status"],
                "observed_mesh_change_lo_S":lo.get("observed_pair_mesh_change_scale_S"),
                "observed_mesh_change_hi_S":hi.get("observed_pair_mesh_change_scale_S"),
                "mesh_near_tie_factor":lo["mesh_near_tie_factor"],
                "phase_A_lo":lo["phase_A"],"phase_B_lo":lo["phase_B"],
                "phase_A_hi":hi["phase_A"],"phase_B_hi":hi["phase_B"],
                "source_A_lo":lo["source_A"],"source_B_lo":lo["source_B"],
                "source_A_hi":hi["source_A"],"source_B_hi":hi["source_B"],
                "fine_source_A_lo":lo.get("fine_source_A"),"fine_source_B_lo":lo.get("fine_source_B"),
                "fine_source_A_hi":hi.get("fine_source_A"),"fine_source_B_hi":hi.get("fine_source_B"),
                "direct_COMSOL_confirmation":"two directly solved signs; no interpolated zero asserted",
                "uniqueness_claim":False,
            }
            fine_lo, fine_hi = num(row["fine_D_lo_S"]), num(row["fine_D_hi_S"])
            if fine_lo is None or fine_hi is None:
                row["bracket_mesh_status"]="awaiting_both_roles_at_both_endpoints"
            elif "near_tie_at_observed_mesh_change" in {row["mesh_status_lo"],row["mesh_status_hi"]}:
                row["bracket_mesh_status"]="near_tie_requires_wider_or_better_resolved_bracket"
            elif fine_lo*fine_hi<0 and row["mesh_status_lo"]==row["mesh_status_hi"]=="same_sign_on_finer_mesh":
                row["bracket_mesh_status"]="sign_bracket_retained_on_finer_mesh"
            else:
                row["bracket_mesh_status"]="current_bracket_not_confirmed_on_finer_mesh"
            # Preserve the raw adjacent sign change, but report a wider bracket
            # when a nearby endpoint cannot resolve ranking at the observed
            # coarse/fine change scale. Never synthesize another target value.
            eligible_lo=[p for p in part if p["sigma_f_S_m"]<=sl and
                         p["mesh_rank_status"]=="same_sign_on_finer_mesh" and
                         p["D_point_sign"]==lo["D_point_sign"]]
            eligible_hi=[p for p in part if p["sigma_f_S_m"]>=sh and
                         p["mesh_rank_status"]=="same_sign_on_finer_mesh" and
                         p["D_point_sign"]==hi["D_point_sign"]]
            if eligible_lo and eligible_hi:
                rl=max(eligible_lo,key=lambda p:p["sigma_f_S_m"])
                rh=min(eligible_hi,key=lambda p:p["sigma_f_S_m"])
                rsl,rsh=rl["sigma_f_S_m"],rh["sigma_f_S_m"]
                expanded=rsl<sl or rsh>sh
                row.update(
                    report_sigma_lo_S_m=rsl,report_sigma_hi_S_m=rsh,
                    report_bracket_width_S_m=rsh-rsl,
                    report_bracket_width_over_midpoint=(rsh-rsl)/((rsh+rsl)/2),
                    report_width_within_one_percent=(rsh-rsl)/((rsh+rsl)/2)<=.01,
                    report_D_lo_S=rl["D_B_minus_A_S"],report_D_hi_S=rh["D_B_minus_A_S"],
                    report_fine_D_lo_S=rl["fine_D_B_minus_A_S"],report_fine_D_hi_S=rh["fine_D_B_minus_A_S"],
                    report_observed_mesh_change_lo_S=rl["observed_pair_mesh_change_scale_S"],
                    report_observed_mesh_change_hi_S=rh["observed_pair_mesh_change_scale_S"],
                    report_fine_source_A_lo=rl.get("fine_source_A"),report_fine_source_B_lo=rl.get("fine_source_B"),
                    report_fine_source_A_hi=rh.get("fine_source_A"),report_fine_source_B_hi=rh.get("fine_source_B"),
                    report_bracket_expanded=expanded,
                    report_bracket_mesh_status="expanded_to_resolved_outer_fine_endpoints" if expanded else "sign_bracket_retained_on_finer_mesh",
                    report_bracket_scope="Nearest enclosing directly solved endpoints with same coarse/fine signs and min absolute D above the explicit mesh-change diagnostic. This is not a rigorous error enclosure or uniqueness proof.",
                )
            else:
                row["report_bracket_mesh_status"]="awaiting_resolved_enclosing_fine_endpoints"
            brackets.append(row)
        for row in brackets:
            row.update(sampled_sigma_min_S_m=part[0]["sigma_f_S_m"],sampled_sigma_max_S_m=part[-1]["sigma_f_S_m"],paired_states_in_layout=len(part),observed_sign_change_brackets=len(brackets))
        output.extend(brackets)
        ties = [r for r in part if abs(r["D_B_minus_A_S"])<=tie_floor]
        for r in ties:
            output.append({"layout":layout,"entry_type":"sampled_numerical_tie_not_exact_root","sigma_lo_S_m":r["sigma_f_S_m"],"sigma_hi_S_m":r["sigma_f_S_m"],"D_lo_S":r["D_B_minus_A_S"],"fine_D_lo_S":r.get("fine_D_B_minus_A_S"),"mesh_status_lo":r["mesh_rank_status"],"numerical_tie_threshold_S":tie_floor,"uniqueness_claim":False})
        if not brackets and not ties:
            output.append({"layout":layout,"entry_type":"no_crossing_observed_in_sampled_range","sampled_sigma_min_S_m":part[0]["sigma_f_S_m"],"sampled_sigma_max_S_m":part[-1]["sigma_f_S_m"],"paired_states_in_layout":len(part),"observed_sign_change_brackets":0,"signs_observed":";".join(sorted({r["D_point_sign"] for r in part})),"scope":"No crossing observed at sampled points; does not prove absence between points or outside this range.","uniqueness_claim":False})
    return output


def prediction_bounds(pred: dict, tie_floor: float) -> tuple[float | None,float | None,str]:
    upper=num(pred.get("predicted_deltaG_S")); lower=num(pred.get("certified_lower_deltaG_S"))
    bound=num(pred.get("residual_error_bound_S"))
    if upper is None or lower is None or bound is None:
        return lower,upper,"unavailable"
    if lower>upper+tie_floor or bound<0 or truthy(pred.get("upper_bound_violation")):
        return lower,upper,"recorded_bounds_invalid"
    return lower,upper,"recorded_discrete_bounds_conditional_on_frozen_method_assumptions"


def comparison_rows(states: list[dict], pairs: list[dict], tie_floor: float) -> list[dict]:
    by_pair: dict[tuple,dict] = defaultdict(dict)
    for row in states:
        by_pair[(row["layout"],sigma_key(row["sigma_f_S_m"]))][row["role"]]=row
    pair_lookup={(r["layout"],sigma_key(r["sigma_f_S_m"])):r for r in pairs}
    output=[]
    for state in states:
        pred=state["_prediction"]; estimate=num(pred.get("predicted_deltaG_S"))
        if estimate is None:
            continue
        actual=state["deltaG_S"]; error=estimate-actual; g0=state["G0_S"]
        lower,upper,bstatus=prediction_bounds(pred,tie_floor)
        phase=str(pred.get("phase") or state["phase"])
        if state["layout"]=="original":
            evaluation="original_seen_geometry_existing_target_reference" if state["round4_reused"] else "original_seen_geometry_new_sigma_reference"
        elif "adaptive" in phase or "location" in phase:
            evaluation="adaptive_location_after_initial_targets; not initial_preselected_batch"
        else:
            evaluation="new_geometry_preselected_batch_as_recorded; timing_not_audited_by_summary"
        row={
            **scalar_fields(pred),"case_id":state["case_id"],"layout":state["layout"],"role":state["role"],
            "sigma_f_S_m":state["sigma_f_S_m"],"phase":phase,"truth_phase":state["phase"],
            "evaluation_scope":evaluation,"round4_reused":state["round4_reused"],
            "truth_G0_S":g0,"truth_G_S":state["G_S"],"truth_deltaG_S":actual,
            "predicted_deltaG_S":estimate,"increment_error_S":error,
            "increment_absolute_error_S":abs(error),
            "increment_relative_error":abs(error)/abs(actual) if abs(actual)>tie_floor else None,
            "increment_error_over_background":abs(error)/abs(g0) if g0 else None,
            "increment_near_zero":abs(actual)<=tie_floor,"near_zero_threshold_S":tie_floor,
            "predicted_lower_deltaG_S":lower,"predicted_upper_deltaG_S":upper,"bound_status":bstatus,
            "truth_deltaG_within_recorded_interval":lower-tie_floor<=actual<=upper+tie_floor if lower is not None and upper is not None and bstatus.startswith("recorded_discrete") else None,
            "truth_source_file":state["truth_source_file"],
            "prediction_source_file":str(Path(state["truth_source_file"]).parent/"prediction.json").replace("\\","/"),
            "truth_solve_seconds":state.get("solve_seconds"),
        }
        ident=(state["layout"],sigma_key(state["sigma_f_S_m"]))
        pair=pair_lookup.get(ident);roles=by_pair[ident]
        if pair:
            row.update({k:v for k,v in pair.items() if k not in {"layout","sigma_f_S_m"}})
            pa,pb=roles["A"]["_prediction"],roles["B"]["_prediction"]
            ea,eb=num(pa.get("predicted_deltaG_S")),num(pb.get("predicted_deltaG_S"))
            la,ua,sa=prediction_bounds(pa,tie_floor);lb,ub,sb=prediction_bounds(pb,tie_floor)
            row["pair_prediction_phase_A"]=pa.get("phase")
            row["pair_prediction_phase_B"]=pb.get("phase")
            row["paired_values_repeated_for_A_and_B"]=True
            if ea is not None and eb is not None:
                pd=eb-ea; row["predicted_D_B_minus_A_S"]=pd
                row["predicted_D_point_sign"]=sign(pd,tie_floor)
                row["D_prediction_error_S"]=pd-pair["D_B_minus_A_S"]
                row["point_prediction_sign_matches_COMSOL"]=sign(pd,tie_floor)==pair["D_point_sign"] if pair["D_point_sign"]!="numerical_tie" else None
                row["point_sign_comparison_scope"]="Current-mesh COMSOL sign; mesh uncertainty and interval ranking are separate columns."
            if sa.startswith("recorded_discrete") and sb.startswith("recorded_discrete"):
                dlo,dhi=lb-ua,ub-la
                interval_sign="B_greater" if dlo>tie_floor else "A_greater" if dhi<-tie_floor else "overlap_zero_cannot_rank"
                row.update(
                    predicted_D_lower_S=dlo,predicted_D_upper_S=dhi,
                    predicted_D_interval_ranking=interval_sign,
                    COMSOL_D_inside_prediction_interval=dlo-tie_floor<=pair["D_B_minus_A_S"]<=dhi+tie_floor,
                    interval_rank_matches_COMSOL=(interval_sign==pair["D_point_sign"]) if interval_sign!="overlap_zero_cannot_rank" and pair["D_point_sign"]!="numerical_tie" else None,
                    pair_bound_scope="[lower_B-upper_A,upper_B-lower_A], conditional on frozen discrete-error assumptions; no mesh or rock-model uncertainty included.",
                )
            else:
                row["predicted_D_interval_ranking"]="both_valid_recorded_bounds_required"
        else:
            row["pair_status"]="counterpart_not_yet_available"
        output.append(row)
    return output


def write_csv(path: Path, rows: list[dict], preferred: list[str]) -> int:
    if not rows:
        return 0
    clean=[{k:v for k,v in row.items() if not k.startswith("_")} for row in rows]
    fields=set().union(*(row.keys() for row in clean))
    names=[k for k in preferred if k in fields]+sorted(fields-set(preferred))
    temporary=path.with_suffix(".csv.tmp")
    try:
        with temporary.open("w",encoding="utf-8-sig",newline="") as handle:
            writer=csv.DictWriter(handle,fieldnames=names);writer.writeheader();writer.writerows(clean)
        os.replace(temporary,path)
    finally:
        if temporary.exists(): temporary.unlink()
    return len(clean)


def summarize(root: Path,tie_floor: float=1e-12,mesh_factor: float=3.) -> dict:
    warnings=[]
    fine=load_refinements(root,warnings)
    states=load_states(root,fine,warnings)
    pairs=pair_states(states,tie_floor,mesh_factor)
    cross=crossover_rows(pairs,tie_floor)
    predictions=comparison_rows(states,pairs,tie_floor)
    preferred=["layout","role","case_id","sigma_f_S_m","phase","G0_S","G_S","deltaG_S","relative_gain","D_B_minus_A_S","D_point_sign","fine_D_B_minus_A_S","mesh_rank_status"]
    outputs={"conductivity_response.csv":states,"crossover_summary.csv":cross,"prediction_vs_comsol.csv":predictions}
    counts={}
    for filename,rows in outputs.items():
        counts[filename]=write_csv(root/filename,rows,preferred)
        if not rows and (root/filename).exists():
            warnings.append(f"No current rows for {filename}; existing derived file was left untouched.")
    return {
        "root":str(root),"written_rows":counts,
        "successful_states":len(states),"complete_AB_pairs":len(pairs),
        "state_predictions_awaiting_truth":sum(1 for branch in ("original","periphery") for p in (root/branch).glob("*/sf_*/prediction.json") if not (p.parent/"SUCCESS.json").exists()),
        "refined_case_sigma_states":len(fine),
        "sign_change_brackets":[{k:r.get(k) for k in ("layout","sigma_lo_S_m","sigma_hi_S_m","D_lo_S","D_hi_S","bracket_width_over_midpoint","bracket_mesh_status","report_sigma_lo_S_m","report_sigma_hi_S_m","report_bracket_width_over_midpoint","report_bracket_mesh_status")} for r in cross if r.get("entry_type")=="observed_sign_change_bracket"],
        "numerical_tie_threshold_S":tie_floor,"mesh_near_tie_factor":mesh_factor,
        "mesh_scope":"Observed coarse/fine increment changes only; not a rigorous FEM error bound.",
        "warnings":warnings,
    }


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1]/"runs"/"round4_crossover")
    parser.add_argument("--numerical-tie-S",type=float,default=1e-12,dest="tie_floor")
    parser.add_argument("--mesh-near-tie-factor",type=float,default=3.,dest="mesh_factor",help="Near-tie diagnostic: min(|D_coarse|,|D_fine|) <= factor*(|deltaA_fine-deltaA|+|deltaB_fine-deltaB|). Not an error certificate.")
    args=parser.parse_args()
    if not args.root.is_dir():parser.error(f"Missing result root: {args.root}")
    if not math.isfinite(args.tie_floor) or args.tie_floor<0:parser.error("Tie threshold must be nonnegative and finite")
    if not math.isfinite(args.mesh_factor) or args.mesh_factor<=0:parser.error("Mesh near-tie factor must be positive and finite")
    print(json.dumps(summarize(args.root.resolve(),args.tie_floor,args.mesh_factor),ensure_ascii=True,indent=2))


if __name__=="__main__":
    main()
