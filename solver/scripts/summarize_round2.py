"""Independent, read-only audit of round-2 COMSOL cases; never starts a solver.

Writes only runs/2026-10-21_round2/{summary,comparisons}.csv,
validation_report.json and README_CN.md. Existing round-1 artifacts are untouched.
"""
from __future__ import annotations

import argparse
import csv
import io
import itertools
import json
import math
from datetime import datetime, timedelta, timezone
from pathlib import Path


LIMITS = {"background_relative": 1e-4, "area_absolute_mm2": 0.01,
          "area_relative_to_ellipse": 1e-4, "conservation_relative": 1e-6,
          "ellipse_relative": 1e-4, "endcap_absolute_mm2": 0.001,
          "endcap_relative": 1e-3, "endcap_high_fraction_min": 0.9,
          "endcap_low_fraction_max": 0.1, "mesh_relative": 1e-4,
          "effect_to_mesh_change_min": 5}
GEOMETRY_KEYS = ("L_mm", "cx_mm", "cz_mm", "a_mm", "b_mm", "theta_deg")
MATERIAL_KEYS = ("sigma_low_S_m", "sigma_high_S_m", "sigma_frac_S_m", "out_of_plane_thickness_m")
BACKGROUND_KEYS = ("background", "L_mm", "h_mm", "z_interfaces_mm", "layer_materials", "interface_amplitudes_mm", "sigma_low_S_m", "sigma_high_S_m", "out_of_plane_thickness_m")
FIXED_ENDCAP_FAMILIES = {"equal_tile_internal_rearrangement", "equal_tiles_wave_inversion_complement"}


def num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def rel(a, b):
    return 2 * abs(a-b) / max(abs(a)+abs(b), 1e-30)


def clean(x):
    if isinstance(x, dict):
        return {k: clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    return None if isinstance(x, float) and not math.isfinite(x) else x


def vector(x):
    return [] if x is None else x if isinstance(x, list) else [x]


def get(d, *keys):
    return next((d[k] for k in keys if k in d and d[k] is not None), None)


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write(path, text):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(text, encoding="utf-8")
    temp.replace(path)


def csvout(path, rows, preferred):
    keys = list(preferred)
    for row in rows:
        keys.extend(k for k in row if k not in keys)
    s = io.StringIO(newline="")
    writer = csv.DictWriter(s, fieldnames=keys)
    writer.writeheader()
    for row in rows:
        writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v for k, v in clean(row).items()})
    write(path, "\ufeff" + s.getvalue())


def sig(p, keys):
    return json.dumps({k: p.get(k) for k in keys}, sort_keys=True)


def active(p):
    return bool(p.get("fracture")) and not p.get("passive_ellipse", False)


def sigma(c):
    # Do not silently replace an unavailable reaction value by the old gradient.
    r = c["r"]
    value = get(r, "sigma_reaction_S_m", "sigma_eff_reaction_S_m")
    if num(value):
        return value
    current = get(r, "reaction_positive_A", "reaction_current_positive_A")
    if num(current):
        p = c["p"]
        return abs(current) / p.get("out_of_plane_thickness_m", 1)
    if r.get("conductivity_method") in ("reaction", "reaction_current", "conservative_reaction"):
        return r.get("sigma_eff_S_m")
    return None


def cap_data(r, side):
    # Canonical names; aliases allow the solver to finish its first-case schema.
    prefix = "capminus" if side == "minus" else "capplus"
    return {"area": get(r, f"{prefix}_area_mm2", f"endcap_{side}_area_mm2"),
            "high": get(r, f"{prefix}_high_area_mm2", f"{prefix}_high_mm2", f"endcap_{side}_high_area_mm2"),
            "low": get(r, f"{prefix}_low_area_mm2", f"{prefix}_low_mm2", f"endcap_{side}_low_area_mm2"),
            "fraction": get(r, f"{prefix}_high_fraction", f"endcap_{side}_high_fraction")}


def audit(c):
    p, r, folder = c["p"], c["r"], c["folder"]
    checks = []

    def check(name, test, value=None, limit=None, note=""):
        checks.append(dict(name=name, status="pending" if test is None else "pass" if test else "fail", value=value, limit=limit, note=note))

    sr = sigma(c)
    check("reaction_conductivity_available_positive", sr > 0 if num(sr) else None, sr)
    check("primary_sigma_uses_reaction_current", rel(sr,r.get("sigma_eff_S_m"))<=1e-10 if num(sr) and num(r.get("sigma_eff_S_m")) else None, r.get("sigma_eff_S_m"),1e-10)
    for name, keys in (("reaction_flux", ("reaction_flux_relative_error", "reaction_flux_balance_relative_error")), ("reaction_energy", ("reaction_power_relative_error", "reaction_energy_relative_error"))):
        v = get(r, *keys)
        check(name, v <= LIMITS["conservation_relative"] if num(v) else None, v, LIMITS["conservation_relative"])
    energy_sigma = get(r, "sigma_energy_S_m", "sigma_eff_energy_S_m", "power_volume_W_per_m")
    check("reaction_matches_energy", rel(sr, energy_sigma) <= 1e-6 if num(sr) and num(energy_sigma) else None, rel(sr, energy_sigma) if num(sr) and num(energy_sigma) else None, 1e-6)
    total = sum(r.get("area_"+phase+"_mm2", 0) for phase in ("low", "high", "fracture"))
    check("area_closure", abs(total-p.get("L_mm",120)**2)<=0.01, total-p.get("L_mm",120)**2, 0.01)
    ideal = math.pi*p["a_mm"]*p["b_mm"]
    if active(p) and p.get("intervention", "none") in ("none", "full"):
        af = r.get("area_fracture_mm2")
        err = abs(af-ideal) if num(af) else None
        check("ellipse_area", err<=0.01 and err/ideal<=1e-4 if num(err) else None, err, min(0.01,1e-4*ideal))
    if p.get("fracture"):
        original = r.get("original_ellipse_area_mm2")
        check("original_ellipse_area", abs(original-ideal)<=min(0.01,1e-4*ideal) if num(original) else None, original,min(0.01,1e-4*ideal))
        length = p.get("endcap_length_mm",3)
        q = 1-length/p["a_mm"]
        capideal = p["a_mm"]*p["b_mm"]*(math.acos(q)-q*math.sqrt(1-q*q)) if -1<=q<=1 else None
        for side in ("minus","plus"):
            cap = cap_data(r,side)
            err = abs(cap["area"]-capideal) if num(cap["area"]) and num(capideal) else None
            check(f"{side}_endcap_geometry", err<=0.001 and err/capideal<=1e-3 if num(err) else None, err, min(0.001,1e-3*capideal) if num(capideal) else None)
            allareas = all(num(cap[k]) for k in ("area","high","low"))
            check(f"{side}_endcap_phase_closure", abs(cap["high"]+cap["low"]-cap["area"])<=0.001 if allareas else None, cap)
            if allareas and cap["area"]>0:
                fraction = cap["high"]/cap["area"]
                check(f"{side}_endcap_high_fraction", 0<=fraction<=1+1e-8 and (not num(cap["fraction"]) or abs(cap["fraction"]-fraction)<1e-7), fraction)
    required = [f"{c['id']}.mph", "parameters.json", "selections.json", "domain_integrals.csv", "mesh_stats.json", "solve_round2_case_snapshot.m", "run_round2_snapshot.m", "fields.csv", "solver_progress.log"]
    check("artifacts_saved", all((folder/f).is_file() for f in required), [f for f in required if not (folder/f).is_file()])
    try:
        sel = read(folder/"selections.json")
        dom = vector(sel.get("all_domains"))
        material = vector(sel.get("domain_background_material"))
        low, high, frac = [set(vector(sel.get(k))) for k in ("low_domains","high_domains","fracture_domains")]
        check("material_domain_coverage", not(low&high or low&frac or high&frac) and low|high|frac==set(dom))
        with (folder/"domain_integrals.csv").open(encoding="utf-8-sig",newline="") as f:
            domains = {int(float(row["domain_id"])):row for row in csv.DictReader(f)}
        check("domain_integral_coverage", set(domains)==set(dom))
        check("domain_area_sum", abs(sum(float(row["area_mm2"]) for row in domains.values())-total)<=0.01)
        if len(material)==len(dom):
            phase_map = dict(zip(dom,material))
            for phase, name in ((0,"low"),(1,"high")):
                replaced = sum(float(domains[i]["area_mm2"]) for i in frac if phase_map[i]==phase)
                reported = r.get(f"replaced_{name}_mm2")
                check(f"replacement_by_actual_background_{name}", abs(replaced-reported)<=0.01 if num(reported) else None, replaced-reported if num(reported) else None, 0.01, "Uses actual rearranged material sequence, never layer parity.")
                background_area=sum(float(domains[i]["area_mm2"]) for i in dom if phase_map[i]==phase)
                saved_bg=r.get(f"background_{name}_area_mm2")
                check(f"background_actual_{name}_area",abs(background_area-saved_bg)<=0.01 if num(saved_bg) else None,background_area)
                original_ids=set(vector(sel.get("inclusion_geometry_domains")))
                original_replaced=sum(float(domains[i]["area_mm2"]) for i in original_ids if phase_map[i]==phase)
                saved_original=r.get(f"original_ellipse_replaced_{name}_mm2")
                if p.get("fracture"):
                    check(f"original_ellipse_replaced_{name}",abs(original_replaced-saved_original)<=0.01 if num(saved_original) else None,original_replaced)
            if p.get("intervention") in ("cut","side"):
                phase = p.get("intervention_phase")
                removed, comp = set(vector(sel.get("removed_domains"))),set(vector(sel.get("compensation_domains")))
                check("intervention_restores_correct_phase", bool(removed) and all(phase_map[i]==phase for i in removed) and removed<= (low if phase==0 else high), sorted(removed))
                check("compensation_replaces_matching_phase", bool(comp) and all(phase_map[i]==phase for i in comp) and comp<=frac, sorted(comp))
        else:
            check("actual_background_material_labels", None, material)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        check("domain_audit_inputs", None, note=str(exc))
    if p.get("intervention") in ("cut","side"):
        removed, comp = r.get("removed_area_mm2"),r.get("compensation_area_mm2")
        err = abs(removed-comp) if num(removed) and num(comp) else None
        check("removed_equals_compensated", err<=0.01 and err/ideal<=1e-4 if num(err) else None, err,min(0.01,1e-4*ideal))
    if p.get("intervention","none")!="none":
        window=r.get("window_area_mm2")
        check("common_local_window_area",abs(window-36)<=0.001 if num(window) else None,window,0.001)
        check("local_mechanism_fields_saved",all((folder/f).is_file() for f in ("local_window_fields.csv","section_current.csv")))
    rawmax=r.get("field_coordinate_max_raw")
    check("raw_field_coordinates_are_geometry_mm",abs(rawmax-p.get("L_mm",120))<1e-5 if num(rawmax) else None,rawmax,1e-5)
    # For flat rearranged sequences the analytic result depends on actual thicknesses.
    z, material = p.get("z_interfaces_mm"),p.get("layer_materials")
    amplitudes=p.get("interface_amplitudes_mm",[])
    if not active(p) and isinstance(z,list) and isinstance(material,list) and len(material)==len(z)+1 and not any(amplitudes) and not p.get("A_mm",0):
        edges=[0]+z+[p.get("L_mm",120)]
        weights=[(b-a)/p.get("L_mm",120) for a,b in zip(edges,edges[1:])]
        conductivities=[p["sigma_high_S_m"] if m==1 else p["sigma_low_S_m"] for m in material]
        expected=1/sum(w/s for w,s in zip(weights,conductivities)) if p.get("load","z")=="z" else sum(w*s for w,s in zip(weights,conductivities))
        check("flat_sequence_analytic", rel(sr,expected)<=1e-5 if num(sr) else None,dict(observed=sr,expected=expected),1e-5)
    status="fail" if any(t["status"]=="fail" for t in checks) else "pending" if any(t["status"]=="pending" for t in checks) else "pass"
    return dict(case_id=c["id"],status=status,checks=checks)


def background(c, axis, byid):
    explicit=c["p"].get("background_case_"+axis)
    if explicit:
        return byid.get(explicit)
    candidates=[a for a in byid.values() if not active(a["p"]) and not a["p"].get("passive_ellipse") and a["p"].get("load","z")==axis and sig(a["p"],BACKGROUND_KEYS)==sig(c["p"],BACKGROUND_KEYS)]
    candidates.sort(key=lambda a:abs(a["p"].get("mesh_level",1)-c["p"].get("mesh_level",1)))
    return candidates[0] if candidates else None


def physical_profile(p):
    """Remove same-material partition edges before assessing physical equivalence."""
    z,amps,layers=p["z_interfaces_mm"],p["interface_amplitudes_mm"],p["layer_materials"]
    kept=[i for i in range(len(z)) if layers[i]!=layers[i+1]]
    return [z[i] for i in kept],[amps[i] for i in kept],[layers[0]]+[layers[i+1] for i in kept]


def transformed_equivalence(a,b):
    """Exact parameter-family D2 audit; 90-degree maps do not preserve z electrodes."""
    pa,pb=a["p"],b["p"]
    required=("z_interfaces_mm","layer_materials","interface_amplitudes_mm")
    if not all(k in p for p in (pa,pb) for k in required):
        return dict(status="pending",reason="missing complete interface family")
    if pa.get("load","z")!=pb.get("load","z") or sig(pa,MATERIAL_KEYS)!=sig(pb,MATERIAL_KEYS):
        return dict(status="pending",reason="loads or material values differ")
    L=pa.get("L_mm",120)
    z,amps,layers=physical_profile(pa)
    zb,ab,lb=physical_profile(pb)
    transforms={"identity":(z,amps,layers,pa["cx_mm"],pa["cz_mm"],pa["theta_deg"]),
                "horizontal_reflection":(z,[-x for x in amps],layers,L-pa["cx_mm"],pa["cz_mm"],180-pa["theta_deg"]),
                "vertical_reflection":([L-x for x in z[::-1]],[-x for x in amps[::-1]],layers[::-1],pa["cx_mm"],L-pa["cz_mm"],-pa["theta_deg"]),
                "rotation_180":([L-x for x in z[::-1]],amps[::-1],layers[::-1],L-pa["cx_mm"],L-pa["cz_mm"],pa["theta_deg"]+180)}
    out=[]
    def same(u,v):
        return len(u)==len(v) and all(abs(x-y)<1e-8 for x,y in zip(u,v))
    for name,(zt,at,lt,x,y,theta) in transforms.items():
        bg=same(zt,zb) and same(at,ab) and lt==lb
        ellipse=abs(x-pb["cx_mm"])<1e-8 and abs(y-pb["cz_mm"])<1e-8 and abs(math.sin(math.radians(theta-pb["theta_deg"])))<1e-10 and pa["a_mm"]==pb["a_mm"] and pa["b_mm"]==pb["b_mm"]
        out.append(dict(transform=name,background_equal=bg,ellipse_equal=ellipse,full_equal=bg and ellipse))
    return dict(status="equivalent" if any(x["full_equal"] for x in out) else "non_equivalent_within_exact_sine_layer_family",transforms=out,excluded_maps="90/270-degree rotations and diagonal reflections exchange z electrodes with insulated x sides")


def centerline_phases(p):
    """Geometric diagnostic only; 4001 fixed samples, independent of FE solution."""
    if not all(k in p for k in ("z_interfaces_mm","layer_materials","interface_amplitudes_mm")):
        return None
    theta=math.radians(p["theta_deg"]); phases=[]
    for i in range(4001):
        s=-p["a_mm"]+2*p["a_mm"]*i/4000
        x=p["cx_mm"]+s*math.cos(theta);z=p["cz_mm"]+s*math.sin(theta)
        layer=sum(z>h+amp*math.sin(2*math.pi*x/p.get("L_mm",120)) for h,amp in zip(p["z_interfaces_mm"],p["interface_amplitudes_mm"]))
        phase=p["layer_materials"][layer]
        if not phases or phase!=phases[-1]:phases.append(phase)
    return phases


def same_buffer(a,b,depth):
    """Exact sufficient equality inside horizontal strips for this sine family."""
    def description(p,bottom):
        z,amps,layers=physical_profile(p)
        edges=[(base,amp,layers[i],layers[i+1]) for i,(base,amp) in enumerate(zip(z,amps)) if (base-abs(amp)<depth if bottom else base+abs(amp)>p.get("L_mm",120)-depth)]
        return dict(outer_phase=layers[0] if bottom else layers[-1],interfaces=edges)
    try:return description(a,True)==description(b,True) and description(a,False)==description(b,False)
    except (KeyError,IndexError,TypeError):return None


def comparisons(cases,audits):
    byid={c["id"]:c for c in cases}; rows=[]; pair_groups={}
    for c in cases:
        p,r=c["p"],c["r"]
        if not active(p):
            continue
        bg=background(c,p.get("load","z"),byid)
        s,s0=sigma(c),sigma(bg) if bg else None
        rows.append(dict(kind="gain",case_id=c["id"],reference_id=bg["id"] if bg else None,status="available" if num(s) and num(s0) else "pending",sigma_reaction_S_m=s,sigma0_reaction_S_m=s0,delta_sigma_S_m=s-s0 if num(s) and num(s0) else None,relative_gain=s/s0-1 if num(s) and num(s0) else None))
        if bg and p.get("intervention","none")=="none" and num(s) and num(s0):
            bound=bg["r"].get("prediction_upper_bound_delta_sigma_S_m")
            dv=bg["r"].get("background_neartip_dV_V")
            if num(bound) or num(dv):
                gain=s-s0
                rows.append(dict(kind="background_predictor",case_id=c["id"],reference_id=bg["id"],status="available",actual_gain_S_m=gain,unfitted_energy_upper_bound_S_m=bound,bound_over_actual_gain=bound/gain if num(bound) and gain>0 else None,actual_gain_below_bound=gain<=bound if num(bound) else None,background_neartip_dV_V=dv,background_neartip_dV2_V2=dv*dv if num(dv) else None,interpretation="unfitted_bound_and_local_driving_indicator_not_a_fitted_law"))
        if p.get("pair_id") and p.get("configuration_role") and p.get("intervention","none")=="none" and p.get("load","z")=="z":
            key=(p["pair_id"],p.get("mesh_level",1),p.get("A_mm",0),p.get("a_mm"),p.get("h_mm"))
            pair_groups.setdefault(key,[]).append(c)
        parent=byid.get(p.get("parent_full_case"))
        if parent and p.get("intervention") in ("cut","side"):
            sp=sigma(parent);ideal=math.pi*p["a_mm"]*p["b_mm"]
            differences={f"{phase}_area_difference_mm2":r.get(f"area_{phase}_mm2",0)-parent["r"].get(f"area_{phase}_mm2",0) for phase in ("low","high","fracture")}
            rows.append(dict(kind="intervention_vs_full",case_id=c["id"],reference_id=parent["id"],status="available" if num(s) and num(sp) else "pending",intervention=p["intervention"],intervention_phase=p.get("intervention_phase"),same_background=sig(p,BACKGROUND_KEYS)==sig(parent["p"],BACKGROUND_KEYS),area_match=all(abs(v)<=min(0.01,1e-4*ideal) for v in differences.values()),conductivity_difference_S_m=s-sp if num(s) and num(sp) else None,conductivity_loss_relative_to_full=1-s/sp if num(s) and num(sp) else None,**differences))
    for key,group in pair_groups.items():
        for a,b in itertools.combinations(sorted(group,key=lambda c:c["p"]["configuration_role"]),2):
            if a["p"]["configuration_role"]==b["p"]["configuration_role"]:
                continue
            pa,pb,ra,rb=a["p"],b["p"],a["r"],b["r"]
            astat,bstat=audits[a["id"]]["status"],audits[b["id"]]["status"]
            controls={"same_fracture_geometry":sig(pa,GEOMETRY_KEYS)==sig(pb,GEOMETRY_KEYS),"same_materials":sig(pa,MATERIAL_KEYS)==sig(pb,MATERIAL_KEYS),"case_numerics":False if "fail" in (astat,bstat) else None if "pending" in (astat,bstat) else True}
            row=dict(kind="configuration_pair",pair_id=key[0],mesh_level=key[1],A_mm=key[2],case_id=a["id"],reference_id=b["id"],role_a=pa["configuration_role"],role_b=pb["configuration_role"],full_length_to_layer_thickness=2*pa["a_mm"]/pa["h_mm"])
            ideal=math.pi*pa["a_mm"]*pa["b_mm"]
            for phase in ("low","high"):
                aa=get(ra,f"original_ellipse_replaced_{phase}_mm2",f"replaced_{phase}_mm2")
                ab=get(rb,f"original_ellipse_replaced_{phase}_mm2",f"replaced_{phase}_mm2")
                d=aa-ab if num(aa) and num(ab) else None
                row[f"replaced_{phase}_difference_mm2"]=d
                row[f"replaced_{phase}_difference_relative_to_ellipse"]=abs(d)/ideal if num(d) else None
                controls["replacement_"+phase]=abs(d)<=0.01 and abs(d)/ideal<=1e-4 if num(d) else None
            for phase in ("low","high","fracture"):
                d=ra[f"area_{phase}_mm2"]-rb[f"area_{phase}_mm2"]
                row[f"final_{phase}_area_difference_mm2"]=d
                controls["final_area_"+phase]=abs(d)<=0.01
            ba,bz=None,None
            for axis in ("x","z"):
                bg_a,bg_b=background(a,axis,byid),background(b,axis,byid)
                sa,sb=sigma(bg_a) if bg_a else None,sigma(bg_b) if bg_b else None
                row[f"background_{axis}_a"]=bg_a["id"] if bg_a else None
                row[f"background_{axis}_b"]=bg_b["id"] if bg_b else None
                row[f"background_{axis}_a_mesh_level"]=bg_a["p"].get("mesh_level") if bg_a else None
                row[f"background_{axis}_b_mesh_level"]=bg_b["p"].get("mesh_level") if bg_b else None
                error=rel(sa,sb) if num(sa) and num(sb) else None
                row[f"background_{axis}_relative_difference"]=error
                controls["background_"+axis]=error<=1e-4 if num(error) else None
                if bg_a and bg_b:
                    for phase in ("low","high"):
                        d=bg_a["r"][f"area_{phase}_mm2"]-bg_b["r"][f"area_{phase}_mm2"]
                        row[f"background_{axis}_{phase}_area_difference_mm2"]=d
                        controls[f"background_{axis}_{phase}_area"]=abs(d)<=0.01
                if axis=="z": ba,bz=bg_a,bg_b
            equivalence=transformed_equivalence(a,b)
            controls["non_equivalent_full_config"]=None if equivalence["status"]=="pending" else equivalence["status"]!="equivalent"
            row["equivalence_audit"]=equivalence
            row["same_bottom_top_10mm_buffer"]=same_buffer(pa,pb,10)
            row["same_bottom_top_20mm_buffer"]=same_buffer(pa,pb,20)
            controls["same_electrode_near_10mm_buffer"]=row["same_bottom_top_10mm_buffer"]
            for name,c in (("a",a),("b",b)):
                fractions=[]
                for side in ("minus","plus"):
                    cap=cap_data(c["r"],side)
                    f=cap["high"]/cap["area"] if num(cap["high"]) and num(cap["area"]) and cap["area"]>0 else cap["fraction"]
                    row[f"{name}_{side}_high_fraction"]=f;fractions.append(f)
                row[f"{name}_both_caps_majority_high"]=all(f>=0.9 for f in fractions) if all(num(f) for f in fractions) else None
                row[f"{name}_both_caps_majority_low"]=all(f<=0.1 for f in fractions) if all(num(f) for f in fractions) else None
                sequence=centerline_phases(c["p"])
                bridges=sum(sequence[i:i+3]==[1,0,1] for i in range(len(sequence)-2)) if sequence else None
                row[f"{name}_centerline_phase_sequence"]=sequence
                row[f"{name}_geometric_high_low_high_connections"]=bridges
                role=c["p"].get("configuration_role","").upper()
                # The fixed-endcap designs deliberately keep one high and one low
                # cap. They must not inherit the central-slab B high/high target.
                if c["p"].get("geometry_family") in FIXED_ENDCAP_FAMILIES:
                    expected=[c["p"].get("expected_endcap_high_fraction_lower"),c["p"].get("expected_endcap_high_fraction_upper")]
                    controls[f"{role}_finite_caps_match_prescribed_fixed_phases"]=all(abs(f-e)<=LIMITS["endcap_relative"] for f,e in zip(fractions,expected)) if all(num(v) for v in fractions+expected) else None
                    expected_bridges=c["p"].get("expected_completely_crossed_low_material_runs")
                    controls[f"{role}_prescribed_internal_low_layer_crossings"]=bridges==expected_bridges if num(bridges) and num(expected_bridges) else None
                else:
                    if role=="B":controls["B_finite_high_caps_and_low_layer_crossing"]=row[f"{name}_both_caps_majority_high"] and bridges>0 if bridges is not None and row[f"{name}_both_caps_majority_high"] is not None else None
                    if role=="A":controls["A_low_caps_or_no_high_low_high_connection"]=row[f"{name}_both_caps_majority_low"] or bridges==0 if bridges is not None and row[f"{name}_both_caps_majority_low"] is not None else None
            fixed_caps=all(p.get("geometry_family") in FIXED_ENDCAP_FAMILIES for p in (pa,pb))
            row["endcap_design_criterion"]="prescribed_fixed_endcaps_and_internal_rearrangement" if fixed_caps else "central_slab_end_access_targets"
            if fixed_caps:
                for side in ("minus","plus"):
                    ac,bc=cap_data(ra,side),cap_data(rb,side)
                    values=[ac[k] for k in ("area","high","low","fraction")]+[bc[k] for k in ("area","high","low","fraction")]
                    controls[f"fixed_{side}_cap_phase_areas_and_fraction_match"]=all(abs(ac[k]-bc[k])<=LIMITS["endcap_absolute_mm2"] for k in ("area","high","low")) and rel(ac["area"],bc["area"])<=LIMITS["endcap_relative"] and abs(ac["fraction"]-bc["fraction"])<=LIMITS["endcap_relative"] if all(num(v) for v in values) else None
            sa,sb=sigma(a),sigma(b)
            s0a,s0b=sigma(ba) if ba else None,sigma(bz) if bz else None
            row.update(sigma_a_S_m=sa,sigma_b_S_m=sb,sigma_relative_difference=rel(sa,sb) if num(sa) and num(sb) else None,sigma_b_minus_a_S_m=sb-sa if num(sa) and num(sb) else None)
            if all(num(x) for x in (sa,sb,s0a,s0b)):
                row.update(gain_a_S_m=sa-s0a,gain_b_S_m=sb-s0b,gain_b_minus_a_S_m=(sb-s0b)-(sa-s0a),relative_gain_a=sa/s0a-1,relative_gain_b=sb/s0b-1)
            row["controls"]=controls
            row["status"]="exploratory_controls_failed" if any(x is False for x in controls.values()) else "pending_controls" if any(x is None for x in controls.values()) else "comparable_not_yet_mechanistically_proven"
            rows.append(row)
    interventions={}
    for c in cases:
        p=c["p"]
        if p.get("intervention") in ("cut","side"):
            key=(p.get("parent_full_case"),p.get("mesh_level"),p.get("intervention_phase"),p.get("cut_u_mm"),p.get("cut_width_mm"))
            interventions.setdefault(key,{})[p["intervention"]]=c
    for key,group in interventions.items():
        if not all(k in group for k in ("cut","side")):continue
        cut,side=group["cut"],group["side"];cp,sp,cr,sr=cut["p"],side["p"],cut["r"],side["r"]
        ss,sc=sigma(side),sigma(cut);ideal=math.pi*cp["a_mm"]*cp["b_mm"]
        deltas={f"{field}_difference":cr[field]-sr[field] for field in ("area_low_mm2","area_high_mm2","area_fracture_mm2","removed_area_mm2","compensation_area_mm2") if num(cr.get(field)) and num(sr.get(field))}
        compkeys=("comp_x_mm","comp_z_mm","comp_radius_mm")
        samecomp=sig(cp,compkeys)==sig(sp,compkeys)
        samebg=sig(cp,BACKGROUND_KEYS)==sig(sp,BACKGROUND_KEYS)
        matched=len(deltas)==5 and all(abs(v)<=min(0.01,1e-4*ideal) for v in deltas.values())
        row=dict(kind="cut_vs_side",case_id=cut["id"],reference_id=side["id"],parent_full_case=key[0],mesh_level=key[1],intervention_phase=key[2],cut_u_mm=key[3],same_background=samebg,same_compensation_geometry=samecomp,matched_phase_and_removed_areas=matched,status="controlled_intervention" if samebg and samecomp and matched else "exploratory_controls_failed",sigma_cut_minus_side_S_m=sc-ss if num(sc) and num(ss) else None,**deltas)
        for field in ("window_dV_V","window_dissipation_W_per_m","window_low_dissipation_W_per_m","window_high_dissipation_W_per_m","window_fracture_dissipation_W_per_m","section_current_signed_A_per_m","section_outside_original_ellipse_A_per_m","section_restored_matrix_A_per_m","section_original_ellipse_current_A_per_m","section_background_bypass_current_A_per_m"):
            if num(cr.get(field)) and num(sr.get(field)):
                row[field+"_cut_minus_side"]=cr[field]-sr[field]
        rows.append(row)
    # Repeat comparisons over three mesh levels, retaining both strong and weak effects.
    families={}
    for c in cases:
        p=c["p"]
        family=(p.get("family_id"),sig(p,GEOMETRY_KEYS+MATERIAL_KEYS+BACKGROUND_KEYS+("load","fracture","passive_ellipse","intervention","cut_u_mm","cut_width_mm","side_u1_mm","side_u2_mm","side_width_mm","comp_radius_mm","comp_x_mm","comp_z_mm")))
        families.setdefault(family,{})[p.get("mesh_level",1)]=c
    for family,levels in families.items():
        ordered=[levels[k] for k in sorted(levels)]
        if len(ordered)<2: continue
        fine=ordered[-1]
        for coarse in ordered[:-1]:
            sf,sc=sigma(fine),sigma(coarse)
            error=rel(sf,sc) if num(sf) and num(sc) else None
            rows.append(dict(kind="mesh_conductivity",case_id=coarse["id"],reference_id=fine["id"],status="available" if num(error) else "pending",available_mesh_levels=sorted(levels),three_mesh_levels=len(levels)>=3,sigma_relative_difference=error,below_1e_4=error<=1e-4 if num(error) else None))
    pair_trends={}
    for row in rows:
        if row["kind"]=="configuration_pair":
            pair_trends.setdefault((row["pair_id"],row["A_mm"],row["role_a"],row["role_b"]),{})[row["mesh_level"]]=row
    for key,levels in pair_trends.items():
        if len(levels)<3: continue
        ordered=[levels[k] for k in sorted(levels)][-3:]
        values=[r.get("gain_b_minus_a_S_m") for r in ordered]
        if not all(num(v) for v in values): continue
        change=max(abs(values[2]-values[1]),abs(values[1]-values[0]));same_sign=all(v>0 for v in values) or all(v<0 for v in values)
        ratio=abs(values[-1])/change if change else None
        rows.append(dict(kind="mesh_pair_effect",pair_id=key[0],case_id=ordered[-1]["case_id"],reference_id=ordered[-1]["reference_id"],status="available",gain_differences_S_m=values,conservative_observed_mesh_change_S_m=change,effect_to_mesh_change_ratio=ratio,same_sign_three_levels=same_sign,numerically_resolved=same_sign and (change==0 or ratio>=5),interpretation="observed_mesh_change_not_rigorous_error_bound"))
    return rows


def fmt(x):
    return f"{x:.8g}" if num(x) else "待验证"


def readme(report,rows,comparison_rows):
    lines=["# 二维第二轮结果与独立审计", "", f"北京时间 {report['snapshot_time_beijing']}。本汇总仅处理第二轮目录，不修改首轮67例及原汇总。", "", f"成功保存 {len(rows)} 例；审计通过 {report['audit_counts']['pass']}，未通过 {report['audit_counts']['fail']}，必要量待验证 {report['audit_counts']['pending']}。尚未开始与进行中的算例不算失败。", "", "主比较使用守恒反力电流电导，体耗散独立验证；原梯度通量只作为诊断。背景x/z相对差≤1e-4，替换面积差每相同时≤0.01 mm²且≤1e-4πab，阈值先于结果确定。", "", "端部为长轴两端各3 mm的真实椭圆端帽，报告实际分相面积，不能用尖端单点标签代替。A/B完整配置按保留电极的反射/180°旋转与材料标签核对，完整等价配置仅作自检。", "", "h在核心设计中表示B中央关键低导层厚；A中央高导段q由面积匹配决定，远处层厚也调整。因此不是全层等厚扫描，2a/h不是唯一相似参数，实际层厚及解释限制逐例保留。", "", "## 已有A/B/C比较", "", "| 配对/网格 | 状态 | σB−σA (S/m) | 电导对称相对差 | 增益差 (S/m) |", "|---|---|---:|---:|---:|"]
    pairs=[r for r in comparison_rows if r["kind"]=="configuration_pair"]
    for r in pairs:
        lines.append(f"| {r['pair_id']} / {r['mesh_level']} / {r['role_a']}-{r['role_b']} | {r['status']} | {fmt(r.get('sigma_b_minus_a_S_m'))} | {fmt(r.get('sigma_relative_difference'))} | {fmt(r.get('gain_b_minus_a_S_m'))} |")
    if not pairs: lines += ["", "尚无可汇总的完整配置配对；不提前声称现象成立或不存在。"]
    lines += ["", "端帽验收按设计家族区分：core检验A低/低、B高/高及关键低导段跨越；strict与tile_wave检验两角色预定的一端高、一端低保持不变，并逐端核对实际面积/分相面积/高相比例及内部跨层序列。不能把core的B高/高目标施加给固定端帽对照。背景电性通过1e-4阈值只表示数值容差内可比，不自动表示精确相同。"]
    for r in pairs:
        if r["mesh_level"]==3 and r["pair_id"].startswith("strict_h10_L60_wave"):
            lines += ["", f"{r['pair_id']}：背景x/z相对差 {fmt(r.get('background_x_relative_difference'))} / {fmt(r.get('background_z_relative_difference'))}；" + ("在预设阈值内近似匹配，非精确相同。" if r["status"]=="comparable_not_yet_mechanistically_proven" else "存在超阈值控制，只作为探索。")]
    lines += ["", "配对通过控制仍不等于桥接是唯一原因。端帽是否同时连接关键低导层两侧、缓冲逐点相同、局部路径和干预相别需连同结果审核。完整检查和缺项见 validation_report.json；失败案例保留 FAILURE.txt，缺字段标待验证。", "", "## 网格与机制", "", f"已产生 {sum(r['kind']=='mesh_pair_effect' for r in comparison_rows)} 项三档配对效应趋势和 {sum(r['kind']=='intervention_vs_full' for r in comparison_rows)} 项对完整裂缝的机制比较。三档检查效应本身及总电导，不按任何目标百分比挑选结果。"]
    effects=[r for r in comparison_rows if r["kind"]=="mesh_pair_effect"]
    if effects:
        lines += ["", "| 配对 | 细档增益差 (S/m) | 观察到的网格变化量 (S/m) | 三档同号且效应/变化≥5 |", "|---|---:|---:|---|"]
        for r in effects:lines.append(f"| {r['pair_id']} | {fmt(r['gain_differences_S_m'][-1])} | {fmt(r['conservative_observed_mesh_change_S_m'])} | {r['numerically_resolved']} |")
    direct=[r for r in comparison_rows if r["kind"]=="cut_vs_side"]
    if direct:
        lines += ["", "| 切断 / 侧削 | 网格 | 原相(0低/1高) | 控制 | σcut−σside (S/m) | 局部ΔV差 (V) |", "|---|---:|---:|---|---:|---:|"]
        for r in direct:lines.append(f"| {r['case_id']} / {r['reference_id']} | {r['mesh_level']} | {r['intervention_phase']} | {r['status']} | {fmt(r.get('sigma_cut_minus_side_S_m'))} | {fmt(r.get('window_dV_V_cut_minus_side'))} |")
        lines += ["", "低/高切口分别匹配其原基质相，并保持同一补回几何。局部有限截面只用于路径诊断；其离散分段积分的闭合残差另列，不能用作总电流。高导域形状和界面仍改变，不称纯拓扑效应。"]
    predictors=[r for r in comparison_rows if r["kind"]=="background_predictor"]
    if predictors:
        lines += ["", f"已保存 {len(predictors)} 项无拟合背景指标检查；comparisons.csv列出原背景近端电势差、能量增量上界及上界/真实增益。上界即便满足不等式，也不等于有准确预测能力；排序和跨配置适用范围须另行分析。"]
    lines += ["", "## 解释边界与下一步", "", "需共同回答首轮约0.83%为何有限、何种配置更强/更弱、剩余含量/界面/边界混杂、不同2a/h下的适用范围，以及第一组现象能否独立成立。未完成的项目不以结论代替。候选预测只用无裂缝场，不用含裂缝真值拟合；不自动进入三维。", "", "重跑：`python scripts/summarize_round2.py`。只读原始结果并写本目录汇总，不启动COMSOL或MATLAB，不作图。", ""]
    return "\n".join(lines)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project",type=Path,default=Path(__file__).resolve().parents[1])
    args=parser.parse_args();run=args.project/"runs"/"2026-10-21_round2"
    if not run.is_dir(): parser.error(f"Directory not found: {run}")
    expected={};manifest_errors=[]
    for path in sorted((run/"manifests").glob("*.json")):
        try:
            values=read(path);values=values.get("cases",[values]) if isinstance(values,dict) else values
            for item in values:
                if isinstance(item,dict) and item.get("id"):
                    expected.setdefault(item["id"],dict(item));expected[item["id"]].setdefault("group",path.stem)
        except (OSError,ValueError,TypeError) as exc:
            manifest_errors.append(dict(file=str(path),error=str(exc)))
    cases=[];states=[];read_errors=[]
    folders={p.name:p for p in run.iterdir() if p.is_dir() and (p/"parameters.json").exists()}
    for cid in sorted(set(expected)|set(folders)):
        folder=run/cid;state="not_started"
        if (folder/"SUCCESS.json").exists():
            try:
                p=read(folder/"parameters.json");r=read(folder/"SUCCESS.json")
                # Manifest metadata fills missing bookkeeping only, never overrides saved physics.
                p={**expected.get(cid,{}),**p}
                cases.append(dict(id=cid,p=p,r=r,folder=folder));state="success"
            except (OSError,ValueError,TypeError) as exc:
                state="pending_read_retry";read_errors.append(dict(case_id=cid,error=str(exc)))
        elif (folder/"FAILURE.txt").exists():state="failed_attempt"
        elif (folder/"parameters.json").exists():state="in_progress_or_interrupted"
        states.append(dict(case_id=cid,status=state,group=expected.get(cid,{}).get("group")))
    audits=[audit(c) for c in cases];audit_byid={a["case_id"]:a for a in audits}
    comparison_rows=comparisons(cases,audit_byid);rows=[]
    for c in cases:
        row={**c["r"],**{"param_"+k:v for k,v in c["p"].items()}}
        row.update(case_id=c["id"],audit_status=audit_byid[c["id"]]["status"],sigma_reaction_audited_S_m=sigma(c),full_length_to_layer_thickness=2*c["p"]["a_mm"]/c["p"]["h_mm"])
        edges=[0]+c["p"].get("z_interfaces_mm",[])+[c["p"].get("L_mm",120)]
        thicknesses=[b-a for a,b in zip(edges,edges[1:])]
        row.update(actual_layer_thicknesses_mm=thicknesses,actual_layer_count=len(thicknesses),all_layers_equal_nominal_h=all(abs(h-c["p"]["h_mm"])<1e-8 for h in thicknesses))
        r=c["r"]
        if all(num(r.get(k)) for k in ("window_dissipation_W_per_m","window_low_dissipation_W_per_m","window_fracture_dissipation_W_per_m")):
            row["window_high_dissipation_W_per_m_derived"]=r["window_dissipation_W_per_m"]-r["window_low_dissipation_W_per_m"]-r["window_fracture_dissipation_W_per_m"]
        if all(num(r.get(k)) for k in ("section_current_signed_A_per_m","section_original_ellipse_current_A_per_m","section_background_bypass_current_A_per_m")):
            residual=r["section_current_signed_A_per_m"]-r["section_original_ellipse_current_A_per_m"]-r["section_background_bypass_current_A_per_m"]
            row["sampled_section_partition_residual_A_per_m"]=residual
            row["sampled_section_partition_relative_residual"]=abs(residual)/max(abs(r["section_current_signed_A_per_m"]),1e-30)
        rows.append(row)
    report=dict(snapshot_time_beijing=datetime.now(timezone(timedelta(hours=8))).isoformat(timespec="seconds"),thresholds=LIMITS,successful_cases=len(cases),expected_manifest_cases=len(expected),audit_counts={s:sum(a["status"]==s for a in audits) for s in ("pass","pending","fail")},case_states=states,cases=audits,manifest_read_errors=manifest_errors,case_read_errors=read_errors,comparison_counts={k:sum(r["kind"]==k for r in comparison_rows) for k in sorted({r["kind"] for r in comparison_rows})},notes=["No missing group is a failed solve.","Missing controls remain pending, never zero-filled to pass.","Thresholds are predeclared and do not target an effect magnitude."])
    csvout(run/"summary.csv",rows,("case_id","audit_status","sigma_reaction_audited_S_m"))
    csvout(run/"comparisons.csv",comparison_rows,("kind","case_id","reference_id","status"))
    write(run/"validation_report.json",json.dumps(clean(report),ensure_ascii=False,indent=2)+"\n")
    write(run/"README_CN.md",readme(report,rows,comparison_rows))
    print(json.dumps(dict(success=len(cases),audit=report["audit_counts"],comparisons=report["comparison_counts"]),ensure_ascii=False))


if __name__=="__main__":main()
