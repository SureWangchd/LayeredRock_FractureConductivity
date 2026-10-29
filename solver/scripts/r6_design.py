# -*- coding: utf-8 -*-
"""Resolve the fixed Round-6 geometry design; never reads target responses.

The supplied JSON prescribes 16 geometries and 32 conductivity states. Existing
Round-5 material maps, including real electrode buffers, remain authoritative.
All reported areas here are analytical geometry, not COMSOL integration results.
"""
from __future__ import annotations

import copy
import csv
import json
import math
from pathlib import Path

from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/round6_applicability"
R5 = ROOT / "runs/round5_modal_closure"
PACK = ROOT / "round6_taskpack"


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_csv(path, rows):
    with Path(path).open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sigma_tag(value):
    return f"{value:.10g}".replace(".", "p")


def phase(c, z):
    for lo, hi, material in zip([0, *c["z_interfaces_mm"]], [*c["z_interfaces_mm"], 120], c["layer_materials"]):
        if lo <= z <= hi:
            return material
    raise ValueError(z)


def physical_runs(c):
    runs = []
    for lo, hi, material in zip([0, *c["z_interfaces_mm"]], [*c["z_interfaces_mm"], 120], c["layer_materials"]):
        if runs and runs[-1][2] == material:
            runs[-1][1] = hi
        else:
            runs.append([lo, hi, material])
    return runs


def projected_cdf(z, c):
    """Exact cumulative ellipse area below horizontal z (mm)."""
    theta = math.radians(c["theta_deg"])
    radius = math.hypot(c["a_mm"] * math.sin(theta), c["b_mm"] * math.cos(theta))
    q = max(-1.0, min(1.0, (z-c["cz_mm"])/radius))
    return c["a_mm"]*c["b_mm"]*(math.pi/2+math.asin(q)+q*math.sqrt(max(0,1-q*q)))


def cap_strip_area(c, u0, u1, z0, z1):
    """Quadrature of a finite axial cap intersected with a horizontal strip.

    Splits at all ellipse/interface intersections so integrand changes are
    explicitly located; no point-tip classification and no geometry raster.
    """
    a, b = c["a_mm"], c["b_mm"]
    st, ct = math.sin(math.radians(c["theta_deg"])), math.cos(math.radians(c["theta_deg"]))
    assert ct > 0
    points = [u0, u1]
    for z in [z0, z1]:
        d = z-c["cz_mm"]
        aa, bb, cc = st*st+(ct*b/a)**2, -2*d*st, d*d-(ct*b)**2
        discr = bb*bb-4*aa*cc
        if discr >= 0:
            for u in [(-bb-math.sqrt(discr))/(2*aa), (-bb+math.sqrt(discr))/(2*aa)]:
                if u0 < u < u1:
                    points.append(u)
    points = sorted(set(points))

    def width(u):
        v = b*math.sqrt(max(0, 1-(u/a)**2))
        low = max(-v, (z0-c["cz_mm"]-u*st)/ct)
        high = min(v, (z1-c["cz_mm"]-u*st)/ct)
        return max(0, high-low)

    result, error = 0.0, 0.0
    for lo, hi in zip(points[:-1], points[1:]):
        value, err = quad(width, lo, hi, epsabs=2e-11, epsrel=2e-12, limit=100)
        result += value
        error += err
    return result, error


def analyze_geometry(c):
    a, b = c["a_mm"], c["b_mm"]
    theta = math.radians(c["theta_deg"])
    rz, rx = math.hypot(a*math.sin(theta), b*math.cos(theta)), math.hypot(a*math.cos(theta), b*math.sin(theta))
    runs = physical_runs(c)
    ellipse_total = math.pi*a*b
    high = sum(projected_cdf(hi, c)-projected_cdf(lo, c) for lo, hi, p in runs if p)
    low = ellipse_total-high
    primitive = lambda u: b*(u*math.sqrt(max(0, 1-(u/a)**2))+a*math.asin(u/a))
    cap_length = c["endcap_length_mm"]
    cap_area = primitive(a)-primitive(a-cap_length)
    caps, max_error = {}, 0.0
    for name, u0, u1 in [("minus", -a, -a+cap_length), ("plus", a-cap_length, a)]:
        values = {0: 0.0, 1: 0.0}
        for z0, z1, p in runs:
            value, err = cap_strip_area(c, u0, u1, z0, z1)
            values[p] += value
            max_error = max(max_error, err)
        assert abs(sum(values.values())-cap_area) < 2e-9
        fraction = values[1]/cap_area
        caps[name] = {"axial_interval_mm": [u0, u1], "axial_length_mm": cap_length,
                      "area_mm2": cap_area, "high_area_mm2": values[1], "low_area_mm2": values[0],
                      "high_fraction": fraction,
                      "material_class": "high" if fraction > 1-1e-9 else ("low" if fraction < 1e-9 else "mixed")}
    layers = [{"z_interval_mm": [lo, hi], "material": "high" if p else "low",
               "ellipse_overlap_mm2": projected_cdf(hi,c)-projected_cdf(lo,c)}
              for lo, hi, p in runs if hi > c["cz_mm"]-rz and lo < c["cz_mm"]+rz]
    return {
        "scope": "Analytical geometry only; actual COMSOL phase integrals saved separately by solver",
        "ellipse_area_analytic_mm2": ellipse_total,
        "ellipse_replaced_high_analytic_mm2": high,
        "ellipse_replaced_low_analytic_mm2": low,
        "ellipse_vertical_envelope_mm": [c["cz_mm"]-rz,c["cz_mm"]+rz],
        "ellipse_horizontal_envelope_mm": [c["cx_mm"]-rx,c["cx_mm"]+rx],
        "clearance_to_fixed_core_interface_mm": min(c["cz_mm"]-rz-30,90-c["cz_mm"]-rz),
        "clearance_to_external_boundary_mm": min(c["cz_mm"]-rz,120-c["cz_mm"]-rz,c["cx_mm"]-rx,120-c["cx_mm"]-rx),
        "background_high_area_analytic_mm2": 120*sum(hi-lo for lo,hi,p in runs if p),
        "background_low_area_analytic_mm2": 120*sum(hi-lo for lo,hi,p in runs if not p),
        "final_high_area_analytic_mm2": 7200-high,
        "final_low_area_analytic_mm2": 7200-low,
        "final_fracture_area_analytic_mm2": ellipse_total,
        "finite_endcaps": caps,
        "overlapped_physical_material_runs": layers,
        "cap_quadrature_max_estimated_abs_error_mm2": max_error,
        "G0z_analytic_S": .12/sum((hi-lo)*.001/(.1 if p else .01) for lo,hi,p in runs),
        "sigma_eff_x_analytic_S_m": .055,
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    prescribed = read_json(PACK/"round6_design.json")
    old = read_json(R5/"design.json")
    source = {(c["configuration_role"],c["layout"]): c for c in old["cases"]}
    cases, states, sources, geometry_rows = [], [], [], []
    by = {}
    for state in prescribed["states"]:
        length, role, layout = int(state["full_length_mm"]), state["role"], state["layout"]
        cid = f"r6_L{length}_{role}_{layout}"
        source_key = f"L{length}_{role}"
        if cid not in by:
            original = source[(role, layout)]
            c = copy.deepcopy(original)
            for stale in ["geometry_only_checks", "expected_completely_crossed_low_material_runs", "area_matching_analytic_abs_error_mm2"]:
                c.pop(stale, None)
            c.update({"id":cid, "group":"round6_applicability", "phase":"prospective_new_source_geometry",
                      "family_id":f"strict_h10_L{length}_fixed_core_peripheral_band_position", "pair_id":f"r6_L{length}_{layout}",
                      "full_length_mm":length,"full_width_mm":3,"a_mm":length/2,"b_mm":1.5,
                      "source_key":source_key,"reference_case_id":f"r6_L{length}_{role}_P77",
                      "source_parameters":"runs/round5_modal_closure/design.json", "source_case_id":original["id"],
                      "parent_full_case":original["id"], "existing_geometry_case_id":"", "existing_case_directory":"",
                      "background_case_z":cid+"_bgz", "background_case_x":cid+"_bgx",
                      "solver_input_guard":"Only backgrounds/reference predictions allowed until all 32 prospective predictions and intervals are saved.",
                      "design_basis":"Lengths40/64 and fourlayouts fixed by supplied Round6 JSON before new target solutions; no response-based screening.",
                      "interpretation_limit":"New source lengths alter core layer coverage and finite endcap phases; do not transfer original60mm matched-endcap claim.",
                      "geometry_only_area_note":"Analytical and quadrature values; independently export actualCOMSOL integrations.",
                      "coordinate_convention":"x horizontal,z vertical; counterclockwise long axis60deg from+x; COMSOL2D y denotes physicalz; geometry mm, internalSI",
                      "target_sigma_f_S_m":[.18,10],"external_voltage_V":1})
            checks = analyze_geometry(c)
            c["geometry_only_checks"] = checks
            c["expected_replaced_high_area_mm2"] = checks["ellipse_replaced_high_analytic_mm2"]
            c["expected_replaced_low_area_mm2"] = checks["ellipse_replaced_low_analytic_mm2"]
            c["endcap_area_each_mm2"] = checks["finite_endcaps"]["minus"]["area_mm2"]
            c["expected_endcap_high_fraction_lower"] = checks["finite_endcaps"]["minus"]["high_fraction"]
            c["expected_endcap_high_fraction_upper"] = checks["finite_endcaps"]["plus"]["high_fraction"]
            assert checks["background_high_area_analytic_mm2"] == checks["background_low_area_analytic_mm2"] == 7200
            assert checks["clearance_to_fixed_core_interface_mm"] > 0
            assert abs(checks["clearance_to_fixed_core_interface_mm"]-state["clearance_to_fixed_core_interface_mm"]) < 1e-12
            for z in [5,15]: assert phase(c,z)==0
            for z in [105,115]: assert phase(c,z)==1
            assert c["z_interfaces_mm"] == original["z_interfaces_mm"] and c["layer_materials"] == original["layer_materials"]
            by[cid] = c
            cases.append(c)
            geometry_rows.append({"case_id":cid,"source_key":source_key,"role":role,"full_length_mm":length,"layout":layout,
              "distance_lower_mm":c["distance_lower_mm"],"distance_upper_mm":c["distance_upper_mm"],
              **{key:checks[key] for key in ["ellipse_area_analytic_mm2","ellipse_replaced_high_analytic_mm2","ellipse_replaced_low_analytic_mm2",
                   "clearance_to_fixed_core_interface_mm","background_high_area_analytic_mm2","background_low_area_analytic_mm2",
                   "final_high_area_analytic_mm2","final_low_area_analytic_mm2","final_fracture_area_analytic_mm2","G0z_analytic_S"]},
              "endcap_length_mm":c["endcap_length_mm"],"endcap_area_each_mm2":c["endcap_area_each_mm2"],
              "capminus_high_fraction":c["expected_endcap_high_fraction_lower"],"capplus_high_fraction":c["expected_endcap_high_fraction_upper"],
              "capminus_material":checks["finite_endcaps"]["minus"]["material_class"],"capplus_material":checks["finite_endcaps"]["plus"]["material_class"],
              "is_COMSOL_result":False})
        sf = state["sigma_f_S_m"]
        states.append({"state_id":cid+"_sf_"+sigma_tag(sf),"taskpack_state_id":state["case_id"],"case_id":cid,"source_key":source_key,
                       "role":role,"full_length_mm":length,"layout":layout,"sigma_f_S_m":sf,"phase":"prospective_new_source_geometry",
                       "target_relative_folder":f"truth/{cid}/sf_{sigma_tag(sf)}",
                       "prediction_relative_folder":f"predictions/{cid}/sf_{sigma_tag(sf)}",
                       "prediction_file":f"predictions/{cid}/sf_{sigma_tag(sf)}/prediction.json",
                       "target_status_at_design":"not_computed"})
    for length in [40,64]:
        for role in ["A","B"]:
            c = by[f"r6_L{length}_{role}_P77"]
            sources.append({"source_key":f"L{length}_{role}","role":role,"full_length_mm":length,"reference_case_id":c["id"],
                            "reference_relative_folder":f"reference/L{length}_{role}_base",
                            "geometry_source_kind":"new passiveellipse P77 COMSOL background, never a60mm target field",
                            "clearance_to_fixed_core_interface_mm":c["geometry_only_checks"]["clearance_to_fixed_core_interface_mm"]})
    assert len(cases)==16 and len(states)==32 and len(sources)==4
    for length in [40,64]:
        for role in ["A","B"]:
            group = [c for c in cases if c["full_length_mm"]==length and c["configuration_role"]==role]
            for c in group[1:]:
                for key in ["ellipse_replaced_high_analytic_mm2","ellipse_replaced_low_analytic_mm2","finite_endcaps"]:
                    assert c["geometry_only_checks"][key] == group[0]["geometry_only_checks"][key]
    layout_specs = [copy.deepcopy(s) for s in old["layout_specs"] if s["layout_id"] in ["P11","P77","P17","P71"]]
    for s in layout_specs: s["phase"]="round6_fixed_periphery_new_source"
    design = {"schema":"round6_applicability_geometry_v1","status":"prospective_design_only_no_new_target_answers",
      "design_date_Asia_Shanghai":"2026-10-26","root":ROOT.as_posix(),"source_taskpack_design":"round6_taskpack/round6_design.json",
      "actual_geometry_authority":"ExistingRound5 strictA/B materialmaps; prescribednewlengths40/64; unchangedreal20mm electrodebuffers",
      "common_z_interfaces_mm":old["common_z_interfaces_mm"],"core_z_mm":[30,90],"out_of_plane_thickness_m":1,
      "sample_width_mm":120,"sample_height_mm":120,"M0_equivalent_tensor":old["M0_equivalent_tensor"],
      "layout_specs":layout_specs,"cases":cases,"new_states":states,"source_groups":sources,
      "new_geometry_count":16,"new_initial_target_count":32,"modes":prescribed["modes"],"reference_order_checks":[32,64],
      "primary_gain_relative_tolerance":.02,"secondary_gain_tolerances":[.01,.05],"single_side_relative_tolerance":.1,
      "maximum_planned_target_background_refinement_pairs":6,"prediction_before_truth":True,
      "all_32_predictions_and_intervals_required_before_first_new_target":True,
      "high_mode_tail":"All unretained modes keep common anisotropicM0 exterior; never insulated",
      "no_new_target_field_read_by_design":True,"units":{"geometry":"mm","conductivity":"S/m","conductance":"S","voltage":"V","out_of_plane_thickness":"m"}}
    write_json(OUT/"design.json",design)
    write_csv(OUT/"geometry_candidates.csv",geometry_rows)
    with (R5/"port_coefficients.csv").open(encoding="utf-8-sig",newline="") as f:
        port_rows = [row for row in csv.DictReader(f) if row["layout"] in ["P11","P77","P17","P71"]]
    write_csv(OUT/"port_coefficients.csv",port_rows)
    write_json(OUT/"geometry_design_ready.json",{"status":"success","new_geometry_count":16,"new_target_count":32,
         "source_groups":sources,"scope":"Geometry plus analyticalcoefficients only, notnewCOMSOL results",
         "all_new_target_status":"not_computed_by_this_script"})
    rows = []
    for c in cases:
        if c["layout"]!="P77":continue
        q=c["geometry_only_checks"]
        rows.append(f"| {c['full_length_mm']} | {c['configuration_role']} | {q['clearance_to_fixed_core_interface_mm']:.9f} | {q['ellipse_area_analytic_mm2']:.9f} | {q['ellipse_replaced_high_analytic_mm2']:.9f} | {q['ellipse_replaced_low_analytic_mm2']:.9f} | {c['expected_endcap_high_fraction_lower']:.9f} | {c['expected_endcap_high_fraction_upper']:.9f} |")
    note="""# 第六轮新裂缝长度的几何记录

2026-10-26。本文件由设计脚本生成，只有解析几何与有限端帽积分，不是COMSOL求解结果。16种新几何、32种电导率状态完全按任务包指定；不读目标场，不按效果筛选。

沿用第五轮真实层序：下0–20 mm低导电极缓冲，上100–120 mm高导电极缓冲。核心z=30–90 mm固定；外围只有原来的2 mm薄带在d=1或7 mm处。每个长度/角色的四种外围布局保持核心、椭圆以及替换相面积相同。A/B的核心原本不同。

椭圆长轴全长40或64 mm、短轴全长3 mm，中心(60,60) mm，长轴从+x逆时针60°。端帽采用既有求解器的有限轴向长度3 mm并包括完整横截面，不采用数学尖端单点。端帽材料已重新计算，没有复制60 mm的端帽标签。

| 全长mm | 角色 | 到核心界面净距mm | 椭圆面积mm² | 替换高导mm² | 替换低导mm² | 下端帽高导比例 | 上端帽高导比例 |
|---:|:---:|---:|---:|---:|---:|---:|---:|
"""+"\n".join(rows)+"""

两种长度的A/B仍具有相同总替换高/低导面积，但40 mm时B的端帽类别已改变。因此本轮不能沿用“端部完全匹配”的60 mm强命题，也不能把长度变化等同于只改变净距；它同时改变覆盖层段及激发的空间分布。

整椭圆的相面积使用水平投影解析累积式；有限端帽与各水平材料带交集用局部长轴坐标一维积分，显式分割所有椭圆/界面交点。数值积分误差记录在design.json，不代替COMSOL域积分。真实模型须分别导出背景相面积、原椭圆替换面积和有限端帽积分。

外围系数直接继承已确认真实缓冲层的第五轮系数，四种布局与长度无关。下Y0=0.438596491228 S/m²，上Y0=0.980392156863 S/m²。M0外围仍为原等效张量；n>K高阶响应不设为绝缘。

design.json：cases为16个完整求解器参数字典；new_states为32个预先规定状态；source_groups为4个长度/角色的M0新参考入口，均从对应P77被动椭圆背景建立。target_relative_folder和prediction_file记录统一断点续算目录。尚未完成全部32例预测与区间前，不得求解新完整目标。

复算入口：`D:/Users/PC/anaconda3/envs/efield310/python.exe scripts/r6_design.py`。只写本轮设计文件，不覆盖旧轮记录，不生成散列审计。
"""
    (OUT/"geometry_note_CN.md").write_text(note,encoding="utf-8")
    print(json.dumps({"cases":len(cases),"states":len(states),"source_groups":sources,
                      "design":str(OUT/"design.json"),"geometry_summary":rows},ensure_ascii=False))


if __name__ == "__main__":
    main()
