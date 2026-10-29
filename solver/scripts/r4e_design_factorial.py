"""Prepare four geometry-only controls for the round-4 peripheral factorial.

No COMSOL process is started. Existing round-4 files are read only.
"""
from __future__ import annotations

import copy
import csv
import json
from pathlib import Path

from r4_design_periphery import (
    exact_changed_measure_mm,
    integrated_high_mm,
    phase,
    physical_interfaces,
    potential_at,
    read_json,
)


ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "runs" / "round4_extension_20261024"
OLD = ROOT / "runs" / "round4_crossover"
SIGMA = [0.1, 0.159375, 0.3]
LAYOUTS = {
    "NF": {"lower": [25, 27], "upper": [95, 97], "lower_position": "near", "upper_position": "far"},
    "FN": {"lower": [23, 25], "upper": [93, 95], "lower_position": "far", "upper_position": "near"},
}


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    existing = read_json(OLD / "periphery_manifest.json")["cases"]
    old_by = {(c["configuration_role"], c["layout"]): c for c in existing}
    cases, rows, references = [], [], []
    for role in ["A", "B"]:
        nn = old_by[(role, "swap_near")]
        ff = old_by[(role, "swap_far")]
        original = read_json(ROOT / nn["source_parameters"])
        for label, old_case in [("NN", nn), ("FF", ff)]:
            references.append({
                "factorial_layout": label,
                "configuration_role": role,
                "case_id": old_case["id"],
                "source_folder": (OLD / "periphery" / old_case["id"]).as_posix(),
                "sigma_f_S_m": SIGMA,
                "reuse_policy": "Reuse exact existing solved states; do not overwrite or recompute by default.",
            })
        for label, spec in LAYOUTS.items():
            c = copy.deepcopy(nn)
            c.update({
                "id": f"r4e_strict_{role}_{label}",
                "group": "round4_extension_peripheral_factorial",
                "family_id": "strict_h10_L60_peripheral_position_factorial",
                "geometry_family": "strict_fixed_core_independent_upper_lower_strip_positions",
                "pair_id": f"r4e_strict_{label}",
                "layout": label,
                "factorial_layout": label,
                "lower_position": spec["lower_position"],
                "upper_position": spec["upper_position"],
                "passive_ellipse": True,
                "fracture": True,
                "load": "z",
                "mesh_level": 3,
                "candidate_state": "geometry_only_not_COMSOL_solved",
                "design_description_CN": "补齐已有NN/FF的两因子位置对照；只移动固定宽度外围条带，不再增加材料界面。",
                "interpretation_limit": "First letter is lower low-band position; second is upper high-band position. Comparisons among NN/NF/FN/FF preserve physical interface count and each peripheral side's phase areas. This is a mechanism experiment designed after known round-4 results, not an independent blind geometry test.",
                "lower_high_to_low_strip_z_mm": spec["lower"],
                "upper_low_to_high_strip_z_mm": spec["upper"],
                "lower_shift_from_virtual_neutral_mm": 1 if spec["lower_position"] == "near" else -1,
                "upper_shift_from_virtual_neutral_mm": -1 if spec["upper_position"] == "near" else 1,
                "cross_periphery_exchange_center_distance_mm": (sum(spec["upper"]) - sum(spec["lower"])) / 2,
                "factorial_reference_NN_case_id": nn["id"],
                "factorial_reference_FF_case_id": ff["id"],
                "no_new_material_interfaces_relative_to_NN_or_FF": True,
                "background_core_potential_relative_to_NN_analytic_V": 0.0,
                "background_core_potential_relative_to_FF_analytic_V": 0.0,
            })
            # Keep the complete original round-4 partition list, including same-phase boundaries.
            edges = [0] + c["z_interfaces_mm"] + [120]
            materials = []
            for lo, hi in zip(edges[:-1], edges[1:]):
                z = (lo + hi) / 2
                p = phase(original, z)
                if spec["lower"][0] < z < spec["lower"][1]:
                    assert p == 1
                    p = 0
                if spec["upper"][0] < z < spec["upper"][1]:
                    assert p == 0
                    p = 1
                materials.append(p)
            c["layer_materials"] = materials
            c["background_case_z"] = c["id"] + "_bgz"
            c["background_case_x"] = c["id"] + "_bgx"
            c["physical_interfaces_z_mm"] = physical_interfaces(c)
            c["physical_interface_count"] = len(c["physical_interfaces_z_mm"])
            c["background_core_potential_offset_analytic_V"] = potential_at(c, 60) - potential_at(original, 60)
            checks = {
                "changed_core_area_relative_original_analytic_mm2": 120 * exact_changed_measure_mm(original, c, 30, 90),
                "changed_bottom_buffer_relative_original_analytic_mm2": 120 * exact_changed_measure_mm(original, c, 0, 20),
                "changed_top_buffer_relative_original_analytic_mm2": 120 * exact_changed_measure_mm(original, c, 100, 120),
                "high_area_total_analytic_mm2": 120 * integrated_high_mm(c),
                "high_area_lower_periphery_analytic_mm2": 120 * integrated_high_mm(c, 20, 30),
                "high_area_upper_periphery_analytic_mm2": 120 * integrated_high_mm(c, 90, 100),
                "high_area_lower_difference_relative_NN_analytic_mm2": 120 * (integrated_high_mm(c, 20, 30) - integrated_high_mm(nn, 20, 30)),
                "high_area_upper_difference_relative_NN_analytic_mm2": 120 * (integrated_high_mm(c, 90, 100) - integrated_high_mm(nn, 90, 100)),
                "changed_area_relative_NN_analytic_mm2": 120 * exact_changed_measure_mm(nn, c, 0, 120),
                "changed_area_relative_FF_analytic_mm2": 120 * exact_changed_measure_mm(ff, c, 0, 120),
                "interface_count_difference_relative_NN": len(physical_interfaces(c)) - len(physical_interfaces(nn)),
                "interface_count_difference_relative_FF": len(physical_interfaces(c)) - len(physical_interfaces(ff)),
                "core_background_potential_max_difference_relative_NN_analytic_V": max(abs(potential_at(c, z) - potential_at(nn, z)) for z in [30, 31, 35, 45, 55, 65, 75, 85, 89, 90]),
                "core_background_potential_max_difference_relative_FF_analytic_V": max(abs(potential_at(c, z) - potential_at(ff, z)) for z in [30, 31, 35, 45, 55, 65, 75, 85, 89, 90]),
            }
            for key in ["changed_core_area_relative_original_analytic_mm2", "changed_bottom_buffer_relative_original_analytic_mm2", "changed_top_buffer_relative_original_analytic_mm2", "high_area_lower_difference_relative_NN_analytic_mm2", "high_area_upper_difference_relative_NN_analytic_mm2", "interface_count_difference_relative_NN", "interface_count_difference_relative_FF"]:
                assert checks[key] == 0, (c["id"], key, checks[key])
            assert checks["high_area_total_analytic_mm2"] == 7200
            assert checks["high_area_lower_periphery_analytic_mm2"] == 960
            assert checks["high_area_upper_periphery_analytic_mm2"] == 240
            assert checks["changed_area_relative_NN_analytic_mm2"] == 480
            assert checks["changed_area_relative_FF_analytic_mm2"] == 480
            assert checks["core_background_potential_max_difference_relative_NN_analytic_V"] < 1e-14
            assert checks["core_background_potential_max_difference_relative_FF_analytic_V"] < 1e-14
            assert abs(c["background_core_potential_offset_analytic_V"] - 3 / 110) < 1e-14
            assert c["z_interfaces_mm"] == nn["z_interfaces_mm"] == ff["z_interfaces_mm"]
            assert len(materials) == len(c["interface_amplitudes_mm"]) + 1
            c["geometry_only_checks"] = checks
            c["geometry_only_area_note"] = "Analytical design quantities only; actual COMSOL domain areas and x/z backgrounds must be recorded separately."
            cases.append(c)
            rows.append({"case_id": c["id"], "configuration_role": role, "factorial_layout": label, "lower_start_mm": spec["lower"][0], "lower_end_mm": spec["lower"][1], "upper_start_mm": spec["upper"][0], "upper_end_mm": spec["upper"][1], "physical_interface_count": c["physical_interface_count"], **checks})
    manifest = {
        "schema": "round4_extension_peripheral_factorial_v1",
        "status": "geometry_only_pre_solve",
        "date": "2026-10-24",
        "case_count": 4,
        "new_target_count_prespecified": 12,
        "sigma_f_S_m": SIGMA,
        "prediction_sigma_points": SIGMA,
        "factorial_letter_order": ["lower_low_band_position", "upper_high_band_position"],
        "known_round4_results_used_to_design_mechanism_experiment": True,
        "independent_blind_geometry_test": False,
        "all_new_predictions_to_be_saved_before_any_new_target": True,
        "frozen_method_version": "r3_background_vector4_residual_rankguard_v1",
        "old_references": references,
        "common_z_interfaces_mm": cases[0]["z_interfaces_mm"],
        "cases": cases,
    }
    (DEST / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    with (DEST / "geometry_candidates.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (DEST / "design_CN.md").write_text(DESIGN_TEXT, encoding="utf-8")
    print(json.dumps({"manifest": str(DEST / "manifest.json"), "new_cases": [c["id"] for c in cases], "prespecified_target_count": 12, "analytical_checks": "passed; no COMSOL solved by this script"}, ensure_ascii=False))


DESIGN_TEXT = """# 第四轮追加机制实验：把上、下外围的位置作用分开

日期：2026-10-24。这里只生成候选几何和求解前设计，没有新的 COMSOL 结果。

已有 NN（near）与 FF（far）同时移动上下两条带，无法区分是哪一侧造成响应。新增 NF、FN 补齐 2×2 因子表。第一个字母始终指**下方低导带**，第二个字母指**上方高导带**。N 为靠近核心，F 为远离核心。

|布局|下方高→低薄带 z（mm）|上方低→高薄带 z（mm）|安排|
|---|---|---|---|
|NN|25–27|93–95|复用第四轮 swap_near|
|NF|25–27|95–97|新增|
|FN|23–25|93–95|新增|
|FF|23–25|95–97|复用第四轮 swap_far|

每一布局均分别保留 strict A、B 两个内部层序，不引入 core 家族。新增 4 个背景，首批固定 σf=0.1、0.159375、0.3 S/m，共 12 个完整含裂缝目标。旧 NN/FF 这三个点全部复用。新背景先计算 x/z 响应及冻结第三轮 v1 预测，所有首批预测保存后才运行新目标；不调整预测器，不把这些依据已知第四轮现象设计的机制实验称为独立几何盲测。

## 几何与背景控制

继续 120×120 mm，椭圆60×3 mm、60°、中心(60,60) mm，σlow=0.01、σhigh=0.1 S/m。保持全宽z=30–90 mm核心、z=0–20与100–120 mm电极缓冲、椭圆内部原相图和有限端帽及其保护邻域逐点不变。

沿用旧四背景相同的17条内部几何分界（包括同相分块），不添加新分界：10,20,23,25,27,30,40,50,60,70,80,90,93,95,97,100,110 mm。全部起伏振幅为0，沿用 m3 网格与选择集构建方式，实际域编号由新 COMSOL 几何提取。

四个布局相互比较时，每侧外围相面积也固定：下方z=20–30 mm高导960、低导240 mm²，上方z=90–100 mm高导240、低导960 mm²；全域两相各7200 mm²。物理界面数/长度相同（A为7条，B为11条；不能把同相分块当物理界面）。单独移动一条2 mm宽带2 mm，改变材料类别的面积合计480 mm²，但两相的总量及各侧总量均不变。

因此 NN↔NF 与 FN↔FF 是只动上带；NN↔FN 与 NF↔FF 是只动下带。原始未交换布局不进入这个两因子差分；原始→交换新增4条材料界面的变化不能混入纯位置作用。

各布局的无裂缝方向电导率解析值为σeff,z=0.0181818181818、σeff,x=0.055 S/m。由于每侧总串联电阻也相同，**核心中连绝对背景电位都相同**，不只是端帽电位差相同。相对于原始未交换布局，核心仍整体偏移3/110 V。核心背景电场、椭圆内原相图和冻结场潜在增益U均相同。以上为解析设计检查；实际背景响应、相面积和端帽量必须另存 COMSOL 积分值。

## 预先规定的比较量

每个角色r独立使用其真实无裂缝背景定义 g[r,L,U]=Gwith−G0，再计算 D[L,U]=g[B,L,U]−g[A,L,U]。

对 gA、gB、D 都报告四个条件效应：

- 只将下带移近：上带近时 YNN−YFN，上带远时 YNF−YFF。
- 只将上带移近：下带近时 YNN−YNF，下带远时 YFN−YFF。
- 交互量 I=YNN−YNF−YFN+YFF，检验两侧的作用是否依赖另一侧位置。
- 下侧平均位置效应 EL=[(YNN−YFN)+(YNF−YFF)]/2；上侧 EU=[(YNN−YNF)+(YFN−YFF)]/2。EL+EU=YNN−YFF 是代数恒等式；有交互时不要把这两项说成唯一的独立物理份额。

对 Sf=dG/dσf 及 M=Sf,B−Sf,A 同样计算条件差，分开回答已累计的优势和继续增大σf的边际响应。Sf与分段场均使用完整目标解，明确属于事后机制诊断。只在共同物理坐标的固定8个椭圆区段比较，不把分段耗散份额直接当因果贡献。

## 可被结果否定的解释

候选一：上侧高导带距离主要控制 A 的上半部响应。如果只动上带对 gA、Sf,A 和上半部场产生稳定变化，而只动下带较小，才支持此局部解释；若相反，则应否定。不能只因已有 NN−FF 的上半部Sf变化较大就先认定上带是原因。

候选二：上下外围可近似独立叠加。若 I 在代表性细网格下明确大于所观察的差分网格变化，说明两侧通过全域电位重分布存在可分辨耦合，应否定简单叠加；若 I 接近数值分辨尺度，只能说尚未分辨交互，不能证明完全独立。

候选三：一个端帽电位差或背景局部电场足以决定增益。本组所有新背景的这些量都相同；只要新目标的增益差达到数值可分辨程度，这种仅依赖局部标量的解释就不充分。但这不自动证明某个唯一绕流通道，也不是新预测算法。

首批不追求所有混合布局的精确交叉根。先用三个公共电导率点识别条件作用；只有混合布局的排序/交互显示明确新问题时，才少量补点，并保留新增点的后验选点属性。必要的细网格只针对新条件差或交互最有代表性的一组，使用各自同网格背景，避免把全局守恒当作局部差值精度证明。

## 输出与复现

生成器：F:/岩层裂缝/scripts/r4e_design_factorial.py。

本目录 manifest.json 包含4个新增完整参数、旧NN/FF的来源、首批电导率和解析保护检查；geometry_candidates.csv 记录候选几何量。生成器不启动 MATLAB/COMSOL，不写入前三轮或第四轮原目录。实际求解另行保存模型、预测、场、相面积及响应，不能用候选解析值替代计算结果。
"""


if __name__ == "__main__":
    main()
