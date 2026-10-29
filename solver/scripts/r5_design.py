# -*- coding: utf-8 -*-
"""Round-5 geometry and analytical exterior coefficients, without target solves.

Existing strict A/B material maps and the actual 20 mm electrode buffers are
authoritative. The taskpack's isolated-band illustration omits those buffers;
this script documents that difference and does not alter the supplied files.
"""
from __future__ import annotations

import copy
import csv
import importlib
import io
import json
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "runs/round5_modal_closure"
PACK = ROOT / "round5_taskpack"
sys.path.insert(0, str(PACK))

from periphery_modal_response import (  # noqa: E402
    Layer, equivalent_modal_admittances, equivalent_tensor,
    input_admittance, modal_admittances, original_transfer_check,
)
from r4_design_periphery import (  # noqa: E402
    exact_changed_measure_mm, integrated_high_mm, phase,
    physical_interfaces, potential_at, read_json,
)
from r4e_peripheral_modes import admittance, layers_from_electrode  # noqa: E402


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_csv(path, rows):
    with Path(path).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sigma_tag(value):
    return f"{value:.10g}".replace(".", "p")


def core_outward_layers(case, side):
    """Merge same-phase geometric partitions and reverse old electrode order."""
    result = []
    for h, sigma in reversed(layers_from_electrode(case, side)):
        if result and result[-1].conductivity_S_m == sigma:
            result[-1] = Layer(result[-1].thickness_m + h, sigma)
        else:
            result.append(Layer(h, sigma))
    return result


def run_reference_tests():
    # Importing the module avoids its __main__ writer in the original taskpack.
    module = importlib.import_module("test_reference")
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(module.TestReference)
    capture = io.StringIO()
    result = unittest.TextTestRunner(stream=capture, verbosity=2).run(suite)
    (DEST / "reference_selftest_output.txt").write_text(capture.getvalue(), encoding="utf-8")
    record = {
        "scope": "Taskpack analytical/algebra self-tests only; not COMSOL results.",
        "execution_date_Asia_Shanghai": "2026-10-26",
        "original_taskpack_files_modified": False,
        "tests_run": result.testsRun,
        "passed": result.testsRun - len(result.failures) - len(result.errors),
        "failures": len(result.failures), "errors": len(result.errors),
        "numpy_version": module.np.__version__,
    }
    write_json(DEST / "reference_selftest.json", record)
    if not result.wasSuccessful():
        raise RuntimeError("Reference algebra tests failed; see saved output")
    return record


def old_development_cases():
    cases = {}
    for c in read_json(ROOT / "runs/round4_crossover/periphery_manifest.json")["cases"]:
        label = "NN" if c["layout"] == "swap_near" else "FF"
        folder = ROOT / "runs/round4_crossover/periphery" / c["id"]
        cases[(c["configuration_role"], label)] = (c, folder)
    for c in read_json(ROOT / "runs/round4_extension_20261024/manifest.json")["cases"]:
        folder = ROOT / "runs/round4_extension_20261024/cases" / c["id"]
        cases[(c["configuration_role"], c["layout"])] = (c, folder)
    return cases


def limited_prior_geometry_check(new_cases):
    """Check only geometry metadata relevant to peripheral positions, no fields."""
    paths = sorted((ROOT / "runs/round3_predict/test").glob("r3_*/parameters.json"))
    prior = [(p.as_posix(), read_json(p)) for p in paths]
    for (_, _), (c, folder) in old_development_cases().items():
        prior.append((folder.as_posix(), c))
    collisions = []
    for c in new_cases:
        for source, old in prior:
            if any(abs(float(x)) > 1e-14 for x in old.get("interface_amplitudes_mm", [0])):
                continue
            keys = ["L_mm", "a_mm", "b_mm", "theta_deg", "cx_mm", "cz_mm"]
            if any(old.get(k) != c[k] for k in keys):
                continue
            if exact_changed_measure_mm(old, c, 0, 120) == 0:
                collisions.append({"new_case": c["id"], "old_geometry_source": source})
    if collisions:
        raise RuntimeError("Seen geometries require prespecified replacements: " + repr(collisions))
    return {
        "scope": "Geometry metadata only for 12 round3 candidate geometries and 8 round4 peripheral geometries; no target responses read.",
        "metadata_files_or_cases_checked": len(prior), "matches": collisions,
        "earlier_round1_round2_basis": "Existing equal-thickness or wave geometry families contain no movable peripheral 2 mm band pair.",
    }


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    template = read_json(PACK / "round5_design_template.json")
    old = old_development_cases()
    layouts = []
    for label, (dl, du) in template["known_development_layouts"].items():
        layouts.append({"layout_id": label, "distance_lower_mm": dl, "distance_upper_mm": du,
                        "lower_low_band_z_mm": [30-dl-2, 30-dl],
                        "upper_high_band_z_mm": [90+du, 92+du], "phase": "development"})
    layouts.extend({**x, "phase": "prospective_initial"} for x in template["new_layouts"])
    interfaces = set(old[("A", "NN")][0]["z_interfaces_mm"])
    for layout in layouts:
        interfaces.update(layout["lower_low_band_z_mm"] + layout["upper_high_band_z_mm"])
    interfaces = sorted(interfaces)
    edges = [0, *interfaces, 120]
    originals = {role: read_json(ROOT / f"runs/2026-10-21_round2/strict_h10_L60_{role}_m3_full/parameters.json") for role in ["A", "B"]}
    cases, geometry_rows, dev_states, new_states = [], [], [], []
    for layout in layouts:
        label = layout["layout_id"]
        for role in ["A", "B"]:
            original = originals[role]
            c = copy.deepcopy(original)
            c.update({
                "id": f"r5_strict_{role}_{label}", "configuration_role": role,
                "layout": label, "layout_id": label,
                "group": "round5_modal_closure", "family_id": "strict_h10_L60_fixed_core_peripheral_band_position",
                "pair_id": f"r5_strict_{label}", "geometry_family": "strict_fixed_core_peripheral_2mm_bands",
                "passive_ellipse": True, "fracture": True, "mesh_level": 3, "load": "z",
                "source_parameters": f"runs/2026-10-21_round2/strict_h10_L60_{role}_m3_full/parameters.json",
                "source_case_id": original["id"], "phase": layout["phase"],
                "distance_lower_mm": layout["distance_lower_mm"], "distance_upper_mm": layout["distance_upper_mm"],
                "lower_low_band_z_mm": layout["lower_low_band_z_mm"], "upper_high_band_z_mm": layout["upper_high_band_z_mm"],
                "lower_high_to_low_strip_z_mm": layout["lower_low_band_z_mm"],
                "upper_low_to_high_strip_z_mm": layout["upper_high_band_z_mm"],
                "z_interfaces_mm": interfaces, "interface_amplitudes_mm": [0] * len(interfaces),
                "actual_partition_thicknesses_mm": [hi-lo for lo, hi in zip(edges[:-1], edges[1:])],
                "protected_core_z_mm": [30, 90], "protected_electrode_buffers_z_mm": [[0, 20], [100, 120]],
                "solver_input_guard": "Background candidate. Full targets only after all 48 prospective predictions are saved.",
                "geometry_only_area_note": "Analytical design values; actual COMSOL integrations must be saved separately.",
                "interpretation_limit": "Fixed strict A/B core with actual low bottom and high top 20 mm buffers; only the two 2 mm peripheral bands move.",
            })
            materials = []
            for z0, z1 in zip(edges[:-1], edges[1:]):
                z = (z0 + z1) / 2
                p = phase(original, z)
                if layout["lower_low_band_z_mm"][0] < z < layout["lower_low_band_z_mm"][1]:
                    assert p == 1
                    p = 0
                if layout["upper_high_band_z_mm"][0] < z < layout["upper_high_band_z_mm"][1]:
                    assert p == 0
                    p = 1
                materials.append(p)
            c["layer_materials"] = materials
            c["background_case_z"], c["background_case_x"] = c["id"] + "_bgz", c["id"] + "_bgx"
            c["physical_interfaces_z_mm"] = physical_interfaces(c)
            c["physical_interface_count"] = len(c["physical_interfaces_z_mm"])
            nn = old[(role, "NN")][0]
            checks = {
                "changed_core_relative_original_mm2": 120 * exact_changed_measure_mm(original, c, 30, 90),
                "changed_lower_buffer_relative_original_mm2": 120 * exact_changed_measure_mm(original, c, 0, 20),
                "changed_upper_buffer_relative_original_mm2": 120 * exact_changed_measure_mm(original, c, 100, 120),
                "background_high_area_analytic_mm2": 120 * integrated_high_mm(c),
                "background_low_area_analytic_mm2": 120 * (120-integrated_high_mm(c)),
                "lower_periphery_high_area_analytic_mm2": 120 * integrated_high_mm(c, 0, 30),
                "upper_periphery_high_area_analytic_mm2": 120 * integrated_high_mm(c, 90, 120),
                "ellipse_area_analytic_mm2": math.pi*c["a_mm"]*c["b_mm"],
                "ellipse_replaced_high_analytic_mm2": c["expected_replaced_high_area_mm2"],
                "ellipse_replaced_low_analytic_mm2": c["expected_replaced_low_area_mm2"],
                "physical_interface_count": c["physical_interface_count"],
                "core_background_V_max_difference_vs_NN_analytic_V": max(abs(potential_at(c, z)-potential_at(nn, z)) for z in [30, 31, 35, 45, 55, 65, 75, 85, 89, 90]),
                "G0z_analytic_S": .12/sum((hi-lo)*.001/(c["sigma_high_S_m"] if p else c["sigma_low_S_m"]) for lo, hi, p in zip(edges[:-1], edges[1:], materials)),
                "sigma_eff_x_analytic_S_m": .055,
            }
            assert all(checks[k] == 0 for k in list(checks)[:3])
            assert checks["background_high_area_analytic_mm2"] == checks["background_low_area_analytic_mm2"] == 7200
            assert checks["lower_periphery_high_area_analytic_mm2"] == 960
            assert checks["upper_periphery_high_area_analytic_mm2"] == 2640
            assert checks["core_background_V_max_difference_vs_NN_analytic_V"] < 1e-14
            assert c["physical_interface_count"] == nn["physical_interface_count"]
            c["geometry_only_checks"] = checks
            if layout["phase"] == "development":
                oldcase, folder = old[(role, label)]
                assert exact_changed_measure_mm(oldcase, c, 0, 120) == 0
                c["existing_geometry_case_id"] = oldcase["id"]
                c["existing_case_directory"] = folder.as_posix()
            else:
                c["existing_geometry_case_id"] = ""
                c["existing_case_directory"] = ""
            cases.append(c)
            geometry_rows.append({"case_id": c["id"], "role": role, "layout": label, "phase": layout["phase"],
                                  "distance_lower_mm": layout["distance_lower_mm"], "distance_upper_mm": layout["distance_upper_mm"], **checks})
            points = template["known_development_sigma_f_S_m"] if layout["phase"] == "development" else template["new_sigma_f_S_m"]
            for sf in points:
                state = {"state_id": c["id"] + "_sf_" + sigma_tag(sf), "case_id": c["id"], "role": role,
                         "layout": label, "sigma_f_S_m": sf, "phase": layout["phase"],
                         "target_relative_folder": f"truth/{c['id']}/sf_{sigma_tag(sf)}"}
                if layout["phase"] == "development":
                    state["existing_target_directory"] = (old[(role, label)][1] / ("sf_"+sigma_tag(sf))).as_posix()
                    dev_states.append(state)
                else:
                    new_states.append(state)
    prior_check = limited_prior_geometry_check([c for c in cases if c["phase"] == "prospective_initial"])
    by = {(c["configuration_role"], c["layout"]): c for c in cases}
    port_rows, layer_rows, comparisons, side_specs = [], [], [], {}
    for layout in layouts:
        label = layout["layout_id"]
        layout["actual_layers_core_to_electrode"] = {}
        layout["Y_actual_S_m2"] = {}
        layout["Y_reference_S_m2"] = {}
        for side in ["lower", "upper"]:
            c = by[("A", label)]
            layers = core_outward_layers(c, side)
            assert layers == core_outward_layers(by[("B", label)], side)
            actual = modal_admittances(layers, .12, 64)
            reference = equivalent_modal_admittances(layers, .12, 64)
            h, sx, sz = equivalent_tensor(layers)
            side_specs[side] = {"height_m": h, "sigma_xx_S_m": sx, "sigma_zz_S_m": sz,
                                "Y0_S_m2": float(actual[0]), "tensor_coordinates": "x,z"}
            serial = [{"thickness_m": x.thickness_m, "conductivity_S_m": x.conductivity_S_m} for x in layers]
            layout["actual_layers_core_to_electrode"][side] = serial
            layout["Y_actual_S_m2"][side] = actual.tolist()
            layout["Y_reference_S_m2"][side] = reference.tolist()
            for index, layer in enumerate(layers, 1):
                layer_rows.append({"layout": label, "side": side, "order_from_core": index,
                                   "thickness_m": layer.thickness_m, "conductivity_S_m": layer.conductivity_S_m})
            for n in range(65):
                yold = admittance(layers_from_electrode(c, side), n)
                comparisons.append(abs(actual[n]/yold - 1))
                if n <= 8:
                    comparisons.append(abs(actual[n]/original_transfer_check(layers, n*math.pi/.12) - 1))
                port_rows.append({"layout": label, "side": side, "n": n, "k_per_m": n*math.pi/.12,
                                  "Y_actual_S_m2": actual[n], "Y_reference_S_m2": reference[n],
                                  "delta_Y_S_m2": actual[n]-reference[n], "Y_units": "S/m^2",
                                  "norm_N_m": .12 if n == 0 else .06, "width_m": .12,
                                  "out_of_plane_thickness_m": 1, "M0_sigma_xx_S_m": sx,
                                  "M0_sigma_zz_S_m": sz, "electrode_V": 0 if side == "lower" else 1,
                                  "layer_order": "core_to_electrode", "layer_source_case_id": c["id"],
                                  "scope": "analytical_from_actual_geometry_not_COMSOL_result"})
    assert max(comparisons) < 1e-13
    for side in ["lower", "upper"]:
        first = layouts[0]["Y_reference_S_m2"][side]
        assert all(max(abs(a/b-1) for a,b in zip(x["Y_reference_S_m2"][side], first)) < 1e-14 for x in layouts)
    selftest = run_reference_tests()
    validation = {
        "scope": "Analytical exterior and candidate geometry checks only, no new COMSOL solve or target field read.",
        "nmax": 64, "max_relative_difference_old_impedance_and_independent_transfer": max(comparisons),
        "taskpack_reference_selftest": selftest,
        "actual_Y0_lower_S_m2": side_specs["lower"]["Y0_S_m2"],
        "actual_Y0_upper_S_m2": side_specs["upper"]["Y0_S_m2"],
        "taskpack_illustrative_Y0_lower_S_m2": 1/(.028/.1+.002/.01),
        "taskpack_illustrative_Y0_upper_S_m2": 1/(.028/.01+.002/.1),
        "taskpack_discrepancy": "Illustration treats all 30 mm as the host plus band. Actual models retain bottom 20 mm LOW and top 20 mm HIGH buffers.",
        "prior_geometry_check": prior_check,
    }
    design = {
        "schema": "round5_modal_closure_geometry_v1", "status": "design_and_analytical_coefficients_only",
        "design_date_Asia_Shanghai": "2026-10-26", "root": ROOT.as_posix(),
        "source_taskpack_design": "round5_taskpack/round5_design_template.json",
        "actual_geometry_authority": "Existing strict_h10_L60_A/B and round4 NN/NF/FN/FF parameters, including actual electrode buffers.",
        "common_z_interfaces_mm": interfaces, "core_z_mm": [30, 90],
        "out_of_plane_thickness_m": 1, "sample_width_mm": 120, "sample_height_mm": 120,
        "M0_equivalent_tensor": side_specs, "layout_specs": layouts,
        "cases": cases, "development_states": dev_states, "new_states": new_states,
        "development_state_count": len(dev_states), "new_initial_target_count": len(new_states),
        "mode_order_candidates": template["mode_order_candidates"], "reference_order_checks": [32, 64],
        "high_mode_tail": template["high_mode_tail"], "prediction_before_truth": True,
        "all_48_predictions_required_before_first_new_target": True,
        "limited_crossover_validation_layouts": template["limited_crossover_validation_layouts"],
        "crossover_sigma_interval_S_m": template["crossover_sigma_interval_S_m"],
        "prior_geometry_check": prior_check,
    }
    assert len(dev_states) == 24 and len(new_states) == 48
    write_json(DEST / "design.json", design)
    write_json(DEST / "analytical_validation.json", validation)
    write_csv(DEST / "port_coefficients.csv", port_rows)
    write_csv(DEST / "port_layers.csv", layer_rows)
    write_csv(DEST / "geometry_candidates.csv", geometry_rows)
    (DEST / "geometry_note_CN.md").write_text(GEOMETRY_NOTE, encoding="utf-8")
    print(json.dumps({"design": str(DEST/"design.json"), "cases": len(cases), "development_states": len(dev_states),
                      "prospective_states": len(new_states), "coefficient_rows": len(port_rows),
                      "M0": side_specs, "analytic_checks": "passed", "reference_tests": selftest}, ensure_ascii=False))


GEOMETRY_NOTE = """# 第五轮真实外围几何与任务包示例的区别

日期：2026-10-26。本文件与设计表、解析系数不是新的COMSOL结果；生成脚本不读取目标场、不启动求解、不改动第三轮或旧模型。

以本机已算strict A/B、NN/NF/FN/FF参数为准：下方0–20 mm始终为低导缓冲，上方100–120 mm始终为高导缓冲。只有20–30 mm高导宿主内2 mm低导带、90–100 mm低导宿主内2 mm高导带移动。任务书要求保留这些缓冲，因此不能采用数学说明中将全部30 mm当成同一宿主的示例数值。

每侧仍为30 mm：下外围高导8 mm、低导22 mm；上外围高导22 mm、低导8 mm。正确均匀响应为：

- 下侧Y0 = 1/(0.008/0.1 + 0.022/0.01) = **0.438596491228 S/m²**。
- 上侧Y0 = 1/(0.022/0.1 + 0.008/0.01) = **0.980392156863 S/m²**。
- 下侧M0张量：σx=0.034 S/m、σz=0.013157894737 S/m。
- 上侧M0张量：σx=0.076 S/m、σz=0.029411764706 S/m。

任务包的2.083333与0.354610 S/m²适用于其简化示例，不适用于本机这组保持缓冲的完整外围。参考函数的递推仍适用，只需传入真实层序。原任务包保持原样，区别在本文件及analytical_validation.json明确记录。

design.json有20个完整背景候选参数（4旧布局×A/B、6新布局×A/B）、24个开发状态与48个新初批状态。cases可交现有solve_round2_case建几何，passive_ellipse默认true；完整含裂缝状态由控制脚本在全批预测保存后另行建立。development_states指向旧结果目录，供开发对照；new_states没有目标结果。layout_specs内包含真实核心向电极层序，以及n=0…64的真实/参考系数数组。外围系数A/B相同，核心材料图保持各自原样。

所有候选共用完整几何分界的并集，补入新带端点；同相分界不算物理材料界面。核心30–90 mm原分界不变。物理界面数保持A为7、B为11。每侧相面积、总相面积、椭圆替换相面积、核心与缓冲不变性按解析几何检查；实际COMSOL积分必须由求解另存，不能用本表替代。

port_coefficients.csv按layout、side、n记录真实Y、共同M0参考Y、带符号差、归一长度N、单位及层序来源。n=0被完整保留。n>K及n>64必须由实际M0参考外围承担，不得照CSV截断为绝缘。port_layers.csv层序方向始终为核心→电极；下侧因此与全局z递增方向相反。

已重新运行任务包16项解析/代数自测，输出在reference_selftest_output.txt与reference_selftest.json；另对真实缓冲几何系数检查了旧第四轮阻抗递推和独立低阶传递矩阵。以上只验证解析系数，不证明新组合模型或目标COMSOL计算成功。

复算入口（现有含NumPy环境）：`D:/Users/PC/anaconda3/envs/efield310/python.exe scripts/r5_design.py`。旧几何仅检查与新布局相关的参数元数据，未读旧目标场或重算旧实验。新布局均未与这部分既有几何重复。
"""


if __name__ == "__main__":
    main()
