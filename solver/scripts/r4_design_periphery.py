"""Geometry-only round-4 peripheral interchange; never calls COMSOL.

Preserves each original strict A/B material map throughout z=30..90 mm,
and preserves the original 20 mm electrode buffers. Analytical quantities
below are candidate-design checks, not measured COMSOL areas or results.
"""
from __future__ import annotations

import bisect
import copy
import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "runs" / "round4_crossover"
ADDED_PARTITIONS = [23, 25, 27, 93, 95, 97]
LAYOUTS = {
    "swap_near": {"lower": [25, 27], "upper": [93, 95], "neutral_shift_mm": 1},
    "swap_far": {"lower": [23, 25], "upper": [95, 97], "neutral_shift_mm": -1},
}


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def phase(c, z):
    return c["layer_materials"][bisect.bisect_right(c["z_interfaces_mm"], z)]


def physical_interfaces(c):
    p = c["layer_materials"]
    return [z for i, z in enumerate(c["z_interfaces_mm"]) if p[i] != p[i + 1]]


def integrated_high_mm(c, lower=0, upper=120):
    edges = [0] + c["z_interfaces_mm"] + [120]
    return sum(max(0, min(hi, upper) - max(lo, lower)) * p
               for lo, hi, p in zip(edges[:-1], edges[1:], c["layer_materials"]))


def exact_changed_measure_mm(original, candidate, lower, upper):
    edges = sorted({lower, upper, *original["z_interfaces_mm"], *candidate["z_interfaces_mm"]})
    edges = [z for z in edges if lower <= z <= upper]
    return sum(hi - lo for lo, hi in zip(edges[:-1], edges[1:])
               if phase(original, (lo + hi) / 2) != phase(candidate, (lo + hi) / 2))


def potential_at(c, z):
    # Top=1 V, bottom=0 V; exact one-dimensional background potential.
    edges = [0] + c["z_interfaces_mm"] + [120]
    resistance = sum((hi - lo) / (0.1 if p else 0.01)
                     for lo, hi, p in zip(edges[:-1], edges[1:], c["layer_materials"]))
    below = sum(max(0, min(hi, z) - lo) / (0.1 if p else 0.01)
                for lo, hi, p in zip(edges[:-1], edges[1:], c["layer_materials"]) if lo < z)
    return below / resistance


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    angle = math.radians(60)
    rz = math.hypot(30 * math.sin(angle), 1.5 * math.cos(angle))
    envelope = [60 - rz, 60 + rz]
    cap_upper = 60 - 27 * math.sin(angle) + 1.5 * math.sqrt(1 - .9 ** 2) * math.cos(angle)
    cap_area = 30 * 1.5 * (math.pi / 2 - .9 * math.sqrt(1 - .9 ** 2) - math.asin(.9))
    ellipse_area = math.pi * 30 * 1.5
    cases, geometry_rows = [], []
    for role in ["A", "B"]:
        source = ROOT / "runs" / "2026-10-21_round2" / f"strict_h10_L60_{role}_m3_full" / "parameters.json"
        original = read_json(source)
        assert original["L_mm"] == 120 and original["a_mm"] == 30 and original["b_mm"] == 1.5
        assert original["theta_deg"] == 60 and original["cx_mm"] == 60 and original["cz_mm"] == 60
        assert integrated_high_mm(original) == 60
        assert phase(original, 25) == 1 and phase(original, 95) == 0
        common_interfaces = sorted(set(original["z_interfaces_mm"] + ADDED_PARTITIONS))
        for layout, spec in LAYOUTS.items():
            c = copy.deepcopy(original)
            c["id"] = f"r4_strict_{role}_{layout}"
            c["source_parameters"] = source.relative_to(ROOT).as_posix()
            c["source_case_id"] = original["id"]
            c["parent_full_case"] = original["id"]
            c["group"] = "round4_periphery_candidates"
            c["family_id"] = "strict_h10_L60_peripheral_interchange"
            c["geometry_family"] = "strict_internal_map_fixed_upper_lower_peripheral_exchange"
            c["pair_id"] = f"r4_strict_{layout}"
            c["layout"] = layout
            c["configuration_role"] = role
            c["passive_ellipse"] = True
            c["fracture"] = True
            c["load"] = "z"
            c["mesh_level"] = 3
            c["h_is_nominal"] = True
            c["A_mm"] = 0
            c["z_interfaces_mm"] = common_interfaces
            c["interface_amplitudes_mm"] = [0] * len(common_interfaces)
            edges = [0] + common_interfaces + [120]
            phases = []
            for lo, hi in zip(edges[:-1], edges[1:]):
                mid = (lo + hi) / 2
                p = phase(original, mid)
                if spec["lower"][0] < mid < spec["lower"][1]:
                    assert p == 1
                    p = 0
                if spec["upper"][0] < mid < spec["upper"][1]:
                    assert p == 0
                    p = 1
                phases.append(p)
            c["layer_materials"] = phases
            c["actual_partition_thicknesses_mm"] = [hi - lo for lo, hi in zip(edges[:-1], edges[1:])]
            c["background_case_z"] = c["id"] + "_bgz"
            c["background_case_x"] = c["id"] + "_bgx"
            c["candidate_state"] = "geometry_only_not_COMSOL_solved"
            c["solver_input_guard"] = "passive_ellipse=true; full targets require separate post-prediction entry"
            c["design_description_CN"] = "上下外围等面积交换；近/远仅比较两个已交换布局的成对位置，原始布局没有孤立可移高导带。"
            c["interpretation_limit"] = (
                "Original-to-candidate exchanges high/low material between disjoint peripheral regions and adds four physical interfaces. "
                "Near-to-far translates the lower low strip and upper high strip by 2 mm away from the protected core. "
                "It is not a 1-2 mm translation of an isolated band in the original geometry."
            )
            c["protected_core_z_mm"] = [30, 90]
            c["protected_electrode_buffers_z_mm"] = [[0, 20], [100, 120]]
            c["ellipse_vertical_envelope_mm"] = envelope
            c["protected_endcap_neighborhood_radius_mm"] = 3
            c["lower_high_to_low_strip_z_mm"] = spec["lower"]
            c["upper_low_to_high_strip_z_mm"] = spec["upper"]
            c["strip_x_mm"] = [0, 120]
            c["strip_width_z_mm"] = 2
            c["exchange_area_each_analytic_mm2"] = 240
            c["background_high_area_analytic_mm2"] = 7200
            c["background_low_area_analytic_mm2"] = 7200
            c["ellipse_area_analytic_mm2"] = ellipse_area
            c["ellipse_replaced_high_area_analytic_mm2"] = ellipse_area / 2
            c["ellipse_replaced_low_area_analytic_mm2"] = ellipse_area / 2
            c["finite_endcap_area_each_analytic_mm2"] = cap_area
            c["finite_endcap_z_ranges_mm"] = [[envelope[0], cap_upper], [120 - cap_upper, envelope[1]]]
            c["physical_interfaces_z_mm"] = physical_interfaces(c)
            c["physical_interface_count"] = len(c["physical_interfaces_z_mm"])
            c["physical_interfaces_added_relative_to_original"] = 4
            c["additional_material_interface_length_analytic_mm"] = 480
            c["virtual_neutral_lower_strip_z_mm"] = [24, 26]
            c["virtual_neutral_upper_strip_z_mm"] = [94, 96]
            c["lower_shift_from_virtual_neutral_mm"] = spec["neutral_shift_mm"]
            c["upper_shift_from_virtual_neutral_mm"] = -spec["neutral_shift_mm"]
            c["virtual_neutral_is_not_original_and_is_not_scheduled"] = True
            c["cross_periphery_exchange_center_distance_mm"] = (sum(spec["upper"]) - sum(spec["lower"])) / 2
            c["background_sigma_eff_z_analytic_S_m"] = 1 / (.5 / .01 + .5 / .1)
            c["background_sigma_eff_x_analytic_S_m"] = .5 * .01 + .5 * .1
            c["background_core_potential_offset_analytic_V"] = potential_at(c, 60) - potential_at(original, 60)
            c["background_core_field_and_endcap_difference_analytic_unchanged"] = True
            checks = {
                "changed_core_area_analytic_mm2": 120 * exact_changed_measure_mm(original, c, 30, 90),
                "changed_bottom_buffer_area_analytic_mm2": 120 * exact_changed_measure_mm(original, c, 0, 20),
                "changed_top_buffer_area_analytic_mm2": 120 * exact_changed_measure_mm(original, c, 100, 120),
                "background_high_area_difference_analytic_mm2": 120 * (integrated_high_mm(c) - integrated_high_mm(original)),
                "changed_material_total_area_analytic_mm2": 120 * exact_changed_measure_mm(original, c, 0, 120),
                "min_distance_changed_strip_to_protected_core_mm": min(30 - spec["lower"][1], spec["upper"][0] - 90),
                "min_vertical_clearance_changed_strip_to_ellipse_mm": min(envelope[0] - spec["lower"][1], spec["upper"][0] - envelope[1]),
                "min_distance_changed_strip_to_electrode_buffer_mm": min(spec["lower"][0] - 20, 100 - spec["upper"][1]),
            }
            assert all(checks[k] == 0 for k in list(checks)[:4])
            assert checks["changed_material_total_area_analytic_mm2"] == 480
            assert len(c["layer_materials"]) == len(c["z_interfaces_mm"]) + 1
            assert len(c["interface_amplitudes_mm"]) == len(c["z_interfaces_mm"])
            assert min(c["actual_partition_thicknesses_mm"]) >= 2
            assert c["physical_interface_count"] - len(physical_interfaces(original)) == 4
            assert abs(c["background_core_potential_offset_analytic_V"] - 3 / 110) < 1e-14
            for z in [31, 35, 45, 55, 65, 75, 85, 89]:
                assert abs(potential_at(c, z) - potential_at(original, z) - 3 / 110) < 1e-14
            c["geometry_only_checks"] = checks
            c["geometry_only_area_note"] = "All areas here are analytic design quantities; actual COMSOL integrals must be saved separately."
            cases.append(c)
            geometry_rows.append({
                "case_id": c["id"], "configuration_role": role, "layout": layout,
                "quantity_status": "analytic_candidate_geometry_not_COMSOL_measurement",
                "lower_high_to_low_z_start_mm": spec["lower"][0], "lower_high_to_low_z_end_mm": spec["lower"][1],
                "upper_low_to_high_z_start_mm": spec["upper"][0], "upper_low_to_high_z_end_mm": spec["upper"][1],
                "exchange_area_each_analytic_mm2": 240,
                "background_high_area_analytic_mm2": 7200, "background_low_area_analytic_mm2": 7200,
                "ellipse_replaced_high_area_analytic_mm2": ellipse_area / 2,
                "ellipse_replaced_low_area_analytic_mm2": ellipse_area / 2,
                "physical_interface_count": c["physical_interface_count"],
                "added_material_interface_length_analytic_mm": 480,
                "cross_periphery_exchange_center_distance_mm": c["cross_periphery_exchange_center_distance_mm"],
                "background_core_potential_offset_analytic_V": c["background_core_potential_offset_analytic_V"],
                **checks,
            })
    manifest = {
        "schema": "round4_peripheral_exchange_candidates_v1",
        "status": "geometry_only_pre_solve",
        "design_date": "2026-10-24",
        "case_count": len(cases),
        "original_reference_case_ids": [f"strict_h10_L60_{r}_m3_full" for r in ["A", "B"]],
        "original_reference_is_separate_from_four_new_candidates": True,
        "physical_design_basis": "Pure-phase same-side translation is impossible; task-authorized upper/lower equal-area peripheral exchange is feasible.",
        "no_target_conductivity_results_used_for_design": True,
        "same_geometry_partitions_added_to_all_candidates_mm": ADDED_PARTITIONS,
        "prediction_sigma_points": "To be selected after original-crossover localization and saved before new target solves; no new targets selected here.",
        "cases": cases,
    }
    (DEST / "periphery_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    with (DEST / "periphery_geometry_candidates.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(geometry_rows[0]))
        writer.writeheader()
        writer.writerows(geometry_rows)
    print(json.dumps({"cases": [c["id"] for c in cases], "geometry_checks": "passed_analytic_only",
                      "ellipse_vertical_envelope_mm": envelope, "finite_endcap_area_each_analytic_mm2": cap_area}, ensure_ascii=False))


if __name__ == "__main__":
    main()
