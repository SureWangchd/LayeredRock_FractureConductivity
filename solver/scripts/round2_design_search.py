"""Analytic 2D geometry design only; never launches COMSOL or selects by response.

Examples:
  python scripts/round2_design_search.py --audit-equal-tiles
  python scripts/round2_design_search.py --core-json PATH --strict-json PATH
Output files are written only when their paths are explicitly requested.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path


def vertical_radius(a: float, b: float = 1.5) -> float:
    return math.hypot(a * math.sin(math.pi / 3), b * math.cos(math.pi / 3))


def ellipse_area_below(z: float, a: float, b: float = 1.5) -> float:
    """Exact area in mm^2 below horizontal z, ellipse center z=60 mm."""
    radius = vertical_radius(a, b)
    u = max(-1.0, min(1.0, (z - 60.0) / radius))
    return a * b * (math.asin(u) + u * math.sqrt(max(0.0, 1 - u * u)) + math.pi / 2)


def central_area(halfwidth: float, a: float, b: float = 1.5) -> float:
    return ellipse_area_below(60 + halfwidth, a, b) - ellipse_area_below(60 - halfwidth, a, b)


def phase_metrics(interfaces: list[float], materials: list[int], a: float) -> dict:
    bounds = [0.0, *interfaces, 120.0]
    assert len(bounds) == len(materials) + 1
    thicknesses = [hi - lo for lo, hi in zip(bounds, bounds[1:])]
    assert min(thicknesses) > 0
    high_thickness = sum(t * m for t, m in zip(thicknesses, materials))
    replaced_high = sum(
        ellipse_area_below(hi, a) - ellipse_area_below(lo, a)
        for lo, hi, mat in zip(bounds, bounds[1:], materials)
        if mat == 1
    )
    # Merge equal-material adjacent tiles: these are not physical interfaces.
    runs = []
    for lo, hi, mat in zip(bounds, bounds[1:], materials):
        if runs and runs[-1]["material"] == mat:
            runs[-1]["z_max_mm"] = hi
        else:
            runs.append({"z_min_mm": lo, "z_max_mm": hi, "material": mat})
    return {
        "actual_partition_thicknesses_mm": thicknesses,
        "physical_material_runs": runs,
        "expected_background_high_area_mm2": 120 * high_thickness,
        "expected_replaced_high_area_mm2": replaced_high,
        "expected_replaced_low_area_mm2": math.pi * a * 1.5 - replaced_high,
    }


def endcap_geometry(a: float, b: float = 1.5) -> dict:
    """Fixed last 3 mm of the long-axis coordinate, full transverse extent."""
    u = 1 - 3.0 / a
    minimum_abs_z = (a - 3) * math.sin(math.pi / 3) - b * math.cos(math.pi / 3) * math.sqrt(1 - u * u)
    return {
        "endcap_definition": "Each finite endcap is s in [a-3 mm,a] or [-a,-a+3 mm], including the full ellipse transverse section.",
        "endcap_axial_length_mm": 3,
        "endcap_area_each_mm2": a * b * (math.acos(u) - u * math.sqrt(1 - u * u)),
        "endcap_min_abs_z_offset_mm": minimum_abs_z,
    }


def pure_cap_states(interfaces: list[float], materials: list[int], a: float) -> tuple:
    """0/1 if an entire finite cap has that material; 'mixed' otherwise."""
    bounds = [0.0, *interfaces, 120.0]
    v = endcap_geometry(a)["endcap_min_abs_z_offset_mm"]
    radius = vertical_radius(a)
    result = []
    for cap_lo, cap_hi in [(60 - radius, 60 - v), (60 + v, 60 + radius)]:
        seen = {
            mat for lo, hi, mat in zip(bounds, bounds[1:], materials)
            if min(cap_hi, hi) - max(cap_lo, lo) > 1e-10
        }
        result.append(next(iter(seen)) if len(seen) == 1 else "mixed")
    return tuple(result)


def common_case(case_id: str, pair: str, h: int, length: int, role: str) -> dict:
    return dict(
        id=case_id, pair_id=pair, configuration_role=role,
        h_mm=h, h_is_nominal=True, a_mm=length / 2, b_mm=1.5,
        theta_deg=60, cx_mm=60, cz_mm=60, L_mm=120,
        full_length_mm=length, A_mm=0,
        design_basis="Prescribed geometry and phase areas only; no conductivity-response target or screening.",
    )


def core_candidates() -> list[dict]:
    cases = []
    for h, length in [(10, 20), (15, 30), (20, 40), (10, 40), (10, 60), (20, 60)]:
        a = length / 2
        radius = vertical_radius(a)
        assert h < 2 * radius
        envelope = radius + 1
        target_high = math.pi * a * 1.5 - central_area(h / 2, a)
        lo, hi = 0.0, radius
        for _ in range(90):
            mid = (lo + hi) / 2
            if central_area(mid, a) < target_high:
                lo = mid
            else:
                hi = mid
        q = (lo + hi) / 2
        pair = f"core_h{h}_L{length}"
        for role in ["A", "B"]:
            if role == "A":
                extent, remote_end = q, 40 - q
                materials = [0, 1, 0, 0, 1, 0, 0, 1, 0]
                description, expected_cap = "central_high_endcaps_low", 0
            else:
                extent, remote_end = h / 2, 40 - envelope + h / 2
                materials = [0, 1, 0, 1, 0, 1, 0, 1, 0]
                description, expected_cap = "central_low_endcaps_high", 1
            interfaces = [10, remote_end, 60 - envelope, 60 - extent,
                          60 + extent, 60 + envelope, 120 - remote_end, 110]
            case = common_case(pair + "_" + role, pair, h, length, role)
            case.update(phase_metrics(interfaces, materials, a))
            case.update(endcap_geometry(a))
            case.update(
                z_interfaces_mm=interfaces, layer_materials=materials,
                configuration_description=description,
                central_actual_thickness_mm=2 * extent,
                central_thickness_matches_nominal_h=role == "B",
                buffer_bottom_mm=10, buffer_top_mm=10, buffer_material=0,
                expected_endcap_high_fraction_lower=expected_cap,
                expected_endcap_high_fraction_upper=expected_cap,
                ellipse_area_analytic_mm2=math.pi * a * 1.5,
                ellipse_vertical_halfheight_mm=radius,
                local_envelope_halfheight_mm=envelope,
                area_matching_analytic_abs_error_mm2=abs(case["expected_replaced_high_area_mm2"] - target_high),
                interpretation_limit="h is the exact central low-layer thickness in B. A central high layer and remote compensating layers differ; this is not a uniform-layer-thickness sweep.",
            )
            assert abs(case["expected_background_high_area_mm2"] - 7200) < 1e-8
            assert case["area_matching_analytic_abs_error_mm2"] < 1e-8
            assert pure_cap_states(interfaces, materials, a) == (expected_cap, expected_cap)
            cases.append(case)
    return cases


def strict_candidates() -> list[dict]:
    sequences = {
        10: {"A": [0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1],
             "B": [0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1]},
        20: {"A": [0, 1, 1, 0, 0, 1], "B": [0, 1, 0, 1, 0, 1]},
    }
    cases = []
    for h, roles in sequences.items():
        for role, materials in roles.items():
            interfaces = list(range(h, 120, h))
            pair = f"strict_h{h}_L60"
            case = common_case(pair + "_" + role, pair, h, 60, role)
            case.update(phase_metrics(interfaces, materials, 30))
            case.update(endcap_geometry(30))
            case.update(
                z_interfaces_mm=interfaces, layer_materials=materials,
                buffer_bottom_mm=20, buffer_top_mm=20,
                buffer_bottom_material=0, buffer_top_material=1,
                expected_endcap_high_fraction_lower=1,
                expected_endcap_high_fraction_upper=0,
                expected_completely_crossed_low_material_runs=0 if role == "A" else (2 if h == 10 else 1),
                interpretation_limit="Adjacent equal-material tiles have no physical conductivity interface. h10_A and h20_A are the same material distribution with different partitions, not independent thickness evidence.",
            )
            assert abs(case["expected_background_high_area_mm2"] - 7200) < 1e-8
            assert abs(case["expected_replaced_high_area_mm2"] - math.pi * 30 * 1.5 / 2) < 1e-8
            assert pure_cap_states(interfaces, materials, 30) == (1, 0)
            cases.append(case)
    return cases


def equal_tile_audit() -> list[dict]:
    """Exhaustive finite search, fixed bottom-low/top-high buffers >=20 mm."""
    results = []
    for h in [10, 15, 20]:
        n = 120 // h
        buffered_tiles = math.ceil(20 / h)
        interior = list(range(buffered_tiles, n - buffered_tiles))
        high_needed = n // 2 - buffered_tiles
        for length in [20, 30, 40, 60]:
            a = length / 2
            interfaces = list(range(h, 120, h))
            patterns = []
            for selected in itertools.combinations(interior, high_needed):
                high = set(selected) | set(range(n - buffered_tiles, n))
                sequence = tuple(int(i in high) for i in range(n))
                area = phase_metrics(interfaces, list(sequence), a)["expected_replaced_high_area_mm2"]
                patterns.append((sequence, area, pure_cap_states(interfaces, list(sequence), a)))
            strict = ordered_cap_change = pure_hh_ll = 0
            for first, second in itertools.combinations(patterns, 2):
                if abs(first[1] - second[1]) > 1e-8:
                    continue
                if second[0] == tuple(reversed(first[0])):
                    continue
                strict += 1
                ordered_cap_change += first[2] != second[2]
                pure_hh_ll += {first[2], second[2]} == {(0, 0), (1, 1)}
            hh = [p[1] for p in patterns if p[2] == (1, 1)]
            ll = [p[1] for p in patterns if p[2] == (0, 0)]
            min_gap = min((abs(x - y) for x in hh for y in ll), default=None)
            radius = vertical_radius(a)
            results.append(dict(
                h_mm=h, full_length_mm=length, pattern_count=len(patterns),
                buffer_each_mm=buffered_tiles * h,
                strict_area_tolerance_mm2=1e-8,
                strict_nonmirror_pair_count=strict,
                strict_pairs_with_changed_ordered_cap_states=ordered_cap_change,
                strict_pure_HH_vs_LL_pair_count=pure_hh_ll,
                min_HH_vs_LL_area_gap_mm2=min_gap,
                ellipse_only_in_two_center_tiles=radius < h,
                note="This exhaustive result applies only to these equal tiles and fixed buffers; it is not a general impossibility theorem for arbitrary interfaces.",
            ))
    return results


def write_cases(path: str, cases: list[dict]) -> None:
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    payload = dict(status="geometry_candidates_not_COMSOL_results", cases=cases)
    dest.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-json")
    parser.add_argument("--strict-json")
    parser.add_argument("--audit-equal-tiles", action="store_true")
    args = parser.parse_args()
    if args.core_json:
        write_cases(args.core_json, core_candidates())
    if args.strict_json:
        write_cases(args.strict_json, strict_candidates())
    if args.audit_equal_tiles:
        print(json.dumps(equal_tile_audit(), ensure_ascii=False, indent=2))
    if not any((args.core_json, args.strict_json, args.audit_equal_tiles)):
        print(json.dumps(dict(core_cases=core_candidates(), strict_cases=strict_candidates()), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
