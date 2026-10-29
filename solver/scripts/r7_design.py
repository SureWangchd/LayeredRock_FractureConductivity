# -*- coding: utf-8 -*-
"""Resolve Round 7 finite-3D geometry from the actual Round 5 layer maps.

This script reads geometry/configuration only.  It never reads any target
conductance, field, or fitting result.  Volumes below are analytical CAD-design
values, and must be distinguished from the subsequently integrated COMSOL ones.
"""
from __future__ import annotations

import copy
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/round7_3d_ct"
R5 = ROOT / "runs/round5_modal_closure/design.json"
LAYOUTS = ("P11", "P17", "P71", "P77")
SIGMAS = (0.18, 10.0)


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def physical_runs(c):
    runs = []
    for lo, hi, phase in zip([0.0, *c["z_interfaces_mm"]], [*c["z_interfaces_mm"], c["L_mm"]], c["layer_materials"]):
        if runs and runs[-1]["phase"] == phase:
            runs[-1]["z_max_mm"] = hi
        else:
            runs.append({"z_min_mm": lo, "z_max_mm": hi, "phase": phase})
    return runs


def volume_cdf_z(z, c):
    """Exact volume of a rotated ellipsoid lying below z, in mm^3."""
    radius = math.hypot(c["a_mm"] * c["basis_u"][2], c["c_mm"] * c["basis_n"][2])
    t = min(1.0, max(-1.0, (z-c["cz_mm"])/radius))
    volume = 4.0*math.pi*c["a_mm"]*c["b_mm"]*c["c_mm"]/3.0
    return volume*(0.5+0.75*t-0.25*t*t*t)


def cap_geometry(c, sign):
    """Finite u-endcaps; extrema prove phase labels for this prescribed design."""
    a, b, thin = c["a_mm"], c["b_mm"], c["c_mm"]
    q = 1-c["endcap_axial_length_mm"]/a
    t0, t1 = ((q, 1.0) if sign > 0 else (-1.0, -q))
    uz, nz = c["basis_u"][2], c["basis_n"][2]
    hz = math.hypot(a*uz, thin*nz)
    ts = [t0, t1]
    for critical in [a*uz/hz, -a*uz/hz]:
        if t0 <= critical <= t1:
            ts.append(critical)
    zs = [c["cz_mm"]+a*uz*t+s*thin*nz*math.sqrt(max(0.0, 1-t*t)) for t in ts for s in [-1, 1]]
    zlo, zhi = min(zs), max(zs)
    phases = {r["phase"] for r in physical_runs(c) if max(r["z_min_mm"], zlo) < min(r["z_max_mm"], zhi)}
    if len(phases) != 1:
        raise ValueError("Prescribed endcap crosses a phase boundary; use finite-volume integration, not tip labels")
    phase = next(iter(phases))
    # Integration of pi*b*c*(1-u^2/a^2) du over a finite axial cap.
    volume = math.pi*a*b*thin*((t1-t1**3/3)-(t0-t0**3/3))
    return {
        "u_interval_mm": [a*t0, a*t1],
        "z_envelope_mm": [zlo, zhi],
        "volume_mm3_analytic": volume,
        "high_fraction_analytic": float(phase),
        "low_fraction_analytic": float(1-phase),
        "phase_label": "high" if phase else "low",
        "classification_basis": "The entire finite 3-mm axial cap lies within one actual material run; not a point-tip classification.",
    }


def geometry_checks(c):
    L = c["L_mm"]
    runs = physical_runs(c)
    phase_background = {phase: L*L*sum(r["z_max_mm"]-r["z_min_mm"] for r in runs if r["phase"] == phase) for phase in (0, 1)}
    phase_replaced = {phase: sum(volume_cdf_z(r["z_max_mm"], c)-volume_cdf_z(r["z_min_mm"], c) for r in runs if r["phase"] == phase) for phase in (0, 1)}
    ellipse_volume = 4*math.pi*c["a_mm"]*c["b_mm"]*c["c_mm"]/3
    half_spans = [math.sqrt(sum((semi*direction[axis])**2 for semi, direction in zip([c["a_mm"], c["b_mm"], c["c_mm"]], [c["basis_u"], c["basis_v"], c["basis_n"]]))) for axis in range(3)]
    centers = [c["cx_mm"], c["cy_mm"], c["cz_mm"]]
    assert all(center-span > 0 and center+span < L for center, span in zip(centers, half_spans))
    assert centers[2]-half_spans[2] > 30 and centers[2]+half_spans[2] < 90
    assert abs(sum(phase_background.values())-L**3) < 1e-8
    assert abs(sum(phase_replaced.values())-ellipse_volume) < 1e-8
    assert all(abs(phase_background[p]-L**3/2) < 1e-8 for p in (0, 1))
    series_resistance_per_area = sum((r["z_max_mm"]-r["z_min_mm"])*1e-3/(c["sigma_high_S_m"] if r["phase"] else c["sigma_low_S_m"]) for r in runs)
    sigma_z = (L*1e-3)/series_resistance_per_area
    sigma_xy = sum((r["z_max_mm"]-r["z_min_mm"])*(c["sigma_high_S_m"] if r["phase"] else c["sigma_low_S_m"]) for r in runs)/L
    return {
        "kind": "analytical geometry only; COMSOL volume integrations remain required",
        "sample_volume_mm3": L**3,
        "ellipsoid_volume_mm3": ellipse_volume,
        "ellipsoid_half_spans_xyz_mm": half_spans,
        "ellipsoid_bounds_xyz_mm": [[center-span, center+span] for center, span in zip(centers, half_spans)],
        "minimum_external_clearance_mm": min(min(center-span, L-center-span) for center, span in zip(centers, half_spans)),
        "minimum_core_interface_clearance_mm": min(centers[2]-half_spans[2]-30, 90-centers[2]-half_spans[2]),
        "background_low_volume_mm3": phase_background[0],
        "background_high_volume_mm3": phase_background[1],
        "replaced_low_volume_mm3": phase_replaced[0],
        "replaced_high_volume_mm3": phase_replaced[1],
        "final_low_matrix_volume_mm3": phase_background[0]-phase_replaced[0],
        "final_high_matrix_volume_mm3": phase_background[1]-phase_replaced[1],
        "endcap_minus": cap_geometry(c, -1),
        "endcap_plus": cap_geometry(c, 1),
        "sigma_z_background_analytic_S_m": sigma_z,
        "sigma_x_background_analytic_S_m": sigma_xy,
        "sigma_y_background_analytic_S_m": sigma_xy,
        "Gz_background_analytic_S": sigma_z*(L*1e-3),
        "Gx_background_analytic_S": sigma_xy*(L*1e-3),
        "Gy_background_analytic_S": sigma_xy*(L*1e-3),
        "AB_replacement_match_analytic": abs(phase_replaced[0]-phase_replaced[1]) < 1e-8,
    }


def sigma_tag(sigma):
    return f"{sigma:g}".replace(".", "p")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    old = read_json(R5)
    cases = []
    for role in ("A", "B"):
        for layout in LAYOUTS:
            matches = [v for v in old["cases"] if v["configuration_role"] == role and v["layout"] == layout]
            assert len(matches) == 1
            source = matches[0]
            # Copy actual material map, not earlier cartoon layer thicknesses or
            # any response-dependent inferred configuration.
            c = {key: copy.deepcopy(source[key]) for key in [
                "z_interfaces_mm", "layer_materials", "sigma_low_S_m", "sigma_high_S_m",
                "distance_lower_mm", "distance_upper_mm", "lower_low_band_z_mm", "upper_high_band_z_mm",
                "buffer_bottom_mm", "buffer_top_mm", "buffer_bottom_material", "buffer_top_material",
                "protected_core_z_mm", "protected_electrode_buffers_z_mm",
            ]}
            c.update({
                "id": f"r7_strict_{role}_{layout}", "role": role, "configuration_role": role,
                "layout": layout, "layout_id": layout, "kind": "synthetic3d_explicit",
                "dimension": 3, "geometry_input_unit": "mm", "solver_unit": "SI",
                "L_mm": 120.0, "Lx_mm": 120.0, "Ly_mm": 120.0, "Lz_mm": 120.0,
                "cx_mm": 60.0, "cy_mm": 60.0, "cz_mm": 60.0,
                "a_mm": 30.0, "b_mm": 30.0, "c_mm": 1.5, "theta_deg": 60.0,
                "basis_u": [0.5, 0.0, math.sqrt(3)/2],
                "basis_v": [0.0, 1.0, 0.0],
                "basis_n": [-math.sqrt(3)/2, 0.0, 0.5],
                "rotation_convention": "r=center+u*local_u+v*local_v+n*local_n; u rises 60 degrees above +x toward +z; COMSOL Ellipsoid axes [30,30,1.5] mm with Rot y=-60 degrees.",
                "ellipse_equation_local": "(local_u/30[mm])^2+(local_v/30[mm])^2+(local_n/1.5[mm])^2<=1",
                "endcap_axial_length_mm": 3.0, "endcap_length_mm": 3.0,
                "source_design": str(R5.relative_to(ROOT)).replace("\\", "/"),
                "source_case_id": source["id"],
                "load": "z", "deltaV_V": 1.0,
                "electrodes": {"top_z_mm": 120.0, "top_voltage_V": 1.0, "bottom_z_mm": 0.0, "bottom_voltage_V": 0.0, "other_external_faces": "electric insulation"},
                "same_role_peripheral_comparison": "Exactly identical core, ellipsoid and electrode buffers; only two remote 2-mm bands move.",
                "mesh_level": "base", "sigma_f_values_S_m": list(SIGMAS),
            })
            c["physical_material_runs"] = physical_runs(c)
            c["physical_z_interfaces_mm"] = [r["z_max_mm"] for r in c["physical_material_runs"][:-1]]
            c["geometry_checks"] = geometry_checks(c)
            cases.append(c)

    lower = old["M0_equivalent_tensor"]["lower"]
    upper = old["M0_equivalent_tensor"]["upper"]
    m0_tensors = {
        side: {
            "z_range_mm": [0.0, 30.0] if side == "lower" else [90.0, 120.0],
            "sigma_xx_S_m": src["sigma_xx_S_m"],
            "sigma_yy_S_m": src["sigma_xx_S_m"],
            "sigma_zz_S_m": src["sigma_zz_S_m"],
            "off_diagonal_S_m": 0.0,
            "includes_actual_20mm_electrode_buffer": True,
        } for side, src in [("lower", lower), ("upper", upper)]
    }
    m0_cases = []
    for role in ("A", "B"):
        c = copy.deepcopy(next(c for c in cases if c["role"] == role and c["layout"] == "P77"))
        c.update({
            "id": f"r7_strict_{role}_M0", "layout": "M0", "layout_id": "M0", "kind": "synthetic3d_M0",
            "M0_equivalent_tensor": copy.deepcopy(m0_tensors),
            "homogenize_z_ranges_mm": [[0.0, 30.0], [90.0, 120.0]],
            "phase_volume_interpretation": "Explicit original layer map retained for provenance; homogeneous anisotropic exterior no longer has literal separate low/high phase domains.",
            "M0_reference_only": True,
        })
        m0_cases.append(c)

    def target(c, sigma):
        return {
            "state_id": f"{c['id']}_sf_{sigma_tag(sigma)}", "case_id": c["id"],
            "role": c["role"], "layout": c["layout"], "kind": c["kind"],
            "sigma_f_S_m": sigma,
            "target_relative_folder": f"synthetic3d/{c['id']}/sf_{sigma_tag(sigma)}",
            "background_case_id": c["id"]+"_bg_z",
        }
    targets = [target(c, sigma) for c in cases for sigma in SIGMAS]
    m0_targets = [target(c, sigma) for c in m0_cases for sigma in SIGMAS]
    background_states = [
        {"state_id": c["id"]+"_bg_"+axis, "case_id": c["id"], "role": c["role"], "layout": c["layout"], "load": axis, "passive_ellipsoid": True}
        for c in cases for axis in ("z", "x")
    ]
    background_states.append({"state_id": cases[0]["id"]+"_bg_y", "case_id": cases[0]["id"], "role": cases[0]["role"], "layout": cases[0]["layout"], "load": "y", "passive_ellipsoid": True, "purpose": "Single explicit y-direction check of analytically y-invariant background; not a full tensor claim."})
    m0_background_states = [{"state_id": c["id"]+"_bg_z", "case_id": c["id"], "role": c["role"], "layout": "M0", "load": "z", "passive_ellipsoid": True} for c in m0_cases]
    refinement_targets = []
    for role, sigma in [("A", 0.18), ("B", 10.0)]:
        for layout in ("P11", "P77"):
            c = next(c for c in cases if c["role"] == role and c["layout"] == layout)
            entry = target(c, sigma)
            entry["reason"] = "Preselected before any Round-7 target solve. Paired near/near versus far/far targets test the controlled joint peripheral change; two pairs cover A/B and both conductivities. Corresponding true background is refined too."
            refinement_targets.append(entry)

    design = {
        "schema": "round7_finite3d_geometry_v1", "status": "prescribed_design_and_analytical_geometry_only",
        "source_task": "round7_3d_ct_plan/01_Codex_三维与CT补充验证计划.md",
        "source_actual_layer_design": str(R5.relative_to(ROOT)).replace("\\", "/"),
        "case_count": len(cases), "target_count": len(targets), "M0_case_count": len(m0_cases), "M0_target_count": len(m0_targets),
        "geometry_convention": "Finite ellipsoid (a,b,c)=(30,30,1.5) mm inside a 120-mm cube; no extrusion or periodic y direction.",
        "G_definition": "I/DeltaV, S; sigma_z_app=G*Lz/(Lx*Ly), S/m; g=Gf-G0; eta=g/G0. No out-of-plane thickness parameter.",
        "cases": cases, "targets": targets, "new_states": targets,
        "m0_cases": m0_cases, "m0_targets": m0_targets,
        "background_states": background_states, "m0_background_states": m0_background_states,
        "refinement_targets": refinement_targets, "max_refined_full_targets": 4,
        "M0_equivalent_tensor": m0_tensors,
        "provenance_warning": "Analytical volumes and directions are design checks. Actual COMSOL integrations, mesh, current/energy balance and solved models must be saved independently.",
        "scope_stop": "16 explicit targets, four M0 references and at most four refined explicit target/background pairs; no parameter scan, crossing refinement, new modal solver or added physical fields.",
    }
    assert len(cases) == 8 and len(targets) == 16 and len(m0_targets) == 4 and len(refinement_targets) == 4
    write_json(OUT/"design.json", design)
    write_json(OUT/"refinement_targets.json", {"selection_basis": "Fixed geometry/conductivity coverage before target results, not effect-size selection.", "targets": refinement_targets})
    rows = []
    for c in cases:
        g = c["geometry_checks"]
        rows.append({
            "case_id": c["id"], "role": c["role"], "layout": c["layout"],
            "Lx_mm": c["Lx_mm"], "Ly_mm": c["Ly_mm"], "Lz_mm": c["Lz_mm"],
            "a_mm": c["a_mm"], "b_mm": c["b_mm"], "c_mm": c["c_mm"],
            "ellipsoid_volume_mm3_analytic": g["ellipsoid_volume_mm3"],
            "background_low_volume_mm3_analytic": g["background_low_volume_mm3"],
            "background_high_volume_mm3_analytic": g["background_high_volume_mm3"],
            "replaced_low_volume_mm3_analytic": g["replaced_low_volume_mm3"],
            "replaced_high_volume_mm3_analytic": g["replaced_high_volume_mm3"],
            "cap_minus_high_fraction_analytic": g["endcap_minus"]["high_fraction_analytic"],
            "cap_plus_high_fraction_analytic": g["endcap_plus"]["high_fraction_analytic"],
            "cap_volume_each_mm3_analytic": g["endcap_minus"]["volume_mm3_analytic"],
            "minimum_core_clearance_mm": g["minimum_core_interface_clearance_mm"],
            "lower_band_z0_mm": c["lower_low_band_z_mm"][0], "lower_band_z1_mm": c["lower_low_band_z_mm"][1],
            "upper_band_z0_mm": c["upper_high_band_z_mm"][0], "upper_band_z1_mm": c["upper_high_band_z_mm"][1],
            "background_Gz_S_analytic": g["Gz_background_analytic_S"],
            "background_sigma_z_S_m_analytic": g["sigma_z_background_analytic_S_m"],
            "background_sigma_xy_S_m_analytic": g["sigma_x_background_analytic_S_m"],
            "value_kind": "analytical prescribed geometry, not COMSOL result",
        })
    with (OUT/"geometry_candidates.csv").open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"output": str(OUT/"design.json"), "explicit_geometries": len(cases), "explicit_targets": len(targets), "M0_targets": len(m0_targets), "explicit_background_direction_states": len(background_states), "preselected_refinements": len(refinement_targets), "ellipsoid_volume_mm3": cases[0]["geometry_checks"]["ellipsoid_volume_mm3"], "minimum_core_clearance_mm": cases[0]["geometry_checks"]["minimum_core_interface_clearance_mm"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
