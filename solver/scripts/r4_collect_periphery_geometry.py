"""Join four new peripheral candidate designs with their actual COMSOL backgrounds.

Reads only round-4 files. No solve, old-case audit, or replacement of candidates.
The output explicitly separates analytical design quantities from FE integrals.
"""
import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "runs" / "round4_crossover"


def read_json(p):
    return json.loads(p.read_text(encoding="utf-8-sig"))


def main():
    manifest = read_json(BASE / "periphery_manifest.json")
    candidates = {r["case_id"]: r for r in csv.DictReader(
        (BASE / "periphery_geometry_candidates.csv").open(encoding="utf-8-sig"))}
    measurements = {}
    for c in manifest["cases"]:
        folder = BASE / "periphery" / c["id"]
        m = read_json(folder / "background_matching.json")
        s = read_json(folder / "source_background" / "SUCCESS.json")
        assert s["status"] == "success" and m["case_id"] == c["id"]
        measurements[c["id"]] = (m, s)
    rows = []
    for c in manifest["cases"]:
        cid = c["id"]
        m, s = measurements[cid]
        row = dict(candidates[cid])
        row["quantity_status"] = "analytic_design_and_actual_COMSOL_background_side_by_side"
        row.update({
            "matching_record_pass": m["strict_matching_pass"],
            "background_response_tolerance_relative": m["background_response_tolerance_relative"],
            "phase_area_tolerance_mm2": m["phase_area_tolerance_mm2"],
            "background_high_area_actual_mm2": m["background_high_area_mm2"],
            "background_low_area_actual_mm2": m["background_low_area_mm2"],
            "background_high_area_actual_minus_original_design_mm2": m["background_high_area_mm2"] - c["background_high_area_analytic_mm2"],
            "background_low_area_actual_minus_original_design_mm2": m["background_low_area_mm2"] - c["background_low_area_analytic_mm2"],
            "ellipse_area_analytic_mm2": c["ellipse_area_analytic_mm2"],
            "ellipse_area_actual_mm2": m["ellipse_area_mm2"],
            "ellipse_area_actual_minus_analytic_mm2": m["ellipse_area_mm2"] - c["ellipse_area_analytic_mm2"],
            "ellipse_replaced_high_area_actual_mm2": m["replaced_high_area_mm2"],
            "ellipse_replaced_low_area_actual_mm2": m["replaced_low_area_mm2"],
            "ellipse_replaced_high_actual_minus_original_design_mm2": m["replaced_high_area_mm2"] - c["ellipse_replaced_high_area_analytic_mm2"],
            "ellipse_replaced_low_actual_minus_original_design_mm2": m["replaced_low_area_mm2"] - c["ellipse_replaced_low_area_analytic_mm2"],
            "G0z_analytic_S": c["background_sigma_eff_z_analytic_S_m"],
            "G0z_actual_S": m["G0z_S"],
            "G0z_original_recorded_S": m["G0z_original_S"],
            "G0z_actual_minus_original_S": m["G0z_S"] - m["G0z_original_S"],
            "G0z_relative_difference_from_original": m["G0z_relative_to_original"],
            "G0z_actual_minus_analytic_S": m["G0z_S"] - c["background_sigma_eff_z_analytic_S_m"],
            "G0x_analytic_S": m["G0x_analytic_S"],
            "G0x_actual_S": m["G0x_S"],
            "G0x_actual_minus_analytic_S": m["G0x_S"] - m["G0x_analytic_S"],
            "G0x_relative_difference_from_analytic": m["G0x_relative_to_analytic"],
            "G0x_original_measured_reference_available": False,
            "G0_convention": "t=1 m and Lx=Lz; conductance in S equals effective-conductivity numerical value but not unit",
            "capminus_area_analytic_mm2": c["finite_endcap_area_each_analytic_mm2"],
            "capplus_area_analytic_mm2": c["finite_endcap_area_each_analytic_mm2"],
            "capminus_area_actual_mm2": s["capminus_area_mm2"],
            "capplus_area_actual_mm2": s["capplus_area_mm2"],
            "capminus_high_fraction_actual": m["capminus_high_fraction"],
            "capplus_high_fraction_actual": m["capplus_high_fraction"],
            "cap_dV_actual_V": m["cap_dV_V"],
            "cap_dV_original_recorded_V": m["original_cap_dV_V"],
            "cap_dV_actual_minus_original_V": m["cap_dV_V"] - m["original_cap_dV_V"],
            "triangles": s["triangles"],
            "min_mesh_quality": s["min_mesh_quality"],
            "z_reaction_flux_relative_error": s["reaction_flux_relative_error"],
            "z_reaction_power_relative_error": s["reaction_power_relative_error"],
            "x_flux_relative_error": m["x_flux_relative_error"],
            "x_energy_relative_error": m["x_energy_relative_error"],
            "measurement_source": "Actual COMSOL domain/boundary quadrature and solved background; no candidate area substituted",
            "original_area_comparison_basis": "Original analytical design from inherited candidate record; no new read or audit of old measured areas",
            "matching_source": f"periphery/{cid}/background_matching.json",
            "actual_domain_source": f"periphery/{cid}/source_background/SUCCESS.json",
        })
        other_layout = "swap_far" if c["layout"] == "swap_near" else "swap_near"
        other_id = f"r4_strict_{c['configuration_role']}_{other_layout}"
        other, _ = measurements[other_id]
        row["other_layout_case_id"] = other_id
        for key in ["G0z_S", "G0x_S", "cap_dV_V", "background_high_area_mm2", "background_low_area_mm2",
                    "replaced_high_area_mm2", "replaced_low_area_mm2"]:
            row["current_minus_other_layout_" + key] = m[key] - other[key]
        rows.append(row)
    output = BASE / "periphery_geometry.csv"
    with output.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"rows": len(rows), "matching_pass_count": sum(r["matching_record_pass"] for r in rows),
                      "output": str(output),
                      "max_G0z_relative_to_original": max(r["G0z_relative_difference_from_original"] for r in rows),
                      "max_G0x_relative_to_analytic": max(r["G0x_relative_difference_from_analytic"] for r in rows),
                      "max_cap_dV_change_V": max(abs(r["cap_dV_actual_minus_original_V"]) for r in rows)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
