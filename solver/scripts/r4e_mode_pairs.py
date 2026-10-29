"""Posthoc two-field DtN identity from saved COMSOL modal traces.

No solver call, fit, completion of missing modes, or predictor modification.
Requires r4e_modes_run to finish first. All old state files are read only.
"""
from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/round4_extension_20261024"
PAIRS = [("FF", "NN", "both_move_near"), ("FF", "NF", "lower_only_near"),
         ("FN", "NN", "lower_only_near"), ("FF", "FN", "upper_only_near"),
         ("NF", "NN", "upper_only_near")]


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def target_path(role, layout):
    if layout in ("NN", "FF"):
        old_layout = "swap_near" if layout == "NN" else "swap_far"
        folder = ROOT / "runs/round4_crossover/periphery" / f"r4_strict_{role}_{old_layout}"
    else:
        folder = OUT / "cases" / f"r4e_strict_{role}_{layout}"
    return folder / "sf_0p159375/SUCCESS.json"


def load_state(role, layout):
    folder = OUT / "modes" / f"{role}_{layout}"
    summary = read_json(folder / "mode_summary.json")
    with (folder / "mode_coefficients.csv").open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    data = {(r["side"], int(r["n"])): {"a": float(r["delta_V_cosine_V"]), "Y": float(r["Y_S_per_m2"])} for r in rows}
    assert len(data) >= 66
    truth = read_json(target_path(role, layout))
    assert math.isclose(truth["sigma_f_S_m"], .159375, rel_tol=1e-13)
    return summary, data, truth


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    summaries, terms, closures = [], [], []
    for role in ["A", "B"]:
        states = {layout: load_state(role, layout) for layout in ["NN", "NF", "FN", "FF"]}
        for layout, (state_summary, _, _) in states.items():
            for side in state_summary["sides"]:
                closures.append({"role": role, "layout": layout, "side": side["side"],
                                 "mode0_power_W": side["mode0_power_W"],
                                 "n1to32_power_W": side["n1to32_power_W"],
                                 "direct_outer_perturbation_power_W": side["direct_outer_perturbation_power_W"],
                                 "mode32_minus_direct_W": side["mode32_minus_direct_W"],
                                 "mode32_direct_relative_difference": side["mode32_direct_relative_difference"],
                                 "mode16_to32_increment_relative_to_direct": side["mode16_to32_increment_relative_to_direct"],
                                 "mode32_fraction_from_nge1": side["mode32_fraction_from_nge1"]})
        for p, q, comparison in PAIRS:
            sp, dp, tp = states[p]
            sq, dq, tq = states[q]
            for field in ["width_m", "thickness_m", "DV_V"]:
                assert math.isclose(sp[field], sq[field], rel_tol=1e-13)
            factor = sp["thickness_m"] * sp["width_m"] / sp["DV_V"]**2
            actual = tq["deltaG_S"] - tp["deltaG_S"]
            sums = {side: [] for side in ["lower", "upper"]}
            for side in sums:
                for n in range(33):
                    ap, aq = dp[side, n]["a"], dq[side, n]["a"]
                    change_y = dq[side, n]["Y"] - dp[side, n]["Y"]
                    term = factor * (1 if n == 0 else .5) * change_y * ap * aq
                    sums[side].append(term)
                    terms.append({"role": role, "layout_p": p, "layout_q": q, "comparison": comparison,
                                  "side": side, "n": n, "delta_Y_S_per_m2": change_y,
                                  "deltaV_p_V": ap, "deltaV_q_V": aq, "bilinear_increment_S": term})
            row = {"role": role, "layout_p": p, "layout_q": q, "comparison": comparison,
                   "sigma_f_S_m": .159375, "COMSOL_deltaG_q_minus_p_S": actual,
                   "G0_q_minus_p_S": tq["G0_S"] - tp["G0_S"], "U_q_minus_p_S": tq["U_S"] - tp["U_S"],
                   "mode0_term_S": sums["lower"][0] + sums["upper"][0],
                   "lower_modes0to32_S": math.fsum(sums["lower"]),
                   "upper_modes0to32_S": math.fsum(sums["upper"])}
            for nmax in [4, 8, 16, 32]:
                total = math.fsum(math.fsum(vals[:nmax + 1]) for vals in sums.values())
                row[f"bilinear_modes0to{nmax}_S"] = total
                row[f"modes{nmax}_minus_COMSOL_S"] = total - actual
                row[f"modes{nmax}_absolute_relative_error"] = abs(total - actual) / abs(actual) if actual else None
            row["mode16to32_change_S"] = row["bilinear_modes0to32_S"] - row["bilinear_modes0to16_S"]
            summaries.append(row)
    dest = OUT / "modes"
    write_csv(dest / "two_field_identity.csv", summaries)
    write_csv(dest / "two_field_mode_terms.csv", terms)
    write_csv(dest / "modal_energy_closure.csv", closures)
    one_side = [r for r in summaries if r["comparison"] != "both_move_near"]
    result = {
        "status": "posthoc_two_target_field_identity_evaluated", "pair_count": len(summaries),
        "formula": "deltaG_q-deltaG_p = t*L/DV^2*[deltaY0*a0p*a0q + .5*sum(n>=1,deltaYn*anp*anq)] summed over two peripheral interfaces",
        "requirements": "Same core/background gradient, ellipse contrast, and electrode loads; n0 peripheral admittance identical. Both actual target fields are used.",
        "max_abs_relative_residual_n16": max(r["modes16_absolute_relative_error"] or 0 for r in summaries),
        "max_abs_relative_residual_n32": max(r["modes32_absolute_relative_error"] or 0 for r in summaries),
        "max_mode16to32_abs_change_S": max(abs(r["mode16to32_change_S"]) for r in summaries),
        "max_abs_residual_n32_S": max(abs(r["modes32_minus_COMSOL_S"]) for r in summaries),
        "single_side_pair_count": len(one_side),
        "single_side_n32_relative_residual_min": min(r["modes32_absolute_relative_error"] for r in one_side),
        "single_side_n32_relative_residual_median": statistics.median(r["modes32_absolute_relative_error"] for r in one_side),
        "single_side_n32_relative_residual_max": max(r["modes32_absolute_relative_error"] for r in one_side),
        "max_abs_mode0_difference_term_S": max(abs(r["mode0_term_S"]) for r in summaries),
        "max_single_periphery_n32_volume_closure_relative_error": max(r["mode32_direct_relative_difference"] for r in closures),
        "interpretation": "A continuum two-field identity, truncated to saved modes and evaluated on finite-element potentials. Residuals are retained, not corrected. Not an advance prediction or unique causal split when both sides move."
    }
    (dest / "two_field_identity_summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()
