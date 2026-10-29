"""Analytic, background-only DtN admittances of the flat peripheral rectangles.

Does not run COMSOL, read target fields, or change the frozen predictor.
All calculations use SI units; manifest input geometry is in mm.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "runs" / "round4_extension_20261024" / "peripheral_modes"


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def layers_from_electrode(case, side):
    """Return (thickness_m, sigma_S_m), starting at grounded-delta-V electrode."""
    if any(abs(float(a)) > 1e-14 for a in case["interface_amplitudes_mm"]):
        raise ValueError("Separated cosine modes require flat peripheral layers")
    edges = [0.0, *case["z_interfaces_mm"], case["L_mm"]]
    lo, hi = (0.0, 30.0) if side == "lower" else (90.0, 120.0)
    layers = []
    for z0, z1, phase in zip(edges[:-1], edges[1:], case["layer_materials"]):
        h = min(z1, hi) - max(z0, lo)
        if h > 0:
            sigma = case["sigma_high_S_m"] if phase else case["sigma_low_S_m"]
            assert sigma > 0
            layers.append((h * 1e-3, float(sigma)))
    if side == "upper":
        layers.reverse()
    assert math.isclose(sum(h for h, _ in layers), 0.03, abs_tol=1e-14)
    return layers


def admittance(layers, n, width_m=0.12):
    """Input DtN coefficient Y_n [S/m^2] at the core boundary.

    Propagate impedance from zero impedance (delta V = 0) at the electrode.
    tanh recursion avoids overflow of direct cosh/sinh transfer matrices.
    """
    assert int(n) == n and n >= 0 and width_m > 0
    if n == 0:
        return 1.0 / math.fsum(h / sigma for h, sigma in layers)
    k = n * math.pi / width_m
    impedance = 0.0
    for h, sigma in layers:
        beta = sigma * k
        q = math.tanh(k * h)
        impedance = (impedance + q / beta) / (1.0 + beta * impedance * q)
    return 1.0 / impedance


def validations():
    errors = []
    for sigma in [0.01, 0.1, 3.0]:
        for n in [0, 1, 4, 16, 10000]:
            layers = [(0.007, sigma), (0.008, sigma), (0.015, sigma)]
            target = sigma / 0.03 if n == 0 else sigma * (n * math.pi / 0.12) / math.tanh(n * math.pi / 0.12 * 0.03)
            got = admittance(layers, n)
            errors.append(abs(got / target - 1.0))
    assert max(errors) < 2e-14
    a = [(0.02, 0.01), (0.01, 0.1)]
    b = list(reversed(a))
    assert admittance(a, 0) == admittance(b, 0)
    assert not math.isclose(admittance(a, 1), admittance(b, 1), rel_tol=1e-4)
    return {"uniform_layer_limit_max_relative_error": max(errors),
            "uniform_partition_invariance": "passed",
            "zero_mode_series_resistance_and_order_invariance": "passed",
            "nonzero_mode_can_depend_on_order": "passed",
            "high_mode_no_overflow_n": 10000}


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    cases = []
    for c in read_json(ROOT / "runs/round4_crossover/periphery_manifest.json")["cases"]:
        cases.append(("NN" if c["layout"] == "swap_near" else "FF", c))
    for c in read_json(ROOT / "runs/round4_extension_20261024/manifest.json")["cases"]:
        cases.append((c["layout"], c))
    rows, layer_rows = [], []
    for layout, c in cases:
        for side in ["lower", "upper"]:
            layers = layers_from_electrode(c, side)
            for i, (h, sigma) in enumerate(layers, 1):
                layer_rows.append({"case_id": c["id"], "role": c["configuration_role"], "layout": layout,
                                   "side": side, "order_from_electrode": i, "thickness_m": h, "sigma_S_m": sigma})
            for n in range(17):
                rows.append({"case_id": c["id"], "role": c["configuration_role"], "layout": layout,
                             "side": side, "n": n, "k_per_m": n * math.pi / .12,
                             "Y_S_per_m2": admittance(layers, n), "width_m": .12,
                             "out_of_plane_thickness_m": c["out_of_plane_thickness_m"]})
    for name, data in [("admittances.csv", rows), ("layers_from_electrode.csv", layer_rows)]:
        with (DEST / name).open("w", newline="", encoding="utf-8-sig") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    check = validations()
    check.update({"status": "analytic_background_coefficients_only_not_target_solution", "n_max": 16,
                  "DtN_units": "S/m^2", "fourier_coefficients_units": "V",
                  "perturbation_energy_formula": "t*L*(Y0*a0^2+0.5*sum(n>=1,Yn*an^2)) [W]",
                  "restriction": "Flat x-uniform peripheral layers, insulated x sides and fixed electrode delta V=0; does not apply to wavy peripherals.",
                  "interpretation": "Yn depends only on background; mode amplitudes from solved target are posthoc. This is not a replacement predictor or a proof of the sign of gain differences."})
    (DEST / "analytic_validation.json").write_text(json.dumps(check, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"rows": len(rows), "output": str(DEST), "validation": check}, ensure_ascii=False))


if __name__ == "__main__":
    main()
