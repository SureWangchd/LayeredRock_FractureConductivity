"""Render the matched-configuration phenomenon from archived COMSOL results.

Geometry panels are exact parameter-based material diagrams, not solved fields.
All numerical marks are read from archived CSV values. No solver is invoked.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Ellipse, Patch, Rectangle

from style import COLORS, apply_style, export_figure, panel_label


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "runs" / "2026-10-21_round2"
DEFAULT_OUT = ROOT / "outputs" / "JGRSE_20261028_LabelRevision" / "pilot02_phenomenon"
FAMILIES = ("core_h20_L60", "strict_h10_L60")
ROLES = ("A", "B")
WIDTH_MM = 183
HEIGHT_MM = 135


def load_data():
    summary = pd.read_csv(SOURCE / "summary.csv").set_index("case_id")
    configs, controls, gains, layers, outlines = {}, [], [], [], []
    for family in FAMILIES:
        for role in ROLES:
            case = f"{family}_{role}_m3_full"
            params = json.loads((SOURCE / case / "parameters.json").read_text(encoding="utf-8-sig"))
            row = summary.loc[case]
            bgx, bgz = summary.loc[params["background_case_x"]], summary.loc[params["background_case_z"]]
            sigma = float(row["sigma_reaction_audited_S_m"])
            sigma0_recorded = float(bgz["sigma_reaction_audited_S_m"])
            sigma0 = 1.0 / 55.0  # exact flat-layer harmonic mean, S/m
            configs[(family, role)] = params
            controls.append({
                "family": family, "role": role, "case_id": case,
                "sigma0_x_S_m": float(bgx["sigma_reaction_audited_S_m"]),
                "sigma0_z_S_m": sigma0_recorded,
                "background_high_area_mm2": float(row["background_high_area_mm2"]),
                "background_low_area_mm2": float(row["background_low_area_mm2"]),
                "ellipse_area_mm2": float(row["original_ellipse_area_mm2"]),
                "replaced_high_area_mm2": float(row["original_ellipse_replaced_high_mm2"]),
                "replaced_low_area_mm2": float(row["original_ellipse_replaced_low_mm2"]),
                "final_high_area_mm2": float(row["area_high_mm2"]),
                "final_low_area_mm2": float(row["area_low_mm2"]),
                "final_fracture_area_mm2": float(row["area_fracture_mm2"]),
                "lower_cap_high_fraction": float(row["capminus_high_fraction"]),
                "upper_cap_high_fraction": float(row["capplus_high_fraction"]),
                "background_case_x": params["background_case_x"],
                "background_case_z": params["background_case_z"],
                "background_mesh_x": int(bgx["mesh_level"]),
                "background_mesh_z": int(bgz["mesh_level"]),
                "source_csv": str(SOURCE / "summary.csv"),
            })
            gains.append({
                "family": family, "role": role, "case_id": case,
                "sigma_with_S_m": sigma, "sigma0_z_S_m": sigma0,
                "sigma0_z_recorded_M1_S_m": sigma0_recorded,
                "relative_gain_recorded_M1_percent": 100 * (sigma - sigma0_recorded) / sigma0_recorded,
                "baseline_definition": "Exact flat-layer harmonic mean, 1/55 S/m",
                "delta_sigma_S_m": sigma - sigma0,
                "relative_gain_percent": 100 * (sigma - sigma0) / sigma0,
                "mesh_level": int(row["mesh_level"]), "triangles": int(row["triangles"]),
                "out_of_plane_thickness_m": params["out_of_plane_thickness_m"],
                "source_csv": str(SOURCE / "summary.csv"),
            })
            bounds = [0.0] + params["z_interfaces_mm"] + [params["L_mm"]]
            materials = params["layer_materials"]
            for index, material in enumerate(materials):
                layers.append({
                    "family": family, "role": role, "case_id": case,
                    "partition_index_bottom_to_top": index + 1,
                    "x_min_mm": 0, "x_max_mm": params["L_mm"],
                    "z_min_mm": bounds[index], "z_max_mm": bounds[index + 1],
                    "material": "high" if material else "low",
                    "sigma_S_m": params["sigma_high_S_m"] if material else params["sigma_low_S_m"],
                    "top_is_physical_internal_interface": bool(index < len(materials) - 1 and material != materials[index + 1]),
                    "source_parameters": str(SOURCE / case / "parameters.json"),
                })
            theta = np.deg2rad(params["theta_deg"])
            for t in np.linspace(0, 2 * np.pi, 721):
                u, v = params["a_mm"] * np.cos(t), params["b_mm"] * np.sin(t)
                outlines.append({"family": family, "role": role, "angle_parameter_rad": t,
                                 "x_mm": params["cx_mm"] + u * np.cos(theta) - v * np.sin(theta),
                                 "z_mm": params["cz_mm"] + u * np.sin(theta) + v * np.cos(theta)})
    return configs, pd.DataFrame(controls), pd.DataFrame(gains), pd.DataFrame(layers), pd.DataFrame(outlines)


def physical_runs(params):
    """Merge artificial same-material partitions for a physical material map."""
    bounds = [0.0] + params["z_interfaces_mm"] + [params["L_mm"]]
    runs = []
    for i, material in enumerate(params["layer_materials"]):
        if runs and runs[-1][2] == material:
            runs[-1] = (runs[-1][0], bounds[i + 1], material)
        else:
            runs.append((bounds[i], bounds[i + 1], material))
    return runs


def geometry_axis(ax, params, role, show_y, angle=False):
    length = params["L_mm"]
    for lo, hi, material in physical_runs(params):
        ax.add_patch(Rectangle((0, lo), length, hi - lo,
                               facecolor=COLORS["high" if material else "low"], edgecolor="none"))
    for lo, hi, _ in physical_runs(params)[:-1]:
        ax.plot([0, length], [hi, hi], color="black", lw=0.4)
    for lo, hi in ((0, params["buffer_bottom_mm"]),
                   (length - params["buffer_top_mm"], length)):
        ax.add_patch(Rectangle((0, lo), length, hi - lo, facecolor="none", edgecolor="black",
                               hatch="///", linewidth=0.0, alpha=0.24, zorder=3))
    ax.add_patch(Ellipse((params["cx_mm"], params["cz_mm"]), 2 * params["a_mm"],
                         2 * params["b_mm"], angle=params["theta_deg"],
                         facecolor=COLORS["fracture"], edgecolor="black", lw=0.55, zorder=5))
    # The endcaps have a finite axial length of 3 mm, rather than point tips.
    theta = np.deg2rad(params["theta_deg"])
    for sign in (-1, 1):
        u = sign * (params["a_mm"] - params["endcap_length_mm"])
        vmax = params["b_mm"] * np.sqrt(1 - (u / params["a_mm"]) ** 2)
        points = np.array([[params["cx_mm"] + u * np.cos(theta) - v * np.sin(theta),
                            params["cz_mm"] + u * np.sin(theta) + v * np.cos(theta)]
                           for v in (-vmax, vmax)])
        ax.plot(points[:, 0], points[:, 1], color="black", lw=0.5, zorder=6)
    if angle:
        cx, cz = params["cx_mm"], params["cz_mm"]
        ax.plot([cx, cx + 24], [cz, cz], color="black", lw=0.5, zorder=7)
        ax.add_patch(Arc((cx, cz), 34, 34, theta1=0, theta2=params["theta_deg"],
                         color="black", lw=0.5, zorder=7))
        ax.text(cx + 26, cz + 11, "60°", fontsize=8, ha="center", va="center", zorder=8, color="black")
    ax.plot([0, length], [0, 0], color="black", lw=0.8, clip_on=False, zorder=9)
    ax.plot([0, length], [length, length], color="black", lw=0.8, clip_on=False, zorder=9)
    ax.text(0.5, 1.035, "1 V", transform=ax.transAxes, ha="center", va="bottom", fontsize=8, color="black")
    ax.text(0.98, 0.025, "0 V", transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="black",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.5, "alpha": 0.8}, zorder=10)
    ax.text(0.035, 0.94, role, transform=ax.transAxes, fontsize=9, color="black",
            va="top", fontstyle="italic", zorder=10)
    ax.set(xlim=(0, length), ylim=(0, length), aspect="equal", xticks=[0, 60, 120], yticks=[0, 60, 120])
    ax.set_xlabel(r"$x$ (mm)", labelpad=1.2)
    ax.tick_params(axis="both", labelsize=8, length=2.2, pad=2, direction="in", top=True, right=True)
    if show_y:
        ax.set_ylabel(r"$z$ (mm)", labelpad=1.2)
    else:
        ax.tick_params(axis="y", labelleft=False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.6)


CONTROL_ROWS = [
    ("sigma0_x_S_m", r"$\sigma_{0,x}$ (S m$^{-1}$)", ".5f"),
    ("sigma0_z_S_m", r"$\sigma_{0,z}$ (S m$^{-1}$)", ".7f"),
    ("background_high_area_mm2", "High area (mm²)", ".3f"),
    ("background_low_area_mm2", "Low area (mm²)", ".3f"),
    ("ellipse_area_mm2", "Ellipse area (mm²)", ".3f"),
    ("replaced_high_area_mm2", "Replaced high (mm²)", ".3f"),
    ("replaced_low_area_mm2", "Replaced low (mm²)", ".3f"),
]


def control_table(ax, controls):
    ax.set_axis_off()
    ax.set(xlim=(0, 1), ylim=(0, 1))
    panel_label(ax, "c", x=-0.095, y=1.08)
    ax.text(0, 1.08, "Matched quantities", va="top", fontsize=9, color="black")
    ax.text(0.66, 0.95, "Core A/B", ha="center", va="top", fontsize=8.5, color="black")
    ax.text(0.90, 0.95, "Strict A/B", ha="center", va="top", fontsize=8.5, color="black")
    ax.plot([0, 1], [0.865, 0.865], color="black", lw=0.6)
    display = []
    for i, (key, label, fmt) in enumerate(CONTROL_ROWS):
        y = 0.80 - i * 0.107
        ax.text(0, y, label, ha="left", va="center", fontsize=8.2, color="black")
        row = {"quantity": key, "label": label}
        for family, x in zip(FAMILIES, [0.66, 0.90]):
            vals = controls.loc[controls.family.eq(family), key].to_numpy(float)
            rounded = [format(v, fmt) for v in vals]
            if len(set(rounded)) != 1:
                raise ValueError(f"Control table would hide visible A/B difference: {family}, {key}, {rounded}")
            ax.text(x, y, rounded[0], ha="center", va="center", fontsize=8.2, color="black")
            row[family] = rounded[0]
        display.append(row)
    ax.plot([0, 1], [0.085, 0.085], color="black", lw=0.6)
    return pd.DataFrame(display)


def gain_panel(ax, gains):
    panel_label(ax, "d", x=-0.257, y=1.08)
    centers = np.array([0, 1.2])
    for center, family in zip(centers, FAMILIES):
        pair = gains.loc[gains.family.eq(family), "relative_gain_percent"].to_numpy(float)
        ax.plot([center, center], [pair.min(), pair.max()], color="black", lw=0.55, zorder=1)
    for role, marker in zip(ROLES, ["o", "^"]):
        values = [float(gains.loc[gains.family.eq(f) & gains.role.eq(role), "relative_gain_percent"].iloc[0]) for f in FAMILIES]
        ax.plot(centers, values, linestyle="none", marker=marker, ms=5.5, markerfacecolor="white",
                markeredgecolor=COLORS[role], markeredgewidth=1.1, label=role, zorder=3)
        for center, value in zip(centers, values):
            ax.text(center + 0.075, value, f"{value:.2f}", ha="left", va="center", fontsize=8.5, color="black")
    ax.set_xticks(centers, ["Core", "Strict"])
    ax.set_xlim(-0.33, 1.68)
    ax.set_ylim(0, 34)
    ax.set_yticks([0, 10, 20, 30])
    ax.set_ylabel(r"$\Delta\sigma_{\mathrm{eff}}/\sigma_{0,z}$ (%)", labelpad=5)
    ax.grid(False)
    ax.legend(ncol=2, loc="upper left", borderaxespad=0.3,
              columnspacing=1.0, handlelength=1.0, handletextpad=0.3, fontsize=8.5)
    ax.tick_params(axis="both", direction="in", top=True, right=True, pad=3)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.6)


def background_analytic_check(controls):
    """Describe only the four plotted models' archived flat-layer references."""
    sigma_x = 0.5 * 0.01 + 0.5 * 0.1
    sigma_z = 1.0 / (0.5 / 0.01 + 0.5 / 0.1)
    rows = []
    for row in controls.to_dict(orient="records"):
        for direction, expected in (("x", sigma_x), ("z", sigma_z)):
            actual = float(row[f"sigma0_{direction}_S_m"])
            rows.append({
                "case_id": row["case_id"], "family": row["family"], "role": row["role"],
                "background_case_id": row[f"background_case_{direction}"],
                "background_mesh_level": int(row[f"background_mesh_{direction}"]),
                "direction": direction, "sigma0_recorded_S_m": actual,
                "sigma0_analytic_S_m": expected,
                "absolute_relative_deviation": abs(actual - expected) / abs(expected),
            })
    return pd.DataFrame(rows)


def write_caption(outdir):
    caption = (ROOT / "source_data/Fig01_configuration_effect/caption.txt").read_text(encoding="utf-8").strip()
    (outdir / "caption.txt").write_text(caption + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    apply_style()
    configs, controls, gains, layers, outlines = load_data()
    controls.to_csv(outdir / "panel_c_controls_all_states.csv", index=False)
    analytic_check = background_analytic_check(controls)
    analytic_check.to_csv(outdir / "panel_c_background_analytic_comparison.csv", index=False)
    gains.to_csv(outdir / "panel_d_conductivity_gain.csv", index=False)
    layers.to_csv(outdir / "panel_ab_layer_geometry.csv", index=False)
    outlines.to_csv(outdir / "panel_ab_ellipse_outline.csv", index=False)
    pd.DataFrame([{k: p[k] for k in ["id", "L_mm", "cx_mm", "cz_mm", "a_mm", "b_mm", "theta_deg",
                                                   "endcap_length_mm", "buffer_bottom_mm", "buffer_top_mm",
                                                   "sigma_low_S_m", "sigma_high_S_m", "sigma_frac_S_m"]}
                  for p in configs.values()]).to_csv(outdir / "panel_ab_ellipse_parameters.csv", index=False)

    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    geom_width = 0.182
    geom_height = geom_width * WIDTH_MM / HEIGHT_MM
    xs = [0.075, 0.303, 0.565, 0.793]
    axes = []
    for i, ((family, role), params) in enumerate(configs.items()):
        ax = fig.add_axes([xs[i], 0.64, geom_width, geom_height])
        geometry_axis(ax, params, role, show_y=(role == "A"), angle=(i == 0))
        axes.append(ax)
    panel_label(axes[0], "a", x=-0.22, y=1.23)
    panel_label(axes[2], "b", x=-0.22, y=1.23)
    fig.text(0.278, 0.974, "Core", ha="center", va="top", fontsize=9, color="black")
    fig.text(0.278, 0.945, "Different end contacts", ha="center", va="top", fontsize=7.3, color="black")
    fig.text(0.768, 0.974, "Strict", ha="center", va="top", fontsize=9, color="black")
    fig.text(0.768, 0.945, "Matched end contacts", ha="center", va="top", fontsize=7.3, color="black")
    handles = [Patch(facecolor=COLORS["low"], edgecolor="black", linewidth=0.4, label="Low"),
               Patch(facecolor=COLORS["high"], edgecolor="black", linewidth=0.4, label="High"),
               Patch(facecolor=COLORS["fracture"], edgecolor="black", linewidth=0.4, label="Fracture"),
               Patch(facecolor="white", edgecolor="black", hatch="///", linewidth=0.4, label="Buffer")]
    fig.legend(handles=handles, loc="center", bbox_to_anchor=(0.52, 0.541), ncol=4, columnspacing=1.8,
               handlelength=1.6, handleheight=0.8, handletextpad=0.5, fontsize=8.5)

    table_ax = fig.add_axes([0.075, 0.075, 0.415, 0.36])
    display = control_table(table_ax, controls)
    display.to_csv(outdir / "panel_c_display_values.csv", index=False)
    gain_ax = fig.add_axes([0.615, 0.105, 0.35, 0.33])
    gain_panel(gain_ax, gains)
    export_figure(fig, outdir, "figure02_phenomenon")
    plt.close(fig)
    write_caption(outdir)

    note_path = outdir / "figure_notes.json"
    notes = json.loads(note_path.read_text(encoding="utf-8")) if note_path.exists() else {}
    mismatch_rows = []
    for family in FAMILIES:
        pair = controls[controls.family.eq(family)].set_index("role")
        for key, _, _ in CONTROL_ROWS:
            va, vb = float(pair.loc["A", key]), float(pair.loc["B", key])
            mismatch_rows.append({"family": family, "quantity": key, "value_A": va, "value_B": vb,
                                  "B_minus_A": vb - va, "symmetric_relative_abs_difference": abs(vb - va) / ((abs(va) + abs(vb)) / 2)})
    pd.DataFrame(mismatch_rows).to_csv(outdir / "panel_c_pairwise_mismatches.csv", index=False)
    notes.update({
        "script": str(Path(__file__).resolve()), "source_case_count": 4,
        "source_case_ids": gains.case_id.tolist(), "dimensions_mm": [WIDTH_MM, HEIGHT_MM],
        "gain_definition": "100 * (sigma_with - sigma0_z) / sigma0_z; sigma0_z = 1/55 S/m, the exact flat-layer harmonic mean",
        "statistics": "four deterministic selected configurations; no replicates, p values or statistical error bars",
        "background_mesh_note": "Gains use the exact analytical background 1/55 S/m. Archived M1 background values remain in the control and original-source columns; the full-ellipse target solutions are M3.",
        "background_reference_analytic_check": {
            "formula_x": "0.5 * sigma_low + 0.5 * sigma_high = 0.055 S/m",
            "formula_z": "1 / (0.5 / sigma_low + 0.5 / sigma_high) = 0.01818181818181818 S/m",
            "max_relative_deviation_x": float(analytic_check.loc[analytic_check.direction.eq("x"), "absolute_relative_deviation"].max()),
            "max_relative_deviation_z": float(analytic_check.loc[analytic_check.direction.eq("z"), "absolute_relative_deviation"].max()),
            "scope": "Comparison of the four plotted states' archived M1 background values to flat-layer analytical conductivities. No background or fractured target was re-solved. This does not bound M3 target-mesh error; normalization is not a same-mesh subtraction.",
        },
        "actual_gain_percent": gains[["family", "role", "relative_gain_percent"]].to_dict(orient="records"),
        "control_display_rule": "A/B displayed values must agree after the declared rounding; the script raises if they do not. Full-precision values and all pairwise mismatches are separately exported.",
        "qa_status": "rendered; visual inspection pending",
        "label_revision": "Short end-contact explanations under Core/Strict; table columns explicitly identify Core A/B and Strict A/B. Geometry, values, colors and frame positions unchanged.",
        "end_contact_label_scope": "Labels refer to finite endcap material contacts within each pair; no claim of matching the entire internal layer sequence or background field.",
    })
    note_path.write_text(json.dumps(notes, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(outdir / "figure02_phenomenon.png")


if __name__ == "__main__":
    main()
