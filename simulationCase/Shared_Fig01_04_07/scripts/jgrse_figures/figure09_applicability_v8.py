"""Plot prospective adequacy and output-specific errors from saved results.

No COMSOL solve, fitted point, data sampling or hash is performed. The modal
curves are algebraic sums of existing prospective M0 interface potentials and
existing analytical response coefficients at every integer mode 1 through 64.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
import numpy as np
import pandas as pd

from style import COLORS, apply_style, export_figure, panel_label

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "runs" / "round6_applicability"
OUT = ROOT / "outputs" / "JGRSE_20261028_PointMarkers" / "pilot09_applicability"


def csv(name):
    return pd.read_csv(SOURCE / name)


def save_csv(df, out, name):
    df.to_csv(out / name, index=False, encoding="utf-8-sig", float_format="%.16g")


def prepare(out):
    validation = csv("applicability_validation.csv")
    pred = csv("applicability_predictions.csv")
    gain = csv("closure_new_source_comparison.csv")
    side = csv("higher_order_side_comparison.csv")
    valid = validation.loc[(validation.reference_K == 0) & (validation.target_N == 64)
        & np.isclose(validation.tolerance, .02)].copy()
    before = pred.loc[(pred.reference_K == 0) & (pred.target_N == 64)].copy()
    valid = valid.merge(before[["state_id", "relative_reference_error_bound", "target_field_used",
        "target_gain_used_in_criterion", "Q_S", "B_S", "reference_gain_S", "background_difference_S"]],
        on="state_id", validate="one_to_one")
    assert len(valid) == 32 and valid.state_id.nunique() == 32
    assert not valid.target_field_used.any() and not valid.target_gain_used_in_criterion.any()
    valid["bound_percent"] = valid.relative_reference_error_bound * 100
    valid["actual_gain_error_percent"] = valid.COMSOL_relative_gain_error * 100
    assert (valid.bound_percent > 0).all() and (valid.actual_gain_error_percent > 0).all()
    assert int(valid.prospective_accepted.sum()) == 19
    assert int((~valid.prospective_accepted & valid.COMSOL_pass).sum()) == 5
    assert not (valid.prospective_accepted & ~valid.COMSOL_pass).any()
    valid = valid.sort_values(["full_length_mm", "role", "sigma_f_S_m", "layout"]).reset_index(drop=True)
    save_csv(valid, out, "all_32_states_source_data.csv")
    for letter, length in [("a", 40), ("b", 64)]:
        selected = valid.loc[valid.full_length_mm == length]
        assert len(selected) == 16
        save_csv(selected, out, f"panel_{letter}_source_data.csv")

    # Read pretarget interface traces. M0 includes the fracture-bearing core.
    coefficients = csv("port_coefficients.csv")
    coefficients = coefficients.loc[(coefficients.layout == "P11") & coefficients.n.between(1, 64)].copy()
    assert len(coefficients) == 128
    save_csv(coefficients, out, "panel_c_peripheral_coefficients.csv")
    modal_rows, cumulative, files = [], [], []
    old = csv("source_scale_response_pretruth.csv")
    checks = []
    for length in [40, 64]:
        for role in ["A", "B"]:
            path = SOURCE / "predictions" / f"r6_L{length}_{role}_P11" / "sf_0p18" / "prediction.json"
            record = json.loads(path.read_text(encoding="utf-8"))
            assert not record["target_truth_used"]
            method = next(m for m in record["methods"] if m["K"] == 0)
            assert not method["target_truth_used"]
            trace = np.asarray(method["port_v"], dtype=float)
            assert trace.shape == (2, 65)
            assert np.array_equal(np.asarray(method["mode_numbers"]), np.arange(65))
            absolute_sum = np.zeros(64)
            for index, name in enumerate(["lower", "upper"]):
                coeff = coefficients.loc[coefficients.side == name].sort_values("n")
                for _, row in coeff.iterrows():
                    n = int(row.n)
                    q = float(method["thickness_m"] * row.norm_N_m * row.delta_Y_S_m2
                        * trace[index, n]**2 / method["DV_V"]**2)
                    absolute_sum[n-1] += abs(q)
                    modal_rows.append({"full_length_mm": length, "role": role, "layout": "P11",
                        "sigma_f_S_m": .18, "reference_K": 0, "side": name, "n": n,
                        "v_M0_V": trace[index, n], "delta_Y_S_m2": row.delta_Y_S_m2,
                        "norm_N_m": row.norm_N_m, "thickness_m": method["thickness_m"],
                        "DV_V": method["DV_V"], "Q_S": q, "absolute_Q_S": abs(q),
                        "target_field_used": False, "prediction_source": str(path)})
            assert absolute_sum.sum() > 0
            cumul = np.cumsum(absolute_sum) / absolute_sum.sum()
            for n, value in enumerate(cumul, start=1):
                cumulative.append({"full_length_mm": length, "role": role, "sigma_f_S_m": .18,
                    "layout": "P11", "reference_K": 0, "n": n,
                    "sum_both_interfaces_absolute_Q_n_S": absolute_sum[n-1],
                    "sum_modes_1_64_absolute_Q_S": absolute_sum.sum(),
                    "cumulative_weight_fraction": value, "target_field_used": False})
            prior = old.loc[(old.full_length_mm == length) & (old.role == role)
                & np.isclose(old.sigma_f_S_m, .18) & (old.layout == "P11") & (old.reference_K == 0)]
            assert len(prior) == 2
            tail_columns = ["absolute_Q_n9_16_S", "absolute_Q_n17_32_S", "absolute_Q_n33_64_S"]
            prior_tail = prior[tail_columns].to_numpy().sum() / prior.absolute_Q_side_sum_S.sum()
            np.testing.assert_allclose(1-cumul[7], prior_tail, atol=2e-13, rtol=2e-10)
            checks.append({"length_mm": length, "role": role, "cumulative_at_n8": cumul[7],
                "tail_from_saved_aggregate": prior_tail, "absolute_reconstruction_difference": abs(1-cumul[7]-prior_tail)})
            files.append(str(path))
    modes = pd.DataFrame(cumulative)
    assert len(modes) == 256 and len(modal_rows) == 512
    save_csv(modes, out, "all_modal_curves_source_data.csv")
    for letter, length in [("c", 40), ("d", 64)]:
        save_csv(modes.loc[modes.full_length_mm == length], out, f"panel_{letter}_source_data.csv")
    save_csv(pd.DataFrame(modal_rows), out, "interface_mode_source_data.csv")

    # All four actually computed truncation orders, without fitted intermediates.
    raw, maxima = [], []
    gain_part = gain.loc[(gain.role == "B") & np.isclose(gain.sigma_f_S_m, .18)]
    side_part = side.loc[(side.role == "B") & np.isclose(side.sigma_f_S_m, .18)
        & side.contrast.isin(["upper_7_to_1_at_lower_1", "upper_7_to_1_at_lower_7"])]
    assert len(gain_part) == 32 and len(side_part) == 16
    for length in [40, 64]:
        for k in [0, 8, 32, 64]:
            for output, frame, error, key, prediction, target in [
                ("gain", gain_part, "gain_relative_error", "state_id", "gain_prediction_S", "gain_COMSOL_S"),
                ("upper_side", side_part, "relative_amplitude_error", "contrast", "predicted_side_effect_S", "COMSOL_side_effect_S")]:
                rows = frame.loc[(frame.K == k) & (frame.full_length_mm == length)]
                assert len(rows) == (4 if output == "gain" else 2)
                for _, row in rows.iterrows():
                    raw.append({"output": output, "K": k, "state_or_contrast": row[key],
                        "prediction_S": row[prediction], "COMSOL_S": row[target],
                        "error_percent": 100*row[error], "full_length_mm": length, "role": "B", "sigma_f_S_m": .18})
                maxima.append({"output": output, "K": k, "maximum_error_percent": 100*rows[error].max(),
                    "full_length_mm": length, "constituent_count": len(rows),
                    "aggregation": "Maximum over all designated constituents"})
    errors = pd.DataFrame(maxima)
    assert (errors.maximum_error_percent > 0).all()
    save_csv(pd.DataFrame(raw), out, "error_constituent_source_data.csv")
    for letter, length in [("e", 40), ("f", 64)]:
        save_csv(errors.loc[errors.full_length_mm == length], out, f"panel_{letter}_source_data.csv")
    return valid, modes, errors, files, checks


def decision_panel(ax, all_states, length, legend=False):
    frame = all_states.loc[all_states.full_length_mm == length]
    for role in ["A", "B"]:
        for sigma in [.18, 10.]:
            subset = frame.loc[(frame.role == role) & np.isclose(frame.sigma_f_S_m, sigma)]
            ax.scatter(subset.bound_percent, subset.actual_gain_error_percent,
                marker="o" if role == "A" else "^", s=15,
                facecolor="white" if sigma == .18 else COLORS[role],
                edgecolor=COLORS[role], linewidth=.8, zorder=3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set(xlim=(.05, 30), ylim=(.001, 30))
    ax.xaxis.set_major_locator(FixedLocator([.1, 1, 2, 10]))
    ax.xaxis.set_major_formatter(FixedFormatter(["0.1", "1", "2", "10"]))
    ax.yaxis.set_major_locator(FixedLocator([.001, .1, 2, 10]))
    ax.yaxis.set_major_formatter(FixedFormatter([r"$10^{-3}$", "0.1", "2", "10"]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())
    ax.axvline(2, color="black", ls="--", lw=.65, zorder=0)
    ax.axhline(2, color="black", ls="--", lw=.65, zorder=0)
    ax.set_xlabel(r"Prospective bound, $\eta$ (%)", labelpad=2)
    ax.set_ylabel(r"Gain error, $E_g$ (%)", labelpad=3)
    ax.text(.42, .97, fr"$L={length}$ mm", transform=ax.transAxes,
        ha="center", va="top", fontsize=7.3)
    if legend:
        handles = []
        for role in ["A", "B"]:
            for sigma in [.18, 10.]:
                handles.append(Line2D([], [], ls="none", marker="o" if role == "A" else "^", ms=3,
                    markeredgecolor=COLORS[role], markerfacecolor="white" if sigma == .18 else COLORS[role],
                    label=f"{role}, {sigma:g}"))
        ax.legend(handles=handles, loc="lower right", ncol=2, fontsize=6.3,
            columnspacing=.55, handlelength=.75, handletextpad=.3, labelspacing=.1, borderaxespad=.25)


def modal_panel(ax, modes, length):
    for role in ["A", "B"]:
        row = modes.loc[(modes.role == role) & (modes.full_length_mm == length)].sort_values("n")
        ax.plot(row.n, row.cumulative_weight_fraction, color=COLORS[role], ls="-", lw=1.15,
            marker="o" if role == "A" else "^", markevery=[0, 1, 3, 7, 15, 63],
            ms=2.8, markerfacecolor="white", markeredgewidth=.6, label=role)
    ax.set_xscale("log", base=2)
    ax.set(xlim=(.95, 68), ylim=(0, 1.05))
    ax.xaxis.set_major_locator(FixedLocator([1, 2, 4, 8, 16, 32, 64]))
    ax.xaxis.set_major_formatter(FixedFormatter(["1", "2", "4", "8", "16", "32", "64"]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_yticks([0, .5, 1])
    ax.axvline(8, color="black", ls=":", lw=.6, zorder=0)
    ax.set_xlabel(r"$n$", labelpad=2)
    ax.set_ylabel("Cumulative\n" + r"response weight, $C_n$", labelpad=3)
    ax.text(.97, .08, r"P11, $\sigma_f=0.18$ S/m", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=6.8)


def error_panel(ax, errors, length, legend=False):
    errors = errors.loc[errors.full_length_mm == length]
    for output, color, marker, line, label in [
        ("gain", COLORS["A"], "o", "-", "Gain"),
        ("upper_side", COLORS["B"], "^", "--", "Upper side")]:
        values = errors.loc[errors.output == output].sort_values("K")
        ax.plot(values.K, values.maximum_error_percent, color=color, marker=marker, ls=line,
            ms=3, lw=1.15, markerfacecolor="white", markeredgewidth=.7, label=label)
    ax.set_yscale("log")
    ax.set(xlim=(-2, 66), ylim=(1e-7, 300))
    ax.set_xticks([0, 8, 32, 64])
    ax.yaxis.set_major_locator(FixedLocator([1e-6, .01, 100]))
    ax.yaxis.set_major_formatter(FixedFormatter([r"$10^{-6}$", r"$10^{-2}$", r"$10^2$"]))
    ax.yaxis.set_minor_locator(NullLocator())
    if length == 64:
        # Label the two actual K=8 values without tolerance guide lines.
        gain_value = float(errors.loc[(errors.output == "gain") & (errors.K == 8), "maximum_error_percent"].iloc[0])
        side_value = float(errors.loc[(errors.output == "upper_side") & (errors.K == 8), "maximum_error_percent"].iloc[0])
        ax.annotate(f"{side_value:.2f}%", xy=(8, side_value), xytext=(7, 0),
            textcoords="offset points", ha="left", va="bottom", fontsize=6.2, color="black")
        ax.annotate(f"{gain_value:.3f}%", xy=(8, gain_value), xytext=(0, -6),
            textcoords="offset points", ha="center", va="top", fontsize=6.2, color="black")
    ax.set_xlabel(r"$K$", labelpad=2)
    ax.set_ylabel(r"$E_{\max}$ (%)", labelpad=3)
    ax.text(.97, .14, r"B, $\sigma_f=0.18$ S/m", transform=ax.transAxes,
        ha="right", va="center", fontsize=6.7)
    if legend:
        ax.legend(loc="upper right", ncol=2, fontsize=6.2, handlelength=1.7,
            columnspacing=.8, labelspacing=.15, borderaxespad=.25)


CAPTION = (ROOT / "source_data/Fig07_applicability/caption.txt").read_text(encoding="utf-8").strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=OUT)
    args = parser.parse_args()
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    notes_path = out / "figure_notes.json"
    notes = json.loads(notes_path.read_text(encoding="utf-8")) if notes_path.exists() else {}
    notes.update({"status": "data_preparation_before_rendering", "planned_size_mm": [183, 120],
        "core_conclusion": "Prospective adequacy is output dependent; a gain tolerance does not certify controlled small peripheral effects.",
        "panel_contract": {
            "a": "All 16 L40 states: M0 finite-M64 prospective bound versus later full-COMSOL gain error.",
            "b": "All 16 L64 states with identical axes and encodings; 32 new states in a+b without deletion.",
            "c": "L40, A/B, common-sigma prospective M0 response curves at actual modes n=1..64.",
            "d": "L64, A/B, common-sigma prospective M0 response curves at actual modes n=1..64.",
            "e": "K=0,8,32,64 for all four gain states and both upper-side contrasts of L40/B/0.18.",
            "f": "K=0,8,32,64 for all four gain states and both upper-side contrasts of L64/B/0.18."
        }})
    notes_path.write_text(json.dumps(notes, ensure_ascii=False, indent=2), encoding="utf-8")
    states, modes, errors, fields, checks = prepare(out)
    apply_style()
    plt.rcParams.update({"font.size": 7.3, "axes.labelsize": 7.5,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5})
    fig = plt.figure(figsize=(183/25.4, 120/25.4))
    frame_width, frame_height = .36, .18
    axes = [fig.add_axes([x, y, frame_width, frame_height])
        for x, y in [(.105, .75), (.615, .75), (.105, .435), (.615, .435), (.105, .12), (.615, .12)]]
    decision_panel(axes[0], states, 40, legend=True)
    decision_panel(axes[1], states, 64)
    modal_panel(axes[2], modes, 40)
    modal_panel(axes[3], modes, 64)
    error_panel(axes[4], errors, 40, legend=True)
    error_panel(axes[5], errors, 64)
    for ax, letter in zip(axes, "abcdef"):
        panel_label(ax, letter, x=-.21, y=1.105).set_fontsize(8)
    preview = export_figure(fig, out, "figure07_applicability")
    plt.close(fig)
    # Crop the exported full figure without resizing or altering panel geometry.
    from PIL import Image
    with Image.open(preview) as full:
        width, height = full.size
        crop_box = (round(.525*width), round(.65*height), width, round(.975*height))
        full.crop(crop_box).save(out / "figure07_panel_f_preview.png")
    (out / "caption.txt").write_text(CAPTION, encoding="utf-8")
    notes.update({"status": "rendered_pending_visual_check", "preview": str(preview),
        "caption_words": len(CAPTION.split()), "backend": "python",
        "version_change": "Relative to v7, only the existing n=2 data points in panels c/d receive open-circle/open-triangle markers. Marker modes are 1,2,4,8,16,64. All data, lines, labels, axes and layout remain unchanged; prior versions are retained.",
        "panel_f_preview": str(out / "figure07_panel_f_preview.png"),
        "panel_f_preview_crop": {"normalized_full_image_box": [.525, .65, 1, .975], "resized": False},
        "export": {"width_mm": 183, "height_mm": 120, "TIFF_dpi": 600, "PNG_preview_dpi": 220},
        "source_files": [str(SOURCE / name) for name in ["applicability_validation.csv",
            "applicability_predictions.csv", "closure_new_source_comparison.csv",
            "higher_order_side_comparison.csv", "port_coefficients.csv",
            "source_scale_response_pretruth.csv"]] + fields,
        "source_field_definition": "Saved prospective M0 combination-model interface potentials. M0 retains the fracture-bearing core; these are not the unfractured background potentials.",
        "modal_formula": "Q_s,n=t*N_n*deltaY_s,n*v_M0,s,n^2/DeltaV^2; C(n)=sum_{s,m=1..n}|Q_s,m|/sum_{s,m=1..64}|Q_s,m|. No full-target fields used.",
        "modal_reconstruction_checks": checks,
        "numerical_scope": "Conditional finite-M64 intervals in a common FE space and electrode drive, with stated energy-domination conditions; complete COMSOL is external validation. Target/core meshes and truncation remain separate checks; six representative targets were refined.",
        "data_accounting": {"plotted_new_target_states": 32, "states_each_length_panel": 16,
            "deleted_states": 0, "modal_curves": 4, "actual_modes_per_curve": 64,
            "interface_mode_source_rows": 512, "computed_K": [0, 8, 32, 64],
            "panels_e_f_constituents": 48},
        "plot_integrity": "No target-category lines, dummy data, fitted equalities, 1:1 prediction line, smoothed model curves, hidden failed states, extra solves or hashes. Line segments in e/f follow the actual ordered truncation parameter.",
        "style": "Times New Roman, black text, pure blue/red, thin full frames, inward ticks. L01-F11 provides six short-wide axes in a 3x2 arrangement, with columns for fixed length conditions and response/error rows.",
        "style_source_review": {
            "primary_template": "L01 Figure 11, six short-wide quantitative panels in a 3x2 arrangement",
            "primary_template_path": "F:/Ling_Skill/Ling_Skill_V0.2/assets/L01/figure_11.png",
            "independently_viewed": ["L05 Figure 7", "L01 Figure 9", "L01 Figure 10", "L01 Figure 11"],
            "secondary_influences": "The former 2x2 L05 arrangement remains in v3 and is not the template claimed for this label-revision version.",
            "scientific_adaptation": "Axes retain the present project's real bounds, integer modes and computed K orders; no parameter axes or trajectories are invented to imitate reference shapes.",
            "frame_dimensions_mm": [183*frame_width, 120*frame_height]
        },
        "interpretation_limits": [
            "Not confirmed does not mean failed.",
            "A finite-model bound is not a confidence limit for continuous or complete-COMSOL truth.",
            "Only a specified common 0.18 S/m P11 family is shown in response panels c/d, matching the controlled error families in e/f, while all32 target states remain in a+b.",
            "Cumulative response weights are neither physical energy fractions nor unique causal shares.",
            "K counts corrected nonzero response modes per interface, not total model degrees of freedom."]})
    notes_path.write_text(json.dumps(notes, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"preview": str(preview), "caption_words": len(CAPTION.split()),
        "states": len(states), "real_modal_rows": len(modes)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
