"""Render the peripheral-mechanism pilot from saved COMSOL results only."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from style import apply_style, panel_label, export_figure, COLORS

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / 'runs/round4_extension_20261024'
OUT = ROOT / 'outputs/JGRSE_20261028_LabelRevision/pilot06_mechanism'
SIGMA = .159375


def load_data():
    fine = pd.read_csv(RUN / 'fine_factorial_comparison.csv')
    fine = fine[np.isclose(fine.sigma_f_S_m, SIGMA)].copy()
    selected = fine[fine.quantity.isin(['gA', 'gB']) & fine.effect.isin(
        ['lower_average', 'upper_average', 'joint_NN_minus_FF'])].copy()
    selected['effect_uS'] = selected.fine_value * 1e6
    selected['diagnostic_halfwidth_uS'] = (selected.diagnostic_factor *
        selected.observed_weighted_sum_underlying_mesh_changes * 1e6)
    coeffs = []
    for role in ['A', 'B']:
        for layout in ['NN', 'NF', 'FN', 'FF']:
            d = pd.read_csv(RUN / f'modes/{role}_{layout}/mode_coefficients.csv')
            d['role'] = role
            d['layout'] = layout
            coeffs.append(d)
    coeffs = pd.concat(coeffs, ignore_index=True)
    spectral = []
    for side in ['lower', 'upper']:
        near = coeffs.query("role == 'A' and layout == 'NN' and side == @side").set_index('n')
        far = coeffs.query("role == 'A' and layout == 'FF' and side == @side").set_index('n')
        for n in near.index:
            spectral.append({'side': side, 'n': int(n),
                'Y_near_S_m2': near.loc[n, 'Y_S_per_m2'],
                'Y_far_S_m2': far.loc[n, 'Y_S_per_m2'],
                'delta_Y_S_m2': near.loc[n, 'Y_S_per_m2'] - far.loc[n, 'Y_S_per_m2']})
    spectral = pd.DataFrame(spectral)
    # These curves are interpolation of stored mode amplitudes, not a new solve.
    traces = []
    x_mm = np.linspace(0, 120, 481)
    for role in ['A', 'B']:
        for side in ['lower', 'upper']:
            d = coeffs.query("role == @role and layout == 'NN' and side == @side and n > 0")
            v = np.cos(np.outer(x_mm / 120 * np.pi, d.n.to_numpy())) @ d.delta_V_cosine_V.to_numpy()
            for x, voltage in zip(x_mm, v):
                traces.append({'role': role, 'side': side, 'layout': 'NN',
                    'x_mm': x, 'zero_mean_deltaV_mV': 1000 * voltage, 'n_min': 1, 'n_max': 32})
    traces = pd.DataFrame(traces)
    energies = pd.read_csv(RUN / 'modes/modal_energy_closure.csv')
    energies = energies[energies.layout == 'NN'].copy()
    energies['nonuniform_energy_uW'] = 1e6 * energies.n1to32_power_W
    cases = json.loads((ROOT / 'runs/round4_crossover/periphery_manifest.json').read_text(encoding='utf-8-sig'))['cases']
    cases += json.loads((RUN / 'manifest.json').read_text(encoding='utf-8-sig'))['cases']
    geometry = []
    for case in cases:
        if case['configuration_role'] != 'A':
            continue
        label = {'swap_near': 'NN', 'swap_far': 'FF'}.get(case['layout'], case['layout'])
        if label not in ['NN', 'NF', 'FN', 'FF']:
            continue
        for side, key in [('lower', 'lower_high_to_low_strip_z_mm'),
                          ('upper', 'upper_low_to_high_strip_z_mm')]:
            lo, hi = case[key]
            geometry.append({'layout': label, 'side': side, 'z_start_mm': lo, 'z_end_mm': hi})
    geometry = pd.DataFrame(geometry).drop_duplicates()
    assert len(selected) == 6 and len(coeffs) == 528 and len(geometry) == 8
    assert np.allclose(spectral.query('n == 0').delta_Y_S_m2, 0, atol=1e-12)
    return fine, selected, coeffs, spectral, traces, energies, geometry


def geometry_panel(ax, geometry):
    ax.set_axis_off()
    panel_label(ax, 'a', x=-.13, y=1.10)
    ax.set_title('Peripheral layers', pad=5)
    # True z coordinates within each window; the intervening core is omitted.
    for side, bounds, y in [('upper', (90, 100), .53), ('lower', (20, 30), .07)]:
        sub = ax.inset_axes([.12, y, .86, .32])
        for i, label in enumerate(['NN', 'NF', 'FN', 'FF']):
            base = COLORS['low'] if side == 'upper' else COLORS['high']
            band = COLORS['high'] if side == 'upper' else COLORS['low']
            sub.add_patch(Rectangle((i - .33, bounds[0]), .66, 10,
                                   facecolor=base, edgecolor='black', linewidth=.45))
            row = geometry.query('layout == @label and side == @side').iloc[0]
            sub.add_patch(Rectangle((i - .33, row.z_start_mm), .66,
                                   row.z_end_mm - row.z_start_mm,
                                   facecolor=band, edgecolor='black', linewidth=.5))
        sub.set_xlim(-.6, 3.6)
        sub.set_ylim(*bounds)
        sub.set_yticks([bounds[0], bounds[1]])
        sub.set_ylabel('$z$ (mm)', labelpad=2)
        sub.spines[['top', 'right', 'bottom']].set_visible(False)
        sub.tick_params(axis='both', top=False, right=False)
        sub.tick_params(axis='x', length=0)
        sub.set_xticks(range(4), ['NN', 'NF', 'FN', 'FF'] if side == 'lower' else [])
        sub.text(.01, 1.09, side.capitalize(), transform=sub.transAxes,
                 fontsize=7.5)
        boundary = 90 if side == 'upper' else 30
        sub.axhline(boundary, color='black', linewidth=.5, linestyle='--')
    ax.text(.56, .48, 'Fixed core', ha='center', va='center',
            transform=ax.transAxes, fontsize=7.5, color='black')


def effects_panel(ax, data):
    panel_label(ax, 'b', x=-.18, y=1.10)
    ax.set_title('Change in gain', pad=5)
    specs = [('upper_average', 'Upper'), ('lower_average', 'Lower'),
             ('joint_NN_minus_FF', 'Joint')]
    for role, dx, hatch in [('A', -.18, None), ('B', .18, '///')]:
        for j, (effect, label) in enumerate(specs):
            row = data.loc[(data.quantity == 'g' + role) & (data.effect == effect)].iloc[0]
            x = j + dx
            value = row.effect_uS
            ax.bar(x, value, width=.30, facecolor='white', edgecolor=COLORS[role],
                   hatch=hatch, linewidth=1.0, label=role if j == 0 else None, zorder=2)
            ax.errorbar(x, value, yerr=row.diagnostic_halfwidth_uS, fmt='none',
                        color='black', ecolor='black', capsize=2.2, elinewidth=.7, zorder=4)
            if effect == 'joint_NN_minus_FF' and role == 'B':
                ax.annotate('Unresolved', xy=(x, value - row.diagnostic_halfwidth_uS),
                            xytext=(2.05, -1.22), fontsize=7, color='black',
                            ha='center', va='center',
                            arrowprops={'arrowstyle': '->', 'color': 'black', 'lw': .5})
    ax.axhline(0, color='black', lw=.6, zorder=1)
    ax.set_xlim(-.6, 2.6)
    ax.set_ylim(-1.75, 6.0)
    ax.set_yticks([-1, 0, 2, 4, 6])
    ax.set_xticks(range(3), [s[1] for s in specs])
    ax.set_ylabel(r'$\delta g$ ($\mu$S)', labelpad=2)
    ax.legend(loc='upper right', ncol=2, handlelength=1.2, columnspacing=1.1)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    apply_style()
    fine, selected, coeffs, spectral, traces, energies, geometry = load_data()
    for name, data in [('all_fine_factorial_contrasts', fine), ('panel_b_effects', selected),
                       ('all_mode_coefficients', coeffs), ('panel_c_admittances', spectral),
                       ('panel_d_voltage_reconstruction', traces), ('panel_e_energies', energies),
                       ('panel_a_geometry', geometry)]:
        data.to_csv(OUT / (name + '.csv'), index=False, encoding='utf-8-sig')
    fig = plt.figure(figsize=(183 / 25.4, 143 / 25.4))
    a = fig.add_axes([.07, .58, .39, .34])
    b = fig.add_axes([.61, .58, .36, .34])
    c = fig.add_axes([.08, .12, .235, .32])
    d = fig.add_axes([.407, .12, .235, .32])
    e = fig.add_axes([.735, .12, .235, .32])
    geometry_panel(a, geometry)
    effects_panel(b, selected)

    panel_label(c, 'c', x=-.24, y=1.13)
    c.set_title('Outer response')
    for side, ls in [('upper', '-'), ('lower', '--')]:
        s = spectral[spectral.side == side]
        c.plot(s.n, s.delta_Y_S_m2, color=COLORS[side],
               linestyle=ls, label=side.capitalize())
    c.axhline(0, color='black', lw=.5, ls=':')
    c.plot([0], [0], 'o', color='black', markerfacecolor='white', ms=3)
    c.set_xlim(0, 32)
    c.set_xticks([0, 8, 16, 24, 32])
    c.set_xlabel('$n$')
    c.set_ylabel(r'$\Delta Y_n$ (S m$^{-2}$)')
    c.legend(loc='lower right', fontsize=7)

    panel_label(d, 'd', x=-.24, y=1.13)
    d.set_title('Interface potential (NN)')
    for role in ['A', 'B']:
        for side, ls in [('upper', '-'), ('lower', '--')]:
            tr = traces.query('role == @role and side == @side')
            d.plot(tr.x_mm, tr.zero_mean_deltaV_mV, color=COLORS[role],
                   linestyle=ls, label=f'{role}, {side}')
    d.axhline(0, color='black', lw=.5, ls=':')
    d.set_xlim(0, 120)
    d.set_xticks([0, 60, 120])
    d.set_xlabel('$x$ (mm)')
    d.set_ylabel(r'$\delta V$ (mV)')
    role_legend = d.legend(handles=[Line2D([0], [0], color=COLORS['A'], label='A'),
                      Line2D([0], [0], color=COLORS['B'], label='B')],
             loc='lower left', fontsize=7, handlelength=1.5, labelspacing=.3)
    d.add_artist(role_legend)
    d.legend(handles=[Line2D([0], [0], color='black', linestyle='-', label='Upper'),
                      Line2D([0], [0], color='black', linestyle='--', label='Lower')],
             loc='lower right', fontsize=7, handlelength=1.5, labelspacing=.3)

    panel_label(e, 'e', x=-.24, y=1.13)
    e.set_title('Peripheral energy (NN)')
    for i, role in enumerate(['A', 'B']):
        for side, dx, marker in [('lower', -.18, 'o'), ('upper', .18, '^')]:
            value = energies.query('role == @role and side == @side').iloc[0].nonuniform_energy_uW
            e.vlines(i + dx, 0, value, color=COLORS[role], linewidth=.7)
            e.plot(i + dx, value, marker=marker, color=COLORS[role],
                   markerfacecolor='white', markersize=4.5, markeredgewidth=.8)
            e.text(i + dx, value+1.2, f'{value:.2f}', ha='center', fontsize=7, color='black')
    e.set_xticks([0, 1], ['A', 'B'])
    e.set_xlim(-.6, 1.6)
    e.set_ylim(-2, 46)
    e.set_yticks([0, 20, 40])
    e.set_ylabel(r'$P^{\prime}$ ($\mu$W)')
    e.legend(handles=[Line2D([], [], color='black', markerfacecolor='white', marker='o', linestyle='none', label='Lower'),
                      Line2D([], [], color='black', markerfacecolor='white', marker='^', linestyle='none', label='Upper')],
             loc='upper right', fontsize=7, handletextpad=.25)
    export_figure(fig, OUT, 'pilot06_mechanism')
    plt.close(fig)
    # Standalone view of the revised panel requested by the user.
    panel_fig, panel_ax = plt.subplots(figsize=(96 / 25.4, 72 / 25.4))
    panel_fig.subplots_adjust(left=.19, bottom=.17, right=.96, top=.83)
    effects_panel(panel_ax, selected)
    export_figure(panel_fig, OUT, 'panel_b_effects')
    plt.close(panel_fig)
    caption = (ROOT / "source_data/Fig04_mechanism/caption.txt").read_text(encoding="utf-8").strip()
    (OUT / 'caption.txt').write_text(caption, encoding='utf-8')
    notes = {
        'original_blueprint_number': 6, 'planned_main_figure': 'standalone peripheral-mechanism figure',
        'status': 'rendered_from_existing_results', 'new_solves': 0,
        'style_revision': 'v4: only a black solid/dashed Upper/Lower legend is added to panel (d); the existing A/B legend, curves, line widths, axes, layout, and other panels are retained from v3',
        'dimensions_mm': [183, 143],
        'figure_contract': 'figure_contract.json',
        'source_files': [str(RUN / 'fine_factorial_comparison.csv'), str(RUN / 'manifest.json'),
            str(ROOT / 'runs/round4_crossover/periphery_manifest.json'),
            str(RUN / 'modes/modal_energy_closure.csv')],
        'selection_counts': {'fine_contrast_rows_exported': len(fine), 'displayed_contrasts': len(selected),
             'mode_rows_exported': len(coeffs), 'reconstructed_trace_rows': len(traces)},
        'statistics': 'Deterministic numerical states; no sampled replicates, p-values, or statistical confidence intervals.',
        'not_shown': 'Other conductivities and the small factorial interaction remain in existing data; the complete fine-grid contrast table is exported alongside the plot.',
        'field_processing': 'Only cosine reconstruction of saved coefficients; the mean mode is deliberately omitted from voltage disturbance curves to show nonuniform excitation.',
        'limitations': ['No direct claim that nonuniform perturbation energies equal causal effects.',
            'Core A and B are different; each is fixed only across its own peripheral layouts.',
            'The 3× observed mesh scale is not a bound on all discretization error.',
            'Ranking-crossover evidence is arranged separately, outside this standalone mechanism figure.'],
    }
    (OUT / 'figure_notes.json').write_text(json.dumps(notes, ensure_ascii=False, indent=2), encoding='utf-8')
    print(OUT / 'pilot06_mechanism.png')


if __name__ == '__main__':
    main()
