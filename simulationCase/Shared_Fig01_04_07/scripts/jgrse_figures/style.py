"""Shared publication rendering settings for the layered-rock figures."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

COLORS = {
    'A': '#0000FF', 'B': '#FF0000', 'low': '#FFFFFF',
    'high': '#C0C0C0', 'fracture': '#008000', 'neutral': '#000000',
    'M0': '#000000', 'M8': '#0000FF', 'M32': '#FF0000', 'M64': '#008000',
    'upper': '#0000FF', 'lower': '#FF0000', 'both': '#000000',
}
WIDTH_MM = 183

def apply_style():
    plt.rcParams.update({
        'font.family': 'serif', 'font.serif': ['Times New Roman'],
        'font.size': 8, 'axes.labelsize': 8.5, 'axes.titlesize': 8.5,
        'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 7.5,
        'axes.titleweight': 'normal', 'axes.titlepad': 4,
        'axes.spines.top': True, 'axes.spines.right': True,
        'axes.linewidth': .6, 'lines.linewidth': 1.0,
        'axes.edgecolor': '#000000', 'axes.labelcolor': '#000000',
        'text.color': '#000000', 'xtick.color': '#000000', 'ytick.color': '#000000',
        'xtick.direction': 'in', 'ytick.direction': 'in',
        'xtick.top': True, 'ytick.right': True,
        'xtick.major.width': .6, 'ytick.major.width': .6,
        'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
        'legend.frameon': False, 'legend.handlelength': 1.6,
        'svg.fonttype': 'none', 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'figure.facecolor': 'white', 'savefig.facecolor': 'white',
        'mathtext.fontset': 'custom', 'mathtext.rm': 'Times New Roman',
        'mathtext.it': 'Times New Roman:italic',
        'mathtext.bf': 'Times New Roman:bold', 'mathtext.fallback': 'stix',
        'text.usetex': False,
        'axes.unicode_minus': True,
    })

def panel_label(ax, label, x=-.14, y=1.08):
    return ax.text(x, y, f'({label})', transform=ax.transAxes, fontsize=9,
                   fontstyle='italic', fontweight='normal', color='black',
                   va='top', ha='left', clip_on=False)

def export_figure(fig, outdir, stem):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f'{stem}.svg')
    fig.savefig(outdir / f'{stem}.pdf')
    fig.savefig(outdir / f'{stem}.png', dpi=220)
    fig.savefig(outdir / f'{stem}.tiff', dpi=600,
                pil_kwargs={'compression': 'tiff_lzw'})
    return outdir / f'{stem}.png'
