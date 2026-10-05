"""Shared layout and Table 1/2 data for Offner vs SPISEA v2.5 plots."""
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


_OFFNER_COLOR = '#8b3a2a'
_TABLE1_COLOR = '#2f6db3'
_BD_SHADE = '#e8d5b5'
_FIGSIZE = (11.2, 4.6)
_DPI = 160
_BD_XLIM = (0.012, 0.20)
_FULL_XLIM = (0.015, 40.0)

# Table 1 γ_trunc (1–100 au / 1–102 au) with 1σ. Masses are geometric
# means of the tabulated M1 intervals. The L/early-T γ=2.5 text value
# is an interpolation knot only and is not plotted as a Table 1 point.
_TABLE1_GAMMA = [
    # name, M_lo, M_hi, gamma, gamma_err
    ('Fontanive+2018', 0.019, 0.058, 4.8, 2.2),
    ('Close+2003', 0.080, 0.095, 3.3, 1.2),
    ('Allen+2007', 0.06, 0.15, 1.7, 0.5),
    ('Winters mid-M', 0.15, 0.30, 0.7, 0.5),
    ('Winters early-M', 0.3, 0.6, 0.1, 0.4),
    ('Raghavan+2010', 0.75, 1.25, 0.2, 0.4),
    ('De Rosa A', 1.6, 2.4, -1.3, 0.4),
    ('MDS 3-5', 3.0, 5.0, -1.0, 0.5),
    ('MDS 5-8', 5.0, 8.0, -1.7, 0.5),
    ('MDS 8-17', 8.0, 17.0, -1.6, 0.5),
    ('Sana O', 17.0, 50.0, -1.4, 0.4),
]

# Table 1 ã_all (au) with 1σ.
_TABLE1_A_ALL = [
    ('Fontanive+2018', 0.019, 0.058, 2.9, 1.1),
    ('Close+2003', 0.080, 0.095, 3.7, 1.3),
    ('Allen+2007', 0.06, 0.15, 6.9, 1.4),
    ('Winters late-M', 0.075, 0.15, 3.9, 1.2),
    ('Winters mid-M', 0.15, 0.30, 10.0, 3.0),
    ('Winters early-M', 0.3, 0.6, 26.0, 4.0),
    ('Raghavan+2010', 0.75, 1.25, 49.0, 6.0),
    ('Tokovinin 2014b', 0.85, 1.5, 31.0, 5.0),
    ('Moe & Kratter', 1.6, 2.4, 32.0, 8.0),
    ('MDS 3-5', 3.0, 5.0, 28.0, 7.0),
    ('MDS 5-8', 5.0, 8.0, 25.0, 7.0),
    ('MDS 8-17', 8.0, 17.0, 23.0, 7.0),
    ('Sana O', 17.0, 50.0, 19.0, 6.0),
]

# Table 2 lognormal μ (au) at the three published bins.
_TABLE2_MU = [
    ('late-M', 0.075, 0.15, 4.0),
    ('early-M', 0.3, 0.6, 25.0),
    ('FGK', 0.75, 1.25, 40.0),
]


def geom(lo, hi):
    return float(np.sqrt(lo * hi))


def table_xy(rows, y_idx=3, e_idx=4):
    m = np.array([geom(r[1], r[2]) for r in rows])
    y = np.array([r[y_idx] for r in rows], dtype=float)
    err = np.array([r[e_idx] for r in rows], dtype=float)
    return m, y, err


def mean_q_from_gamma(gamma, q_min=0.01):
    """⟨q⟩ for P(q) ∝ q^γ on [q_min, 1]."""
    g = np.asarray(gamma, dtype=float)
    qmin = float(q_min)
    g_flat = np.atleast_1d(g).astype(float)
    out_flat = np.empty(g_flat.shape, dtype=float)
    near_m1 = np.abs(g_flat + 1.0) < 1e-12
    near_m2 = np.abs(g_flat + 2.0) < 1e-12
    ok = ~near_m1 & ~near_m2
    if np.any(near_m1):
        out_flat[near_m1] = (1.0 - qmin) / (-np.log(qmin))
    if np.any(near_m2):
        out_flat[near_m2] = -np.log(qmin) / (1.0 / qmin - 1.0)
    if np.any(ok):
        gp = g_flat[ok]
        num = (1.0 - np.power(qmin, gp + 2.0)) / (gp + 2.0)
        den = (1.0 - np.power(qmin, gp + 1.0)) / (gp + 1.0)
        out_flat[ok] = num / den
    out = out_flat.reshape(np.shape(g))
    return float(out) if np.isscalar(gamma) else out


def gamma_step_masses():
    """Dense sampling so the 0.08 Msun γ step renders as a vertical jump."""
    return np.concatenate([
        np.logspace(np.log10(0.012), np.log10(0.07999), 300),
        np.array([0.08, 0.08001]),
        np.logspace(np.log10(0.081), np.log10(40.0), 300),
    ])


def bd_shade(ax, xlim):
    ax.axvspan(xlim[0], 0.08, color=_BD_SHADE, alpha=0.55, zorder=0)
    ax.axvline(0.08, color='#c4a574', ls='--', lw=1.2, zorder=1)


def finish_panel(ax, xlim, ylim, title, ylabel, ylog=False):
    ax.set_xscale('log')
    if ylog:
        ax.set_yscale('log')
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(r'Primary mass $M_1$ ($M_\odot$)')
    ax.set_ylabel(ylabel)
    ax.tick_params(which='both', direction='in', top=True, right=True)


def two_axes(suptitle):
    fig, axes = plt.subplots(1, 2, figsize=_FIGSIZE, gridspec_kw={'wspace': 0.28})
    fig.suptitle(suptitle, fontsize=13, y=1.02)
    return fig, axes


def save(fig, filename):
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)
    fig.savefig(out, dpi=_DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('Wrote', out)
    return out


# Re-export legend artists so plot scripts need one import.
Patch = Patch
Line2D = Line2D
