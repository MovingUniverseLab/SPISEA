#!/usr/bin/env python
"""
Generate docs/figures/csf_offner_vs_spisea2.5.png

Two-panel comparison of companion star fraction vs primary mass:
SPISEA v2.5 ``MultiplicityUnresolved.companion_star_fraction``,
Offner et al. 2023 logistic in log-mass, and Offner Table 1 CF
points.

Run from the repository root::

    python docs/figures/plot_csf_offner_vs_spisea2.5.py
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# Allow running without installing the package.
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          '..', '..'))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from spisea.imf import multiplicity


# Offner et al. 2023 Table 1: (M_lo, M_hi, CF)
_TABLE1_CF = [
    (0.019, 0.058, 0.08),
    (0.05, 0.08, 0.16),
    (0.080, 0.095, 0.19),
    (0.06, 0.15, 0.20),
    (0.075, 0.15, 0.21),
    (0.15, 0.30, 0.27),
    (0.3, 0.6, 0.38),
    (0.75, 1.25, 0.60),
    (0.85, 1.5, 0.62),
    (1.6, 2.4, 0.99),
    (3.0, 5.0, 1.28),
    (5.0, 8.0, 1.55),
    (8.0, 17.0, 1.80),
    (17.0, 50.0, 2.10),
]


def _table1_xy():
    m = np.array([np.sqrt(lo * hi) for lo, hi, _ in _TABLE1_CF])
    cf = np.array([row[2] for row in _TABLE1_CF])
    return m, cf


def _style_panel(ax, m_off, csf_off, m_lu, csf_lu, m_tab, cf_tab,
                 xlim, ylim, title):
    ax.axvspan(xlim[0], 0.08, color='#e8d5b5', alpha=0.55, zorder=0)
    ax.axvline(0.08, color='#c4a574', ls='--', lw=1.2, zorder=1)
    ax.plot(m_lu, csf_lu, color='0.25', ls='--', lw=1.6, zorder=3,
            label=r'SPISEA v2.5 $0.50\,M^{0.45}$')
    ax.plot(m_off, csf_off, color='#8b3a2a', ls='-', lw=2.4, zorder=4,
            label='Offner logistic in log M')
    ax.plot(m_tab, cf_tab, 'o', color='#2f6db3', ms=5.5, mfc='white',
            mew=1.3, zorder=5, label='Offner Table 1 CF')
    ax.set_xscale('log')
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(r'Primary mass $M_1$ ($M_\odot$)')
    ax.set_ylabel('Companion star fraction')
    ax.tick_params(which='both', direction='in', top=True, right=True)


def main():
    offner = multiplicity.MultiplicityUnresolvedOffner2023()
    lu = multiplicity.MultiplicityUnresolved()

    m_wide = np.logspace(np.log10(0.012), np.log10(40.0), 800)
    csf_off = offner.companion_star_fraction(m_wide)
    csf_lu = lu.companion_star_fraction(m_wide)
    m_tab, cf_tab = _table1_xy()

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6),
                             gridspec_kw={'wspace': 0.28})
    fig.suptitle('Offner 2023 vs SPISEA v2.5: companion star fraction',
                 fontsize=13, y=1.02)

    _style_panel(
        axes[0], m_wide, csf_off, m_wide, csf_lu, m_tab, cf_tab,
        xlim=(0.012, 0.20), ylim=(0.0, 0.50),
        title='Brown-dwarf regime')
    _style_panel(
        axes[1], m_wide, csf_off, m_wide, csf_lu, m_tab, cf_tab,
        xlim=(0.015, 20.0), ylim=(0.0, 3.0),
        title='BD through early B')

    legend_handles = [
        Line2D([0], [0], color='0.25', ls='--', lw=1.6,
               label=r'SPISEA v2.5 $0.50\,M^{0.45}$'),
        Line2D([0], [0], color='#8b3a2a', ls='-', lw=2.4,
               label='Offner logistic in log M'),
        Line2D([0], [0], marker='o', color='#2f6db3', ls='none',
               mfc='white', mew=1.3, ms=6, label='Offner Table 1 CF'),
        Patch(facecolor='#e8d5b5', edgecolor='none', alpha=0.8,
              label=r'BD ($M\leq 0.08$)'),
    ]
    axes[0].legend(handles=legend_handles, loc='upper left', fontsize=8,
                   frameon=True, fancybox=False, edgecolor='0.7')

    out = os.path.join(os.path.dirname(__file__),
                       'csf_offner_vs_spisea2.5.png')
    fig.savefig(out, dpi=160, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('Wrote', out)


if __name__ == '__main__':
    main()
