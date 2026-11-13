#!/usr/bin/env python
"""
Generate docs/figures/sep_offner_vs_spisea2.5.png

Offner ``a_mean`` vs SPISEA v2.5 ``MultiplicityResolvedDK.a_mean``.

Run from the repository root::

    python docs/figures/plot_sep_offner_vs_spisea2.5.py
"""
import os
import sys

import numpy as np

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
_FIGDIR = os.path.dirname(os.path.abspath(__file__))
if _FIGDIR not in sys.path:
    sys.path.insert(0, _FIGDIR)

from spisea.imf import multiplicity
from _offner_vs_spisea25 import (
    Line2D, Patch, _BD_SHADE, _BD_XLIM, _FULL_XLIM, _OFFNER_COLOR,
    _TABLE1_A_ALL, _TABLE1_COLOR, _TABLE2_MU, bd_shade, finish_panel, geom,
    save, table_xy, two_axes,
)


def main():
    resolved = multiplicity.MultiplicityResolvedOffner2023()
    dk = multiplicity.MultiplicityResolvedDK()

    m_wide = np.logspace(np.log10(0.012), np.log10(40.0), 800)
    a_off = resolved.a_mean(m_wide)
    a_lu = dk.a_mean(m_wide)
    m_tab, a_tab, err_tab = table_xy(_TABLE1_A_ALL)
    m_t2 = np.array([geom(r[1], r[2]) for r in _TABLE2_MU])
    a_t2 = np.array([r[3] for r in _TABLE2_MU], dtype=float)

    fig, axes = two_axes(
        r'Offner 2023 vs SPISEA v2.5: characteristic separation')
    ylabel = r'$\mu(a)$ (AU)'
    for ax, xlim, title in (
        (axes[0], _BD_XLIM, 'Brown-dwarf regime'),
        (axes[1], _FULL_XLIM, 'BD through O'),
    ):
        bd_shade(ax, xlim)
        ax.plot(m_wide, a_lu, color='0.25', ls='--', lw=1.6, zorder=3)
        ax.plot(m_wide, a_off, color=_OFFNER_COLOR, ls='-', lw=2.4, zorder=4)
        ax.errorbar(m_tab, a_tab, yerr=err_tab, fmt='o', color=_TABLE1_COLOR,
                    ms=5.5, mfc='white', mew=1.3, elinewidth=1.1, capsize=2.5,
                    zorder=5)
        ax.plot(m_t2, a_t2, 's', color=_OFFNER_COLOR, mfc=_OFFNER_COLOR,
                ms=7, zorder=6, mew=0.6)
        finish_panel(ax, xlim, (1.0, 400.0), title, ylabel, ylog=True)

    legend_handles = [
        Line2D([0], [0], color='0.25', ls='--', lw=1.6,
               label=r'SPISEA v2.5 mean $a$'),
        Line2D([0], [0], color=_OFFNER_COLOR, ls='-', lw=2.4,
               label=r'Offner smooth broken PL ($s=0.1$ dex)'),
        Line2D([0], [0], marker='o', color=_TABLE1_COLOR, ls='none',
               mfc='white', mew=1.3, ms=6, label=r'Table 1 $\tilde{a}_\mathrm{all}$'),
        Line2D([0], [0], marker='s', color=_OFFNER_COLOR, ls='none',
               mfc=_OFFNER_COLOR, ms=7, label=r'Table 2 $\mu$ knots'),
        Patch(facecolor=_BD_SHADE, edgecolor='none', alpha=0.8,
              label=r'BD ($M\leq 0.08$)'),
    ]
    axes[0].legend(handles=legend_handles, loc='upper left', fontsize=8,
                   frameon=True, fancybox=False, edgecolor='0.7')
    save(fig, 'sep_offner_vs_spisea2.5.png')


if __name__ == '__main__':
    main()
