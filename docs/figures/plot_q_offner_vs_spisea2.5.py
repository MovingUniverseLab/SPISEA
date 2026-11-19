#!/usr/bin/env python
"""
Generate docs/figures/q_offner_vs_spisea2.5.png

Offner ``q_power_at_mass`` vs SPISEA v2.5 ``MultiplicityUnresolved``.

Run from the repository root::

    python docs/figures/plot_q_offner_vs_spisea2.5.py
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
    _TABLE1_COLOR, _TABLE1_GAMMA, bd_shade, finish_panel, gamma_step_masses,
    save, table_xy, two_axes,
)


def main():
    offner = multiplicity.MultiplicityUnresolvedOffner2023()
    lu = multiplicity.MultiplicityUnresolved()

    m_wide = np.logspace(np.log10(0.012), np.log10(40.0), 800)
    g_off = offner.q_power_at_mass(m_wide)
    m_step = gamma_step_masses()
    g_lu = lu.q_power_at_mass(m_step)
    m_tab, g_tab, err_tab = table_xy(_TABLE1_GAMMA)

    fig, axes = two_axes(r'Offner 2023 vs SPISEA v2.5: mass-ratio index $\gamma$')
    ylabel = r'$\gamma$  ($P(q)\propto q^{\gamma}$)'
    for ax, xlim, title in (
        (axes[0], _BD_XLIM, 'Brown-dwarf regime'),
        (axes[1], _FULL_XLIM, 'BD through O'),
    ):
        bd_shade(ax, xlim)
        ax.axhline(0.0, color='0.55', ls=':', lw=1.1, zorder=2)
        ax.plot(m_step, g_lu, color='0.25', ls='--', lw=1.6, zorder=3)
        ax.plot(m_wide, g_off, color=_OFFNER_COLOR, ls='-', lw=2.4, zorder=4)
        ax.errorbar(m_tab, g_tab, yerr=err_tab, fmt='o', color=_TABLE1_COLOR,
                    ms=5.5, mfc='white', mew=1.3, elinewidth=1.1, capsize=2.5,
                    zorder=5)
        finish_panel(ax, xlim, (-2.2, 7.0), title, ylabel)

    legend_handles = [
        Line2D([0], [0], color='0.25', ls='--', lw=1.6,
               label=r'SPISEA v2.5 ($\gamma=6.1$ / $-0.4$)'),
        Line2D([0], [0], color=_OFFNER_COLOR, ls='-', lw=2.4,
               label=r'Offner err-wt logistic in log $M$'),
        Line2D([0], [0], marker='o', color=_TABLE1_COLOR, ls='none',
               mfc='white', mew=1.3, ms=6, label=r'Table 1 $\gamma_\mathrm{trunc}$'),
        Patch(facecolor=_BD_SHADE, edgecolor='none', alpha=0.8,
              label=r'BD ($M\leq 0.08$)'),
    ]
    axes[1].legend(handles=legend_handles, loc='upper right', fontsize=8,
                   frameon=True, fancybox=False, edgecolor='0.7')
    save(fig, 'q_offner_vs_spisea2.5.png')


if __name__ == '__main__':
    main()
