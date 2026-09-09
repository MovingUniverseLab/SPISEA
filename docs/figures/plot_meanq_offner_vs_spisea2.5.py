#!/usr/bin/env python
"""
Generate docs/figures/meanq_offner_vs_spisea2.5.png

Mean q implied by ``q_power_at_mass`` on Offner vs SPISEA v2.5
``MultiplicityUnresolved``.

Run from the repository root::

    python docs/figures/plot_meanq_offner_vs_spisea2.5.py
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
    bd_shade, finish_panel, gamma_step_masses, mean_q_from_gamma, save,
    two_axes,
)


def main():
    offner = multiplicity.MultiplicityUnresolvedOffner2023()
    lu = multiplicity.MultiplicityUnresolved()

    m_wide = np.logspace(np.log10(0.012), np.log10(40.0), 800)
    q_off = mean_q_from_gamma(offner.q_power_at_mass(m_wide),
                             q_min=offner.q_min)
    m_step = gamma_step_masses()
    q_lu = mean_q_from_gamma(lu.q_power_at_mass(m_step), q_min=lu.q_min)

    fig, axes = two_axes(
        r'Offner 2023 vs SPISEA v2.5: mean mass ratio $\langle q\rangle$')
    ylabel = r'$\langle q\rangle$  on $[0.01,\,1]$'
    for ax, xlim, title in (
        (axes[0], _BD_XLIM, 'Brown-dwarf regime'),
        (axes[1], _FULL_XLIM, 'BD through O'),
    ):
        bd_shade(ax, xlim)
        ax.plot(m_step, q_lu, color='0.25', ls='--', lw=1.6, zorder=3)
        ax.plot(m_wide, q_off, color=_OFFNER_COLOR, ls='-', lw=2.4, zorder=4)
        finish_panel(ax, xlim, (0.0, 1.0), title, ylabel)

    legend_handles = [
        Line2D([0], [0], color='0.25', ls='--', lw=1.6,
               label=r'SPISEA v2.5 from $\gamma$ step'),
        Line2D([0], [0], color=_OFFNER_COLOR, ls='-', lw=2.4,
               label=r'Offner  from $\gamma(M)$ logistic'),
        Patch(facecolor=_BD_SHADE, edgecolor='none', alpha=0.8,
              label=r'BD ($M\leq 0.08$)'),
    ]
    axes[1].legend(handles=legend_handles, loc='upper right', fontsize=8,
                   frameon=True, fancybox=False, edgecolor='0.7')
    save(fig, 'meanq_offner_vs_spisea2.5.png')


if __name__ == '__main__':
    main()
