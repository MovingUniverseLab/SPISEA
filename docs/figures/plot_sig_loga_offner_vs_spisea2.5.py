#!/usr/bin/env python
"""
Generate docs/figures/sig_loga_offner_vs_spisea2.5.png

Offner ``sigma_log_a`` vs SPISEA v2.5
``MultiplicityResolvedDK.sigma_log_a``.

Run from the repository root::

    python docs/figures/plot_sig_loga_offner_vs_spisea2.5.py
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
    bd_shade, finish_panel, save, two_axes,
)


def main():
    resolved = multiplicity.MultiplicityResolvedOffner2023()
    dk = multiplicity.MultiplicityResolvedDK()

    m_wide = np.logspace(np.log10(0.012), np.log10(40.0), 800)
    sig_off = resolved.sigma_log_a(m_wide)
    sig_lu = dk.sigma_log_a(m_wide)
    m_t2 = np.array(resolved.sep_sig_mass, dtype=float)
    sig_t2 = np.array(resolved.sep_sig, dtype=float)

    fig, axes = two_axes(
        r'Offner 2023 vs SPISEA v2.5: $\sigma(\log_{10} a)$')
    ylabel = r'$\sigma(\log_{10} a)$'
    for ax, xlim, title in (
        (axes[0], _BD_XLIM, 'Brown-dwarf regime'),
        (axes[1], _FULL_XLIM, 'BD through O'),
    ):
        bd_shade(ax, xlim)
        ax.plot(m_wide, sig_lu, color='0.25', ls='--', lw=1.6, zorder=3)
        ax.plot(m_wide, sig_off, color=_OFFNER_COLOR, ls='-', lw=2.4, zorder=4)
        ax.plot(m_t2, sig_t2, 's', color=_OFFNER_COLOR, mfc=_OFFNER_COLOR,
                ms=7, zorder=6, mew=0.6)
        finish_panel(ax, xlim, (0.0, 2.05), title, ylabel)

    legend_handles = [
        Line2D([0], [0], color='0.25', ls='--', lw=1.6,
               label=r'SPISEA v2.5 DK $\sigma_{\log a}$'),
        Line2D([0], [0], color=_OFFNER_COLOR, ls='-', lw=2.4,
               label=r'Offner 2-param logistic $\sigma$'),
        Line2D([0], [0], marker='s', color=_OFFNER_COLOR, ls='none',
               mfc=_OFFNER_COLOR, ms=7, label=r'Table 2 $\sigma$ knots'),
        Patch(facecolor=_BD_SHADE, edgecolor='none', alpha=0.8,
              label=r'BD ($M\leq 0.08$)'),
    ]
    axes[1].legend(handles=legend_handles, loc='upper left', fontsize=8,
                   frameon=True, fancybox=False, edgecolor='0.7')
    save(fig, 'sig_loga_offner_vs_spisea2.5.png')


if __name__ == '__main__':
    main()
