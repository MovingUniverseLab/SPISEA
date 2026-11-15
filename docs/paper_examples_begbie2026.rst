.. _paper_examples_begbie2026:

==========================================
Begbie et al. (2026) paper examples
==========================================
These notebooks recreate figures from
`Begbie et al. (2026)
<https://arxiv.org/abs/2607.14292>`_, which adds brown-dwarf
evolution and atmosphere support to SPISEA.

They use ``get_merged_atmosphere_v1`` on purpose, because that is what
the published figures used (Figure 5 also relies on the IsochronePhot
default atmosphere, which is v1). For new work, prefer
``get_merged_atmosphere_v3`` (or ``get_merged_atmosphere_v3bb`` for
compact objects), as in the tutorials.

.. toctree::
   :maxdepth: 1
   :titlesonly:

   Figure 2: Brown-dwarf IMFs <notebooks/papers/begbie2026/fig02_bd_imfs>
   Figure 3: Atmosphere interpolation <notebooks/papers/begbie2026/fig03_atmosphere_interpolation>
   Figure 4: Evolution interpolation <notebooks/papers/begbie2026/fig04_evolution_interpolation>
   Figure 5: HRD ages <notebooks/papers/begbie2026/fig05_hrd_ages>
   Figure 6: Pleiades (Gaia) <notebooks/papers/begbie2026/fig06_pleiades_gaia>
   Figure 7: Pleiades (UKIDSS) <notebooks/papers/begbie2026/fig07_pleiades_ukidss>
   Figure 8: Upper Sco <notebooks/papers/begbie2026/fig08_upper_sco>
   Figure 9: M44 <notebooks/papers/begbie2026/fig09_m44>
   Figure 10: Mass fraction <notebooks/papers/begbie2026/fig10_mass_fraction>
   Figure 11: Companion counts <notebooks/papers/begbie2026/fig11_companion_counts>
   Figure 12: Brown-dwarf separations <notebooks/papers/begbie2026/fig12_bd_separations>
