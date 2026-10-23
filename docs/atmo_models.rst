.. _atmo_models:

========================================
Atmosphere Model Object
========================================
Stellar atmosphere models are defined as functions in
``spisea/atmospheres.py``. These can be called by::

 from spisea import atmospheres
 atmo = atmospheres.<function_name>

To call an atmosphere for a particular star, user must define the
metallicity ([Z]), temperature (in K), and gravity (in cgs)::

  spectrum = atmo(metallicity=0, temperature=5800, gravity=4)
  from astropy import units
  wave = spectrum.waveset.to(units.AA).value   # Wavelength in Angstroms
  flux = spectrum(spectrum.waveset)  # Flux (synphot Quantity; typically FLAM)

The atmosphere function is an input for the :ref:`isochrone_objects`,
and will automatically be used to define the
spectrum of each star in the isochrone model.

SPISEA uses the stsynphot/synphot stack to extract the model atmosphere,
and the output spectrum is a synphot ``SourceSpectrum`` built from the CDBS grid via ``stsynphot.grid_to_spec``.

Below is a table of atmosphere model grids currently supported by
SPISEA. Note that the resolution column reports the original
resolution of the atmosphere model grid. These are available in the
`spisea_cdbs_highres.tar.gz` file on the installation page. However, the default grid SPISEA
uses has degraded the resolution of all atmosphere grids to R = 250
(the `spisea_cdbs.tar.gz` file). 

.. list-table:: Atmosphere Models
   :header-rows: 2
   :widths: 25 15 12 15 12 18 20
 
   * - Model Name
     - T\ :sub:`eff` Range (K)
     - log *g* Range (cgs)
     - Metallicity Range [Fe/H]
     - λ Range (μm)
     - Resolution\ :sup:`a` λ/Δλ
     - Ref
   * - ``get_merged_atmosphere``
     - 250 – 50000
     - \ :sup:`b`
     - \ :sup:`b`
     - \ :sup:`b`
     - \ :sup:`b`
     - Appendix B; Hosek et al. (2020)
   * - ``get_castelli_atmosphere``
     - 3500 – 50000
     - 0 – 5.0
     - -2.5 – 0.2
     - 0.1 – 10
     - ~250
     - Castelli & Kurucz (2004)
   * - ``get_phoenixv16_atmosphere``
     - 2300 – 12000
     - 0.0 – 6.0
     - -4.0 – +1.0
     - 0.05 – 5.5
     - 100,000 – 500,000
     - Husser et al. (2013)
   * - ``get_BTSettl_2015_atmosphere``
     - 1200 – 7000
     - 2.5 – 5.5
     - 0
     - 0.01 – 30
     - 2000 – 700,000
     - Baraffe et al. (2015)
   * - ``get_BTSettl_atmosphere``\ :sup:`d`
     - 2600 – 7000
     - 4.5 – 5.5
     - -2.5 – 0.5
     - 0.1 – 6.9
     - 20,000 – 250,000
     - Allard et al. (2012b,a)
   * - ``get_kurucz_atmosphere``
     - 3000 – 50000
     - 0 – 5.0
     - -5.0 – 1.0
     - 0.1 – 10
     - ~250
     - \ :sup:`c`
   * - ``get_phoenix_atmosphere``
     - 2100 – 69000
     -
     - -4.0 – 0.5
     - 0.001 – 995
     - ~280
     - Allard et al. (2003, 2007)
   * - ``get_Phillips2020_atmosphere``
     - 200 - 3000
     - 2.5 - 5.5
     - 0
     - 0.2 - ~1980.2
     - 0.5 - 5000
     - Phillips et al. (2020)
   * - ``get_Meisner2023_atmosphere``
     - 250 - 1200
     - 2.5 - 5.5
     - -1.0 - 0.3
     - 0.2 - 30
     - ~3000
     - Meisner et al. (2023)
   * - ``get_wd_atmosphere``\ :sup:`e`
     - –
     - –
     - –
     - 0.1 – 3.0
     - 200 – 500,000
     - Koester (2010)
   * - ``get_bb_atmosphere``
     - –
     - –
     - –
     - –
     - –
     - Blackbody Spectrum
 
.. rubric:: Footnotes
 
:sup:`a` Resolution column reports the original resolution of the atmosphere model grid.
The default SPISEA grid degrades all atmosphere grids to R = 250 (``spisea_cdbs.tar.gz``).
 
:sup:`b` See Appendix B; values depend on the underlying model selected by ``get_merged_atmosphere``.
 
:sup:`c` Kurucz (1993); see CDBS documentation.
 
:sup:`d` Solar metallicity only for BTSettl 2015.
 
:sup:`e` White dwarf atmospheres only.

Sonora, PHOENIX and ATLAS suite
-----------------------------------------
``get_merged_atmosphere_atlas_phoenix_sonora`` adds the C/O = 1 equilibrium
Flame Skimmer grid to the stellar atmosphere suite. Below 2500 K, requests
with log(g) >= 2 use Flame through 2200 K and the precomputed Flame/PHOENIX
flux blend between 2200 and 2500 K. Lower-gravity requests use PHOENIX.
The log(g) division follows the lower Flame grid boundary; it is a coverage
choice, not a general physical dwarf/giant classification.

Cool high-gravity requests below [M/H] = -1 explicitly use a -1 Sonora
atmosphere proxy and warn. PHOENIX requests retain their requested metallicity.
At 2500-5000 K the suite uses PHOENIX, at 5000-5500 K the ATLAS/PHOENIX
blend, and at higher temperatures ATLAS. The separate
``get_merged_atmosphere_atlas_phoenix_sonora_w_bb_supplement`` wrapper adds
white-dwarf and blackbody fallbacks. The matching ``_grid`` function reports
the same gravity-dependent choices, including the proxy coordinates.

Model-specific loaders first try the requested spectrum. On an out-of-bounds
lookup, Flame, Flame/PHOENIX and PHOENIX clip temperature and gravity to their
common interpolation support while preserving an in-range metallicity.
This numerical fallback prevents a quarter-dex metallicity from silently
snapping to a native node. Parameter changes produce warnings; missing or
corrupt files are allowed to raise their original errors.

Clipping is an approximation: PHOENIX cannot supply a genuinely cooler
spectrum below its 2300-K floor, and Flame has temperature-dependent gravity
limits. The Flame/PHOENIX blend has only the weighted PHOENIX contribution
outside Flame's native wavelength support, underestimating blue flux there.
The low-metallicity proxy does not extend native Sonora metallicity coverage.

Model Atmosphere Classes
-------------------------
.. autofunction:: atmospheres.get_merged_atmosphere

.. autofunction:: atmospheres.get_flameskimmer_atmosphere

.. autofunction:: atmospheres.get_flameskimmer_phoenix_atmosphere

.. autofunction:: atmospheres.get_merged_atmosphere_atlas_phoenix_sonora

.. autofunction:: atmospheres.get_merged_atmosphere_atlas_phoenix_sonora_w_bb_supplement

.. autofunction:: atmospheres.get_wd_atmosphere

.. autofunction:: atmospheres.get_bb_atmosphere
   
.. autofunction:: atmospheres.get_castelli_atmosphere

.. autofunction:: atmospheres.get_phoenixv16_atmosphere
		  
.. autofunction:: atmospheres.get_BTSettl_2015_atmosphere
		  
.. autofunction:: atmospheres.get_BTSettl_atmosphere
		 
.. autofunction:: atmospheres.get_wdKoester_atmosphere

.. autofunction:: atmospheres.get_kurucz_atmosphere

.. autofunction:: atmospheres.get_phoenix_atmosphere

.. autofunction:: atmospheres.get_Phillips2020_atmosphere

.. autofunction:: atmospheres.get_Meisner2023_atmosphere

