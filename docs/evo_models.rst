.. _evo_models:

========================================
Evolution Model Object
========================================

Stellar evolution models are defined as classes in
``spisea/evolution.py``. These can be called by::

  from spisea import evolution
  evo = evolution.<model_name>()

The evolution object is an input for the :ref:`isochrone_objects`. 
  
Below is a table of the evolution model grids currently supported by SPISEA.

.. list-table:: Evolution Models
   :header-rows: 2
   :widths: 30 15 18 30 25

   * - Model Name
     - Mass Range
     - log(Age) Range
     - Metallicity Values
     - Ref
   * -
     - M\ :sub:`⊙`
     - Years
     - [Fe/H]
     -
   * - ``MISTv1``
     - 0.10 – 300
     - 5.01 – 10.30
     - -4.0, -3.5, -3.0, -2.5, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5
     - v1.2; Choi et al. (2016)
   * - ``MergedBaraffePisaEkstromParsec``
     - 0.08 – 120
     - 6.00 – 10.09
     - 0
     - Appendix B; Hosek et al. (2020)
   * - ``MergedPhillipsBaraffePisaEkstromParsec``
     - 0.01 – 120
     - 6.00 – 10.00
     - 0
     - Appendix; Begbie et al. (2026)
   * - ``Parsec``
     - 0.10 – 65
     - 6.60 – 10.12
     - 0
     - Bressan et al. (2012)
   * - ``Baraffe15``
     - 0.07 – 1.4
     - 5.70 – 10.0
     - 0
     - Baraffe et al. (2015)
   * - ``Ekstrom12``
     - 0.80 – 300
     - 6.00 – 8.0
     - 0
     - Ekström et al. (2012)
   * - ``Pisa``
     - 0.20 – 7
     - 6.00 – 8.0
     - 0
     - Tognelli et al. (2011)
   * - ``Phillips2020``
     - 0.0005–0.075
     - 6.00 - 10.00
     - 0
     - Phillips et al. (2020)
   * - ``Marley2021``
     - 0.0005–0.083\ :sup:`a`
     - 6.00 - 10.00
     - -0.5, 0, 0.5
     - Marley et al. (2021)
   * - ``COSMIC`` :sup:`b`
     - 0.08 - 150
     - -
     - -2.3 - 0.18
     - Breivik et al. (2020)

.. rubric:: Footnotes
 
:sup:`a` Actual maximum value given by the Sonora models (Marley et al., 2021) relies on the age of the cluster. For example, for log(Age)=6.0, the mass range is limited to 0.0005 - 0.011 M\ :sub:`⊙`.

:sup:`b` COSMIC evolves the stars externally and does not use SPISEA's standard isochrone-grid architecture. Instead, it uses a custom atmosphere grid that is created on the fly. See Breivik et al. (2020) for more details on COSMIC. When COSMIC is used, the IFMR is ignored.

Please note the stellar mass range, age range, and metallicity values of the evolution
model grid you choose:

* Stars generated from the :ref:`imf_objects` (controlled by the
  `massLimits` variable) that fall outside of the mass range of the
  evolution model will be dropped from the final cluster table, unless
  the IFMR is defined (which will catch all the high mass stars and
  get their compact object values)

* SPISEA will throw an error if you request an age outside of the
  evolution model age range

* SPISEA will use the evolution model with the closest available
  metallicity to the requested value

If you require other evolution models or need to
expand the existing grids, please see
:ref:`add_evo_models`. 

Merged Sonora and MIST grid
-----------------------------------------
``MergedFlameSkimmerRedDiamondbackMIST`` reads precomputed FITS isochrones
from ``SPISEA_MODELS/evolution/merged/flame_red_diamondback_mist``.
It selects the closest stored metallicity and log age; it does not construct
connections while loading. Install the associated model files as well as
the Python code. Filenames and available parameters come from
``grid_manifest.csv`` and, when installed, ``low_metallicity_manifest.csv``.

The reference grid contains 742 isochrones: seven metallicities from
[M/H] = -1 to +0.5 in 0.25-dex steps, at 106 ages from 0.1 Myr through 15 Gyr.
The regular log10(age/yr) spacing is 0.05, from 5.00 through 10.15.
Separate exact 9-Gyr and 15-Gyr files represent the family handoff and grid
endpoint; both are selected through the manifest without special filename rules.

Flame Skimmer supplies ages through 9 Gyr. Above 9 Gyr, Red Diamondback
supplies -0.5 <= [M/H] <= +0.5, with linear interpolation in metallicity at
-0.25 and +0.25. Older -1 and -0.75 sequences use flagged Flame age
extrapolation. Intermediate metallicities interpolate logarithmic stellar
fields linearly in [M/H]. The adopted mass connections use limited cubics in
log10(initial mass), with radius derived from luminosity and effective
temperature. Gravity is interpolated independently and is not forced to
equal GM/R squared. Retained source rows remain unchanged.

The optional extension adds 848 saved isochrones at [M/H] = -4, -3.5, -3,
-2.5, -2, -1.75, -1.5 and -1.25, for 1,590 total. Each combines the -1
Sonora sequence with MIST at the saved metallicity. These are explicitly
labeled proxies, not native low-metallicity Sonora calculations. FITS metadata
and row provenance identify these approximations. Distinct retained MIST
evolutionary points can share an initial mass; the extension preserves these
points and their ``mist_eep`` identifiers rather than choosing a branch for
downstream interpolation by mass. Review the grid's consistency and mass-tie
reports when using these approximations in scientific work.

ModelMismatch Error
-----------------------------------------
Since SPISEA v2.1.4, we began tracking the version of the evolution
model grid (e.g., the evolution models stored in
``<SPISEA_MODELS>/evolution``).
Each evolution model class has a required model grid
version assigned to it. If your evolution model grid
does not match or exceed the minimum version required by
your desired evolution model, a ``ModelMismatch``
exception will be raised.

To resolve the ``ModelMismatch`` error, please
re-download the latest version of the evolution
model grid in the installation instructions (:ref:`models`).
     

Base Evolution Model Class
------------------------------------
.. autoclass:: evolution.StellarEvolution
	       :members:
 
Specific Evolution Model Classes
--------------------------------------
.. autoclass:: evolution.MISTv1
	       :show-inheritance:
		  
.. autoclass:: evolution.MergedBaraffePisaEkstromParsec
	       :show-inheritance:
		  		  
.. autoclass:: evolution.Baraffe15
	       :show-inheritance:
		  
.. autoclass:: evolution.Ekstrom12
	       :show-inheritance:
		  
.. autoclass:: evolution.Parsec
	       :show-inheritance:	       

.. autoclass:: evolution.Pisa
	       :show-inheritance:

.. autoclass:: evolution.MergedPhillipsBaraffePisaEkstromParsec
	       :show-inheritance:

.. autoclass:: evolution.MergedFlameSkimmerRedDiamondbackMIST
               :show-inheritance:
   
.. autoclass:: evolution.Phillips2020
	       :show-inheritance:

.. autoclass:: evolution.COSMIC
	       :show-inheritance:
