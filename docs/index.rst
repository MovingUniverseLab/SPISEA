.. sampledoc documentation main file, created by
   sphinx-quickstart on Tue Aug 11 05:04:40 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. figure:: ../SPISEA_logo_final.png
	    :width: 380
            :height: 336
	    :align: center

==========================================================================
SPISEA: Stellar Population Interface for Stellar Evolution and Atmospheres
==========================================================================
SPISEA is an open-source python package that generates single-age, single-metallicity
populations (i.e. star clusters). It gives the user control over many parameters:

* Cluster characteristics (age, metallicity, mass, distance)
* Total extinction, differential extinction, and extinction law
* Stellar evolution and atmosphere models
* Stellar multiplicity and Initial Mass Function
* Initial-Final Mass Relation
* Photometric filters

Here is a brief list of things that SPISEA can do:

* make a cluster isochrone in many filters using different stellar models
  
* make a star cluster at any age with an unusual IMF and unresolved multiplicity
  
* make a spectrum of a star cluster in integrated light

SPISEA can be downloaded from `Github <https://github.com/MovingUniverseLab/spisea>`_. Please cite `Hosek et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020arXiv200606691H/abstract>`_  if you use SPISEA in
your research. 

Getting Started
----------------
.. toctree::
   :maxdepth: 1

   getting_started.rst
   quick_start.rst
   more_examples.rst


Documentation
-------------------------
.. toctree::
   :maxdepth: 1

   make_isochrone.rst
   make_cluster.rst
   evo_models.rst
   atmo_models.rst
   imf.rst
   multiplicity.rst
   ifmr.rst
   extinction.rst
   filters.rst

Advanced Documentation
-----------------------------
.. toctree::
   :maxdepth: 1
	      
   add_evo_models.rst
   add_atmo_models.rst
   add_filters.rst

Contributions
---------------
Have an idea? We'd love to have you join our team!
We welcome and encourage contributions to SPISEA.

For feature additions, we ask that users make their
own fork of the repository, make their changes, and then submit a pull
request to the "dev" branch.

All contributions will be acknowledged on the :ref:`contributors` page
(with permission). Contributors with features used in code
releases will be co-authors in future SPISEA software papers.


Change Log
----------
2.5 (2026-07-17)

* *Major Changes*
	* Addition of brown dwarf models (see `Begbie et al. (2026) <https://arxiv.org/abs/2607.14292>`_)
		* New evolution grids: `Phillips2020`, Marley2021, and mergedPhillipsMergedPhillipsBaraffePisaEkstromParsec
		* New atmosphere grids: Phillips2020, Meisner2023
		* get_merged_atmosphere now uses Meisner2023 for T < 1000 K and interpolated BTSettl/Meisner2023 for 1000 < T <= 1200
  		* Updated evolution and atmosphere grid data sets w/ grid_version=3.0
	* Added option to evolve binary systems using `COSMIC <https://cosmic-popsynth.github.io/COSMIC/>`_ via new COSMIC evolution model. 
		* When COSMIC is used when creating synthetic clusters, SPISEA creates binary systems normally (using the Multiplicity obj) which are then evolved using COSMIC. 
		* New COSMIC evolution class, which is used with new IsochronePhotExternalEvolution obj ('external evolution' referring to fact that stellar evolution is done outside of standard SPISEA isochrone-grid framework)
		* New merged atmosphere model class get_merged_atmosphere_w_bb_supplement, which reverts to blackbody atmospheres for stars with parameters outside of existing grid support.
		* New tutorial for creating SPISEA clusters using COSMIC: docs/Cluster_w_COSMIC.ipynb

* *Minor Changes*
	* Added opt-in Offner et al. 2023 (Protostars and Planets VII; arXiv:2203.10066) multiplicity model, including brown dwarfs. MF/CSF is a logistic in log-mass; :math:`\gamma(M)` is an error-weighted logistic; :math:`\mu(a)` is a smooth broken power law (:math:`s=0.1` dex); :math:`\sigma(\log_{10} a)` is a 2-parameter logistic. Comparison figures vs SPISEA v2.5 are in :ref:`multi_obj`. The SPISEA v2.5 :class:`~imf.multiplicity.MultiplicityUnresolved` / :class:`~imf.multiplicity.MultiplicityResolvedDK` objects remain the default.
	* Companion mass and separation draws (including brown-dwarf q and binaries-only BD systems) now live on the multiplicity object rather than being special-cased in ``imf.py``.
	* The MISTv1.2-synthpop model extension was modified to include denser sampling in the gap between the base MISTv1.2 grids and 0.1Msun.
  	* Modified default for MISTv1.2 isochrones: synthpop_extension will be True by default to keep a consistent lower mass limit of 0.1Msun across all ages and metallicities.
	* Added option to return synthetic photometry in terms of AB or ST mag units in IsochronePhot. Vega mag units remains the default. New meta keyword `MAGSYS` added to output tables to specify magnitude unit system.
  	* Added SODC extinction law
  	* Filter handling improvements
		* Fix bug where DECam "Y" filter was mislabeled "y", and updated DECam filters to latest version 
		* Add new filters for Hipparcos, Kepler, OGLE, TESS, Tycho, Washington, Subaru HSC
		* Enable use of all Gaia filters with warning recommending latest (EDR3)
		* All pysynphot filters can now be used
	* Minor bugs and case handling
		* Debug edge case for no companion stars
		* Restore functionality where final mass = initial for companion stars with lower mass than the isochrone's range
		* MIST evolution allows input metallicities within 0.1 dex of the allowed range, since the nearest grid value is adopted, and this eliminates floating value precision issues.

2.4 (2026-03-20)
  * Added backward compatibility for isochrone file names created
    before v2.3
  * When reading isochrone from existing file, only keep user
    requested filters in resulting table.
  * Changed the global random seed to a random generator within each class,
    but still retaining the reproducibility. Test cluster files in
    test_data also updated accordingly for testing purposes
  * Added filter support for IRTF L-band
  * Added conversion function between ST mags and Vega mags

2.3 (2026-02-10)
  * Achieves faster cluster generation (factor of about 2x) by using
    replacing ragged arrays with masked arrays when calculating
    multiplicity properties
  * Added new test functions (and associated test data files) ensuring
    that the primary mass and companion
    mass distibutions remain the same as generated with SPISEA <= v2.2
  * Added support to Euclid VIS filter
    

2.2 (2026-01-16)
	* Compatibility updates for SPISEA to work with `SynthPop <https://synthpop.readthedocs.io/en/latest/>`_. Updates include:
		* Low mass objects below the isochrone grid can optionally be kept in clusters (off by default) and will have ``current_mass=initial_mass`` and ``phase=98``, with no other evolutionary information or photometry.
	* Evolution model versions are now logged in IsochronePhot files and checked if present.
	* The option ``synthpop_extension`` is now available for MISTv1 version=1.2 evolution. This fills in the missing parameter space down to initial mass 0.1Msun in isochrones where needed. Use of this option will require downloading updated isochrone files.
	* Added support for Euclid filters. 
    
2.1.15 (2025-10-25)
  * Updated Roman filter name from outdated w146 to current f146. From
    now on a STScI CDBS data version from 2024-03-22 or later is
    required to use the Roman filter photometry.

2.1.14 (2025-10-23)
  * Added support for NACO intermediate bands
  * Added metallicity bounds check when calling atmosphere model grid
  * Added check to make sure desired filter columns are present in
    isochrone file when deciding whether a new isochrone needs to be
    calculated
  * Small bugfix: ensures that imf_mass_limit is reset properly when
    generating multiple clusters from same imf object

2.1.13 (2024-12-10)
  * Small bugfix to make iso_dir creation more robust
    
2.1.12 (2024-12-02)
  * Added support for Rubin Observatory filters

2.1.11 (2023-09-15)
  * Small bugfix to RedLawFritz11 so paths work regardless of
    operating system used

2.1.10 (2023-06-01)
  * Added support for Roman Space Telescope filters (via stsynphot)
    
    * Note: this requires the ``cdbs/mtab`` and
      ``cdbs/comp/`` directories to be at least the 2021-07-16 version or later. See :ref:`getting_started` for how to download these files

2.1.9 (2023-01-10)
  * Small bugfix: updating pytest setup.cfg and conftest.py files to
    work with latest python/astropy versions
  * Added additional print statements in verbose condition flags


2.1.8 (2023-01-09)
  * Small bugfix: Updating numpy array definitions to use
    dtype='objects' instead of depreciated dtype='np.object'

2.1.7 (2022-10-26)
  * Made several updates to the extinction law objects, including:
    
    * Added Indebetouw05, Schoedel10, Fritz11, and NoguerasLara20 laws
    * Improved docs and added plotting functions for Nishiyama09 and RomanZuniga07 laws
    * Updated K-band scale wavelength for Fitzpatrick07 and NoguerasLara18 to match published values
    * Updated assorted "under the hood" functions to handle new laws
    * Added more test functions for the extinction laws
    * Updated :ref:`ext_law` page so laws are in chronological order

2.1.6 (2022-09-15)
  * Bugfix in how evolution object identifies what age and metallicity
    model to pull from the grid of available models. Previously, it used
    first grid model where grid model value >
    requested value (numpy.searchsorted routine). Now, is takes grid model
    where the difference between the grid model value and requested
    value is minimized. This is most relevant for metallicity
    selection for MIST isochrones, where the model grid metallicities are
    spaced out. 

2.1.5 (2022-07-20)
  * Added additional parameters to get_bb_atmosphere function to give
    user control over the min and max wavelength values, as well as
    wavelength sampling. See :ref:`atmo_models` for details.

2.1.4 (2022-05-24)

  * Expanded MISTv1 evolution model grid age range to logAge = 5.01 --
    10.30
  * Begun version tracking for evolution model grid; updated grid
    (with new MISTv1 age range) is grid 1.1
  * A new ModelMismatch
    exception will be raised if evolution model object requires a
    higher version of the model grid than the user has installed. See
    :ref:`evo_models` for details
  * Fixed MISTv1 sub-directory names so proper Z_solar value
    of 0.0142 is used (only naming convention, does not effect output
    results)
  * Updated evolution model and atmosphere model doc pages with tables
    describing properties of installed options
  * Small tweaks to installation instructions for clarity
  * Updated test functions and removed outdated dependence on nose package
 
2.1.3 (2022-04-23)

  * Updates to documentation and Quick-Start Guide
  * Fully depreciated Hosek18 extinction law to avoid confusion with Hosek18b

2.1.2 (2022-04-01)

  * New installation methods using Docker containers and Singularity


2.1.1 (2022-03-30)

  * Bugfix: fixed application of IFMR to objects in companions table
    (ResolvedCluster and associated subclasses: applies when
    multiplicity + IFMR are used)

  * Minor additions to docs

2.1.0 (2022-01-24)

  * Added ability to create resolved multiple systems, with orbital
    elements assigned to each companion (see :ref:`multi_obj`)

  * Added new metallicity-dependent IFMR models (see
    :ref:`ifmr_objects`)

  * Added blackbody atmosphere class

  * Additional test functions, filter functions, and documentation
    updates

2.0.0 (2020-07-09)

  * Top-level software name change from PyPopStar to SPISEA (see
    :ref:`version` for instructions to do this). This version of code was released in `Hosek et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020arXiv200606691H/abstract>`_. 

1.0.1 (2020-06-24)

  * Bug fix for photometric column headers for some filters, added new tests regarding total cluster mass, small documentation edits

1.0.0 (2019-12-01)

  * Initial release


Logo by Natasha Abrams
