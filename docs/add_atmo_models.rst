.. _add_atmo_models:

========================================
Adding Atmosphere Models
========================================
Stellar atmosphere grids are implemented in SPISEA using STSci's `stsynphot <https://stsynphot.readthedocs.io/>`_ / `synphot <https://synphot.readthedocs.io/>`_ stack and CDBS infrastructure.
Adding a new atmosphere grid involves saving the models in a new
sub-directory under your ``PYSYN_CDBS`` directory and creating a new
function in ``atmospheres.py`` to access those files.

We'd recommend looking at the existing atmosphere model directories as
a template for how to make a new one. Some general notes:

* The new atmosphere grid you implement should be placed in the
  ``<SPISEA_MODELS>/cdbs/grid/`` directory
 
* Inside the new directory you need two things: a ``catalog.fits``
  file (which tells SPISEA the Teff, logg, and metallicity of the
  atmospheres in your grid) and then the atmosphere model files
  themselves.
  
* The ``catalog.fits`` file MUST follow the format where one column
  labeled 'INDEX', which gives the '<teff>,<metallicity>,<logg>' for
  each atmosphere, and then one column labeled 'FILENAME', which gives the path to the file containing that atmosphere.

* There is some flexibility as to how the atmosphere model
  files are organized. In ``BTSettle_2015_rebin``, you'll see that each
  file is simply stored in the same directory as the ``catalog.fits``
  file, with one model per file. In ``phoenix_v16_rebin``, the files are
  separated into sub-directory by metallicity, and then there is one
  file for each teff and then multiple columns within that file for
  different logg. It doesn't matter what you choose, as long as it is
  reflected in the 'INDEX' column of catalog.fits.
  
* For the atmospheres, use units of angstroms for the wavelength and
  "FLAM" for flux ([erg/s/cm^2/A]; see `synphot docs
  <https://synphot.readthedocs.io/en/latest/synphot/units.html>`_)

More detailed documentation on this is coming soon. In the meantime, let us know on the  Github `issue tracker
<https://github.com/MovingUniverseLab/spisea/issues>`_ if you'd like to
implement a new atmospheric model grid.
