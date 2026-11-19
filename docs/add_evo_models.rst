.. _add_evo_models:

========================================
Adding Evolution Models
========================================
Need an evolution model not listed among those pre-packaged
with SPISEA? No problem! With the additional documentation provided
below, you can implement your own evolution models into the SPISEA
framework.

If you have questions or run into problems, please raise an issue on
our Github `issue tracker <https://github.com/MovingUniverseLab/spisea/issues>`_. If you are willing to
have the new models you add be added to the SPISEA package and made
available to the community, please fork or branch off of the
development repository and then submit merge requests to add your
changes to the package. These contributions will be added to the
main branch (with attributes given!) in the next code update.

.. _setup:

Setup: Codes, Directories, Files
--------------------------------
The codes associated with the evolution
models are in ``evolution.py``. Each evolution model grid is implemented
as a sub-class of the ``StellarEvolution`` object: for example, the ``MISTv1``
sub-class corresponds to the MISTv1 model grid, etc. If you are
planning to add to an existing model grid, you will need to edit the corresponding
sub-class. If you want to add a brand new model grid, you will need to
create a new sub-class for that model.

The evolution model grids are packaged as individual isochrone files
in SPISEA. These files are originally downloaded from the website
for the given evolution model, modified to match the required SPISEA
formatting, and then saved in the proper location in the model
directory tree.

The isochrone files are stored in the
``$SPISEA_MODELS/evolution/`` directory. Each evolution model grid is contained
in its own sub-directory; for example,
``$SPISEA_MODELS/evolution/PARSECV1.2s`` contains the PARSEC
isochrones.
Each of these sub-directories has a similar structure::

  ---> version (if applicable)
              ---> iso
	                  ---> metallicity
			              ---> rotation (if applicable)
				                  ---> individual files

The metallicity directories are named ``z<mass_fraction>``, where
<mass_fraction> is the mass fraction of metals. The individual
isochrone files saved as FITS tables with the name convention
``iso_<logAge>.fits``, where <logAge> is the log(age) of the
population. While the ``<mass_fraction>`` can be any
number of digits you want, SPISEA expects the <logAge> to
always have two places after the decimal.

The individual ``iso_<logAge>.fits`` files are FITS tables (readable
by `Astropy Table <https://docs.astropy.org/en/stable/table/>`_.), with one isochrone
(i.e., one age) per file. The column names must match those expected
by SPISEA, which are generic (col1, col2, col3, etc). The mapping
between these generic names and the detailed names actually used by the
code is defined in the ``<evolution sub_class>.isochrone()``
function. Within this function, there is a section renaming these
generic column names to more descriptive names, such as Z, logAge,
Teff, etc. If you add new  ``iso_<logAge>.fits``, you must make sure
this mapping is correct for the new files!


Adding New Metallicities to An Existing Model Grid
--------------------------------------------------
To general steps required to add a new metallicity to an existing
model grid are:

1. Download the raw isochrones from the evolution model website.
   Ideally the user would use same age sampling as the solar
   metallicity model grid (see  :ref:`new_ages`). 
2. Reformat the raw isochrones into SPISEA format. They should be
   saved in the appropriate metallicity sub-directory
   (see :ref:`setup`). 
3. Edit the evolution object in evolution.py so it knows about the new
   isochrones and where they live. The following variables in the ``__
   init__``  function need to be updated:
     * ``self.z_list``: add the new metallicities (in terms of mass fraction)
     * ``self.z_file_map``: edit dictionary to map new metallicities with
       new metallicity sub-directory names.
4. (optional, but highly recommended if planning to merge local edits into
   SPISEA development branch for community use)
   Add a test function to ``spisea/tests/test_synthetic.py`` to make sure everything is working properly.

 
.. _new_ages:

Adding New Ages to An Existing Model Grid
--------------------------------------------------
The isochrone ages for a given evolution model are defined in the
``__init__`` function for the evolution model object under the
``self.age_list`` variable.  The pre-packaged model grids have a typical age sampling of 6.0 <
logAge < 10.0, in steps of 0.01. Some models don't offer models across
this entire age range and so the available age range is truncated (the
``Ekstrom12`` models only go from 6.0 < logAge < 8.0, for example). We
highly recommend that users use this same age range and sampling when adding
new models (e.g., new metallicities, etc).

If the user needs to change the available age range for an existing
model grid, please let us know via the  Github `issue tracker
<https://github.com/MovingUniverseLab/spisea/issues>`_.


Creating a New Model Grid
-------------------------
To create an entirely new model set, the user needs to define a new
``StellarEvolution`` sub-class corresponding to that object and set up
the appropriate directory structure in the
``$SPISEA_MODELS/evolution`` directory. Detailed documentation on
this is coming soon. In the meantime, let us know on the  Github `issue tracker
<https://github.com/MovingUniverseLab/spisea/issues>`_ if you'd like to
implement a new model set. 
