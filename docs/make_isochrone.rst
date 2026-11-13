.. _isochrone_objects:

================
Isochrone Object
================
The isochrone classes are defined in ``spisea/synthetic.py``. The primary
inputs to a isochrone object are the stellar population age, distance,
total extinction, and metallicity, along with the :ref:`atmo_models`, 
:ref:`evo_models`, and :ref:`ext_law`. 

If the IsochronePhot sub-class is used then synthetic photometry
will be produced. The :ref:`filters` are defined as additional
inputs. The output photometry is in Vega mags by default (and is always
saved to the iso file in Vega mag), but the user
can opt to return the IsochronePhot object in AB or ST mags. Either way,
the magnitude system is indicated in the MAGSYS isochrone table metadata.

An example of making an IsochronePhot object::

  from spisea import synthetic, evolution
  from spisea import atmospheres, reddening
  import numpy as np

  # Define isochrone input parameters
  logAge = np.log10(5*10**6.) # Age in log(years)
  dist = 4000 # distance in parsec
  AKs = 0.8 # extinction in Ks-band mags
  metallicity = 0 # Metallicity in [M/H]
  evo_model = evolution.MISTv1() 
  atm_func = atmospheres.get_merged_atmosphere
  red_law = reddening.RedLawHosek18b()

  # Specify filters for synthetic photometry. Here we 
  # use  the HST WFC3-IR F127M, F139M, and F153M filters
  filt_list = ['wfc3,ir,f127m', 'wfc3,ir,f139m', 
  'wfc3,ir,f153m']

  # Specify the directory we want the output isochrone
  # table saved in
  iso_dir = './isochrones/'

  # Make IsochronePhot object 
  my_iso = synthetic.IsochronePhot(logAge, AKs, dist, 
                            metallicity=0,
                            evo_model=evo_model, 
                            atm_func=atm_func,
                            red_law=red_law, 
                            filters=filt_list,
			    iso_dir=iso_dir)

See `Quick Start Example
<https://github.com/MovingUniverseLab/spisea/blob/main/docs/Quick_Start_Make_Cluster.ipynb>`_
for a detailed example showing how to interact with the isochrone
object output. Here is a table from Hosek et al. 2020 that describes
the values in the isochrone output table:

.. figure:: images/isochrone_table.png
	    :width: 468
            :height: 400
	    :align: center



		    

Tips and Tricks: The IsochronePhot Object
-----------------------------------------

* It takes ~1-3 mins to make an IsochronePhot object for the first
  time. A FITS table is created with the stellar parameters and
  photometry for each star. This table is saved in the specified
  iso_dir, under the filename iso_<age>_<aks>_<dist>_<z>.fits.
  By default, iso_dir is set to the current working directory unless
  otherwise defined. 

  In future calls of IsochronePhot, the code will first check iso_dir
  to see if the appropriate FITS table already exists. If so,
  then it will simply read the table and be done. This saves
  significant amounts of computation time.

  * **WARNING**: When IsochronePhot checks to see if the desired
    isochrone table already exists, it checks all isochrone properties
    (evolution models, atmosphere models, and reddening law are encoded in the table meta-data).
    If any of these parameters do not match, then the isochrone will
    be re-calculated. 

    However, to keep the isochrone filenames reasonable, only the
    age, extinction, distance, and metallicity are encoded in the
    filename itself. So, if the evolution model, atmosphere model, or
    reddening law have changed, the original file will be overwritten
    by the new isochrone.

    *To avoid files from being unintentionally overwritten, we recommend
    that users specify different iso_dir paths when making isochrones
    with different evolution models, atmosphere models, or reddening
    laws.*

* For external evolution models (i.e. COSMIC), you should use 
  IsochronePhotExternalEvolution
  instead of IsochronePhot. This is because those evolution models do not have isochrones but
  instead evolve the stars externally. The first time you run a new AKs, metallicity, or distance,
  this will take ~10-20 mins because it is creating a new atmosphere grid. This table is saved in the 
  specified iso_dir, under the filename atm_<aks>_<dist>_<z>.fits.

Base Isochrone Class
----------------------------
.. autoclass:: synthetic.Isochrone
	       :members: plot_HR_diagram, plot_mass_luminosity



Isochrone Sub-classes
-----------------------

.. autoclass:: synthetic.IsochronePhot
	       :show-inheritance:
		:members: make_photometry, plot_CMD, plot_mass_magnitude

.. autoclass:: synthetic.IsochronePhotExternalEvolution
	       :show-inheritance:
		:members: make_photometry, plot_CMD, plot_mass_magnitude


Photometry Conversion Functions
-----------------------------
.. _phot_conversions:

.. autofunction:: synthetic.calc_ab_vega_filter_conversion

.. autofunction:: synthetic.calc_st_vega_filter_conversion

