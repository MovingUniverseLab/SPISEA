.. _quick_start:

===================
Quick Start Guide
===================

The tutorials are meant to be read in order:

* :doc:`Making a cluster <notebooks/01_make_cluster>`
* :doc:`Brown dwarfs <notebooks/02_brown_dwarfs>`
* :doc:`Resolved multiple systems <notebooks/03_resolved_multiples>`
* :doc:`A cluster evolved with COSMIC <notebooks/04_cluster_with_cosmic>`

:doc:`Making a cluster <notebooks/01_make_cluster>` demonstrates the steps
required to simulate a cluster and interact with the resulting output.
This guide shows you how to:

Make a Cluster Isochrone
----------------------------------
You will construct a cluster isochrone, setting the age, distance,
extinction, and metallicity, as well as the stellar evolution and
atmosphere models and reddening law. Photometry will be calculated in
several HST filters, and you will construct a color-magnitude diagram,
identifying the position of a 1 solar-mass star.

.. figure:: images/cmd_example.png
	    :align: center


Make a Simulated Cluster with Multiplicity
-----------------------------------------------------------
You will use the isochrone you generated to make a simulated
cluster. Additional parameters you will set include the stellar
multiplicity, IMF, and overall cluster mass. You will make a
color-magnitude diagram of the cluster, and compare it to the base
isochrone. This shows the impact of stellar multiplicity on the
photometric properties of the stellar systems.

.. figure:: images/simulated_cluster_example.png
	    :align: center


Examine Compact Remnants Produced by Different Clusters
----------------------------------------------------------------
You will generate two massive clusters, one at 100 Myr and the other at 10
Gyr, and compare the distribution of compact remnants produced given
the default IFMR. You will also test how
changing the IMF of these clusters changes these distributions.

.. figure:: images/compact_remnants_example.png
	    :align: center


Tutorials and further examples
------------------------------
The pages below expand on this quick start. They appear in the left
sidebar under Quick Start Guide.

.. toctree::
   :maxdepth: 1
   :titlesonly:
   :caption: Tutorials

   notebooks/01_make_cluster
   notebooks/02_brown_dwarfs
   notebooks/03_resolved_multiples
   notebooks/04_cluster_with_cosmic
   more_examples

.. toctree::
   :maxdepth: 1
   :titlesonly:
   :caption: Paper examples

   notebooks/papers/hosek2020/fig02_hrd_extinction
   notebooks/papers/hosek2020/fig03_cmds_metallicity
   notebooks/papers/hosek2020/fig04_multiplicity_extinction
   notebooks/papers/hosek2020/fig05_imfs
   notebooks/papers/hosek2020/fig06_wolf_rayet
   notebooks/papers/hosek2020/fig07_unresolved_spectra
   notebooks/papers/begbie2026/fig02_bd_imfs
   notebooks/papers/begbie2026/fig03_atmosphere_interpolation
   notebooks/papers/begbie2026/fig04_evolution_interpolation
   notebooks/papers/begbie2026/fig05_hrd_ages
   notebooks/papers/begbie2026/fig06_pleiades_gaia
   notebooks/papers/begbie2026/fig07_pleiades_ukidss
   notebooks/papers/begbie2026/fig08_upper_sco
   notebooks/papers/begbie2026/fig09_m44
   notebooks/papers/begbie2026/fig10_mass_fraction
   notebooks/papers/begbie2026/fig11_companion_counts
   notebooks/papers/begbie2026/fig12_bd_separations
