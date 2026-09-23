.. _imf_objects:

==========
IMF Object
==========
Initial Mass Functions (IMFs) are defined as classes in
``spisea/imf/imf.py``. These can be defined by::

  from spisea.imf import imf
  imf_obj = imf.<class_name>

Currently, the IMF is implemented as a power-law:
dn / dm ~ m^-alpha. The Broken Power-Law IMF object
gives users the flexibility to define the mass range
and exponents of the IMF. 

The IMF object is an input for the :ref:`cluster_objects`, and will
be used to draw the inital stellar mass distribution for the cluster.
A :ref:`multi_obj` may be passed to the IMF object to
form multiple star systems. The default is the SPISEA v2.5
:class:`~imf.multiplicity.MultiplicityUnresolved` /
:class:`~imf.multiplicity.MultiplicityResolvedDK` objects;
Offner et al. 2023 is opt-in (see :ref:`multi_obj`).

Base IMF Class
--------------
.. autoclass:: imf.imf.IMF
	       :members: generate_cluster, calc_multi

Broken Power-Law IMFs
---------------------
.. autoclass:: imf.imf.IMF_broken_powerlaw
	       :show-inheritance:
		  
.. autoclass:: imf.imf.IMFSalpeter1955
	       :show-inheritance:
		  
.. autoclass:: imf.imf.Miller_Scalo_1979
	       :show-inheritance:
   
.. autoclass:: imf.imf.Kennicutt_1983
	       :show-inheritance:
		  
.. autoclass:: imf.imf.Kroupa_2001
	       :show-inheritance:
		  
.. autoclass:: imf.imf.Weidner_Kroupa_2004
	       :show-inheritance:

.. autoclass:: imf.imf.Salpeter_Kirkpatrick_2024
	       :show-inheritance: