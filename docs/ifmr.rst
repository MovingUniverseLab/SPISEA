.. _ifmr_objects:

=============
IFMR Object
=============
The Initial-Final Mass Relation (IFMR) maps a star’s initial zero-age main sequence (ZAMS)
mass to the type and mass of the compact object it will form. The
SPISEA IFMR classes support white dwarfs (WD), neutron stars (NS),
and black holes (BH). The associated code is found in ``spisea/ifmr.py``.

To define an IFMR object::

  from spisea import ifmr
  ifmr_obj = ifmr.<class_name>

Compact objects are included in the output tables produced by
:ref:`cluster_objects`, and can be identified by the "phase" column:
     
* WD: phase = 101
* NS: phase = 102
* BH: phase = 103

See `Quick Start Example
<https://github.com/MovingUniverseLab/spisea/blob/main/docs/Quick_Start_Make_Cluster.ipynb>`_
for more examples of how to interact with the :ref:`cluster_objects`
output.

IFMR Base Classes
--------------------
.. autoclass:: ifmr.IFMR_Spera15
	       :members: generate_death_mass
	      
.. autoclass:: ifmr.IFMR_Raithel18
	       :members: NS_mass, generate_death_mass
			 
.. autoclass:: ifmr.IFMR_N20_Sukhbold
               :members: NS_mass, generate_death_mass
	     
