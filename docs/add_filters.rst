.. _add_filters:

========================================
Adding Photometric Filters
========================================
If the user wants to add new photometric filters to SPISEA, there are 4 main steps:

  1) Save the filter transmissions as text files in a new sub-directory in the top-level ``filt_func`` directory.
  2) Define a new function in ``filters.py`` that reads a unique filter string the user assigns to the new filters,
     and then read the appropriate data files from the new ``filt_func`` sub-directory. 
  3) Edit the ``get_filter_info()`` function  in ``synthetic.py`` to
     call the new function in ``filters.py`` when the new filter
     string is called.
  4) If needed (if converting between the ``'_'`` and ``','`` separators is not sufficient), 
     edit the ``get_obs_str()`` and ``get_filter_col_name()`` functions in ``synthetic.py`` to
     convert between column name and filter string for the new filters.
  5) Add the filter to the ``filt_list`` in the ``test_filters()`` function in ``tests/test_models.py`` to
     ensure the filter can be loaded properly.

Additional documentation on this is coming soon. In the meantime, let us know on the  Github `issue tracker
<https://github.com/MovingUniverseLab/spisea/issues>`_ if you'd like to
implement new photometric filters.
