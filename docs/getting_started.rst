.. _getting_started:

Installation
############

==========================
Install From Git
==========================
If you are downloading the code from scratch, please follow the
instructions below. If you had already downloaded version 1 of the
code and are switching to version 2, please see :ref:`version`. 

SPISEA is hosted on `GitHub <https://github.com/MovingUniverseLab/spisea>`_.
To begin, clone the git repository in your desired code directory::

   git clone https://github.com/MovingUniverseLab/spisea.git

The ``main`` branch contains the current release,
while the ``dev`` branch is for code development.

.. _Dependencies:

Dependencies
----------------
The basic SPISEA installation requires ~20 GB of memory, which is
primarily for the suite of stellar evolution and atmosphere models.
Other dependencies:

* python (>=3.7, < 3.12)
* astropy
* synphot (>= 1.5)
* stsynphot (>= 1.7.0)
* scipy
* numpy (>= 1.17, < 2.0)
* matplotlib
* `pytest-astropy-header <https://github.com/astropy/pytest-astropy-header>`_
* STScI CDBS data package (See :ref:`set-up-cdbs` below for instructions)

.. _synphot-install:

Installing synphot and stsynphot
--------------------------------
SPISEA reads model atmospheres with ``stsynphot.grid_to_spec``.
The minimum versions are ``synphot`` 1.5 and ``stsynphot`` 1.7.0.
``stsynphot`` 1.7.0 is the release that will contain the merged
grid-spectrum cache from `stsynphot pull request 244
<https://github.com/spacetelescope/stsynphot_refactor/pull/244>`_.
That cache keeps each grid spectrum after the first read, so
isochrone photometry does not re-read the same atmosphere files
for every star. The pull request's changelog entry is already
filed under 1.7.0.

``stsynphot`` 1.7.0 is not on PyPI yet. Install both packages with
pip from their GitHub repositories. ``synphot`` 1.8.0 is the current
release and satisfies the 1.5 minimum::

   pip install "synphot @ git+https://github.com/spacetelescope/synphot_refactor.git@1.8.0"

Until pull request 244 is merged, that cache is on Lingfeng Wei's
branch::

   pip install "stsynphot @ git+https://github.com/wei-lingfeng/stsynphot_refactor.git@cache-grid-spectra"

After the pull request is merged, install ``stsynphot`` from the
STScI repository. Tag ``1.7.0`` is the minimum version::

   pip install "stsynphot @ git+https://github.com/spacetelescope/stsynphot_refactor.git@1.7.0"

Once that release is on PyPI, the same minimum is::

   pip install "synphot>=1.5" "stsynphot>=1.7"

.. _Build:

Build
------
Once the repository is cloned, you need to add SPISEA to your
``PYTHONPATH`` environment variable. For example, if you are using
BASH, then you need to add the following line to your ``bash_profile``::
  
   export PYTHONPATH=$PYTHONPATH:<path_to_SPISEA>

where ``<path_to_SPISEA>`` is the location of the top-level SPISEA
directory. 

.. _set-up-cdbs:

Set Up CDBS Directory
---------------------------------
SPISEA uses the STScI CDBS infrastructure to store
model atmospheres and HST filter throughputs.
In particular, it expects a top-level ``cdbs`` directory
with 3 sub-directories: ``comp``, ``mtab``,
and ``grid``. If you already have these set-up then you may
proceed to the next section. Otherwise, follow the steps below.

#. Create a cdbs directory on your machine::

     mkdir <your_path>/cdbs

#. The CDBS files can be downloaded from this STScI `website
   <https://archive.stsci.edu/hlsp/reference-atlases>`_.
   Download the ``synphot1_throughput-master.tar`` and
   ``synphot5_pheonix-models.tar`` files and place them in your
   cdbs directory.

#. Untar each of the files. Once completed, a new directory named ``grp`` will appear.

#. The ``comp``, ``mtab``, and ``grid`` sub-directories SPISEA needs
   are under ``grp/redcat/trds/``. Move these directories directly to
   ``<your_path>/cdbs``.

#. You may remove the (now empty) ``grp`` directory and
   the tar files if desired.

You should end up with a cdbs directory
``<your_path>/cdbs`` with sub-directories
called ``comp``, ``mtab``, and ``grid``.
You will add additional atmosphere models to ``cdbs/grid`` in
the :ref:`models` section below.

Note: To use the Roman filters, the CDBS data version should be from
2024-03-22 or later.


.. _models:

Stellar Evolution and Atmosphere Models
-------------------------------------------------------
SPISEA requires grids of stellar evolution and atmosphere models in
order to function. The evolution models are
stand-alone and specific to SPISEA. The atmosphere model grids use the
STScI CDBS conventions and should be placed in the ``cdbs/grid`` directory.

You will need to download 2 files:

* `spisea_models.tar.gz
  <https://w.astro.berkeley.edu/~jlu/spisea/spisea_models.tar.gz>`_ (4.1 GB; 23 GB unzipped)

* `spisea_cdbs.tar.gz <https://w.astro.berkeley.edu/~jlu/spisea/spisea_cdbs.tar.gz>`_  (155 MB; 354 MB unzipped)

You may **optionally** download a third file, which contains
higher-resolution stellar atmospheres. Note that this file is quite
large, and is not necessary for most SPISEA use cases:

* `spisea_cdbs_highres.tar.gz <https://w.astro.berkeley.edu/~jlu/spisea/spisea_cdbs_highres.tar.gz>`_ (54 GB; 80 GB unzipped)

SPISEA uses the low-resolution atmospheres (R = 250) in
``spisea_cdbs.tar.gz`` by default, as
it is then much faster for synthetic photometry and
is sufficient in most applications. However, the user can change
this default; see  :ref:`atmo_models` for
more details. **Unless you change this default,**
``spisea_cdbs_highres.tar.gz`` **is not required.**

Once downloaded, ``spisea_cdbs.tar.gz`` (and
``spisea_cdbs_highres.tar.gz``, if desired) should be
expanded in  the ``cdbs`` directory. The output directories
will automatically be placed in the correct locations. 
The ``spisea_models.tar.gz`` file can be expanded
anywhere, though we recommend expanding it in the same location as 
your ``cdbs`` directory. 

For example, if your cdbs directory is ``/<your_path>/models/cdbs/``
then you should put the .gz files in ``/<your_path>/models/``
and then unzip them from there::

   cd /<your_path>/models/
   tar -xvf spisea_cdbs.tar.gz
   tar -xvf spisea_models.tar.gz


``spisea_cdbs.tar.gz`` will put the model atmospheres in
``cdbs/grid``, and ``spisea_models.tar.gz`` will put the evolution
models in a new directory called ``evolution/``

.. _setup-paths:

Set Paths to Models
-------------------------------------

You need to notify python where the evolution and atmosphere models
live. This is done by setting two environment variables, ``PYSYN_CDBS`` and
``SPISEA_MODELS``, to point to the ``cdbs`` and ``models``
directories (i.e. the directory where the ``evolution/`` directory
lives), respectively. For example, in your .bash_profile::
  
  export PYSYN_CDBS=/<path_to_cdbs_directory>
  export SPISEA_MODELS=/<path_to_models_directory>


.. _test-setup:

Testing Your SPISEA Setup
---------------------------------------

If all goes well, you should be able to import any of the SPISEA
functions an import statement like those at the top
of the Quick Start Guide, such as::
    
    from spisea import synthetic

You may get warnings that Extinction or thermal files are missing,
such as::

    UserWarning: Extinction files not found in /u/mwhosek/models/cdbs/extinction
    
    UserWarning: No thermal tables found, no thermal calculations can be performed
    
However, these can be safely ignored since SPISEA doesn't use those functionalities.

To further test your SPISEA install, try running the `Quick Start
notebook
<https://github.com/MovingUniverseLab/spisea/blob/main/docs/Quick_Start_Make_Cluster.ipynb>`_.
It is also located in ``SPISEA/docs``.

To test the full range of
evolution models, atmosphere models, and photometric filters,
run the test functions by going into the ``SPISEA/spisea`` directory and running::

    pytest

Note that this uses the python
`pytest
<https://docs.pytest.org/en/7.1.x/>`_
package. This will trigger the test functions we have implemented. If all is
well, you shouldn't get any errors. Warnings are (generally) fine.

TroubleShooting
-----------------------
If SPISEA is not properly in your ``PYTHONPATH``, then when you try
to import the SPISEA functions you will get an error message
reporting ``No module named spisea``.

If the ``SPISEA_MODELS`` or ``PYSYN_CDBS`` paths are broken, then
SPISEA won't know where to get the stellar models.
When trying to import ``spisea/synthetic.py``, you will get
either or both of the following warnings::

    UserWarning: PYSYN_CDBS is undefined; functionality will be SEVERELY crippled.
    
    UserWarning: SPISEA_MODELS is undefined; functionality will be SEVERELY crippled.
      
==========================
Build and deploy from Docker
==========================

Build your own SPISEA image for Docker Containers. This installation form contains SPISEA deployed in a container and includes the data sets, models and all the necessary paths and code.

Requirements
-----------------------

- Linux, Windows or MacOS with Docker installed.
- At least 2 CPUs and 4 GB of RAM and 16 GB of storage.

Installation
-----------------------

To create the container image, clone this repository and build the container::

    git clone https://github.com/MovingUniverseLab/spisea.git
    cd spisea
    docker build -t spisea .
    
Usage
-----------------------
To open a shell ready play with SPISEA::

    docker run -ti spisea bash

To execute a script you have in your current folder:: 

    docker run -ti -v $PWD:$PWD -w $PWD spisea  python myscript.py
    
==========================
Deploy from DockerHub
==========================

If you don't want to build the image from scratch you can use a pre-build container image from `DockerHub <https://hub.docker.com/r/amigahub/spisea>`_ using the following commands::

    docker pull amigahub/spisea:v1

Then, to open a shell ready to play with SPISEA::

    docker run -ti amigahub/spisea:v1 bash

To execute a script you have in your current folder::

    docker run -ti -v $PWD:$PWD -w $PWD spisea  python myscript.py

==========================
Deploy from Singularity containers
==========================

Download the image from DockerHub and convert it into a ``.sif`` image for Singularity.::

    singularity pull spisea.sif docker://amigahub/spisea:v1
    
After downloading the image, you can use it in singularity by opening a shell on SPISEA image::

    singularity shell spisea.sif 





    

