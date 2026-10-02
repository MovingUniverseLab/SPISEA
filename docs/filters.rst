.. _filters:

========================
Photometric Filters
========================

The user can specify what filters are used for synthetic photometry
when defining the :ref:`isochrone_objects`.  Each filter is
identified by a unique string, and an array of such strings
are passed into the Isochrone call.

For example::

    # Use the HST WFC3-IR F127M and F153M filters, along with NIRC2 Kp
    filt_list = ['wfc3,ir,f127m', 'wfc3,ir,f153m', 'nirc2,Kp']
    my_iso = synthetic.IsochronePhot(logAge, AKs, dist, metallicity=0,
                            evo_model=evo_model, atm_func=atm_func,
                            red_law=red_law, filters=filt_list)

These strings follow the format ``<telescope/filter_set>,<filter>``.
Note that there is no space after the comma, and case matters.

Filter transmissions are stored in specific sub-directories
in ``SPISEA/filt_func``. The functions in ``spisea/filters.py`` tell
the :ref:`isochrone_objects` how to interpret the filter strings
and pull the corresponding transmission functions from the filt_func
directories.

Available filters:

* 2MASS
* Bessell
* CTIO_OSIRIS
* DeCam
* Euclid
* GAIA
* HAWK-I
* Hipparcos
* Hubble Space Telescope
* Kepler
* IRTF
* Johnson-Cousins
* Johnson-Glass
* JWST
* Keck NIRC
* Keck NIRC2
* NACO
* OGLE
* PanStarrs 1
* Roman Space Telescope
* Subaru
* TESS
* Tycho
* UKIRT
* Vera C. Rubin Observatory
* VISTA
* Washington
* ZTF


Filter Sets
------------


**2MASS**

`Two-Micron Sky Survey <https://old.ipac.caltech.edu/2mass/>`_

Filters: J, H, Ks

Example: ``'2mass,H'``

**Bessell**

`Bessel (1990) <https://ui.adsabs.harvard.edu/abs/1990PASP..102.1181B/abstract>`_ Johnson-Cousins UBVRI filters

Filters: U, B, V, R, I

Example: ``'bessell,U'``

**CTIO_OSIRIS**

`OSIRIS imager
<http://www.ctio.noao.edu/soar/content/ohio-state-infrared-imagerspectrograph-osiris>`_
on the CTIO telescope.

Filters: H, K

Example: ``'ctio_osiris,H'``


**DeCam**

`Dark Energy Camera <https://noirlab.edu/science/documents/scidoc0472>`_

Filters: u, g, r, i, z, Y

Example: ``'decam,r'``

**Euclid**

Euclid space telescope `NISP filters <https://doi.org/10.5270/esa-kx8w57c>`_
and `VIS single filter <https://svo2.cab.inta-csic.es/svo/theory/fps/index.php?id=Euclid/VIS.vis>`_

Filters: VIS, Y, J, H

Example: ``'euclid,Y'``

**GAIA**

The `GAIA Space Telescope filters <https://www.cosmos.esa.int/web/gaia/iow_20180316>`_.
Note that four sets are available: the pre-launch passbands used in DR1
(`Jordi+10
<https://ui.adsabs.harvard.edu/abs/2010A%26A...523A..48J/abstract>`_),
the passbands used for the DR2 published photometry,
the *revised* DR2 passbands based on the DR2 data (October 2017),
and the `(E)DR3 passbands <https://www.cosmos.esa.int/web/gaia/edr3-passbands>`_.

Filters: G, Gbp, Grp

Versions: dr1, dr2, dr2_rev, edr3

Example (gaia G filter from (E)DR3 passbands):
``'gaia,edr3,G'``

**HAWK-I**

The `High Acuity Wide Field K-band Imager
<https://www.eso.org/sci/facilities/paranal/instruments/hawki.html>`_
located on the ESO VLT (Unit Telescope 4).

Filters: J, H, Ks

Example: ``'hawki,J'``

**Hipparcos**

`Hipparcos Hp filter <http://ulisse.pd.astro.it/Astro/ADPS/Systems/Sys_135/index_135.html>`_

Filters: Hp

Example: ``'hipparcos,Hp'``

**Hubble Space Telescope**

HST filters are defined by their `stsynphot / legacy OBSMODE strings
<https://stsynphot.readthedocs.io/en/latest/stsynphot/obsmode.html>`_. These
are defined in the ``cdbs/mtab/`` and ``cdbs/comp/`` directories.

Example: ``'wfc3,ir,f125w'``


**Johnson-Cousins**

Johnson-Cousin filters (downloaded from
http://www.aip.de/en/research/facilities/stella/instruments/data/johnson-ubvri-filter-curves).

Note: As of July 2026, the link for the source is broken. We note that the I filter here is that of
the Johnson standard filter, which has a long red tail. In the transmission curve here, it is cut off
at 1.1 microns. We recommend careful selection for a standard I filter. Note that the Bessell (1990)
filters are available as listed above. We also have specific filter profiles for many existing
instruments.

Filters: U, B, V, R, I

Example: ``'ubv,B'``

**Johnson-Glass**

Johnson-Glass passbands taken from `Bessell et al. 1988
<https://ui.adsabs.harvard.edu//#abs/1988PASP..100.1134B/abstract>`_

Filters: J, H, K

Example: ``'jg,K'``

**JWST**

JWST NIRCam filters, downloaded from `NIRCam website <https://jwst-docs.stsci.edu/display/JTI/NIRCam+Filters#NIRCamFilters-filt_trans>`_. The filter functions in the nircam_throughputs/modAB_mean/nrc_plus_ote folder is used.

Filters: F070W, F090W,  F115W, F140M, F150W, F150W2, F162M, F164N, F182M, F187N, F200W, F210M, F212N, F250M, F277W, F300M, F322W2, F323N, F335M, F356W, F360M, F405N, F410M, F430M,  F444W, F460M, F466N, F470N, F480M

Example: ``'jwst,F356W'``


**Keck NIRC**

`NIRC1 filters <https://www2.keck.hawaii.edu/inst/nirc/>`_ on the Keck
Telescope

Filters: H, K

Example: ``'nirc1,H'``


**Keck NIRC2**

`NIRC2 filters <https://www2.keck.hawaii.edu/inst/nirc2/filters.html>`_
on the Keck Telescope ()

Filters: J, H, Hcont, K, Kp, Ks, Kcont, Lp, Ms, Brgamma, FeII

Example: ``'nirc2,Ks'``

**Kepler**

`Kepler Kp filter <http://svo2.cab.inta-csic.es/theory/fps/index.php?mode=browse&gname=Kepler&asttype=>`_

Filters: Kp

Example: ``'kepler,Kp'``


**NACO**

`ESO NACO filters <https://www.eso.org/sci/facilities/paranal/instruments/naco/inst/filters.html>`_

Filters: J, H, Ks, IB_2.00, IB_2.03, IB_2.06, IB_2.24, IB_2.27,
IB_2.30, IB_2.33, IB_2.36

Example: ``'naco,H'``

**OGLE**

OGLE-IV filters for V- and I-band from `SVO <https://svo2.cab.inta-csic.es/theory/fps/index.php?mode=browse&gname=LCO&gname2=OGLE-IV&asttype=>`_

The R-wide filter is new in 2026, profile provided by Andrzej Udalski

Filters: V, I, Rw

Example: ``'ogle,Rw'``


**PanStarrs1**

PanStarrs 1 filters from `Tonry et al. 2012 <https://ui.adsabs.harvard.edu/#abs/arXiv:1203.0297>`_

Filters: g, r, i, z, y, w

Example: ``'ps1, g'``


**Roman Space Telescope**

Roman Space Telescope WFI filters are defined by their stsynphot-style
`OBSMODE strings
<https://stsynphot.readthedocs.io/en/latest/stsynphot/obsmode.html>`_.
These are defined in the ``cdbs/mtab/`` and ``cdbs/comp/``
directories. Note that the 2021-07-16 version of these directories
must be downloaded from from the `STScI reference atlases
<https://archive.stsci.edu/hlsp/reference-atlases>`_ for these to
work. See :ref:`getting_started` for how to download these files.

Filters: F062, F087, F106, F129, F158, W146, F184, F213

Example: ``'roman,wfi,f062'``

**Subaru**

Subaru filters from the `SVO Filter Profile Service <http://svo2.cab.inta-csic.es/theory/fps/index.php?mode=browse&gname=Subaru&gname2=HSC&asttype=>`_.

Filters: g, r, i, z, Y, nb387, nb468, nb515, nb527, nb656, nb718, nb816, nb921, nb926, nb973

Instruments: hsc

Example: ``'subaru,hsc,i'``

**TESS**

TESS filter: a single wide, red-optical `bandpass <https://heasarc.gsfc.nasa.gov/docs/tess/telescope_information.html>`_

Filters: tess

Example: ``'tess,tess'``

**Tycho**

`Tycho filters <http://ulisse.pd.astro.it/Astro/ADPS/Systems/Sys_135/index_135.html>`_

Filters: B, V

Example: ``'tycho,B'``

**UKIRT**

`UKIRT Telescope filters <http://www.ukidss.org/technical/photom/photom.html>`_

Filters: Z, Y, J, H, K

Example: ``'ukirt,K'``

**Vera C. Rubin Observatory**

`Rubin Telescope filters <https://github.com/lsst/throughputs/tree/main/baseline>`_

Filters: u, g, r, i, z, y

Example: ``'rubin,g'``


**VISTA**

`VISTA Telescope <http://casu.ast.cam.ac.uk/surveys-projects/vista/technical/filter-set>`_

Filters: Z, Y, J, H, K

Example: ``'vista,Y'``

**Washington**

Washington filter system from `Bessell et al. 2001 <https://ui.adsabs.harvard.edu/abs/2001PASP..113...66B/abstract>`_

Filters: C, M, T1, T2

Example: ``'washington,C'``

**ZTF**

`ZTF Telescope <https://www.ztf.caltech.edu/page/technical>`_

Filters: g, r, i

Example: ``'ztf,g'``

**IRTF**

`IRTF NSFCam <https://irtfweb.ifa.hawaii.edu/~nsfcam2/Filter_Profiles.html>`_

Filters: L

Example: ``'nsfcam,L'``
