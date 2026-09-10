import logging
import numpy as np
import os
import glob
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, Column
import time
import pdb
import warnings

import stsynphot as stsyn
from stsynphot import exceptions as stsyn_exceptions
from synphot.models import BlackBody1D, Empirical1D
from synphot.spectrum import SourceSpectrum
from synphot.units import convert_flux
from synphot import units as su

from spisea.utils.synphot_bridge import rebin_spec

log = logging.getLogger('atmospheres')


def get_atmosphere_bounds(model_dir, metallicity=0, temperature=20000, gravity=4, verbose=False):
    """
    Given atmosphere model, get temperature and gravity bounds.

    Parameters
    ----------
    model_dir : str
        The name of the atmosphere model grid.
    metallicity : float
        The metallicity of the atmosphere, in terms of [Fe/H]. Solar = 0.
    temperature : float
        The temperature of the atmosphere, in units of K.
    gravity : float
        The gravity of the atmosphere, in units of log(cgs) (e.g. gravity~5 for main-sequence stars).
    verbose : bool (optional)
        Whether to print verbose output.

    Returns
    -------
    temperature_new : float
        The closest temperature to the input temperature.
    gravity_new : float
        The closest gravity to the input gravity.
    metallicity_new : float
        The closest metallicity to the input metallicity. 
    """
    teff_arr, z_arr, logg_arr = get_atmosphere_grid(model_dir)

    # Filter by metallicity. Will chose the closest metallicity to desired input
    metal_list = np.unique(np.array(z_arr))
    metal_idx = np.argmin(np.abs(metal_list - metallicity))
    metallicity_new = metal_list[metal_idx]

    z_filt = np.where(z_arr == metal_list[metal_idx])
    teff_arr = teff_arr[z_filt]
    logg_arr = logg_arr[z_filt]

    # # Now find the closest atmosphere in parameter space to
    # # the one we want. We'll find the match with the lowest
    # # fractional difference
    # teff_diff = (teff_arr - temperature) / temperature
    # logg_diff = (logg_arr - gravity) / gravity
    #
    # diff_tot = abs(teff_diff) + abs(logg_diff)
    # idx_f = np.argmin(diff_tot)
    #
    # temperature_new = teff_arr[idx_f]
    # gravity_new = logg_arr[idx_f]

    # First check if temperature within bounds
    temperature_new = temperature
    if temperature > np.max(teff_arr):
        temperature_new = np.max(teff_arr)
    if temperature < np.min(teff_arr):
        temperature_new = np.min(teff_arr)

    # If temperature within bounds, then check if metallicity within bounds
    teff_diff = np.abs(teff_arr - temperature)
    sorted_min_diffs = np.unique(teff_diff)

    ## Find two closest temperatures
    teff_close_1 = teff_arr[np.where(teff_diff == sorted_min_diffs[0])[0][0]]
    teff_close_2 = teff_arr[np.where(teff_diff == sorted_min_diffs[1])[0][0]]

    logg_arr_1 = logg_arr[np.where(teff_arr == teff_close_1)]
    logg_arr_2 = logg_arr[np.where(teff_arr == teff_close_2)]

    ## Switch to most conservative bound of logg out of two closest temps
    gravity_new = gravity
    if gravity > np.min([np.max(logg_arr_1), np.max(logg_arr_2)]):
        gravity_new = np.min([np.max(logg_arr_1), np.max(logg_arr_2)])
    if gravity < np.max([np.min(logg_arr_1), np.min(logg_arr_2)]):
        gravity_new = np.max([np.min(logg_arr_1), np.min(logg_arr_2)])

    if verbose:
        # Print out changes, if any
        if temperature_new != temperature:
            teff_msg = 'Changing to T={0:6.0f} for met={1:4.2f} T={2:6.0f} logg={3:4.2f}'
            print( teff_msg.format(temperature_new, metallicity, temperature, gravity))

        if gravity_new != gravity:
            logg_msg = 'Changing to logg={0:4.2f} for met={1:4.2f} T={2:6.0f} logg={3:4.2f}'
            print( logg_msg.format(gravity_new, metallicity, temperature, gravity))

    if metallicity_new != metallicity:
        logg_msg = 'Changing to met={0:4.2f} for met={1:4.2f} T={2:6.0f} logg={3:4.2f}'
        print( logg_msg.format(metallicity_new, metallicity, temperature, gravity))

    return (temperature_new, gravity_new, metallicity_new)

def get_atmosphere_grid(model_dir):
    """
    Gets grid of temps, Zs, and loggs of atmosphere grid

    Parameters
    ----------
    model_dir : str
        Model directory

    Returns
    -------
    teff_arr : array-like
        Temperature array

    z_arr : array-like
        Metallicity array
 
    logg_arr : array-like
        Surface gravity array
    """
    # Open catalog fits file and break out row indices
    catalog = Table.read('{0}/grid/{1}/catalog.fits'.format(os.environ['PYSYN_CDBS'], model_dir))
    
    teff_arr = []
    z_arr = []
    logg_arr = []
    for cur_row_index in range(len(catalog)):
        index = catalog['INDEX'][cur_row_index]
        tmp = index.split(',')
        teff_arr.append(float(tmp[0]))
        z_arr.append(float(tmp[1]))
        logg_arr.append(float(tmp[2]))
    teff_arr = np.array(teff_arr)
    z_arr = np.array(z_arr)
    logg_arr = np.array(logg_arr)

    return teff_arr, z_arr, logg_arr

def get_kurucz_atmosphere(metallicity=0, temperature=20000, gravity=4, rebin=False):
    """
    Return atmosphere from the Kurucz pysnphot grid
    (`Kurucz 1993 <http://www.stsci.edu/hst/observatory/crds/k93models.html>`_).

    Grid Range:

    * Teff: 3000 - 50000 K
    * gravity: 0 - 5 cgs
    * metallicity: -5.0 - 1.0

    Parameters
    ----------
    metallicity: float
        The stellar metallicity, in terms of [Fe/H]. Solar = 0.

    temperature: float
        The stellar temperature, in units of K

    gravity: float
        The stellar gravity, in cgs units

    rebin: boolean
        Always false for this particular function
    """
        
    try:
        sp = stsyn.grid_to_spec('k93models', temperature, metallicity, gravity)
    except:
        # Check atmosphere catalog bounds
        (temperature, gravity, metallicity) = get_atmosphere_bounds('k93models',
                                                   metallicity=metallicity,
                                                   temperature=temperature,
                                                   gravity=gravity)

        sp = stsyn.grid_to_spec('k93models', temperature, metallicity, gravity)

    # Do some error checking for 0 fluxes.
    idx = np.where(sp(sp.waveset) != 0)[0]
    if len(idx) == 0:
        print( 'Could not find Kurucz 1993 atmosphere model for')
        print( '  temperature = %d' % temperature)
        print( '  metallicity = %.1f' % metallicity)
        print( '  log gravity = %.1f' % gravity)

    return sp

def get_kurucz_atmosphere_grid():
    """
    Return atmosphere grid from the Kurucz pysnphot grid 
    (`Kurucz 1993 <http://www.stsci.edu/hst/observatory/crds/k93models.html>`_).

    Grid Range:

    * Teff: 3000 - 50000 K
    * gravity: 0 - 5 cgs
    * metallicity: -5.0 - 1.0
    """
    teff_arr, z_arr, logg_arr = get_atmosphere_grid('k93models')
    return teff_arr, z_arr, logg_arr

def get_castelli_atmosphere(metallicity=0, temperature=20000, gravity=4, rebin=False):
    """
    Return atmospheres from the synphot ATLAS9 atlas
    (`Castelli & Kurucz 2004 <http://www.stsci.edu/hst/observatory/crds/castelli_kurucz_atlas.html>`_).

    Grid Range:

    * Teff: 3500 - 50000 K
    * gravity: 0 - 5.0 cgs
    * [M/H]: -2.5 - 0.2

    Parameters
    ----------
    metallicity: float
        The stellar metallicity, in terms of [Z]

    temperature: float
        The stellar temperature, in units of K

    gravity: float
        The stellar gravity, in cgs units

    rebin: boolean
        If true, rebins the atmospheres so that they are the same
        resolution as the Castelli+04 atmospheres. Default is False,
        which is often sufficient synthetic photometry in most cases.

    verbose: boolean
        True for verbose output
    """
    try:
        sp = stsyn.grid_to_spec('ck04models', temperature, metallicity, gravity)
    except:
        # Check atmosphere catalog bounds
        (temperature, gravity, metallicity) = get_atmosphere_bounds('ck04models',
                                                   metallicity=metallicity,
                                                   temperature=temperature,
                                                   gravity=gravity)

        sp = stsyn.grid_to_spec('ck04models', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp(sp.waveset) != 0)[0]
    if len(idx) == 0:
        print( 'Could not find Castelli and Kurucz 2004 atmosphere model for')
        print( '  temperature = %d' % temperature)
        print( '  metallicity = %.1f' % metallicity)
        print( '  log gravity = %.1f' % gravity)

    return sp

def get_castelli_atmosphere_grid():
    """
    Return atmosphere grid from the synphot ATLAS9 atlas 
    (`Castelli & Kurucz 2004 <http://www.stsci.edu/hst/observatory/crds/castelli_kurucz_atlas.html>`_).

    Grid Range: 

    * Teff: 3500 - 50000 K
    * gravity: 0 - 5.0 cgs
    * [M/H]: -2.5 - 0.2
    """
    
    teff_arr, z_arr, logg_arr = get_atmosphere_grid('ck04models')
    return teff_arr, z_arr, logg_arr

def get_nextgen_atmosphere(metallicity=0, temperature=5000, gravity=4, rebin=False):
    """
    metallicity = [M/H] (def = 0)
    temperature = Kelvin (def = 5000)
    gravity = log gravity (def = 4.0)
    """
    try:
        sp = stsyn.grid_to_spec('nextgen', temperature, metallicity, gravity)
    except:
        # Check atmosphere catalog bounds
        (temperature, gravity, metallicity) = get_atmosphere_bounds('nextgen',
                                                   metallicity=metallicity,
                                                   temperature=temperature,
                                                   gravity=gravity)

        sp = stsyn.grid_to_spec('nextgen', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp(sp.waveset) != 0)[0]
    if len(idx) == 0:
        print( 'Could not find NextGen atmosphere model for')
        print( '  temperature = %d' % temperature)
        print( '  metallicity = %.1f' % metallicity)
        print( '  log gravity = %.1f' % gravity)

    return sp

def get_nextgen_atmosphere_grid():
    """
    metallicity = [M/H] (def = 0)
    temperature = Kelvin (def = 5000)
    gravity = log gravity (def = 4.0)
    """
    teff_arr, z_arr, logg_arr = get_atmosphere_grid('nextgen')
    return teff_arr, z_arr, logg_arr

def get_amesdusty_atmosphere(metallicity=0, temperature=5000, gravity=4, rebin=False):
    """
    metallicity = [M/H] (def = 0)
    temperature = Kelvin (def = 5000)
    gravity = log gravity (def = 4.0)
    """
    sp = stsyn.grid_to_spec('AMESdusty', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp(sp.waveset) != 0)[0]
    if len(idx) == 0:
        print( 'Could not find AMESdusty Allard+ 2000 atmosphere model for')
        print( '  temperature = %d' % temperature)
        print( '  metallicity = %.1f' % metallicity)
        print( '  log gravity = %.1f' % gravity)

    return sp

def get_amesdusty_atmosphere_grid():
    """
    metallicity = [M/H] (def = 0)
    temperature = Kelvin (def = 5000)
    gravity = log gravity (def = 4.0)
    """
    teff_arr, z_arr, logg_arr = get_atmosphere_grid('AMESdusty')
    return teff_arr, z_arr, logg_arr

def get_phoenix_atmosphere(metallicity=0, temperature=5000, gravity=4,
                               rebin=False):
    """
    Return atmosphere from the synphot
    `PHOENIX atlas <http://www.stsci.edu/hst/observatory/crds/SIfileInfo/synphottables/index_phoenix_models_html>`_.

    Parameters
    ----------
    metallicity: float
        The stellar metallicity, in terms of [Z]

    temperature: float
        The stellar temperature, in units of K

    gravity: float
        The stellar gravity, in cgs units

    rebin: boolean
        If true, rebins the atmospheres so that they are the same
        resolution as the Castelli+04 atmospheres. Default is False,
        which is often sufficient synthetic photometry in most cases.
    """
        
    try:
        sp = stsyn.grid_to_spec('phoenix', temperature, metallicity, gravity)
    except:
        # Check atmosphere catalog bounds
        (temperature, gravity, metallicity) = get_atmosphere_bounds('phoenix',
                                                   metallicity=metallicity,
                                                   temperature=temperature,
                                                   gravity=gravity)

        sp = stsyn.grid_to_spec('phoenix', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp(sp.waveset) != 0)[0]
    if len(idx) == 0:
        print( 'Could not find PHOENIX BT-Settl (Allard+ 2011 atmosphere model for')
        print( '  temperature = %d' % temperature)
        print( '  metallicity = %.1f' % metallicity)
        print( '  log gravity = %.1f' % gravity)

    return sp

def get_phoenix_atmosphere_grid():
    """
    Return atmosphere grid from the synphot 
    `PHOENIX atlas <http://www.stsci.edu/hst/observatory/crds/SIfileInfo/pysynphottables/index_phoenix_models_html>`_.

    """
    teff_arr, z_arr, logg_arr = get_atmosphere_grid('phoenix')
    return teff_arr, z_arr, logg_arr

def get_cmfgenRot_atmosphere(metallicity=0, temperature=24000, gravity=4.3, rebin=True, verbose=False):
    """
    metallicity = [M/H] (def = 0)
    temperature = Kelvin (def = 24000)
    gravity = log gravity (def = 4.3)

    rebin=True: pull from atmospheres at ck04model resolution.
    """
        
    # Take care of atmospheres outside the catalog boundaries
    logg_msg = 'Changing to logg={0:3.1f} for T={1:6.0f} logg={2:4.2f}'
    if gravity > 4.3:
        if verbose:
            print( logg_msg.format(4.3, temperature, gravity))
        gravity = 4.3

    if rebin:
        sp = stsyn.grid_to_spec('cmfgen_rot_rebin', temperature, metallicity, gravity)
    else:
        sp = stsyn.grid_to_spec('cmfgen_rot', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp(sp.waveset) != 0)[0]
    if len(idx) == 0:
        print( 'Could not find CMFGEN rotating atmosphere model (Fierro+15) for')
        print( '  temperature = %d' % temperature)
        print( '  metallicity = %.1f' % metallicity)
        print( '  log gravity = %.1f' % gravity)

    return sp

def get_cmfgenRot_atmosphere_grid(rebin=True):
    """
    metallicity = [M/H] (def = 0)
    temperature = Kelvin (def = 24000)
    gravity = log gravity (def = 4.3)

    rebin=True: pull from atmospheres at ck04model resolution.
    """
    if rebin:
        teff_arr, z_arr, logg_arr = get_atmosphere_grid('cmfgen_rot_rebin')
    else:
        teff_arr, z_arr, logg_arr = get_atmosphere_grid('cmfgen_rot')
    return teff_arr, z_arr, logg_arr


def get_cmfgenRot_atmosphere_closest(metallicity=0, temperature=24000, gravity=4.3, rebin=True,
                                         verbose=False):
    """
    For a given stellar atmosphere, get extract the closest possible match in
    Teff/logg space. Note that this is different from the normal routine
    which interpolates along the input grid to get final spectrum. We can't
    do this here because the Fierro+15 atmosphere grid is so sparse

    rebin=True: pull from atmospheres at ck04model resolution.

    If verbose, print out the parameters of the match
    """
    # Set up the proper root directory
    if rebin == True:
        root_dir = os.environ['PYSYN_CDBS'] + '/cmfgen_rot_rebin/'
    else:
        root_dir = os.environ['PYSYN_CDBS'] + '/cmfgen_rot/'

    # Read in catalog, extract atmosphere info
    cat = Table.read('{0}/catalog.fits'.format(root_dir), format='fits')
    teff_arr = []
    z_arr = []
    logg_arr = []
    for ii in range(len(cat)):
        index = cat['INDEX'][ii]
        tmp = index.split(',')
        teff_arr.append(float(tmp[0]))
        z_arr.append(float(tmp[1]))
        logg_arr.append(float(tmp[2]))
    teff_arr = np.array(teff_arr)
    z_arr = np.array(z_arr)
    logg_arr = np.array(logg_arr)

    # Now find the closest atmosphere in parameter space to
    # the one we want. We'll find the match with the lowest
    # fractional difference
    teff_diff = (teff_arr - temperature) / temperature
    logg_diff = (logg_arr - gravity) / gravity

    diff_tot = abs(teff_diff) + abs(logg_diff)
    idx_f = np.where(diff_tot == min(diff_tot))[0][0]

    # Extract the filename of the best-match model and read as
    # synphot object
    infile = cat[idx_f]['FILENAME'].split('.')
    spec = Table.read('{0}/{1}.fits'.format(root_dir, infile[0]))

    # Now, the CMFGEN atmospheres assume a distance of 1 kpc, while the the
    # ATLAS models are in FLAM at the surface. So, we need to multiply the
    # CMFGEN atmospheres by (1000/R)**2. in order to convert to FLAM on surface.
    # We'll calculate radius from Teff and logL, which is given in the Table_*.txt file
    t = Table.read('{0}/Table_rot.txt'.format(root_dir), format='ascii')
    tmp = np.where(t['col1'] == infile[0])

    lum = t['col3'][tmp] * (3.839*10**33) # cgs
    sigma = 5.6704 * 10**-5 # cgs
    teff = teff_arr[idx_f] # cgs

    radius = np.sqrt( lum / (4.0 * np.pi * teff**4. * sigma) ) # in cm
    radius /= 3.08*10**18 # in pc

    # Make the synphot spectrum
    w = spec['Wavelength']
    f = spec['Flux'] * (1000 / radius)**2.
    sp = SourceSpectrum(
        Empirical1D,
        points=np.asarray(w, dtype=float) * u.AA,
        lookup_table=np.asarray(f, dtype=float) * su.FLAM,
    )

    # Print out parameters of match, if desired
    if verbose:
        print('Teff match: Input: {0}, Output: {1}'.format(temperature, teff_arr[idx_f]))
        print('logg match: Input: {0}, Output: {1}'.format(gravity, logg_arr[idx_f]))

    return sp


def get_cmfgenNoRot_atmosphere(metallicity=0, temperature=22500, gravity=3.98, rebin=True):
    """
    metallicity = [M/H] (def = 0)
    temperature = Kelvin (def = 24000)
    gravity = log gravity (def = 4.3)

    rebin=True: pull from atmospheres at ck04model resolution.
    """
        
    if rebin:
        sp = stsyn.grid_to_spec('cmfgen_norot_rebin', temperature, metallicity, gravity)
    else:
        sp = stsyn.grid_to_spec('cmfgen_norot', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp(sp.waveset) != 0)[0]
    if len(idx) == 0:
        print( 'Could not find CMFGEN rotating atmosphere model (Fierro+15) for')
        print( '  temperature = %d' % temperature)
        print( '  metallicity = %.1f' % metallicity)
        print( '  log gravity = %.1f' % gravity)

    return sp

def get_cmfgenNoRot_atmosphere_grid(rebin=True):
    """
    metallicity = [M/H] (def = 0)
    temperature = Kelvin (def = 24000)
    gravity = log gravity (def = 4.3)

    rebin=True: pull from atmospheres at ck04model resolution.
    """
    if rebin:
        teff_arr, z_arr, logg_arr = get_atmosphere_grid('cmfgen_norot_rebin')
    else:
        teff_arr, z_arr, logg_arr = get_atmosphere_grid('cmfgen_norot')
    return teff_arr, z_arr, logg_arr

def get_cmfgenNoRot_atmosphere(metallicity=0, temperature=30000, gravity=4.14):
    """
    metallicity = [M/H] (def = 0)
    temperature = Kelvin (def = 30000)
    gravity = log gravity (def = 4.14)
    """
    sp = stsyn.grid_to_spec('cmfgenF15_noRot', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp(sp.waveset) != 0)[0]
    if len(idx) == 0:
        print( 'Could not find CMFGEN non-rotating atmosphere model (Fierro+15) for')
        print( '  temperature = %d' % temperature)
        print( '  metallicity = %.1f' % metallicity)
        print( '  log gravity = %.1f' % gravity)

    return sp

def get_cmfgenNoRot_atmosphere_grid():
    """
    metallicity = [M/H] (def = 0)
    temperature = Kelvin (def = 30000)
    gravity = log gravity (def = 4.14)
    """
    teff_arr, z_arr, logg_arr = get_atmosphere_grid('cmfgenF15_noRot')
    return teff_arr, z_arr, logg_arr

def get_phoenixv16_atmosphere(metallicity=0, temperature=4000, gravity=4, rebin=True):
    """
    Return PHOENIX v16 atmospheres from
    `Husser et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013A%26A...553A...6H/abstract>`_.

    Models originally downloaded via `ftp <http://phoenix.astro.physik.uni-goettingen.de/?page_id=15>`_.
    Solar metallicity and [alpha/Fe] is used.

    Grid Range:

    * Teff: 2300 - 7000 K, steps of 100 K; 7000 - 12000 in steps of 200 K
    * gravity: 0.0 - 6.0 cgs, steps of 0.5
    * [M/H]: -4.0 - 1.0

    Parameters
    ----------
    metallicity: float
        The stellar metallicity, in terms of [Z]

    temperature: float
        The stellar temperature, in units of K

    gravity: float
        The stellar gravity, in cgs units

    rebin: boolean
        If true, rebins the atmospheres so that they are the same
        resolution as the Castelli+04 atmospheres. Default is False,
        which is often sufficient synthetic photometry in most cases.
    """
        
    atm_model_name = 'phoenix_v16'
    if rebin == True:
        atm_model_name = 'phoenix_v16_rebin'


    # Extract atmosphere. If that fails, then check bounds and try again
    try:
        sp = stsyn.grid_to_spec(atm_model_name, temperature, metallicity, gravity)
    except:
        # Check atmosphere catalog bounds
        (temperature, gravity, metallicity) = get_atmosphere_bounds(atm_model_name,
                                                   metallicity=metallicity,
                                                   temperature=temperature,
                                                   gravity=gravity)

        sp = stsyn.grid_to_spec(atm_model_name, temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp(sp.waveset) != 0)[0]
    if len(idx) == 0:
        print( 'Could not find PHOENIXv16 (Husser+13) atmosphere model for')
        print( '  temperature = %d' % temperature)
        print( '  metallicity = %.1f' % metallicity)
        print( '  log gravity = %.1f' % gravity)

    return sp

def get_phoenixv16_atmosphere_grid(rebin=True):
    """
    Return PHOENIX v16 atmosphere grid from  
    `Husser et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013A%26A...553A...6H/abstract>`_. 
    
    Models originally downloaded via `ftp <http://phoenix.astro.physik.uni-goettingen.de/?page_id=15>`_.
    Solar metallicity and [alpha/Fe] is used.

    Grid Range:

    * Teff: 2300 - 7000 K, steps of 100 K; 7000 - 12000 in steps of 200 K
    * gravity: 0.0 - 6.0 cgs, steps of 0.5
    * [M/H]: -4.0 - 1.0
        
    rebin: boolean
        If true, rebins the atmospheres so that they are the same
        resolution as the Castelli+04 atmospheres. Default is True
    """
        
    atm_model_name = 'phoenix_v16'
    if rebin == True:
        atm_model_name = 'phoenix_v16_rebin'
    
    teff_arr, z_arr, logg_arr = get_atmosphere_grid(atm_model_name)
    return teff_arr, z_arr, logg_arr


def get_BTSettl_2015_atmosphere(metallicity=0, temperature=2500, gravity=4, rebin=True):
    """
    Return atmosphere from CIFIST2011_2015 grid
    (`Allard et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012RSPTA.370.2765A/abstract>`_,
    `Baraffe et al. 2015 <https://ui.adsabs.harvard.edu/abs/2015A%26A...577A..42B/abstract>`_ )

    Grid originally downloaded from `website <https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/FITS/>`_.

    Grid Range:

    * Teff: 1200 - 7000 K
    * gravity: 2.5 - 5.5 cgs
    * [M/H] = 0

    Parameters
    ----------
    metallicity: float
        The stellar metallicity, in terms of [Z]

    temperature: float
        The stellar temperature, in units of K

    gravity: float
        The stellar gravity, in cgs units

    rebin: boolean
        If true, rebins the atmospheres so that they are the same
        resolution as the Castelli+04 atmospheres. Default is False,
        which is often sufficient synthetic photometry in most cases.
    """
    if rebin == True:
        atm_name = 'BTSettl_2015_rebin'
    else:
        atm_name = 'BTSettl_2015'

    try:
        sp = stsyn.grid_to_spec(atm_name, temperature, metallicity, gravity)
    except:
        # Check atmosphere catalog bounds
        (temperature, gravity, metallicity) = get_atmosphere_bounds(atm_name,
                                                   metallicity=metallicity,
                                                   temperature=temperature,
                                                   gravity=gravity)

        sp = stsyn.grid_to_spec(atm_name, temperature, metallicity, gravity)


    # Do some error checking
    idx = np.where(sp(sp.waveset) != 0)[0]
    if len(idx) == 0:
        print( 'Could not find BTSettl_2015 atmosphere model for')
        print( '  temperature = %d' % temperature)
        print( '  metallicity = %.1f' % metallicity)
        print( '  log gravity = %.1f' % gravity)

    return sp

def get_BTSettl_2015_atmosphere_grid(rebin=True):
    """
    Return atmosphere grid from CIFIST2011_2015 grid 
    (`Allard et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012RSPTA.370.2765A/abstract>`_, 
    `Baraffe et al. 2015 <https://ui.adsabs.harvard.edu/abs/2015A%26A...577A..42B/abstract>`_ )

    Grid originally downloaded from `website <https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/FITS/>`_.

    Grid Range:
    
    * Teff: 1200 - 7000 K
    * gravity: 2.5 - 5.5 cgs
    * [M/H] = 0
        
    rebin: boolean
        If true, rebins the atmospheres so that they are the same
        resolution as the Castelli+04 atmospheres. Default is True  
    """
    if rebin == True:
        atm_name = 'BTSettl_2015_rebin'
    else:
        atm_name = 'BTSettl_2015'
        
    teff_arr, z_arr, logg_arr = get_atmosphere_grid(atm_name)
    return teff_arr, z_arr, logg_arr


def get_BTSettl_atmosphere(metallicity=0, temperature=2500, gravity=4.5, rebin=True):
    """
    Return atmosphere from CIFIST2011 grid
    (`Allard et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012RSPTA.370.2765A/abstract>`_)

    Grid originally downloaded `here <https://phoenix.ens-lyon.fr/Grids/BT-Settl/>`_

    Notes
    ------
    Grid Range:

    * [M/H] = -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5

    Teff and gravity ranges depend on metallicity:

    [M/H] = -2.5

    * Teff: 2600 - 4600 K
    * gravity: 4.5 - 5.5

    [M/H] = -2.0

    * Teff: 2600 - 7000
    * gravity: 4.5 - 5.5

    [M/H] = -1.5

    * Teff: 2600 - 7000
    * gravity: 4.5 - 5.5

    [M/H] = -1.0

    * Teff: 2600 - 7000
    * gravity: Teff < 3200 --> 4.5 - 5.5; Teff > 3200 --> 2.5 - 5.5

    [M/H] = -0.5

    * Teff: 1000 -7000
    * gravity: Teff < 3000 --> 4.5 - 5.5; Teff > 3000 --> 3.0 - 6.0

    [M/H] = 0

    * Teff: 750 - 7000
    * gravity: Teff < 2500 --> 3.5 - 5.5; Teff > 2500 --> 0 - 5.5

    [M/H] = 0.5

    * Teff: 1000 - 5000
    * gravity: 3.5 - 5.0


    Alpha enhancement:

    * [M/H]= -0.0, +0.5 no anhancement
    * [M/H]= -0.5 with [alpha/H]=+0.2
    * [M/H]= -1.0, -1.5, -2.0, -2.5 with [alpha/H]=+0.4

    Parameters
    ----------
    metallicity: float
        The stellar metallicity, in terms of [Z]

    temperature: float
        The stellar temperature, in units of K

    gravity: float
        The stellar gravity, in cgs units

    rebin: boolean
        If true, rebins the atmospheres so that they are the same
        resolution as the Castelli+04 atmospheres. Default is False,
        which is often sufficient synthetic photometry in most cases.

        **PRINT STATEMENTS TO DEBUG
        **check get_atmosphere_bounds
        **comment out try/except clause and check break
    """
    if rebin == True:
        atm_name = 'BTSettl_rebin'
    else:
        atm_name = 'BTSettl'

    try:
        sp = stsyn.grid_to_spec(atm_name, temperature, metallicity, gravity)
    except:
        # Check atmosphere catalog bounds
        (temperature, gravity, metallicity) = get_atmosphere_bounds(atm_name,
                                                   metallicity=metallicity,
                                                   temperature=temperature,
                                                   gravity=gravity)

        sp = stsyn.grid_to_spec(atm_name, temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp(sp.waveset) != 0)[0]
    if len(idx) == 0:
        print( 'Could not find BTSettl_2015 atmosphere model for')
        print( '  temperature = %d' % temperature)
        print( '  metallicity = %.1f' % metallicity)
        print( '  log gravity = %.1f' % gravity)

    return sp

def get_Meisner2023_atmosphere(metallicity=0, temperature=1000, gravity=4.5, rebin=True):
    """
    Return atmosphere from Meisner2023 grid
    (`Meisner et al. 2023 <https://ui.adsabs.harvard.edu/abs/2023AJ....166...57M/abstract>`_)

    Grid originally downloaded `here <https://noctis.erc-atmo.eu/fsdownload/Q7MUSoCLR/meisner2023>`_

    Grid range:
    * Teff = 250 - 1200 K
    * gravity: 2.5 - 5.5 cgs (in steps of 0.5)
    * [M/H] = -1.0, -0.5, +0, +0.3

    """
    if rebin == True:
        atm_name = 'Meisner2023_rebin'
    else:
        atm_name = 'Meisner2023'

    try:
        sp = stsyn.grid_to_spec(atm_name, temperature, metallicity, gravity)
    except:
        # Check atmosphere catalog bounds
        (temperature, gravity, metallicity) = get_atmosphere_bounds(atm_name,
                                                   metallicity=metallicity,
                                                   temperature=temperature,
                                                   gravity=gravity)

        sp = stsyn.grid_to_spec(atm_name, temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp(sp.waveset) != 0)[0]
    if len(idx) == 0:
        print( 'Could not find Meisner2023 atmosphere model for')
        print( '  temperature = %d' % temperature)
        print( '  metallicity = %.1f' % metallicity)
        print( '  log gravity = %.1f' % gravity)

    return sp

def get_Meisner2023_atmosphere_grid(rebin=True):
    """
    Return atmosphere grid from Meisner2023.
    """
    if rebin == True:
        atm_name = 'Meisner2023_rebin'
    else:
        atm_name = 'Meisner2023'

    teff_arr, z_arr, logg_arr = get_atmosphere_grid(atm_name)
    return teff_arr, z_arr, logg_arr

def get_Phillips2020_atmosphere(metallicity=0, temperature=1000, gravity=4.5, rebin=True):
    """
        Return atmosphere from Phillips et al., 2020 using ATMO model
        (`Phillips et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020A%26A...637A..38P/abstract>`_)

        Grid originally downloaded `here <https://noctis.erc-atmo.eu/fsdownload/zyU96xA6o/phillips2020>`_

        Grid Range:
        * Teff: 200 - 3000 K
        * gravity: 2.5 - 5.5 cgs
        * [M/H] = 0
    """
    if rebin == True:
        atm_name = 'Phillips2020_rebin'
    else:
        atm_name = 'Phillips2020'

    try:
        sp = stsyn.grid_to_spec(atm_name, temperature, metallicity, gravity)
    except:
        # Check atmosphere catalog bounds
        (temperature, gravity, metallicity) = get_atmosphere_bounds(atm_name,
                                                   metallicity=metallicity,
                                                   temperature=temperature,
                                                   gravity=gravity)

        sp = stsyn.grid_to_spec(atm_name, temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp(sp.waveset) != 0)[0]
    if len(idx) == 0:
        print( 'Could not find Phillips2020 atmosphere model for')
        print( '  temperature = %d' % temperature)
        print( '  metallicity = %.1f' % metallicity)
        print( '  log gravity = %.1f' % gravity)

    return sp

def get_BTSettl_atmosphere_grid(rebin=True):
    """
    Return atmosphere grid from CIFIST2011 grid 
    (`Allard et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012RSPTA.370.2765A/abstract>`_)

    Grid originally downloaded `here <https://phoenix.ens-lyon.fr/Grids/BT-Settl/>`_

    Notes
    ------
    Grid Range:
    
    * [M/H] = -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5
    
    Teff and gravity ranges depend on metallicity:

    [M/H] = -2.5

    * Teff: 2600 - 4600 K
    * gravity: 4.5 - 5.5
    
    [M/H] = -2.0

    * Teff: 2600 - 7000
    * gravity: 4.5 - 5.5

    [M/H] = -1.5

    * Teff: 2600 - 7000
    * gravity: 4.5 - 5.5

    [M/H] = -1.0

    * Teff: 2600 - 7000
    * gravity: Teff < 3200 --> 4.5 - 5.5; Teff > 3200 --> 2.5 - 5.5 

    [M/H] = -0.5

    * Teff: 1000 -7000
    * gravity: Teff < 3000 --> 4.5 - 5.5; Teff > 3000 --> 3.0 - 6.0

    [M/H] = 0

    * Teff: 750 - 7000
    * gravity: Teff < 2500 --> 3.5 - 5.5; Teff > 2500 --> 0 - 5.5

    [M/H] = 0.5

    * Teff: 1000 - 5000
    * gravity: 3.5 - 5.0


    Alpha enhancement:

    * [M/H]= -0.0, +0.5 no anhancement
    * [M/H]= -0.5 with [alpha/H]=+0.2
    * [M/H]= -1.0, -1.5, -2.0, -2.5 with [alpha/H]=+0.4
        
    rebin: boolean
        If true, rebins the atmospheres so that they are the same
        resolution as the Castelli+04 atmospheres. Default is True.
    """
    if rebin == True:
        atm_name = 'BTSettl_rebin'
    else:
        atm_name = 'BTSettl'

    teff_arr, z_arr, logg_arr = get_atmosphere_grid(atm_name)
    return teff_arr, z_arr, logg_arr


def get_wdKoester_atmosphere(metallicity=0, temperature=20000, gravity=7):
    """
    Return white dwarf atmospheres from
    `Koester et al. 2010 <https://ui.adsabs.harvard.edu/abs/2010MmSAI..81..921K/abstract>`_

    Parameters
    ----------
    metallicity: float
        The stellar metallicity, in terms of [Z]

    temperature: float
        The stellar temperature, in units of K

    gravity: float
        The stellar gravity, in cgs units

    rebin: boolean
        If true, rebins the atmospheres so that they are the same
        resolution as the Castelli+04 atmospheres. Default is False,
        which is often sufficient synthetic photometry in most cases.
    """
    sp = stsyn.grid_to_spec('wdKoester', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp(sp.waveset) != 0)[0]
    if len(idx) == 0:
        print( 'Could not find WD Koester (Koester+ 2010 atmosphere model for')
        print( '  temperature = %d' % temperature)
        print( '  metallicity = %.1f' % metallicity)
        print( '  log gravity = %.1f' % gravity)

    return sp

def get_wdKoester_atmosphere_grid():
    """
    Return white dwarf grid atmospheres from  
    `Koester et al. 2010 <https://ui.adsabs.harvard.edu/abs/2010MmSAI..81..921K/abstract>`_
    """
    teff_arr, z_arr, logg_arr = get_atmosphere_grid('wdKoester')
    return teff_arr, z_arr, logg_arr

def get_atlas_phoenix_atmosphere(metallicity=0, temperature=5250, gravity=4):
    """
    Return atmosphere that is a linear merge of atlas ck04 model and phoenixV16.

    Only valid for temps between 5000 - 5500K, gravity from 0 = 5.0
    """
        
    try:
        sp = stsyn.grid_to_spec('merged_atlas_phoenix', temperature, metallicity, gravity)
    except:
        # Check atmosphere catalog bounds
        (temperature, gravity, metallicity) = get_atmosphere_bounds('merged_atlas_phoenix',
                                                   metallicity=metallicity,
                                                   temperature=temperature,
                                                   gravity=gravity)

        sp = stsyn.grid_to_spec('merged_atlas_phoenix', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp(sp.waveset) != 0)[0]
    if len(idx) == 0:
        print( 'Could not find ATLAS-PHOENIX merge atmosphere model for')
        print( '  temperature = %d' % temperature)
        print( '  metallicity = %.1f' % metallicity)
        print( '  log gravity = %.1f' % gravity)

    return sp

def get_atlas_phoenix_atmosphere_grid():
    """
    Return atmosphere that is a linear merge of atlas ck04 model and phoenixV16.

    Only valid for temps between 5000 - 5500K, gravity from 0 = 5.0 
    """
    teff_arr, z_arr, logg_arr = get_atmosphere_grid('merged_atlas_phoenix')
    return teff_arr, z_arr, logg_arr
        
def get_BTSettl_phoenix_atmosphere(metallicity=0, temperature=5250, gravity=4):
    """
    Return atmosphere that is a linear merge of BTSettl_CITFITS2011_2015 model
    and phoenixV16.

    Only valid for temps between 3200 - 3800K, gravity from 2.5 - 5.5
    """
    try:
        sp = stsyn.grid_to_spec('merged_BTSettl_phoenix', temperature, metallicity, gravity)
    except:
        # Check atmosphere catalog bounds
        (temperature, gravity, metallicity) = get_atmosphere_bounds('merged_BTSettl_phoenix',
                                                   metallicity=metallicity,
                                                   temperature=temperature,
                                                   gravity=gravity)

        sp = stsyn.grid_to_spec('merged_BTSettl_phoenix', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp(sp.waveset) != 0)[0]
    if len(idx) == 0:
        print( 'Could not find ATLAS-PHOENIX merge atmosphere model for')
        print( '  temperature = %d' % temperature)
        print( '  metallicity = %.1f' % metallicity)
        print( '  log gravity = %.1f' % gravity)

    return sp


def get_BTSettl_phoenix_atmosphere_grid():
    """
    Return atmosphere grid that is a linear merge of BTSettl_CITFITS2011_2015 model
    and phoenixV16.

    Only valid for temps between 3200 - 3800K, gravity from 2.5 - 5.5 
    """
    teff_arr, z_arr, logg_arr = get_atmosphere_grid('merged_BTSettl_phoenix')
    return teff_arr, z_arr, logg_arr
        
    
def get_BTSettl_meisner_atmosphere(metallicity=0, temperature=5250, gravity=4):
    """
    Return atmosphere that is a linear merge of BTSettl_CITFITS2011_2015 model
    and Meisner2023.

    Only valid for temps between 1000 - 1200K, gravity from 3.5 - 5.5
    """
    try:
        sp = stsyn.grid_to_spec('merged_BTSettl_meisner', temperature, metallicity, gravity)
    except:
        # Check atmosphere catalog bounds
        (temperature, gravity, metallicity) = get_atmosphere_bounds('merged_BTSettl_meisner',
                                                   metallicity=metallicity,
                                                   temperature=temperature,
                                                   gravity=gravity)

        sp = stsyn.grid_to_spec('merged_BTSettl_meisner', temperature, metallicity, gravity)

    # Do some error checking
    idx = np.where(sp(sp.waveset) != 0)[0]
    if len(idx) == 0:
        print( 'Could not find BTSettl-Meisner merge atmosphere model for')
        print( '  temperature = %d' % temperature)
        print( '  metallicity = %.1f' % metallicity)
        print( '  log gravity = %.1f' % gravity)

    return sp

def get_BTSettl_meisner_atmosphere_grid():
    """
    Return atmosphere grid that is a linear merge of BTSettl_CITFITS2011_2015 model
    and Meisner2023.

    Only valid for temps between 1000 - 1200K, gravity from 3.5 - 5.5
    """
    teff_arr, z_arr, logg_arr = get_atmosphere_grid('merged_BTSettl_meisner')
    return teff_arr, z_arr, logg_arr

#---------------------------------------------------------------------#
def get_merged_atmosphere(metallicity=0, temperature=20000, gravity=4.5, verbose=False,
                              rebin=True):
    """
    Return a stellar atmosphere from a suite of different model grids,
    depending  on the input temperature, (all values in K).

    Parameters
    ----------
    metallicity: float
        The stellar metallicity, in terms of [Z]

    temperature: float
        The stellar temperature, in units of K

    gravity: float
        The stellar gravity, in cgs units

    rebin: boolean
        If true, rebins the atmospheres so that they are the same
        resolution as the Castelli+04 atmospheres. Default is False,
        which is often sufficient synthetic photometry in most cases.

    verbose: boolean
        True for verbose output

    Notes
    -----
    The underlying stellar model grid used changes as a function of
    stellar temperature (in K):

    * T > 20,000: ATLAS
    * 5500 <= T < 20,000: ATLAS
    * 5000 <= T < 5500: ATLAS/PHOENIXv16 merge
    * 3800 <= T < 5000: PHOENIXv16

    For T < 3800, there is an additional gravity and metallicity
    dependence:

    If T < 3800 and [M/H] = 0:

    * T < 3800, logg < 2.5: PHOENIX v16
    * 3200 <= T < 3800, logg > 2.5: BTSettl_CIFITS2011_2015/PHOENIXV16 merge
    * 3200 < T <= 1200, logg > 2.5: BTSettl_CIFITS2011_2015
    * 1000 < T <= 1200, logg >= 3.5: BTSettl_CIFITS2011_2015/Meisner2023 merge
    * 250 < T <= 1000, logg > 2.5: Meisner2023

    Otherwise, if T < 3800 and [M/H] != 0:

    * T < 3800: PHOENIX v16

    References:

    * ATLAS: ATLAS9 models (`Castelli & Kurucz 2004 <http://www.stsci.edu/hst/observatory/crds/castelli_kurucz_atlas.html>`_)
    * PHOENIXv16 (`Husser et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013A%26A...553A...6H/abstract>`_)
    * BTSettl_CIFITS2011_2015: Baraffee+15, Allard+ (https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/SPECTRA/)
    * Meisner2023: ATMO 1D models (`Meisner et al. 2023 <https://ui.adsabs.harvard.edu/abs/2023AJ....166...57M/abstract>`_)

    LTE WARNING:

    The ATLAS atmospheres are calculated with LTE, and so they
    are less accurate when non-LTE conditions apply (e.g. T > 20,000
    K). Ultimately we'd like to add a non-LTE atmosphere grid for
    the hottest stars in the future.

    HOW BOUNDARIES BETWEEN MODELS ARE TREATED:

    At the boundary between two models grids a temperature range is defined
    where the resulting atmosphere is a weighted average between the two
    grids. Near one boundary one model
    is weighted more heavily, while at the other boundary the other
    model is weighted more heavily. These are calculated in the
    temperature ranges where we switch between model grids, to
    ensure a smooth transition.
    """
        
    # For T < 3800, atmosphere depends on metallicity + gravity.
    # If solar metallicity, use BTSettl 2015 grid. Only solar metallicity is
    # currently available here, so if non-solar metallicity, just stick with
    # the Phoenix grid
    if (temperature < 1000):
        if verbose:
            print( 'Meisner2023 atmosphere')
        return get_Meisner2023_atmosphere(metallicity=metallicity,
                                                temperature=temperature,
                                                gravity=gravity,
                                                rebin=rebin)

    if (temperature <= 1200) & (temperature >= 1000):
        if (gravity >= 3.5):
            if verbose:
                print( 'BTSettl/Meisner2023 merged atmosphere')
            return get_Meisner2023_atmosphere(metallicity=metallicity,
                                                    temperature=temperature,
                                                    gravity=gravity,
                                                    rebin=rebin)
        if (gravity < 3.5) & (gravity >=2.5):
            if verbose:
                print( 'Meisner2023 atmosphere')
            return get_Meisner2023_atmosphere(metallicity=metallicity,
                                                    temperature=temperature,
                                                    gravity=gravity,
                                                    rebin=rebin)

    if (temperature <= 3800) & (metallicity == 0):
        # High gravity are in BTSettl regime
        if (temperature <= 3200) & (gravity > 2.5):
            if verbose:
                print( 'BTSettl_2015 atmosphere')
            return get_BTSettl_2015_atmosphere(metallicity=metallicity,
                                                temperature=temperature,
                                                gravity=gravity,
                                                rebin=rebin)

        if (temperature >= 3200) & (temperature < 3800) & (gravity > 2.5):
            if verbose:
                print( 'BTSettl/Phoenixv16 merged atmosphere')
            return get_BTSettl_phoenix_atmosphere(metallicity=metallicity,
                                                temperature=temperature,
                                                gravity=gravity)

        # Low gravity is PHOENIX regime
        if gravity <= 2.5:
            if verbose:
                print( 'Phoenixv16 atmosphere')
            return get_phoenixv16_atmosphere(metallicity=metallicity,
                                            temperature=temperature,
                                            gravity=gravity,
                                            rebin=rebin)

    if (temperature <= 3800) & (metallicity != 0):
        if verbose:
            print( 'Phoenixv16 atmosphere')
        return get_phoenixv16_atmosphere(metallicity=metallicity,
                                        temperature=temperature,
                                        gravity=gravity,
                                        rebin=rebin)
    # For T > 3800, no metallicity or gravity dependence
    if (temperature >= 3800) & (temperature < 5000):
        if verbose:
            print( 'Phoenixv16 atmosphere')
        return get_phoenixv16_atmosphere(metallicity=metallicity,
                                      temperature=temperature,
                                      gravity=gravity,
                                      rebin=rebin)

    if (temperature >= 5000) & (temperature < 5500):
        if verbose:
            print( 'ATLAS/Phoenix merged atmosphere')
        return get_atlas_phoenix_atmosphere(metallicity=metallicity,
                                        temperature=temperature,
                                        gravity=gravity)

    if (temperature >= 5500) & (temperature < 20000):
        if verbose:
            print( 'ATLAS merged atmosphere')
        return get_castelli_atmosphere(metallicity=metallicity,
                                      temperature=temperature,
                                      gravity=gravity)

    if temperature >= 20000:
        if verbose:
            print( 'Still ATLAS merged atmosphere')
        return get_castelli_atmosphere(metallicity=metallicity,
                                       temperature=temperature,
                                       gravity=gravity)

        #print('CMFGEN')
        #return get_cmfgenRot_atmosphere_closest(metallicity=metallicity,
        #                               temperature=temperature,
        #                               gravity=gravity)


def get_merged_atmosphere_grid(rebin=True):

    # temp array, metallicity array, logg array
    Meisner2023_atmosphere_arrs = np.array(get_Meisner2023_atmosphere_grid(rebin))

    BTSettl_2015_atmosphere_arrs = np.array(get_BTSettl_2015_atmosphere_grid(rebin))

    BTSettl_phoenix_atmosphere_arrs = np.array(get_BTSettl_phoenix_atmosphere_grid())

    phoenixv16_atmosphere_arrs = np.array(get_phoenixv16_atmosphere_grid(rebin))
    atlas_phoenix_atmosphere_arrs = np.array(get_atlas_phoenix_atmosphere_grid())

    castelli_atmosphere_arrs = np.array(get_castelli_atmosphere_grid())

    Meisner2023_atmosphere_idxs = np.where(Meisner2023_atmosphere_arrs[0] <= 1200)
    Meisner2023_atmosphere_arrs = Meisner2023_atmosphere_arrs[:,Meisner2023_atmosphere_idxs]

    BTSettl_2015_atmosphere_idxs = np.where((BTSettl_2015_atmosphere_arrs[0] > 1200) &\
                                            (BTSettl_2015_atmosphere_arrs[0] <= 3200) &\
                                            (BTSettl_2015_atmosphere_arrs[1] == 0) &\
                                            (BTSettl_2015_atmosphere_arrs[2] > 2.5))
    BTSettl_2015_atmosphere_arrs = BTSettl_2015_atmosphere_arrs[:,BTSettl_2015_atmosphere_idxs]

    BTSettl_phoenix_atmosphere_idxs = np.where((BTSettl_phoenix_atmosphere_arrs[0] >= 3200) &\
                                               (BTSettl_phoenix_atmosphere_arrs[0] < 3800) &\
                                               (BTSettl_phoenix_atmosphere_arrs[1] == 0) &\
                                               (BTSettl_phoenix_atmosphere_arrs[2] > 2.5))
    BTSettl_phoenix_atmosphere_arrs = BTSettl_phoenix_atmosphere_arrs[:,BTSettl_phoenix_atmosphere_idxs]

    phoenixv16_atmosphere_idxs = np.where(((phoenixv16_atmosphere_arrs[0] <= 3800) & (phoenixv16_atmosphere_arrs[1] == 0) & (phoenixv16_atmosphere_arrs[2] <= 2.5)) |\
                                          ((phoenixv16_atmosphere_arrs[0] <= 3800) & (phoenixv16_atmosphere_arrs[1] != 0)) |\
                                          ((phoenixv16_atmosphere_arrs[0] >= 3800) & (phoenixv16_atmosphere_arrs[0] < 5000)))
    phoenixv16_atmosphere_arrs = phoenixv16_atmosphere_arrs[:,phoenixv16_atmosphere_idxs]

    atlas_phoenix_atmosphere_idxs = np.where((atlas_phoenix_atmosphere_arrs[0] >= 5000) & (atlas_phoenix_atmosphere_arrs[0] < 5500))
    atlas_phoenix_atmosphere_arrs = atlas_phoenix_atmosphere_arrs[:,atlas_phoenix_atmosphere_idxs]

    castelli_atmosphere_idxs = np.where(castelli_atmosphere_arrs[0] >= 5500)
    castelli_atmosphere_arrs = castelli_atmosphere_arrs[:,castelli_atmosphere_idxs]

    super_tarr = np.concatenate((Meisner2023_atmosphere_arrs[0][0], BTSettl_2015_atmosphere_arrs[0][0],
                                BTSettl_phoenix_atmosphere_arrs[0][0],
                                phoenixv16_atmosphere_arrs[0][0], atlas_phoenix_atmosphere_arrs[0][0], 
                                castelli_atmosphere_arrs[0][0]))

    super_zarr = np.concatenate((Meisner2023_atmosphere_arrs[1][0], BTSettl_2015_atmosphere_arrs[1][0],
                                BTSettl_phoenix_atmosphere_arrs[1][0],
                                phoenixv16_atmosphere_arrs[1][0], atlas_phoenix_atmosphere_arrs[1][0], 
                                castelli_atmosphere_arrs[1][0]))

    super_loggarr = np.concatenate((Meisner2023_atmosphere_arrs[2][0], BTSettl_2015_atmosphere_arrs[2][0],
                                BTSettl_phoenix_atmosphere_arrs[2][0],
                                phoenixv16_atmosphere_arrs[2][0], atlas_phoenix_atmosphere_arrs[2][0], 
                                castelli_atmosphere_arrs[2][0]))

    return super_tarr, super_zarr, super_loggarr

def get_merged_atmosphere_w_bb_supplement(metallicity=0, temperature=20000, gravity=4.5, verbose=False,
                              rebin=True):
    """
    Return a stellar atmosphere from a suite of different model grids, 
    depending  on the input temperature, (all values in K).
    IF OUTSIDE SUPPORTED GRID, WILL RETURN BB ATMOSPHERE

    Parameters
    ----------
    metallicity: float
        The stellar metallicity, in terms of [Z]

    temperature: float
        The stellar temperature, in units of K

    gravity: float
        The stellar gravity, in cgs units
        
    rebin: boolean
        If true, rebins the atmospheres so that they are the same
        resolution as the Castelli+04 atmospheres. Default is False,
        which is often sufficient synthetic photometry in most cases.

    verbose: boolean
        True for verbose output

    Notes
    -----
    The underlying stellar model grid used changes as a function of 
    stellar temperature (in K):

    * T > 20,000: ATLAS
    * 5500 <= T < 20,000: ATLAS
    * 5000 <= T < 5500: ATLAS/PHOENIXv16 merge
    * 3800 <= T < 5000: PHOENIXv16

    For T < 3800, there is an additional gravity and metallicity
    dependence:

    If T < 3800 and [M/H] = 0: 
    
    * T < 3800, logg < 2.5: PHOENIX v16
    * 3200 <= T < 3800, logg > 2.5: BTSettl_CIFITS2011_2015/PHOENIXV16 merge
    * 1200 < T <= 3200, logg > 2.5: BTSettl_CIFITS2011_2015
    * 1000 <= T <= 1200, logg >= 2.5: Meisner2023
    * 250 <= T < 1000: Meisner2023

    Otherwise, if T < 3800 and [M/H] != 0:
    
    * T < 3800: PHOENIX v16

    References:

    * ATLAS: ATLAS9 models (`Castelli & Kurucz 2004 <http://www.stsci.edu/hst/observatory/crds/castelli_kurucz_atlas.html>`_)
    * PHOENIXv16 (`Husser et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013A%26A...553A...6H/abstract>`_)
    * BTSettl_CIFITS2011_2015: Baraffee+15, Allard+ (https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/SPECTRA/)
    * Meisner2023: ATMO 1D models (`Meisner et al. 2023 <https://ui.adsabs.harvard.edu/abs/2023AJ....166...57M/abstract>`_)

    LTE WARNING: 

    The ATLAS atmospheres are calculated with LTE, and so they
    are less accurate when non-LTE conditions apply (e.g. T > 20,000
    K). Ultimately we'd like to add a non-LTE atmosphere grid for
    the hottest stars in the future.

    HOW BOUNDARIES BETWEEN MODELS ARE TREATED: 

    At the boundary between two models grids a temperature range is defined 
    where the resulting atmosphere is a weighted average between the two 
    grids. Near one boundary one model
    is weighted more heavily, while at the other boundary the other 
    model is weighted more heavily. These are calculated in the 
    temperature ranges where we switch between model grids, to 
    ensure a smooth transition.
    """
    
    if (gravity >= 9.8):
        if verbose:
            print('BB atmosphere')
        return get_bb_atmosphere(temperature=temperature,
                                 metallicity=metallicity,
                                 gravity=gravity,
                                 verbose=verbose)
    if (temperature < 4.6e3) & (gravity >= 6.5):
        if verbose:
            print('BB atmosphere')
        return get_bb_atmosphere(temperature=temperature,
                                 metallicity=metallicity,
                                 gravity=gravity,
                                 verbose=verbose)
        
    if (temperature < 3.5e3) & (gravity  < 6.5) & (gravity > 6):
        if verbose:
            print('BB atmosphere')
        return get_bb_atmosphere(temperature=temperature,
                                 metallicity=metallicity,
                                 gravity=gravity,
                                 verbose=verbose)

    if gravity > 6:
        if verbose:
            print('WD or BB atmosphere')
        return get_wd_atmosphere(metallicity=metallicity,
                                 temperature=temperature,
                                 gravity=gravity,
                                 verbose=verbose)

    return get_merged_atmosphere(metallicity=metallicity,
                                 temperature=temperature,
                                 gravity=gravity,
                                 verbose=verbose,
                                 rebin=rebin)

def get_merged_atmosphere_w_bb_supplement_grid(bb_supplement_tarr='default', bb_supplement_zarr='default', bb_supplement_loggarr='default', rebin=True):
    """
    Return the atmosphere parameter grid used by
    get_merged_atmosphere_w_bb_supplement.

    This is the standard merged atmosphere grid plus the white dwarf grid and
    blackbody supplement points for high-gravity regions.
    """
    super_tarr, super_zarr, super_loggarr = get_merged_atmosphere_grid(rebin=rebin)

    wd_tarr, wd_zarr, wd_loggarr = get_wdKoester_atmosphere_grid()
    super_tarr = np.concatenate((super_tarr, wd_tarr))
    super_zarr = np.concatenate((super_zarr, wd_zarr))
    super_loggarr = np.concatenate((super_loggarr, wd_loggarr))

    if bb_supplement_tarr == 'default':
        X_highg, Y_highg = np.meshgrid(np.logspace(np.log10(2e3), np.log10(2e4), 35),
                                       np.linspace(9.8, 11.6, 8))
        X_cool_highg, Y_cool_highg = np.meshgrid(np.logspace(np.log10(2e3), np.log10(4.6e3), 20),
                                                 np.linspace(6.5, 9.7, 10))
        X_cool_midg, Y_cool_midg = np.meshgrid(np.logspace(np.log10(2e3), np.log10(3.5e3), 15),
                                               np.linspace(6.1, 6.4, 4))

        bb_supplement_tarr = np.concatenate((X_highg.ravel(), X_cool_highg.ravel(), X_cool_midg.ravel()))
        bb_supplement_loggarr = np.concatenate((Y_highg.ravel(), Y_cool_highg.ravel(), Y_cool_midg.ravel()))
        bb_supplement_zarr = np.zeros(len(bb_supplement_tarr))
        

    super_tarr = np.concatenate((super_tarr, bb_supplement_tarr))
    super_zarr = np.concatenate((super_zarr, bb_supplement_zarr))
    super_loggarr = np.concatenate((super_loggarr, bb_supplement_loggarr))

    return super_tarr, super_zarr, super_loggarr
    
def get_wd_atmosphere(metallicity=0, temperature=20000, gravity=4, verbose=False):
    """
    Return the white dwarf atmosphere from
    `Koester et al. 2010 <https://ui.adsabs.harvard.edu/abs/2010MmSAI..81..921K/abstract>`_.
    If desired parameters are
    outside of grid, return a blackbody spectrum instead

    Parameters
    ----------
    metallicity: float
        The stellar metallicity, in terms of [Z]

    temperature: float
        The stellar temperature, in units of K

    gravity: float
        The stellar gravity, in cgs units

    rebin: boolean
        If true, rebins the atmospheres so that they are the same
        resolution as the Castelli+04 atmospheres. Default is False,
        which is often sufficient synthetic photometry in most cases.

    verbose: boolean
        True for verbose output
    """
    try:
        if verbose:
            print('wdKoester atmosphere')

        return get_wdKoester_atmosphere(metallicity=metallicity,
                                            temperature=temperature,
                                            gravity=gravity)

    except stsyn_exceptions.ParameterOutOfBounds:
        # Use a black-body atmosphere.
        bbspec = get_bb_atmosphere(temperature=temperature, verbose=verbose)
        return bbspec


def get_bd_atmosphere(metallicity=0, temperature=1000, gravity=4, verbose=False):
    """
    Return the brown dwarf atmosphere from
    `Meisner et al. 2023 <https://ui.adsabs.harvard.edu/abs/2023AJ....166...57M/abstract>`_.
    If desired parameters are
    outside of grid, return a blackbody spectrum instead

    Parameters
    ----------
    metallicity: float
        The stellar metallicity, in terms of [Z]

    temperature: float
        The stellar temperature, in units of K

    gravity: float
        The stellar gravity, in cgs units

    rebin: boolean
        If true, rebins the atmospheres so that they are the same
        resolution as the Castelli+04 atmospheres. Default is False,
        which is often sufficient synthetic photometry in most cases.

    verbose: boolean
        True for verbose output
    """
    try:
        if verbose:
            print('Meisner2023 atmosphere')

        return get_Meisner2023_atmosphere(metallicity=metallicity,
                                            temperature=temperature,
                                            gravity=gravity)

    except stsyn_exceptions.ParameterOutOfBounds:
        # Use a black-body atmosphere
        bbspec = get_bb_atmosphere(temperature=temperature, verbose=verbose)
        return bbspec

def get_bb_atmosphere(metallicity=None, temperature=20_000, gravity=None,
                      verbose=False, rebin=None,
                      wave_min=500, wave_max=50_000, wave_num=20_000):
    """
    Return a blackbody spectrum

    Parameters
    ----------
    temperature: float, default=20_000
        The stellar temperature, in units of K
    wave_min: float, default=500
        Sets the minimum wavelength (in Angstroms) of the wavelength range
        for the blackbody spectrum
    wave_max: float, default=50_000
        Sets the maximum wavelength (in Angstroms) of the wavelength range
        for the blackbody spectrum
    wave_num: int, default=20_000
        Sets the number of wavelength points in the wavelength range
        Note: the wavelength range is evenly spaced in log space
    """
    if ((metallicity is not None) or (gravity is not None) or
        (rebin is not None)):
        warnings.warn(
            'Only `temperature` keyword is used for black-body atmosphere'
        )

    if verbose:
        print('Black-body atmosphere')

    # Log-spaced wavelength grid (Angstrom); synphot BlackBody in PHOTLAM → FLAM
    w_grid = (
        np.logspace(np.log10(wave_min), np.log10(wave_max), wave_num, dtype=np.float64)
        * u.AA
    )
    bb = SourceSpectrum(BlackBody1D, temperature=temperature)
    y_photlam = bb(w_grid)
    y_flam = convert_flux(w_grid, y_photlam, su.FLAM)
    bbspec = SourceSpectrum(Empirical1D, points=w_grid, lookup_table=y_flam)

    # `BlackBody` spectrum is normalized to solar radius star at 1 kiloparsec.
    # Need to remove this normalization for SPISEA by multiplying bbspec
    # by (1000 * 1 parsec / 1 Rsun)**2 = (1000 * 3.08e18 cm / 6.957e10 cm)**2
    bbspec *= (1000 * 3.086e18 / 6.957e10)**2

    return bbspec


#--------------------------------------#
# Atmosphere formatting functions
#--------------------------------------#

def download_CMFGEN_atmospheres(Table_rot, Table_norot):
    """
    Downloads CMFGEN models from
    https://sites.google.com/site/fluxesandcontinuum/home;
    these contain continuum as well as lines.

    Table_rot, Table_norot are tables with the file prefixes
    and model atmosphere parameters, taken by hand from the
    Fierro+15 paper

    Website addresses are hardcoded

    Puts downloaded models in the current working directory.
    """
    print( 'WARNING: THIS DOES NOT COMPLETELY WORK')
    print( '**********************')
    t_rot = Table.read(Table_rot, format='ascii')
    t_norot = Table.read(Table_norot, format='ascii')

    tables = [t_rot, t_norot]
    filenames = [t_rot['col1'], t_norot['col1']]

    # Hardcoded list of webiste addresses
    web_base1 = 'https://sites.google.com/site/fluxesandcontinuum/home/'
    web_base2 = 'https://sites.google.com/site/modelsobmassivestars/'
    web = [web_base1+'009-solar-masses/',web_base1+'012-solar-masses/',
           web_base1+'015-solar-masses/',web_base1+'020-solar-masses/',
           web_base1+'025-solar-masses/',web_base2+'009-solar-masses-tracks/',
           web_base2+'040-solar-masses/',web_base2+'060-solar-masses/',
           web_base1+'085-solar-masses/',web_base1+'120-solar-masses/']
    # Array of masses that matches the website addresses
    mass_arr = np.array([9.,12.,15.,20.,25.,32.,40.,60.,85.,120.])

    # Loop through rotating and unrotating case. First loop is rot, second unrot
    for i in range(2):
        # Extract masses from filenames
        masses = []
        for j in filenames[i]:
            tmp = j.split('m')
            mass = float(tmp[1][:-1])
            masses.append(mass)

        # Download the models webpage by webpage. A bit tricky because masses
        # change slightly within a particular website. THIS IS WHAT FAILS
        for j in range(len(web)):
            if j == 0:
                good = np.where( (masses <= mass_arr[j]) )
            else:
                g = j - 1
                good = np.where( (masses <= mass_arr[j]) &
                                (masses > mass_arr[g]) )
            # Use wget command to pull down the files, and unzip them
            for k in good[0]:
                full = web[j]+'{1:s}.flx.zip'.format(mass_arr[j],filenames[i][k])
                os.system('wget ' + full)
                os.system('unzip '+ filenames[i][k] + '.flx.zip')

    return

def organize_CMFGEN_atmospheres(path_to_dir):
    """
    Organize CMFGEN grid from Fierro+15
    (http://www.astroscu.unam.mx/atlas/index.html)
    into rot and noRot directories

    path_to_dir is from current working directory to directory
    containing the downloaded models. Assumed that models
    and tables describing parameters are in this directory.

    Tables describing parameters MUST be named Table_rot.txt,
    Table_noRot.txt. Made by hand from Tables 3, 4 in Fierro+15.
    These are located in same original directory as atmosphere files

    Will separate files into 2 subdirectories, one rotating and
    the other non-rotating

    *Can't have any other files starting with "t" in model directory to start!*
    """
    # First, record current working directory to return to later
    start_dir = os.getcwd()

    # Enter atmosphere directory, collect rotating and non-rotating
    # file names (assumed to all start with "t")
    os.chdir(path_to_dir)
    rot_models = glob.glob("t*r.flx*")
    noRot_models = glob.glob("t*n.flx*")

    # Separate into different subdirectories
    if os.path.exists('cmfgenF15_rot'):
        pass
    else:
        os.mkdir('cmfgenF15_rot')
        os.mkdir('cmfgenF15_noRot')

    for mod in rot_models:
        cmd = 'mv {0:s} cmfgenF15_rot'.format(mod)
        os.system(cmd)

    for mod in noRot_models:
        cmd = 'mv {0:s} cmfgenF15_noRot'.format(mod)
        os.system(cmd)

    # Also move Tables with model parameters into correct directory
    os.system('mv Table_rot.txt cmfgenF15_rot')
    os.system('mv Table_noRot.txt cmfgenF15_noRot')

    # Return to original directory
    os.chdir(start_dir)

    return

def make_CMFGEN_catalog(path_to_dir):
    """
    Create cdbs catalog.fits of CMFGEN grid from Fierro+15
    (http://www.astroscu.unam.mx/atlas/index.html).
    THIS IS STEP 2, after organize_CMFGEN_atmospheres has
    been run.

    path_to_dir is from current working directory to directory
    containing the rotating or non-rotating models (i.e. cmfgenF15_rot). Also,
    needs to be a Table*.txt file which contains the parameters for all of the
    original models, since params in filename are not precise enough

    Will create catalog.fits file in atmosphere directory with
    description of each model

    *Can't have any other files starting with "t" in model directory to start!*
    """
    # Record current working directory for later
    start_dir = os.getcwd()

    # Enter atmosphere directory
    os.chdir(path_to_dir)

    # Extract parameters for each atmosphere
    # Note: can't rely on filename for this because not precise enough!!

    #---------OLD: GETTING PARAMS FROM FILENAME-------#
    # Collect file names (assumed to all start with "t")
    #files = glob.glob("t*")
    #for name in files:
    #    tmp = name.split('l')
    #    temp = float(tmp[0][1:]) * 100.0 # In kelvin

    #    lumtmp = tmp[1].split('_')
    #    lum = float(lumtmp[0][:-5]) * 1000.0 # In L_sun

    #    mass = float(lumtmp[0][5:-1]) # In M_sun

        # Need to calculate log g from T and L (cgs)
    #    lum_sun = 3.846 * 10**33 # erg/s
    #    M_sun = 2 * 10**33 # g
    #    G_si = 6.67 * 10**(-8) # cgs
    #    sigma_si = 5.67 * 10**(-5) # cgs

    #    g = (G_si * mass * M_sun * 4 * np.pi * sigma_si * temp**4) / \
    #      (lum * lum_sun)
    #    logg = np.log10(g)
    #---------------------------------------------------#

    # Read table with atmosphere params
    table = glob.glob('Table_*')
    t = Table.read(table[0], format = 'ascii')
    names = t['col1']
    temps = t['col2']
    logg = t['col4']

    # Create catalog.fits file
    index_str = []
    name_str = []
    for i in range(len(names)):
        index = '{0:5.0f},0.0,{1:3.2f}'.format(temps[i], logg[i])

        #---NOTE: THE FOLLOWING DEPENDS ON FINAL LOCATION OF CATALOG FILE---#
        #path = path_to_dir + '/' + names[i]
        path = names[i] + '.fits[Flux]'

        index_str.append(index)
        name_str.append(path)

    catalog = Table([index_str, name_str], names = ('INDEX', 'FILENAME'))

    # Create catalog.fits file in directory with the models
    catalog.write('catalog.fits', format = 'fits')

    # Move back to original directory, create the catalog.fits file
    os.chdir(start_dir)

    return

def cdbs_cmfgen(path_to_dir, path_to_cdbs_dir):
    """
    Code to put cmfgen models into cdbs format and adds proper unit keyword in
    fits header. Save as fits file

    path_to_dir goes from current directory to cmfgen_rot or cmfgen_norot
    directory with the *.flx models. Note that these files have already been
    organized using organize_CMFGEN_atmospheres code.

    path_to_cdbs_dir goes from current directory to cdbs/grid/cmfgen_rot or
    cmfgen_norot directory. Will copy new fits files to this directory.
    This directory must already exist!
    """
    # Save starting directory for later, move into path_to_dir directory
    start_dir = os.getcwd()
    os.chdir(path_to_dir)

    # Collect the filenames, make necessary changes to each one
    files = glob.glob('*.flx')

    # Need to make brand-new fits tables with data we want.
    counter = 0
    for i in files:
        counter += 1
        # Open file, extract useful info
        t = Table.read(i, format='ascii')
        wave = t['col1']
        flux = t['col2'] # Flux is already in erg/cm^2/s/A

        # Need to eliminate duplicate entries (synphot crashes)
        unique = np.unique(wave, return_index=True)
        wave = wave[unique[1]]
        flux = flux[unique[1]]

        # Make fits table from individual columns.
        c0 = fits.Column(name='Wavelength', format='D', array=wave)
        c1 = fits.Column(name='Flux', format='E', array=flux)

        cols = fits.ColDefs([c0, c1])
        tbhdu = fits.BinTableHDU.from_columns(cols)

        #Adding unit keywords
        tbhdu.header['TUNIT1'] = 'ANGSTROM'
        tbhdu.header['TUNIT2'] = 'FLAM'

        prihdu = fits.PrimaryHDU()

        finalhdu = fits.HDUList([prihdu, tbhdu])
        finalhdu.writeto(i[:-4]+'.fits', overwrite=True)

        print( 'Done {0:2.0f} of {1:2.0f}'.format(counter, len(files)))

    # Return to original directory, copy over new .fits files to cdbs directory
    os.chdir(start_dir)
    cmd = 'mv {0:s}/*.fits {1:s}'.format(path_to_dir, path_to_cdbs_dir)
    os.system(cmd)

    return

def rebin_cmfgen(cdbs_path, rot=True):
    """
    Rebin cmfgen_rot and cmfgen_norot models to atlas ck04 resolution;
    this makes spectrophotometry MUCH faster

    cdbs_path: path to cdbs directory
    rot=True for rotating models (cmfgen_rot), False for non-rotating models

    makes new directory in cdbs/grid: cmfgen_rot_rebin or cmfgen_norot_rebin
    """
    # Get an atlas ck04 model, we will use this to set wavelength grid
    sp_atlas = get_castelli_atmosphere()

    # Open a fits table for an existing cmfgen model; we will steal the header.
    # Also define paths to new rebin directories
    if rot == True:
        tmp = cdbs_path+'/grid/cmfgen_rot/t0200l0008m009r.fits'
        path = cdbs_path+'/grid/cmfgen_rot_rebin/'
        orig_path = cdbs_path+'/grid/cmfgen_rot/'
    else:
        tmp = cdbs_path+'/grid/cmfgen_norot/t0200l0007m009n.fits'
        path = cdbs_path+'/grid/cmfgen_norot_rebin/'
        orig_path = cdbs_path+'/grid/cmfgen_norot/'

    cmfgen_hdu = fits.open(tmp)
    header0 = cmfgen_hdu[0].header
    # Create rebin directories if they don't already exist. Copy over
    # catalog.fits file from original directory (will be the same)
    if not os.path.exists(path):
        os.mkdir(path)
        cmd = 'cp {0:s}catalog.fits {1:s}'.format(orig_path, path)
        os.system(cmd)

    # Read in the catalog.fits file
    cat = fits.getdata(orig_path + 'catalog.fits')
    files_all = [cat[ii][1].split('[')[0] for ii in range(len(cat))]

    # First column in new files will be for [atlas] wavelength
    c0 = fits.Column(name='Wavelength', format='D', array=sp_atlas.waveset.to(u.AA).value)

    # For each catalog.fits entry, read the unbinned spectrum and rebin to
    # the atlas resolution. Make a new fits file in rebin directory
    count = 0
    for ff in range(len(files_all)):
        count += 1
        # Extract the temp, Z, logg
        vals = cat[ff][0].split(',')
        temp = float(vals[0])
        metal = float(vals[1])
        grav = float(vals[2])

        # Fetch the spectrum
        if rot == True:
            sp = stsyn.grid_to_spec('cmfgen_rot', temp, metal, grav)
        else:
            sp = stsyn.grid_to_spec('cmfgen_norot', temp, metal, grav)

        # Rebin
        flux_rebin = rebin_spec(sp.waveset, sp, sp_atlas.waveset)
        c1 = fits.Column(name='Flux', format='E', array=flux_rebin)

        # Make the FITS file from the columns with header
        cols = fits.ColDefs([c0,c1])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        prihdu = fits.PrimaryHDU(header=header0)
        tbhdu.header['TUNIT1'] = 'ANGSTROM'
        tbhdu.header['TUNIT2'] = 'FLAM'

        # Write hdu to new directory with same filename
        finalhdu = fits.HDUList([prihdu, tbhdu])
        finalhdu.writeto(path+files_all[ff])

        print( 'Finished file {0} of {1}'.format(count, len(files_all))            )
    return


def organize_PHOENIXv16_atmospheres(path_to_dir, met_str='m00'):
    """
    Construct the Phoenix Husser+13 atmopsheres for each model. Combines the
    fluxes from the *HiRES.fits files and the wavelengths of the
    WAVE_PHONEIX-ACES-AGSS-COND-2011.fits file.

    path_to_dir is the path to the directory containing all of the downloaded
    files

    met_str is the name of the current metallicity

    Creates new fits files for each atmosphere: phoenix<metallicity>_<temp>.fits,
    which contains columns for the log g (column header = g#.#). Puts
    atmospheres in new directory phoenixm00
    """
    # Save current directory for return later, move into working dir
    start_dir = os.getcwd()
    os.chdir(path_to_dir)

    # If it doesn't already exist, create the current metallicity subdirectory
    sub_dir = '../phoenix{0}'.format(met_str)
    if os.path.exists(sub_dir):
        pass
    else:
        os.mkdir(sub_dir)

    # Extract wavelength array, make column for later
    wavefile = fits.open('WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
    wave = wavefile[0].data
    wavefile.close()
    wave_col = Column(wave, name = 'WAVELENGTH')

    # Create temp array for Husser+13 grid (given in paper)
    temp_arr = np.arange(2300, 7001, 100)
    temp_arr = np.append(temp_arr, np.arange(7000, 12001, 200))

    print( 'Looping though all temps')
    # For each temp, build file containing the flux for all gravities
    i = 0
    for temp in temp_arr:
        files = glob.glob('lte{0:05d}-*-HiRes.fits'.format(temp))
        files.sort()
        # Start the table with the wavelength column
        t = Table()
        t.add_column(wave_col)
        for f in files:
            # Extract the logg out of filename
            logg = f[9:13]

            # Extract fluxes from file
            spectrum = fits.open(f)
            flux = spectrum[0].data
            spectrum.close()

            # Make Column object with fluxes, add to table
            col = Column(flux, name = 'g{0:2.1f}'.format(float(logg)))
            t.add_column(col)

        # Now, construct final fits file for the given temp
        outname = 'phoenix{0}_{1:05d}.fits'.format(met_str, temp)
        t.write('{0}/{1}'.format(sub_dir, outname), format = 'fits', overwrite = True)

        # Progress counter for user
        i += 1
        print( 'Done {0:d} of {1:d}'.format(i, len(temp_arr)))

    # Return to original directory
    os.chdir(start_dir)
    return

def make_PHOENIXv16_catalog(path_to_dir, met_str='m00'):
    """
    Makes catalog.fits file for Husser+13 phoenix models. Assumes that
    organize_PHOENIXv16_atmospheres has been run already, and that the models lie
    in subdirectory phoenix[met_str].

    path_to_directory is the path to the directory with the reformatted
    models (i.e. the output from construct_atmospheres, phoenix[met_str])

    Puts catalog.fits file in directory the user starts in
    """
    # Save starting directory for later, move into working directory
    start_dir = os.getcwd()
    os.chdir(path_to_dir)

    # Extract metallicity from metallicity string
    met = float(met_str[1]) + (float(met_str[2]) * 0.1)
    if 'm' in met_str:
        met *= -1.

    # Collect the filenames. Each is a unique temp with many different log g's
    files = glob.glob('phoenix*.fits')
    files.sort()

    # Create the catalog.fits file, row by row
    index_arr = []
    filename_arr = []
    for i in files:
        # Get log g values from the column header the file
        t = Table.read(i, format='fits')
        keys = t.keys()
        logg_vals = keys[1:]

        # Extract temp from filename
        name = i.split('_')
        temp = float(name[1][:-5])
        for j in logg_vals:
            logg = float(j[1:])
            index = '{0:5.0f},{1:2.1f},{2:2.1f}'.format(temp, met, logg)
            filename = path_to_dir + i + '[' + j + ']'
            # Add row to table
            index_arr.append(index)
            filename_arr.append(filename)

    catalog = Table([index_arr, filename_arr], names=('INDEX', 'FILENAME'))

    # Return to starting directory, write catalog
    os.chdir(start_dir)

    if os.path.exists('catalog.fits'):
        from astropy.table import vstack

        prev_catalog = Table.read('catalog.fits', format='fits')
        joined_catalog = vstack([prev_catalog, catalog])

        joined_catalog.write('catalog.fits', format='fits', overwrite=True)
    else:
        catalog.write('catalog.fits', format='fits', overwrite=True)

    return

def cdbs_PHOENIXv16(path_to_cdbs_dir):
    """
    Put the PHOENIXv16 (Husser+13) fits files into cdbs format. This primarily
    consists of adjusting the flux units from [erg/s/cm^2/cm] to [erg/s/cm^2/A]
    and adding the appropriate keywords to the fits header.

    path_to_cdbs_dir goes from current working directory to phoenix[met] directory
    in cdbs/grids/phoenix_v16. Note that these files have already been organized
    using organize_PHOENIXv16_atmospheres code.

    Overwrites original files in directory
    """
    # Save starting directory for later, move into working directory
    start_dir = os.getcwd()
    os.chdir(path_to_cdbs_dir)

    # Collect the filenames, make necessary changes to each one
    files = glob.glob('phoenix*.fits')

    ## Need to sort filenames; glob doesn't always give them in order
    files.sort()

    # Need to make brand-new fits tables with data we want.
    counter = 0
    for i in files:
        counter += 1

        # Read in current FITS table
        cur_table = Table.read(i, format='fits')

        cur_table.columns[0].name = 'Wavelength'

        num_cols = len(cur_table.colnames)

        # Multiplying each flux column by 10^-8 for conversion
        for cur_col_index in range(1, num_cols, 1):
            cur_col_name = cur_table.colnames[cur_col_index]
            cur_table[cur_col_name] = cur_table[cur_col_name] * 10.**-8


        # Construct new FITS file based on old one
        hdu = fits.open(i)
        header_0 = hdu[0].header
        header_1 = hdu[1].header
        sci = hdu[1].data

        tbhdu = fits.table_to_hdu(cur_table)

        # Copying over the older headers, adding unit keywords
        prihdu = fits.PrimaryHDU(header=header_0)
        tbhdu.header['TUNIT1'] = 'ANGSTROM'
        tbhdu.header['TUNIT2'] = 'FLAM'
        tbhdu.header['TUNIT3'] = 'FLAM'
        tbhdu.header['TUNIT4'] = 'FLAM'
        tbhdu.header['TUNIT5'] = 'FLAM'
        tbhdu.header['TUNIT6'] = 'FLAM'
        tbhdu.header['TUNIT7'] = 'FLAM'
        tbhdu.header['TUNIT8'] = 'FLAM'
        tbhdu.header['TUNIT9'] = 'FLAM'
        tbhdu.header['TUNIT10'] = 'FLAM'
        tbhdu.header['TUNIT11'] = 'FLAM'
        tbhdu.header['TUNIT12'] = 'FLAM'
        tbhdu.header['TUNIT13'] = 'FLAM'
        tbhdu.header['TUNIT14'] = 'FLAM'

        # Construct and write out final FITS file
        finalhdu = fits.HDUList([prihdu, tbhdu])
        finalhdu.writeto(i, overwrite=True)

        hdu.close()
        print( 'Done {0:2.0f} of {1:2.0f}'.format(counter, len(files)))

    # Change back to starting directory
    os.chdir(start_dir)

    return

def rebin_phoenixV16(cdbs_path):
    """
    Rebin phoenixV16 models to atlas ck04 resolution; this makes
    spectrophotometry MUCH faster

    makes new directory in cdbs/grid: phoenix_v16_rebin

    cdbs_path: path to cdbs directory
    """
    # Get an atlas ck04 model, we will use this to set wavelength grid
    sp_atlas = get_castelli_atmosphere()

    # Open a fits table for an existing phoenix model; we will steal the header
    ## (This assumes that at least 'm00' metallicity exists)
    tmp = '{0}/grid/phoenix_v16/phoenix{1}/phoenix{1}_02400.fits'.format(cdbs_path, 'm00')
    phoenix_hdu = fits.open(tmp)
    header0 = phoenix_hdu[0].header

    # Create cdbs/grid directory for rebinned models
    path = cdbs_path+'/grid/phoenix_v16_rebin/'
    if not os.path.exists(path):
        os.mkdir(path)


    # Read in the existing catalog.fits file and rebin every spectrum.
    cat = fits.getdata(cdbs_path + '/grid/phoenix_v16/catalog.fits')
    files_all = [cat[ii][1].split('[')[0] for ii in range(len(cat))]
    temp_arr = np.zeros(len(files_all), dtype=float)
    logg_arr = np.zeros(len(files_all), dtype=float)
    metal_arr = np.zeros(len(files_all), dtype=float)

    for ff in range(len(files_all)):
        vals = cat[ff][0].split(',')

        temp_arr[ff] = float(vals[0])
        metal_arr[ff] = float(vals[1])
        logg_arr[ff] = float(vals[2])


    metal_uniq = np.unique(metal_arr)
    temp_uniq = np.unique(temp_arr)

    for mm in range(len(metal_uniq)):
        metal = metal_uniq[mm] # metallicity

        # Construct str for metallicity (for appropriate directory name)
        met_str = str(int(np.abs(metal))) + str(int((metal % 1.0)*10))
        if metal > 0:
            met_str = 'p' + met_str
        else:
            met_str = 'm' + met_str

        # Make directory for current metallicity if it does not exist yet
        if not os.path.exists(path + 'phoenix' + met_str):
            os.mkdir(path + 'phoenix' + met_str)

        for tt in range(len(temp_uniq)):
            temp = temp_uniq[tt] # temperature

            # Pick out the list of gravities for this T, Z combo
            idx = np.where((metal_arr == metal) & (temp_arr == temp))[0]
            logg_exist = logg_arr[idx]

            # All gravities will go in one file. Here is the output
            # file name.
            outfile = path + files_all[idx[0]].split('[')[0]

            ## If the rebinned file already exists, continue
            if os.path.exists(outfile):
                continue

            # Build a columns array. One column for each gravity.
            cols_arr = []

            # Make the wavelength column, which is first in the cols array.
            c0 = fits.Column(name='Wavelength', format='D', array=sp_atlas.waveset.to(u.AA).value)
            cols_arr.append(c0)

            for gg in range(len(logg_exist)):
                grav = logg_exist[gg] # gravity

                # Fetch the spectrum
                sp = stsyn.grid_to_spec('phoenix_v16', temp, metal, grav)
                flux_rebin = rebin_spec(sp.waveset, sp, sp_atlas.waveset)

                # Store the spectrum
                name = 'g{0:3.1f}'.format(grav)
                col = fits.Column(name=name, format='E', array=flux_rebin)
                cols_arr.append(col)


            # Make the FITS file from the columns with header.
            cols = fits.ColDefs(cols_arr)
            tbhdu = fits.BinTableHDU.from_columns(cols)
            prihdu = fits.PrimaryHDU(header=header0)
            tbhdu.header['TUNIT1'] = 'ANGSTROM'
            for gg in range(len(logg_exist)):
                tbhdu.header['TUNIT{0:d}'.format(gg+2)] = 'FLAM'

            # Write hdu
            finalhdu = fits.HDUList([prihdu, tbhdu])
            # don't have overwrite to protect original files.
            finalhdu.writeto(outfile)

            print( 'Finished file ' + outfile + ' with gravities: ', logg_exist)


    return

def organize_BTSettl_2015_atmospheres(path_to_dir):
    """
    Construct cdbs-ready BTSettl_CIFITS_2011_2015 atmospheres for each model.
    Will convert wavelength units to angstroms and flux units to [erg/s/cm^2/A]

    path_to_dir is the path to the directory containing all of the downloaded
    files

    Saves cdbs-ready atmospheres into os.environ['PYSYN_CDBS']/grid/BTSettl_2015
    (assumes this directory exists)
    """
    # Save current directory for return later, move into working dir
    start_dir = os.getcwd()
    os.chdir(path_to_dir)

    # If it doesn't already exist, create the BTSettl subdirectory
    if not os.path.exists('BTSettl_2015'):
        os.mkdir('BTSettl_2015')

    # Process each atmosphere file independently
    print( 'Creating cdbs-ready files')
    files = glob.glob('*.spec.fits')

    for i in files:
        hdu = fits.open(i)
        spec = hdu[1].data
        header_0 = hdu[0].header
        header_1 = hdu[1].header

        wave = spec.field(0)
        flux = spec.field(1)

        # Get units right: convert wave from microns to Angstroms,
        # flux from W /m^2/ micron to erg/s/cm^2/A
        wave_new = wave * 10**4
        flux_new = flux * 10**(-1)

        # Make new fits table
        c0 = fits.Column(name='Wavelength', format='D', array=wave_new)
        c1 = fits.Column(name='Flux', format='E', array=flux_new)

        cols = fits.ColDefs([c0, c1])
        tbhdu = fits.BinTableHDU.from_columns(cols)

        # Copy over headers, update unit keywords
        prihdu = fits.PrimaryHDU(header=header_0)
        tbhdu.header['TUNIT1'] = 'ANGSTROM'
        tbhdu.header['TUNIT2'] = 'FLAM'
        hdu_new = fits.HDUList([prihdu, tbhdu])

        # Write new fits table in cdbs directory
        hdu_new.writeto(os.environ['PYSYN_CDBS']+'grid/BTSettl_2015/'+i, overwrite=True)

        hdu.close()
        hdu_new.close()

    # Return to original directory
    os.chdir(start_dir)
    return

def make_BTSettl_2015_catalog(path_to_dir):
    """
    Create cdbs catalog.fits of BTSettl_CIFITS2011_2015 grid.
    THIS IS STEP 2, after organize_CMFGEN_atmospheres has
    been run.

    path_to_dir is from current working directory to the cdbs directory.
    Will create catalog.fits file in atmosphere directory with
    description of each model
    """
    # Record current working directory for later
    start_dir = os.getcwd()

    # Enter atmosphere directory
    os.chdir(path_to_dir)

    # Extract parameters for each atmosphere from the filename,
    # construct columns for catalog file
    files = glob.glob("*spec.fits")
    index_str = []
    name_str = []
    for name in files:
        tmp = name.split('-')
        temp = float(tmp[0][3:]) * 100.0 # In kelvin
        logg = float(tmp[1])

        index_str.append('{0:5.0f},0.0,{1:3.2f}'.format(temp, logg))
        name_str.append('{0}[Flux]'.format(name))

    # Make catalog
    catalog = Table([index_str, name_str], names = ('INDEX', 'FILENAME'))

    # Create catalog.fits file in directory with the models
    catalog.write('catalog.fits', format = 'fits', overwrite=True)

    # Move back to original directory, create the catalog.fits file
    os.chdir(start_dir)

    return

def rebin_BTSettl_2015(cdbs_path=None):
    """
    Rebin BTSettle_CIFITS2011_2015 models to atlas ck04 resolution; this makes
    spectrophotometry MUCH faster

    makes new directory in cdbs/grid: BTSettl_2015_rebin

    cdbs_path: path to cdbs directory. Defaults to the PYSYN_CDBS
    environment variable.
    """
    # Resolve the default here rather than in the signature, so that merely
    # importing this module does not require PYSYN_CDBS to be set
    if cdbs_path is None:
        cdbs_path = os.environ['PYSYN_CDBS']

    # Get an atlas ck04 model, we will use this to set wavelength grid
    sp_atlas = get_castelli_atmosphere()

    # Open a fits table for an existing phoenix model; we will steal the header
    tmp = cdbs_path+'/grid/phoenix_v16/phoenixm00/phoenixm00_02400.fits'
    phoenix_hdu = fits.open(tmp)
    header0 = phoenix_hdu[0].header
    phoenix_hdu.close()

    # Create cdbs/grid directory for rebinned models
    path = cdbs_path+'/grid/BTSettl_2015_rebin/'
    if not os.path.exists(path):
        os.mkdir(path)

    # Read in the existing catalog.fits file and rebin every spectrum.
    cat = fits.getdata(cdbs_path + '/grid/BTSettl_2015/catalog.fits')
    files_all = [cat[ii][1].split('[')[0] for ii in range(len(cat))]

    print( 'Rebinning BTSettl spectra')
    for ff in range(len(files_all)):
        vals = cat[ff][0].split(',')
        temp = float(vals[0])
        metal = float(vals[1])
        logg = float(vals[2])

        # Fetch the BTSettl spectrum, rebin flux
        sp = stsyn.grid_to_spec('BTSettl_2015', temp, metal, logg)
        flux_rebin = rebin_spec(sp.waveset, sp, sp_atlas.waveset)

        # Make new output
        c0 = fits.Column(name='Wavelength', format='D', array=sp_atlas.waveset.to(u.AA).value)
        c1 = fits.Column(name='Flux', format='E', array=flux_rebin)

        cols = fits.ColDefs([c0, c1])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        prihdu = fits.PrimaryHDU(header=header0)
        tbhdu.header['TUNIT1'] = 'ANGSTROM'
        tbhdu.header['TUNIT2'] = 'FLAM'

        outfile = path + files_all[ff].split('[')[0]
        finalhdu = fits.HDUList([prihdu, tbhdu])
        finalhdu.writeto(outfile, overwrite=True)

    return

def make_wavelength_unique(files, dirname):
    """
    Helper function to go through each BTSettl spectrum and ensure that
    each wavelength point is unique. This is required for rebinning to work.


    files: list of files to run this analysis on
    """
    # Loop through each file, find fix repeated wavelength entries if necessary
    for i in files:
        t = Table.read('{0}/{1}'.format(dirname,i), format='fits')
        test = np.unique(t['Wavelength'], return_index=True)

        if len(t) != len(test[0]):
            t = t[test[1]]

            c0 = fits.Column(name='Wavelength', format='D', array=t['Wavelength'])
            c1 = fits.Column(name='Flux', format='E', array=t['Flux'])
            cols = fits.ColDefs([c0, c1])

            tbhdu = fits.BinTableHDU.from_columns(cols)
            prihdu = fits.PrimaryHDU()
            tbhdu.header['TUNIT1'] = 'ANGSTROM'
            tbhdu.header['TUNIT2'] = 'FLAM'
            finalhdu = fits.HDUList([prihdu, tbhdu])
            finalhdu.writeto('{0}/{1}'.format(dirname,i), overwrite=True)

        # Also make sure wavelength is monotonic. If it is not, then it is
        # a sign that the wavelengths are out of order
        diff = np.diff(t['Wavelength'])
        bad = np.where(diff < 0)
        if len(bad[0]) > 0:
            t.sort('Wavelength')

            c0 = fits.Column(name='Wavelength', format='D', array=t['Wavelength'])
            c1 = fits.Column(name='Flux', format='E', array=t['Flux'])
            cols = fits.ColDefs([c0, c1])

            tbhdu = fits.BinTableHDU.from_columns(cols)
            prihdu = fits.PrimaryHDU()
            tbhdu.header['TUNIT1'] = 'ANGSTROM'
            tbhdu.header['TUNIT2'] = 'FLAM'
            finalhdu = fits.HDUList([prihdu, tbhdu])
            finalhdu.writeto('{0}/{1}'.format(dirname,i), overwrite=True)

        print('Done {0}'.format(i))

    return

def organize_BTSettl_atmospheres():
    """
    Construct cdbs-ready atmospheres for the BTSettl grid (CIFITS2011).
    The code expects tp be run in cdbs/grid/BTSettl, and expects that the
    individual model files have been downloaded from online
    (https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011/SPECTRA/)
    and processed into python-readable ascii files.
    """
    orig_dir = os.getcwd()
    dirs = ['btm25', 'btm20', 'btm15', 'btm10', 'btm05', 'btp00', 'btp05']
    #dirs = ['btm10', 'btm05', 'btp00', 'btp05']


    # Go through each directory, turning each spectrum into a cdbs-ready file.
    # Will convert flux into Ergs/sec/cm**2/A (FLAM) units and save as a fits file,
    # for faster access later
    for ii in dirs:
        print('Starting {0}'.format(ii))
        os.chdir(ii)

        files = glob.glob('*.txt')
        count=0
        for jj in files:
            t = Table.read(jj, format='ascii')
            # First, trim the wavelengths to a more reasonable wavelength range
            good = np.where( (t['col1'] > 1000) & (t['col1']  < 70000) )
            t = t[good]

            # Convert flux units to Flam (Ergs/sec/cm**2/A)
            flux_new = 10**(t['col2'] - 8.0)

            # Save the file as a fits file
            c0 = fits.Column(name='Wavelength', format='D', array=t['col1'])
            c1 = fits.Column(name='Flux', format='E', array=flux_new)

            cols = fits.ColDefs([c0, c1])
            tbhdu = fits.BinTableHDU.from_columns(cols)

            # Add unit keywords
            prihdu = fits.PrimaryHDU()
            tbhdu.header['TUNIT1'] = 'ANGSTROM'
            tbhdu.header['TUNIT2'] = 'FLAM'
            hdu_new = fits.HDUList([prihdu, tbhdu])

            # Write new fits table in cdbs directory
            hdu_new.writeto('{0}.fits'.format(jj[:-4]), overwrite=True)
            hdu_new.close()
            count += 1
            print('Done {0} of {1}'.format(count, len(files)))

        # Now, clean up all the files made when unzipping the spectra
        cmd1 = 'rm *.bz2'
        cmd2 = 'rm *.tmp'
        #cmd3 = 'rm *.txt'
        os.system(cmd1)
        os.system(cmd2)
        #os.system(cmd3)
        print('==============================')
        print('Done {0}'.format(ii))
        print('==============================')

        # Go back to original directory, move to next metallicity directory
        os.chdir(orig_dir)

    return

def make_BTSettl_catalog():
    """
    Create cdbs catalog.fits of BTSettl grid.
    THIS IS STEP 2, after organize_BTSettl_atmospheres has
    been run.

    Code expects to be run in cdbs/grid/BTSettl
    Will create catalog.fits file in atmosphere directory with
    description of each model
    """
    # Record current working directory for later
    start_dir = os.getcwd()
    dirs = ['btm25', 'btm20', 'btm15', 'btm10', 'btm05', 'btp00', 'btp05']
    #dirs = ['btp05']

    # Construct the catalog.fits file input. The input consists of
    # and index string that specifies the stellar paramters, and a
    # name string that points to the file
    # Loop over all the metallicity directories to construct these inputs
    index_str = []
    name_str = []
    for ii in dirs:
        os.chdir(ii)
        files = glob.glob('*.fits')

        # Construct the metallicity val
        if 'm' in ii:
            metal_flag = -1 * float(ii[3:])*0.1
        else:
            metal_flag = float(ii[3:])*0.1

        # Now collect the info from the files
        for jj in files:
            tmp = jj.split('-')

            if metal_flag >= 0:
                temp = float(tmp[0].split('+')[0][3:]) * 100.0 # In kelvin
                try:
                    logg = float(tmp[1])
                except:
                    logg = float(tmp[1].split('+')[0])
            else:
                temp = float(tmp[0][3:]) * 100.0 # In kelvin
                logg = float(tmp[1])

            index_str.append('{0},{1},{2:3.2f}'.format(int(temp), metal_flag, logg))
            name_str.append('{0}/{1}[Flux]'.format(ii, jj))

        # Go back to original directory to move to next metallicity
        print('Done {0}'.format(ii))
        os.chdir(start_dir)

    # Make catalog
    catalog = Table([index_str, name_str], names = ('INDEX', 'FILENAME'))

    # Create catalog.fits file in directory with the models
    catalog.write('catalog.fits', format = 'fits', overwrite=True)

    # Move back to original directory, create the catalog.fits file
    os.chdir(start_dir)

    return

def rebin_BTSettl(make_unique=False):
    """
    Rebin BTSettle models to atlas ck04 resolution; this makes
    spectrophotometry MUCH faster

    makes new directory: BTSettl_rebin

    Code expects to be run in cdbs/grid directory
    """
    # Get an atlas ck04 model, we will use this to set wavelength grid
    sp_atlas = get_castelli_atmosphere()

    # Create cdbs/grid directory for rebinned models
    path = 'BTSettl_rebin/'
    if not os.path.exists(path):
        os.mkdir(path)

    # Read in the existing catalog.fits file and rebin every spectrum.
    cat = fits.getdata('BTSettl/catalog.fits')
    files_all = [cat[ii][1].split('[')[0] for ii in range(len(cat))]

    #==============================#
    #tmp = []
    #for ii in files_all:
    #    if ii.startswith('btp00'):
    #        tmp.append(ii)
    #files_all = tmp
    #=============================#

    print( 'Rebinning BTSettl spectra')
    if make_unique:
        print('Making unique')
        make_wavelength_unique(files_all, 'BTSettl')
        print('Done')

    for ff in range(len(files_all)):
        vals = cat[ff][0].split(',')
        temp = float(vals[0])
        metal = float(vals[1])
        logg = float(vals[2])

        # Fetch the BTSettl spectrum, rebin flux
        try:
            sp = stsyn.grid_to_spec('BTSettl', temp, metal, logg)
            flux_rebin = rebin_spec(sp.waveset, sp, sp_atlas.waveset)

            # Make new output
            c0 = fits.Column(name='Wavelength', format='D', array=sp_atlas.waveset.to(u.AA).value)
            c1 = fits.Column(name='Flux', format='E', array=flux_rebin)

            cols = fits.ColDefs([c0, c1])
            tbhdu = fits.BinTableHDU.from_columns(cols)
            prihdu = fits.PrimaryHDU()
            tbhdu.header['TUNIT1'] = 'ANGSTROM'
            tbhdu.header['TUNIT2'] = 'FLAM'

            outfile = path + files_all[ff].split('[')[0]
            finalhdu = fits.HDUList([prihdu, tbhdu])
            finalhdu.writeto(outfile, overwrite=True)
        except:
            pdb.set_trace()
            orig_file = '{0}/{1}'.format('BTSettl/', files_all[ff].split('[')[0])
            outfile = path + files_all[ff].split('[')[0]
            cmd = 'cp {0} {1}'.format(orig_file, outfile)
            os.system(cmd)

        print('Done {0} of {1}'.format(ff, len(files_all)))

    return

def organize_all_Meisner2023_atmospheres():
    """
    Construct cdbs-ready atmospheres for the Meisner2023 grid.
    The code expects tp be run in cdbs/grid/Meisner2023, and expects that the
    individual model files have been downloaded from online
    and processed into python-readable ascii files.
    """
    orig_dir = os.getcwd()
    dirs = ['mm10', 'mm05', 'mp00', 'mp03']

    # Go through each directory, turning each spectrum into a cdbs-ready file.
    # Save as a fits file, for faster access later
    for ii in dirs:
        print('Starting {0}'.format(ii))
        os.chdir(ii)

        files = glob.glob('*.fits')
        count=0
        for jj in files:
            # Open each .fits file and read the data
            with fits.open(jj) as hdul:
                data = hdul[1].data
                wavelength = data['Wavelength']
                flux = data['Flux']

            # Make flux independent of R&D
            flux_new = flux / 5e-20

            # Create new columns with desired format
            c0 = fits.Column(name='Wavelength', format='D', array=wavelength)
            c1 = fits.Column(name='Flux', format='E', array=flux_new)

            cols = fits.ColDefs([c0, c1])
            tbhdu = fits.BinTableHDU.from_columns(cols)

            # Add unit keywords
            prihdu = fits.PrimaryHDU()
            tbhdu.header['TUNIT1'] = 'ANGSTROM'
            tbhdu.header['TUNIT2'] = 'FLAM'
            hdu_new = fits.HDUList([prihdu, tbhdu])

            # Write the new fits table in the cdbs directory
            output_filename = '{0}.fits'.format(jj[:-5])  # Removing the original .fits extension
            hdu_new.writeto(output_filename, overwrite=True)
            hdu_new.close()
            count += 1
            print('Done {0} of {1}'.format(count, len(files)))

        # Go back to original directory, move to next metallicity directory
        os.chdir(orig_dir)

    return

def make_Meisner2023_catalog():
    """
    Create cdbs catalog.fits of Meisner2023 grid.
    THIS IS STEP 2, after organize_Meisner2023_atmospheres has
    been run.

    Code expects to be run in cdbs/grid/Meisner2023
    Will create catalog.fits file in atmosphere directory with
    description of each model
    """
    # Record current working directory for later
    start_dir = os.getcwd()
    dirs = ['mm10', 'mm05', 'mp00', 'mp03']

    # Construct the catalog.fits file input. The input consists of
    # and index string that specifies the stellar paramters, and a
    # name string that points to the file
    # Loop over all the metallicity directories to construct these inputs
    index_str = []
    name_str = []
    for ii in dirs:
        os.chdir(ii)
        files = glob.glob('spec_jwst_*.fits')

        for jj in files:
            # Parse temperature, log(g), and metallicity from filename
            temp_str = jj.split('_')[2]
            logg_str = jj.split('_')[3]
            metal_str = jj.split('_')[4]

            # Extract temperature, surface gravity, and metallicity
            temp = float(temp_str[1:])         # Temperature in Kelvin
            logg = float(logg_str[1:])         # Surface gravity log(g)

            # Build metallicity value
            if metal_str.startswith('m'):
                metallicity = -1 * float(metal_str[1:])
            else:
                metallicity = float(metal_str[1:])

            # Construct index and filename strings
            index_str.append('{0},{1},{2:3.2f}'.format(int(temp), metallicity, logg))
            name_str.append('{0}/{1}[Flux]'.format(ii, jj))

        print('Processed directory:', ii)
        os.chdir(start_dir)


    # Make catalog
    catalog = Table([index_str, name_str], names = ('INDEX', 'FILENAME'))

    # Create catalog.fits file in directory with the models
    catalog.write('catalog.fits', format = 'fits', overwrite=True)

    # Move back to original directory, create the catalog.fits file
    os.chdir(start_dir)

    return

def rebin_Meisner2023(make_unique=False):
    """
    Rebin Meisner2023 models to atlas ck04 resolution; this makes
    spectrophotometry MUCH faster

    makes new directory: Meisner2023_rebin

    Code expects to be run in cdbs/grid directory
    """
    # Get an atlas ck04 model, we will use this to set wavelength grid
    sp_atlas = get_castelli_atmosphere()

    # Create a directory for rebinned Meisner2023 models
    rebin_path = 'Meisner2023_rebin/'
    if not os.path.exists(rebin_path):
        os.mkdir(rebin_path)

    # Load the catalog.fits file and extract all spectra file paths
    cat = Table.read('Meisner2023/catalog.fits')
    files_all = [cat[ii]['FILENAME'].split('[')[0] for ii in range(len(cat))]

    print('Rebinning Meisner2023 spectra')
    if make_unique:
        print('Making unique')
        make_wavelength_unique(files_all, 'Meisner2023')
        print('Done')

    for ff, file in enumerate(files_all):
        vals = cat[ff]['INDEX'].split(',')
        temp = float(vals[0])
        metal = float(vals[1])
        logg = float(vals[2])

        # Fetch the Meisner2023 spectrum and rebin its flux
        try:
            sp = stsyn.grid_to_spec('Meisner2023', temp, metal, logg)
            flux_rebin = rebin_spec(sp.waveset, sp, sp_atlas.waveset)

            # Create the output FITS file
            c0 = fits.Column(name='Wavelength', format='D', array=sp_atlas.waveset.to(u.AA).value)
            c1 = fits.Column(name='Flux', format='E', array=flux_rebin)

            cols = fits.ColDefs([c0, c1])
            tbhdu = fits.BinTableHDU.from_columns(cols)
            prihdu = fits.PrimaryHDU()
            tbhdu.header['TUNIT1'] = 'ANGSTROM'
            tbhdu.header['TUNIT2'] = 'FLAM'

            # Write the new rebinned file in the Meisner2023_rebin directory
            outfile = os.path.join(rebin_path, os.path.basename(file))
            finalhdu = fits.HDUList([prihdu, tbhdu])
            finalhdu.writeto(outfile, overwrite=True)

        except Exception as e:
            print(f"Error processing {file}: {e}")
            orig_file = os.path.join('Meisner2023', file)
            outfile = os.path.join(rebin_path, os.path.basename(file))
            os.system(f'cp {orig_file} {outfile}')

        print('Done {0} of {1}'.format(ff + 1, len(files_all)))

    return



def organize_WDKoester_atmospheres(path_to_dir):
    """
    Construct cdbs-ready wdKoester WD atmospheres for each model. (from Koester 2010)
    Will convert wavelength units to angstroms and flux units to [erg/s/cm^2/A]

    path_to_dir is the path to the directory containing all of the downloaded
    files

    Saves cdbs-ready atmospheres into os.environ['PYSYN_CDBS']/wdKoeseter
    (assumes this directory exists)
    """
    # Save current directory for return later, move into working dir
    start_dir = os.getcwd()
    os.chdir(path_to_dir)

    # Process each atmosphere file independently
    print( 'Creating cdbs-ready files')
    files = glob.glob('*.dk.dat.txt')

    for i in files:
        data = Table.read(i, format='ascii')

        wave = data['col1']   # angstrom
        flux = data['col2']   # erg/s/cm^2/A

        # Make new fits table
        c0 = fits.Column(name='Wavelength', format='D', array=wave)
        c1 = fits.Column(name='Flux', format='E', array=flux)

        cols = fits.ColDefs([c0, c1])
        tbhdu = fits.BinTableHDU.from_columns(cols)

        # Copy over headers, update unit keywords
        prihdu = fits.PrimaryHDU()
        tbhdu.header['TUNIT1'] = 'ANGSTROM'
        tbhdu.header['TUNIT2'] = 'FLAM'
        hdu_new = fits.HDUList([prihdu, tbhdu])

        # Write new fits table in cdbs directory
        hdu_new.writeto(os.environ['PYSYN_CDBS']+'/grid/wdKoester/'+i.replace('.txt', '.fits'), overwrite=True)

        hdu_new.close()

    # Return to original directory
    os.chdir(start_dir)
    return

def make_WDKoester_catalog(path_to_dir):
    """
    Create cdbs catalog.fits of wdKoester grid.
    THIS IS STEP 2, after organize_WDKoester_atmospheres has
    been run.

    path_to_dir is from current working directory to the cdbs directory.
    Will create catalog.fits file in atmosphere directory with
    description of each model
    """
    # Record current working directory for later
    start_dir = os.getcwd()

    # Enter atmosphere directory
    os.chdir(path_to_dir)

    # Extract parameters for each atmosphere from the filename,
    # construct columns for catalog file
    files = glob.glob("*dk.dat.fits")
    index_str = []
    name_str = []
    for name in files:
        tmp = name.split('.')
        tmp2 = tmp[0].split('_')
        temp = float(tmp2[0][2:]) # Kelvin
        logg = float(tmp2[1]) / 100.0   # log(g)

        index_str.append('{0:5.0f},0.0,{1:3.2f}'.format(temp, logg))
        name_str.append('{0}[Flux]'.format(name))

    # Make catalog
    catalog = Table([index_str, name_str], names = ('INDEX', 'FILENAME'))

    # Create catalog.fits file in directory with the models
    catalog.write('catalog.fits', format = 'fits', overwrite=True)

    # Move back to original directory, create the catalog.fits file
    os.chdir(start_dir)

    return

def rebin_WDKoester(cdbs_path=None):
    """
    Rebin wdKoester models to atlas ck04 resolution; this makes
    spectrophotometry MUCH faster

    makes new directory in cdbs/grid: wdKoester_rebin

    cdbs_path: path to cdbs directory. Defaults to the PYSYN_CDBS
    environment variable.
    """
    # Resolve the default here rather than in the signature, so that merely
    # importing this module does not require PYSYN_CDBS to be set
    if cdbs_path is None:
        cdbs_path = os.environ['PYSYN_CDBS']

    # Get an atlas ck04 model, we will use this to set wavelength grid
    sp_atlas = get_castelli_atmosphere()

    # Open a fits table for an existing model; we will steal the header
    tmp = cdbs_path+'/grid/wdKoester/da70000_800.dk.dat.fits'
    wdkoester_hdu = fits.open(tmp)
    header0 = wdkoester_hdu[0].header
    wdkoester_hdu.close()

    # Create cdbs/grid directory for rebinned models
    path = cdbs_path+'/grid/wdKoester_rebin/'
    if not os.path.exists(path):
        os.mkdir(path)

    # Read in the existing catalog.fits file and rebin every spectrum.
    cat = fits.getdata(cdbs_path + '/grid/wdKoester/catalog.fits')
    files_all = [cat[ii][1].split('[')[0] for ii in range(len(cat))]

    print( 'Rebinning wdKoester spectra')
    for ff in range(len(files_all)):
        vals = cat[ff][0].split(',')
        temp = float(vals[0])
        metal = float(vals[1])
        logg = float(vals[2])

        # Fetch the wdKoester spectrum, rebin flux
        sp = stsyn.grid_to_spec('wdKoester', temp, metal, logg)
        flux_rebin = rebin_spec(sp.waveset, sp, sp_atlas.waveset)

        # Make new output
        c0 = fits.Column(name='Wavelength', format='D', array=sp_atlas.waveset.to(u.AA).value)
        c1 = fits.Column(name='Flux', format='E', array=flux_rebin)

        cols = fits.ColDefs([c0, c1])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        prihdu = fits.PrimaryHDU(header=header0)
        tbhdu.header['TUNIT1'] = 'ANGSTROM'
        tbhdu.header['TUNIT2'] = 'FLAM'

        outfile = path + files_all[ff].split('[')[0]
        finalhdu = fits.HDUList([prihdu, tbhdu])
        finalhdu.writeto(outfile, overwrite=True)

    return


