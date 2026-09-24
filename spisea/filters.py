import numpy as np
from astropy.table import Table
import warnings
from astropy import units as u

from synphot import units as su
from synphot.models import Empirical1D
from synphot.spectrum import SpectralElement

import os
import pdb

# Set path to filter functions
code_dir = os.path.dirname(__file__)
filters_dir = code_dir[:-7]+'/filt_func/'


def get_nirc2_filt(name):
    """
    Define nirc2 filter as a synphot spectrum object
    """
    # Read in filter info
    try:
        t = Table.read('{0}/nirc2/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        pdb.set_trace()
        raise ValueError('Could not find NIRC2 filter file {0}/nirc2/{1}.dat'.format(filters_dir, name))

    wavelength = t[t.keys()[0]]
    transmission = t[t.keys()[1]]

    # Lets fix wavelength array for duplicate values
    diff = np.diff(wavelength)
    idx = np.where(diff <= 0)[0]

    while len(idx) != 0:
        wavelength[idx+1] += 1.0e-8

        diff = np.diff(wavelength)
        idx = np.where(diff <= 0)[0]
        #print( 'Duplicate entry loop' )

    # Get rid of all entries with negative transmission
    idx = np.where(transmission > 1)[0]

    # Convert wavelength to Angstroms, transmission to ratio
    wavelength = wavelength[idx] * 10**4 * u.AA
    transmission = transmission[idx] / 100.0 * su.THROUGHPUT  

    # Make spectrum object
    spectrum = SpectralElement(
        Empirical1D,
        points=wavelength,
        lookup_table=transmission,
        meta={"expr": "NIRC2_{0}".format(name)},
    )

    return spectrum


def get_2mass_filt(name):
    """
    Define the 2mass filters as a synphot spectrum object
    """
    # Read in filter info
    try:
        t = Table.read('{0}/2mass/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find 2MASS filter file {0}/2mass/{1}.dat'.format(filters_dir, name))

    wavelength = t[t.keys()[0]]
    transmission = t[t.keys()[1]]

    # Convert wavelength to Angstroms
    wavelength = wavelength * 10**4 * u.AA
    transmission = transmission * su.THROUGHPUT  

    # Make spectrum object
    spectrum = SpectralElement(Empirical1D, 
                               points=wavelength, 
                               lookup_table=transmission, 
                               meta={"expr": f"2MASS_{name}"})
    return spectrum


def get_vista_filt(name):
    """
    Define vista filter as synphot spectrum object
    """
    # Read in filter info
    try:
        t = Table.read('{0}/vista/VISTA_Filters_at80K_forETC_{1}.dat'.format(filters_dir, name),
                           format='ascii')
    except:
        raise ValueError('Could not find VISTA filter file {0}/vista/VISTA_Filters_at80K_forETC_{1}.dat'.format(filters_dir, name))

   # Wavelength must be in angstroms, transmission in fraction
    wave = t['col1'] * 10
    trans = t['col2'] * 0.01

    # Change any negative numbers to 0, as well as anything shortward
    # of 0.4 microns or longward of 2.9 microns
    # (no VISTA filter transmissions beyond these boundaries)
    bad = np.where( (trans < 0) | (wave < 4000) | (wave > 29000) )
    trans[bad] = 0

    # Now we can define the VISTA filter bandpass objects
    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "VISTA_{0}".format(name)},
    )

    return spectrum

def get_decam_filt(name):
    """
    Define DECAM filter as synphot object
    """
    # Read in filter info
    try:
        t = Table.read('{0}/decam/DECam_filters.txt'.format(filters_dir), format='ascii')

        trans = t[name]
    except:
        if name=='y':
            raise ValueError('DECam has a /"Y/" filter, not /"y/". The /"y/" in SPISEA <v3.0 was a bug.')
        else:
            raise ValueError('Could not find DECAM filter {0} in {1}/decam/DECam_filters.txt'.format(name, filters_dir))

    # Don't allow negative transmission
    trans[trans<0] = 0.0

    # Convert wavelengths from nm to angstroms, while eliminating masked regions
    wave = t['wavelength'] * 10.
    trans = trans

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "decam_{0}".format(name)},
    )

    return spectrum

def get_PS1_filt(name):
    """
    Define PS1 filter as synphot object
    """
    try:
        t = Table.read('{0}/ps1/PS1_filters.txt'.format(filters_dir), format='ascii')
        t.rename_column('col1', 'wave')
        t.rename_column('col2', 'open')
        t.rename_column('col3', 'g')
        t.rename_column('col4', 'r')
        t.rename_column('col5', 'i')
        t.rename_column('col6', 'z')
        t.rename_column('col7', 'y')
        t.rename_column('col8', 'w')

        cols = np.array(t.keys())
        idx = np.where(cols == name)[0][0]

        trans = t[cols[idx]]
    except:
        raise ValueError('Could not find PS1 filter {0} in {1}/ps1'.format(name, filters_dir))

    # Convert wavelengths from nm to angstroms
    wave = t['wave'] * 10.

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "ps1_{0}".format(name)},
    )

    return spectrum

def get_jwst_filt(name):
    """
    Define JWST filter as synphot object
    """
    try:
        t = Table.read('{0}/jwst/{1}.txt'.format(filters_dir, name.upper()), format='ascii')
    except:
        raise ValueError('Could not find JWST filter {0} in {1}/jwst'.format(name, filters_dir))

    # Convert wavelengths to angstroms
    wave = t['microns'] * 10**4. * u.AA
    trans = t['throughput'] * su.THROUGHPUT

    # Change any negative numbers to 0
    bad = np.where(trans < 0)
    trans[bad] = 0

    spectrum = SpectralElement(
        Empirical1D,
        points=wave,
        lookup_table=trans,
        meta={"expr": "jwst_{0}".format(name)},
    )

    return spectrum

def get_Johnson_Glass_filt(name):
    """
    Define Johnson-Glass filters as synphot object
    """
    try:
        t = Table.read('{0}/Johnson_Glass/{1}.txt'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find Johnson-Glass filter {0} in {1}/Johnson_Glass'.format(name, filters_dir))

    # Convert wavelengths to angstroms
    wave = t['col1'] * 10. * u.AA
    trans = t['col2'] * su.THROUGHPUT

    # Change any negative numbers to 0
    bad = np.where(trans < 0)
    trans[bad] = 0

    spectrum = SpectralElement(
        Empirical1D,
        points=wave,
        lookup_table=trans,
        meta={"expr": "jg_{0}".format(name)},
    )

    return spectrum

def get_nirc1_filt(name):
    """
    Define Keck/NIRC filters as synphot object
    """
    try:
        t = Table.read('{0}/nirc1/{1}.txt'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find NIRC1 filter {0} in {1}/nirc1'.format(name, filters_dir))

    # Convert wavelengths to angstroms
    wave = t['col1'] * 10**4
    trans = t['col2']

    # Lets fix wavelength array for duplicate values or negative vals;
    # delete these entries
    diff = np.diff(wave)
    idx = np.where(diff <= 0)[0]

    while(len(idx) != 0):
        bad = idx + 1

        wave = np.delete(wave, bad)
        trans = np.delete(trans, bad)

        diff = np.diff(wave)
        idx = np.where(diff <= 0)[0]

    # Change any negative transmission vals to 0
    bad = np.where(trans < 0)
    trans[bad] = 0

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "nirc1_{0}".format(name)},
    )

    return spectrum

def get_ctio_osiris_filt(name):
    """
    Define CTIO/OSIRIS filters as synphot object
    """
    try:
        t = Table.read('{0}/CTIO_OSIRIS/{1}.txt'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find CTIO/OSIRIS filter {0} in {1}/CTIO_OSIRIS'.format(name, filters_dir))

    # Convert wavelengths to angstroms
    wave = t['col1'] * 10**4
    trans = t['col2']

    # Change any negative numbers to 0
    bad = np.where(trans < 0)
    trans[bad] = 0

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "ctio_osiris_{0}".format(name)},
    )

    return spectrum

def get_naco_filt(name):
    """
    Define VLT NACO filters as synphot object
    """
    try:
        t = Table.read('{0}/naco/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find NACO filter {0} in {1}/naco'.format(name, filters_dir))

    # Convert wavelengths to angstroms
    wave = t['col1'] * 10**4
    trans = t['col2']

    # Change any negative numbers to 0
    bad = np.where(trans < 0)
    trans[bad] = 0

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "naco_{0}".format(name)},
    )

    return spectrum

def get_ubv_filt(name):
    """
    Define ubv (Johnson-Cousin filters) as synphot object
    """
    try:
        t = Table.read('{0}/ubv/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find ubv filter {0} in {1}/ubv'.format(name, filters_dir))

    # Convert wavelength from nm to angstroms
    wave = t[t.keys()[0]] * 10.
    # Convert transmission to ratio (from percent)
    trans = t[t.keys()[1]] / 100.

    # Change any negative numbers to 0
    bad = np.where(trans < 0)
    trans[bad] = 0

    if name=='I':
        warnings.warn("Filter profile ubv,I uses the Johnson filter which extends further to long wavelengths than Cousins. "
            "Here, it is improperly cut off at 1.1 microns where transmission is ~20%. "
            "Consider using the Bessell UBVRI (bessell,I) filter system instead.")

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "ubv_{0}".format(name)},
    )

    return spectrum

def get_bessell_filt(name):
    """
    Define ubv (Johnson-Cousin filters, as defined in Bessell 1990) as synphot object
    """
    try:
        t = Table.read('{0}/bessell/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find bessell filter {0} in {1}/bessell'.format(name, filters_dir))

    wave = t[t.keys()[0]] 
    trans = t[t.keys()[1]]

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "bessell_{0}".format(name)},
    )

    return spectrum

def get_ukirt_filt(name):
    """
    Define UKIRT filters as synphot object
    """
    try:
        t = Table.read('{0}/ukirt/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find ukirt filter {0} in {1}/ukirt'.format(name, filters_dir))

    # Convert wavelengths to angstroms (from microns)
    wave = t[t.keys()[0]] * 10000.
    trans = t[t.keys()[1]]

    # Change any negative numbers to 0
    bad = np.where(trans < 0)
    trans[bad] = 0

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "ukirt_{0}".format(name)},
    )

    return spectrum

def get_keck_osiris_filt(name):
    """
    Define keck osiris filters as synphot object
    """
    try:
        t = Table.read('{0}/keck_osiris/{1}.txt'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find keck_osiris filter {0} in {1}/keck_osiris'.format(name, filters_dir))

    # Convert wavelengths to angstroms (from nm), percentage throughput to fraction
    wave = t['col1'] * 10
    trans = t['col2'] / 100.

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "keck_osiris_{0}".format(name)},
    )

    return spectrum

def get_gaia_filt(version, name):
    """
    Define Gaia filters as synphot object.

    version: specify dr1, dr2, dr2_rev, or edr3
    name: filter name
    """
    # Warn if not using latest version
    if version != 'edr3':
        msg = 'Gaia version {0} not recommended, use edr3 for the latest version'.format(version)
        warnings.warn(msg)

    # Set the filter directory
    if version == 'dr1':
        path = '{0}/gaia/dr1/'.format(filters_dir)
    elif version == 'dr2':
        path = '{0}/gaia/dr2/'.format(filters_dir)
    elif version == 'dr2_rev':
        path = '{0}/gaia/dr2_rev/'.format(filters_dir)
    elif version == 'edr3':
        path = '{0}/gaia/edr3/'.format(filters_dir)
    else:
        raise ValueError('GAIA filter version {0} not understood. Please use dr1, dr2, dr2_rev, or edr3'.format(version))

    # Get the filter info
    try:
        t = Table.read('{0}/Gaia_passbands.txt'.format(path), format='ascii')
        if version == 'dr1':
            t.rename_column('BP', 'Gbp')
            t.rename_column('RP', 'Grp')
        else:
            t.rename_column('col1', 'LAMBDA')
            t.rename_column('col2', 'G')
            t.rename_column('col4', 'Gbp')
            t.rename_column('col6', 'Grp')

        trans = t[name]

        # Change 99 values where filters are undefined into 0, to ensure that
        # it doesn't mess up our flux values
        bad = np.where(trans > 90)
        trans[bad] = 0
    except:
        raise ValueError('Could not find Gaia filter {0} for version {1}'.format(name, version))

    # Convert wavelengths to angstroms (from nm)
    wave = t['LAMBDA'] * 10

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "gaia_{0}_{1}".format(version, name)},
    )

    return spectrum

def get_ztf_filt(name):
    """
    Define ztf filters as synphot object
    """
    try:
        t = Table.read('{0}/ztf/{1}.dat'.format(filters_dir, name),
                       format='ascii')
    except:
        raise ValueError('Could not find ztf filter {0} in {1}/ztf'.format(name, filters_dir))

    wave = t['Wavelength']
    trans = t['Transmission']

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "ztf_{0}".format(name)},
    )

    return spectrum

def get_hawki_filt(name):
    """
    Define the HAWK-I filters as a synphot spectrum object
    """
    # Read in filter info
    try:
        t = Table.read('{0}/hawki/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find HAWK-I filter file {0}/hawki/{1}.dat'.format(filters_dir, name))
    #pdb.set_trace()
    wavelength = t[t.keys()[0]]
    transmission = t[t.keys()[1]]

    # Convert wavelength to Angstroms
    wavelength = wavelength * 10

    # Make spectrum object
    spectrum = SpectralElement(
        Empirical1D,
        points=wavelength * u.AA,
        lookup_table=transmission * su.THROUGHPUT,
        meta={"expr": "hawki_{0}".format(name)},
    )

    return spectrum

def get_rubin_filt(name):
    """
    Define the Rubin Vera C LSST filters as a synphot spectrum object
    """
    # Read in filter info
    try:
        t = Table.read('{0}/rubin/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find Rubin LSST filter file {0}/rubin/{1}.dat'.format(filters_dir, name))

    wavelength = t[t.keys()[0]]
    transmission = t[t.keys()[1]]

    # Convert wavelength to Angstroms
    wavelength = wavelength * 10

    # Make spectrum object
    spectrum = SpectralElement(
        Empirical1D,
        points=wavelength * u.AA,
        lookup_table=transmission * su.THROUGHPUT,
        meta={"expr": "rubin_{0}".format(name)},
    )

    return spectrum

def get_euclid_filt(name):
    """
    Define the Euclid filters as a synphot spectrum object
    """
    # Read in filter info
    try:
        t = Table.read('{0}/euclid/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find Euclid filter file {0}/euclid/{1}.dat'.format(filters_dir, name))

    wavelength = t[t.keys()[0]]
    transmission = t[t.keys()[1]]

    # Convert wavelength to Angstroms
    if name.lower() != 'vis':
        wavelength = wavelength * 10

    # Make spectrum object
    spectrum = SpectralElement(
        Empirical1D,
        points=wavelength * u.AA,
        lookup_table=transmission * su.THROUGHPUT,
        meta={"expr": "euclid_{0}".format(name)},
    )

    return spectrum

def get_nsfcam_filt(name):
    """
    Define irtf nsfcam filters as synphot object
    """
    try:
        t = Table.read('{0}/nsfcam/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find nsfcam filter {0} in {1}/nsfcam'.format(name, filters_dir))

    # Wavelength already in angstrom and and transmission in fraction
    wave = t['col1']
    trans = t['col2']

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "nsfcam_{0}".format(name)},
    )

    return spectrum

def get_tess_filt(name):
    """
    Define the TESS filter as synphot object
    """
    try:
        t = Table.read('{0}/tess/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find tess filter {0} in {1}/tess'.format(name, filters_dir))

    # Wavelength from nanometers to angstroms and and transmission in fraction
    wave = t['col1']*10
    trans = t['col2']

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "tess,{0}".format(name)},
    )

    return spectrum

def get_washington_filt(name):
    """
    Define the Washington filters as synphot object
    """
    try:
        t = Table.read('{0}/washington/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find washington filter {0} in {1}/washington'.format(name, filters_dir))

    # Wavelength from nanometers to angstroms and and transmission in fraction
    wave = t['col1']*10
    trans = t['col2']

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "washington,{0}".format(name)},
    )

    return spectrum

def get_hipparcos_filt(name):
    """
    Define the Hipparcos filter as synphot object
    """
    try:
        t = Table.read('{0}/hipparcos/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find hipparcos filter {0} in {1}/hipparcos'.format(name, filters_dir))

    # Wavelength in angstroms and and transmission in fraction
    wave = t['col1']
    trans = t['col2']

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "hipparcos,{0}".format(name)},
    )

    return spectrum

def get_tycho_filt(name):
    """
    Define the Tycho filters as synphot object
    """
    try:
        t = Table.read('{0}/tycho/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find tycho filter {0} in {1}/tycho'.format(name, filters_dir))

    # Wavelength in angstroms and and transmission in fraction
    wave = t['col1']
    trans = t['col2']

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "tycho,{0}".format(name)},
    )

    return spectrum

def get_kepler_filt(name):
    """
    Define the Kepler filters as synphot object
    """
    try:
        t = Table.read('{0}/kepler/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find kepler filter {0} in {1}/kepler'.format(name, filters_dir))

    # Wavelength in angstroms and and transmission in fraction
    wave = t['col1']
    trans = t['col2']

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "kepler,{0}".format(name)},
    )

    return spectrum

def get_ogle_filt(name):
    """
    Define the OGLE filters as synphot object
    """
    try:
        t = Table.read('{0}/ogle/{1}.dat'.format(filters_dir, name), format='ascii')
    except:
        raise ValueError('Could not find ogle filter {0} in {1}/ogle'.format(name, filters_dir))

    # Wavelength in nm->angstroms and and transmission in percent->fraction
    wave = np.flip(t['col1'])*10
    trans = np.flip(t['col2'])/100

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "ogle,{0}".format(name)},
    )

    return spectrum

def get_subaru_filt(instrument, name):
    """
    Define the subaru filters as synphot object
    """
    try:
        t = Table.read('{0}/subaru/{1}/{2}.dat'.format(filters_dir, instrument, name), format='ascii')
    except:
        raise ValueError('Could not find Subaru filter {0} in {1}/subaru/{2}'.format(name, filters_dir, instrument))

    # Wavelength in nm->angstroms and and transmission in percent->fraction
    wave = t['col1']
    trans = t['col2']

    spectrum = SpectralElement(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=trans * su.THROUGHPUT,
        meta={"expr": "subaru,{0},{1}".format(instrument, name)},
    )

    return spectrum
