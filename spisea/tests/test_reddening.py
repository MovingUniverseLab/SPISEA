import numpy as np
import spisea
from spisea import reddening, synthetic, evolution, atmospheres
import pylab as py
import os
import pdb
from astropy import units as u
import pytest

spisea_path = os.path.dirname(spisea.__file__)


def test_RedLawBrokenPowerLaw(plots=False):
    #===============================#
    # Test 1: 2-segment law
    #===============================#
    alpha1 = 2.23
    alpha2 = 3.0
    lambda_limits = [1.27, 1.63, 2.16]
    alpha_vals = [alpha2, alpha1]
    K_wave = 2.14

    red_law = reddening.RedLawBrokenPowerLaw(lambda_limits, alpha_vals, K_wave)

    # Calculate what the redlaw *should* be across wavelength range,
    # by-hand
    wave_test = np.arange(1.27, 2.15+0.001, 0.01)
    law_test = np.ones(len(wave_test)) * np.nan

    # 2.14 - 1.63
    idx = np.where( (wave_test > 1.63) & (wave_test <= 2.17))
    coeff = 1
    law_test[idx] = coeff * wave_test[idx] ** (-1*alpha1)

    # 1.63 - 1.27
    idx = np.where( (wave_test >= 1.27) & (wave_test <= 1.63))
    coeff  = (1.63 ** (-1*alpha1)) / (1.63 ** (-1*alpha2))
    law_test[idx] = coeff * wave_test[idx] ** (-1*alpha2)
    assert np.sum(np.isnan(law_test)) == 0

    # Put in terms of A_lambda / A_Ks, like the reddening object
    idx_K = np.where(abs(wave_test - K_wave) == min(abs(wave_test - K_wave)))
    law_test /= law_test[idx_K]
    law_test *= u.mag

    # Compare law_test and the output from the redlaw object
    law_output = red_law.broken_powerlaw(wave_test, 1)

    assert len(law_test) == len(law_output)
    assert np.sum(np.isnan(law_output)) == 0

    # Calculate the difference between test calc and code output.
    # Make sure they agree within tolerance
    diff = (law_output - law_test) / u.mag
    assert np.all(diff < 10**-4)

    # Let's also make sure the slope of the extinction law
    # in log-log space matches what it should be
    log_wave = np.log10(wave_test)
    log_output = np.log10(law_output / u.mag)

    idx1 = np.where(abs(wave_test-1.28) == np.min(abs(wave_test-1.28)))
    idx2 = np.where(abs(wave_test-1.53) == np.min(abs(wave_test-1.53)))
    slope = (log_output[idx1] - log_output[idx2]) / (log_wave[idx1] - log_wave[idx2])
    assert abs(slope - (-1.0 * alpha2)) < 10**-4

    idx1 = np.where(abs(wave_test-1.7) == np.min(abs(wave_test-1.7)))
    idx2 = np.where(abs(wave_test-2.1) == np.min(abs(wave_test-2.1)))
    slope = (log_output[idx1] - log_output[idx2]) / (log_wave[idx1] - log_wave[idx2])
    assert abs(slope - (-1.0 * alpha1)) < 10**-4

    # If desired (debug only), make plot to see what law looks like
    if plots:
        # Test plot: these should match nearly exactly
        py.figure()
        py.loglog(wave_test, law_test, 'k.', ms=13, label='Test')
        py.loglog(wave_test, law_output, 'r.', ms=8, label='Code')
        py.loglog(red_law.wave*10**-4, red_law.obscuration, 'r-')
        py.axvline(1.63, linestyle='--', color='blue')
        py.xlabel('Wavelength')
        py.ylabel('A_lambda / A_Ks')
        py.legend()
        py.gca().invert_xaxis()
        py.title('2-segment Extlaw: should match')

    # Additional test here: make sure get_red_law is working
    # for broken power law as intended
    red_law_test = reddening.get_red_law('broken_pl,{0},{1},{2}'.format(lambda_limits, alpha_vals, K_wave))
    assert np.all(red_law_test.obscuration == red_law.obscuration)

    #===============================#
    # Test 1: 4-segment law
    #===============================#
    alpha1 = 3.7
    alpha2 = 2.5
    alpha3 = 3.0
    alpha4 = 2.23
    lambda_limits = [0.4, 0.8, 1.27, 1.63, 2.16]
    alpha_vals = [alpha1, alpha2, alpha3, alpha4]
    K_wave = 2.14

    red_law = reddening.RedLawBrokenPowerLaw(lambda_limits, alpha_vals, K_wave)

    # Calculate what the redlaw *should* be across wavelength range,
    # by-hand
    wave_test = np.arange(0.4, 2.15+0.001, 0.01)
    law_test = np.ones(len(wave_test)) * np.nan

    # 0.4 - 0.8
    idx = np.where( (wave_test >= 0.4) & (wave_test < 0.8))
    coeff = 1
    law_test[idx] = coeff * wave_test[idx] ** (-1.0 * alpha1)

    # 0.8 - 1.27
    idx = np.where( (wave_test >= 0.8) & (wave_test < 1.27))
    coeff  = (0.8 ** (-1.0 * alpha1)) / (0.8 ** (-1.0 * alpha2))
    law_test[idx] = coeff * wave_test[idx] ** (-1.0 * alpha2)

    # 1.27 - 1.63
    idx = np.where( (wave_test >= 1.27) & (wave_test < 1.63))
    coeff1 = (0.8 ** (-1.0 * alpha1)) / (0.8 ** (-1.0 * alpha2))
    coeff2 = (1.27 ** (-1.0 * alpha2)) / (1.27 ** (-1.0 * alpha3))
    coeff_f = coeff1 * coeff2
    law_test[idx] = coeff_f * wave_test[idx] ** (-1.0 * alpha3)

    # 1.63 - 2.14
    idx = np.where( (wave_test >= 1.63) & (wave_test < 2.16))
    coeff1 = (0.80 ** (-1.0 * alpha1)) / (0.80 ** (-1.0 *  alpha2))
    coeff2 = (1.27 ** (-1.0 * alpha2)) / (1.27 ** (-1.0 * alpha3))
    coeff3 = (1.63 ** (-1.0 * alpha3)) / (1.63 ** (-1.0 * alpha4))
    coeff_f = coeff1 * coeff2 * coeff3
    law_test[idx] = coeff_f * wave_test[idx] ** (-1.0 * alpha4)

    assert np.sum(np.isnan(law_test)) == 0

    # Put in terms of A_lambda / A_Ks, like the reddening object
    idx_K = np.where(abs(wave_test - K_wave) == min(abs(wave_test - K_wave)))
    law_test /= law_test[idx_K]
    law_test *= u.mag

    # Compare law_test and the output from the redlaw object
    law_output = red_law.broken_powerlaw(wave_test, 1)

    assert len(law_test) == len(law_output)
    assert np.sum(np.isnan(law_output)) == 0

    # Calculate the difference between test calc and code output.
    # Make sure they agree within tolerance
    diff = (law_output - law_test) / u.mag
    assert np.all(diff < 10**-4) 

    # Let's also make sure the slope of the extinction law
    # in log-log space matches what it should be
    log_wave = np.log10(wave_test)
    log_output = np.log10(law_output / u.mag)

    idx1 = np.argmin(abs(wave_test-0.45))
    idx2 = np.argmin(abs(wave_test-0.7))
    slope = (log_output[idx2] - log_output[idx1]) / (log_wave[idx2] - log_wave[idx1])
    assert abs(slope - (-1.0 * alpha1)) < 10**-4

    idx1 = np.argmin(abs(wave_test-0.9))
    idx2 = np.argmin(abs(wave_test-1.2))
    slope = (log_output[idx2] - log_output[idx1]) / (log_wave[idx2] - log_wave[idx1])
    assert abs(slope - (-1.0 * alpha2)) < 10**-4

    idx1 = np.argmin(abs(wave_test-1.28))
    idx2 = np.argmin(abs(wave_test-1.53))
    slope = (log_output[idx2] - log_output[idx1]) / (log_wave[idx2] - log_wave[idx1])
    assert abs(slope - (-1.0 * alpha3)) < 10**-4

    idx1 = np.argmin(abs(wave_test-1.7))
    idx2 = np.argmin(abs(wave_test-2.1))
    slope = (log_output[idx2] - log_output[idx1]) / (log_wave[idx2] - log_wave[idx1])
    assert abs(slope - (-1.0 * alpha4)) < 10**-4


    # If desired (debug only), make plot to see what law looks like
    if plots:
        # Test plot: these should match nearly exactly
        py.figure()
        py.loglog(wave_test, law_test, 'k.', ms=13, label='Test')
        py.loglog(wave_test, law_output, 'r.', ms=8, label='Code')
        py.loglog(red_law.wave*10**-4, red_law.obscuration, 'r-')
        py.axvline(1.63, linestyle='--', color='blue')
        py.axvline(1.27, linestyle='--', color='blue')
        py.axvline(0.8, linestyle='--', color='blue')
        py.xlabel('Wavelength')
        py.ylabel('A_lambda / A_Ks')
        py.legend()
        py.gca().invert_xaxis()
        py.title('4-segment Extlaw: should match')

    return


@pytest.mark.slow
def test_red_law_IsochronePhot():
    """
    Make sure each reddening law can run with IsochronePhot.

    Isochrone FITS files are written under ``spisea/tests/isochrones``.
    """
    # Define properties of stellar pop to models
    logAge = np.log10(5*10**6.) # Age in log(years)
    dist = 8000 # distance in parsec
    metallicity = 0 # Metallicity in [M/H]
    iso_dir = f'{spisea_path}/tests/isochrones'

    # Define evolution/atmosphere models and extinction law
    evo_model = evolution.MISTv1()
    atm_func = atmospheres.get_merged_atmosphere

    # Also specify filters for synthetic photometry.
    filt_list = ['wfc3,ir,f127m', 'wfc3,ir,f153m', 'nirc2,H', 'nirc2,Kp']

    # Define reddening laws and associated AKs vals
    redlaw_arr = ['F11', 'S10', 'NL20', 'I05', 'N09', 'RZ07',
                      'RL85', 'D16', 'F09,2.5,3', 'S16,1.55,0', 'DM16',
                      'H18b', 'NL18', 'C89,3.1', 'pl,2.12,0.9,2.4',
                      'broken_pl,[0.9,1.63,2.3],[3.0,2.23],2.12']

    aks_arr = [2.62, 2.46, 1.67, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3]
    for ii in range(len(redlaw_arr)):
        redlaw = reddening.get_red_law(redlaw_arr[ii])
        aks = aks_arr[ii]

        # Try to run isochrone phot
        iso_test = synthetic.IsochronePhot(
            logAge,
            aks,
            dist,
            metallicity=0,
            evo_model=evo_model,
            atm_func=atm_func,
            red_law=redlaw,
            filters=filt_list,
            min_mass=0.95,
            max_mass=1.05,
            iso_dir=iso_dir,
            recomp=True
        )
        print('----EL {0} works OK!-----'.format(redlaw_arr[ii]))

    return

def test_all_EL():
    """
    Test to make sure all expected ELs are callable
    and work OK.

    (except for the broken power law, which is tested elsewhere)
    """
    # Will crash if any of these calls crash
    red_law = reddening.RedLawNishiyama09()
    red_law = reddening.RedLawCardelli(3.1)
    red_law = reddening.RedLawRomanZuniga07()
    red_law = reddening.RedLawRiekeLebofsky()
    red_law = reddening.RedLawDamineli16()
    red_law = reddening.RedLawDeMarchi16()
    red_law = reddening.RedLawFitzpatrick09(1.8, 5)
    red_law = reddening.RedLawSchlafly16(1.55, 3.3)
    red_law = reddening.RedLawIndebetouw05()
    red_law = reddening.RedLawFritz11()
    red_law = reddening.RedLawHosek18b()
    red_law = reddening.RedLawSchoedel10()
    red_law = reddening.RedLawNoguerasLara18()
    red_law = reddening.RedLawNoguerasLara20()
    red_law = reddening.RedLawSODC(2.5)

    return


def test_extinction_curve_invalid_AKs():
    """
    ``RedLawBase.extinction_curve`` must accept a real AKs or mag Quantity
    and raise ``synphot.exceptions.SynphotError`` for other types.
    """
    from synphot.exceptions import SynphotError

    red_law = reddening.RedLawNishiyama09()

    # Valid scalar and magnitude Quantity should succeed.
    ext = red_law.extinction_curve(1.0)
    assert ext is not None
    ext = red_law.extinction_curve(1.0 * u.mag)
    assert ext is not None

    # Non-numeric AKs is invalid.
    try:
        red_law.extinction_curve('not-a-number')
    except SynphotError:
        pass
    else:
        raise AssertionError('Expected SynphotError for invalid AKs')

    return
