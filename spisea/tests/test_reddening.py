import numpy as np
from spisea import reddening, synthetic, evolution, atmospheres
import pylab as py
import os
import pdb


def test_RedLawBrokenPowerLaw(plots=False):
    #===============================#
    # Test 1: 2-segment law
    #===============================#
    alpha1 = 2.23
    alpha2 = 3.0
    lambda_limits = [2.16, 1.63, 1.27]
    alpha_vals = [alpha1, alpha2]
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

    # Compare law_test and the output from the redlaw object
    law_output = red_law.broken_powerlaw(wave_test, 1)

    assert len(law_test) == len(law_output)
    assert np.sum(np.isnan(law_output)) == 0

    # Calculate the difference between test calc and code output.
    # Make sure they agree within tolerance
    diff = law_output - law_test
    assert np.all(diff < 10**-4)

    # Let's also make sure the slope of the extinction law
    # in log-log space matches what it should be
    log_wave = np.log10(wave_test)
    log_output = np.log10(law_output)

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
    alpha1 = 2.23
    alpha2 = 3.0
    alpha3 = 2.5
    alpha4 = 3.7
    lambda_limits = [2.16, 1.63, 1.27, 0.8, 0.4]
    alpha_vals = [alpha1, alpha2, alpha3, alpha4]
    K_wave = 2.14

    red_law = reddening.RedLawBrokenPowerLaw(lambda_limits, alpha_vals, K_wave)

    # Calculate what the redlaw *should* be across wavelength range,
    # by-hand
    wave_test = np.arange(0.4, 2.15+0.001, 0.01)
    law_test = np.ones(len(wave_test)) * np.nan

    # 2.14 - 1.63
    idx = np.where( (wave_test >= 1.63) & (wave_test < 2.17))
    coeff = 1
    law_test[idx] = coeff * wave_test[idx] ** (-1*alpha1)

    # 1.63 - 1.27
    idx = np.where( (wave_test >= 1.27) & (wave_test < 1.63))
    coeff  = (1.63 ** (-1*alpha1)) / (1.63 ** (-1*alpha2))
    law_test[idx] = coeff * wave_test[idx] ** (-1*alpha2)

    # 1.27 - 0.8
    idx = np.where( (wave_test >= 0.8) & (wave_test < 1.27))
    coeff1 = (1.63 ** (-1*alpha1)) / (1.63 ** (-1*alpha2))
    coeff2 = (1.27 ** (-1*alpha2)) / (1.27 ** (-1*alpha3))
    coeff_f = coeff1 * coeff2
    law_test[idx] = coeff_f * wave_test[idx] ** (-1*alpha3)

    # 0.8 - 0.4
    idx = np.where( (wave_test >= 0.4) & (wave_test < 0.8))
    coeff1 = (1.63 ** (-1*alpha1)) / (1.63 ** (-1*alpha2))
    coeff2 = (1.27 ** (-1*alpha2)) / (1.27 ** (-1*alpha3))
    coeff3 = (0.8 ** (-1*alpha3)) / (0.8 ** (-1*alpha4))
    coeff_f = coeff1 * coeff2 * coeff3
    law_test[idx] = coeff_f * wave_test[idx] ** (-1*alpha4)

    assert np.sum(np.isnan(law_test)) == 0

    # Put in terms of A_lambda / A_Ks, like the reddening object
    idx_K = np.where(abs(wave_test - K_wave) == min(abs(wave_test - K_wave)))
    law_test /= law_test[idx_K]

    # Compare law_test and the output from the redlaw object
    law_output = red_law.broken_powerlaw(wave_test, 1)

    assert len(law_test) == len(law_output)
    assert np.sum(np.isnan(law_output)) == 0

    # Calculate the difference between test calc and code output.
    # Make sure they agree within tolerance
    diff = law_output - law_test
    assert np.all(diff < 10**-4)

    # Let's also make sure the slope of the extinction law
    # in log-log space matches what it should be
    log_wave = np.log10(wave_test)
    log_output = np.log10(law_output)

    idx1 = np.where(abs(wave_test-1.28) == np.min(abs(wave_test-1.28)))
    idx2 = np.where(abs(wave_test-1.53) == np.min(abs(wave_test-1.53)))
    slope = (log_output[idx1] - log_output[idx2]) / (log_wave[idx1] - log_wave[idx2])
    assert abs(slope - (-1.0 * alpha2)) < 10**-4

    idx1 = np.where(abs(wave_test-1.7) == np.min(abs(wave_test-1.7)))
    idx2 = np.where(abs(wave_test-2.1) == np.min(abs(wave_test-2.1)))
    slope = (log_output[idx1] - log_output[idx2]) / (log_wave[idx1] - log_wave[idx2])
    assert abs(slope - (-1.0 * alpha1)) < 10**-4

    idx1 = np.where(abs(wave_test-0.9) == np.min(abs(wave_test-0.9)))
    idx2 = np.where(abs(wave_test-1.2) == np.min(abs(wave_test-1.2)))
    slope = (log_output[idx1] - log_output[idx2]) / (log_wave[idx1] - log_wave[idx2])
    assert abs(slope - (-1.0 * alpha3)) < 10**-4

    idx1 = np.where(abs(wave_test-0.45) == np.min(abs(wave_test-0.45)))
    idx2 = np.where(abs(wave_test-0.7) == np.min(abs(wave_test-0.7)))
    slope = (log_output[idx1] - log_output[idx2]) / (log_wave[idx1] - log_wave[idx2])
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

def test_red_law_IsochronePhot():
    """
    Make sure each reddening law can run with IsochronePhot
    """
    # Define properties of stellar pop to models
    logAge = np.log10(5*10**6.) # Age in log(years)
    dist = 8000 # distance in parsec
    metallicity = 0 # Metallicity in [M/H]

    # Define evolution/atmosphere models and extinction law
    evo_model = evolution.MISTv1()
    atm_func = atmospheres.get_merged_atmosphere

    # Also specify filters for synthetic photometry.
    filt_list = ['wfc3,ir,f127m', 'wfc3,ir,f153m', 'nirc2,H', 'nirc2,Kp']

    # Define reddening laws and associated AKs vals
    redlaw_arr = ['F11', 'S10', 'NL20', 'I05', 'N09', 'RZ07',
                      'RL85', 'D16', 'F09,2.5,3', 'S16,1.55,0', 'DM16',
                      'H18b', 'NL18', 'C89,3.1', 'pl,2.12,0.9,2.4',
                      'broken_pl,[2.3,1.63,0.9],[2.23, 3.0],2.12']

    aks_arr = [2.62, 2.46, 1.67, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3]
    for ii in range(len(redlaw_arr)):
        redlaw = reddening.get_red_law(redlaw_arr[ii])
        #redlaw_arr[ii]
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
            iso_dir='isochrones/',
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

def test_RedLawDeliwala26():
    """
    Check that adding dust changes RC magnitudes in the measured proportions 
    from RedLawDeliwala26.

    Use the representative RC atmosphere from the 8 Gyr, [M/H] = +0.3
    MIST population, selecting the model nearest the median RC temperature.
    Begin with F11 extinction of 2.5 mag at 2.166 microns, then add
    0.001 mag at that wavelength using Deliwala+26.
    """

    measured_ratios = {
        'f115w': 3.06,
        'f140m': 2.21,
        'f182m': 1.21,
        'f212n': 1.00,
        'f323n': 1.05,
        'f405n': 0.39,
    }

    # Choose an RC model with a typical temperature in the 8 Gyr population.
    isochrone = evolution.MISTv1().isochrone(age=8e9, metallicity=0.3)
    rc_models = isochrone[isochrone['phase'] == 3] 
    median_log_temperature = np.median(rc_models['logT'])
    closest_to_median = np.argmin(abs(rc_models['logT'] - median_log_temperature))
    star = rc_models[closest_to_median]
    rc_spectrum = atmospheres.get_merged_atmosphere(
        temperature=10**star['logT'],
        gravity=star['logg'],
        metallicity=0.3,
    )
    existing_extinction = reddening.RedLawFritz11().reddening(2.5)
    spectrum_before = rc_spectrum * existing_extinction

    # Add a small amount of dust with the Deliwala+26 law.
    law = reddening.RedLawDeliwala26()
    extra_extinction = law.reddening(0.001)
    spectrum_after = spectrum_before * extra_extinction

    # Measure how much fainter the star becomes in each filter.
    magnitude_changes = {}
    for band in measured_ratios:
        bandpass = synthetic.get_filter_info('jwst,' + band)
        magnitude_before = synthetic.mag_in_filter(spectrum_before, bandpass)
        magnitude_after = synthetic.mag_in_filter(spectrum_after, bandpass)
        magnitude_changes[band] = magnitude_after - magnitude_before

    # Compare each change with F212N: e.g., F115W should change 3.06 times as much.
    f212n_change = magnitude_changes['f212n']
    assert np.isfinite(f212n_change) and f212n_change > 0
    for band, measured_ratio in measured_ratios.items():
        recovered_ratio = magnitude_changes[band] / f212n_change
        np.testing.assert_allclose(
            recovered_ratio, measured_ratio, rtol=0.01, atol=0,
            err_msg=f'{band}: Deliwala26 differential ratio',
        )

