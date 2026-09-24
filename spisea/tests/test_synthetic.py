import os
import pdb
import time
from unittest.mock import patch

import pytest
import warnings
import importlib
import spisea
import numpy as np
import pylab as plt
from astropy import units
from astropy.table import Table
from astropy import units as u
from synphot import units as su
from astropy import constants as c
from scipy.spatial import cKDTree as KDTree
from spisea import synthetic as syn
from spisea.imf import imf, multiplicity
from spisea import reddening, evolution, atmospheres, ifmr

spisea_path = os.path.dirname(spisea.__file__)

def test_isochrone(plot=False):
    logAge = 6.7
    AKs = 2.7
    distance = 4000

    startTime = time.time()
    iso = syn.Isochrone(logAge, AKs, distance, verbose=plot)
    print('Test completed in: %d seconds' % (time.time() - startTime))
    # Typically takes 104 - 120 seconds.
    # Limited by stsynphot grid_to_spec / atmosphere lookup in atmospheres.py

    assert iso.points.meta['LOGAGE'] == logAge
    assert iso.points.meta['AKS'] == AKs
    assert iso.points.meta['DISTANCE'] == distance
    assert len(iso.points) > 100

    if plot:
        plt.figure(1)
        iso.plot_HR_diagram()

        plt.figure(2)
        iso.plot_mass_luminosity()

    #return iso

def test_iso_wave():
    """
    Test to make sure isochrones generated have spectra with the proper
    wavelength range, and that the user has control over that wavelength
    range (propagated through IsochronePhot)
    """
    # Define isochrone parameters
    logAge = 6.7 # Age in log(years)
    AKs = 2.7 # extinction in mags
    dist = 4000 # distance in parsec
    metal = 0.0 # metallicity
    iso_dir = f'{spisea_path}/tests/isochrones'

    # Define evolution/atmosphere models and extinction law (optional)
    evo_model = evolution.MergedBaraffePisaEkstromParsec()
    atm_func = atmospheres.get_merged_atmosphere
    red_law = reddening.RedLawHosek18b()

    # Also specify filters for synthetic photometry (optional). Here we use
    # the HST WFC3-IR F127M, F139M, and F153M filters
    filt_list = ['wfc3,ir,f127m']

    # First, let's make sure the vega spectrum has the proper limits
    vega = syn.Vega()

    assert np.min(vega.waveset) == 995 * u.AA
    assert np.max(vega.waveset) == 100200 * u.AA

    # Make Isochrone object. Will use wave_range = [3000,52000].
    # Make sure range matches to resolution of atmosphere.
    wave_range1 = [3000, 52000]
    my_iso = syn.IsochronePhot(
        logAge, AKs, dist,
        metallicity=metal,
        evo_model=evo_model,
        atm_func=atm_func,
        red_law=red_law,
        filters=filt_list,
        mass_sampling=10,
        wave_range=wave_range1,
        recomp=True,
        iso_dir=iso_dir
    )

    test = my_iso.spec_list[0]

    assert np.min(test.waveset) == 3010 * u.AA
    assert np.max(test.waveset) == 51900 * u.AA

    # Now let's try changing the wave range. Is it carried through
    # properly?
    wave_range2 = [1200, 90000]
    my_iso = syn.IsochronePhot(
        logAge,
        AKs,
        dist,
        evo_model=evo_model,
        atm_func=atm_func,
        red_law=red_law,
        filters=filt_list,
        mass_sampling=10,
        wave_range=wave_range2,
        recomp=True,
        iso_dir=iso_dir
    )

    test2 = my_iso.spec_list[0]

    assert np.min(test2.waveset) == 1205 * u.AA
    assert np.max(test2.waveset) == 89800 * u.AA

    # Does the error exception catch the bad wave_range?
    wave_range3 = [1200, 1000000]
    try:
        my_iso = syn.IsochronePhot(
            logAge,
            AKs,
            dist,
            evo_model=evo_model,
            atm_func=atm_func,
            red_law=red_law,
            filters=filt_list,
            mass_sampling=10,
            wave_range=wave_range3,
            recomp=True,
            iso_dir=iso_dir
        )
        print('WAVE TEST FAILED!!! Should have crashed here, wavelength range out of bounds')
        raise ValueError()
    except:
        print('Wavelength out of bound condition passed. Test is good')
        pass
    return

def test_IsochronePhot(plot=False):
    logAge = 6.7
    AKs = 2.7
    distance = 4000
    filt_list = ['wfc3,ir,f127m', 'nirc2,J']
    mass_sampling=1
    iso_dir = f'{spisea_path}/tests/isochrones'

    evo_model = evolution.MISTv1()
    atm_func = atmospheres.get_merged_atmosphere
    redlaw = reddening.RedLawNishiyama09()

    startTime = time.time()
    iso = syn.IsochronePhot(
        logAge,
        AKs,
        distance,
        evo_model=evo_model,
        atm_func=atm_func,
        red_law=redlaw,
        filters=filt_list,
        mass_sampling=mass_sampling,
        iso_dir=iso_dir,
        recomp=True
    )
    endTime = time.time()
    print('IsochronePhot generated in: %d seconds' % (endTime - startTime))
    # Typically takes 40 seconds if file is regenerated.
    # Limited by stsynphot grid_to_spec / atmosphere lookup in atmospheres.py

    assert iso.points.meta['LOGAGE'] == logAge
    assert iso.points.meta['AKS'] == AKs
    assert iso.points.meta['DISTANCE'] == distance
    assert len(iso.points) > 100

    assert 'm_nirc2_J' in iso.points.colnames

    if plot:
        plt.figure(1)
        iso.plot_CMD('mag814w', 'mag160w')

        plt.figure(2)
        iso.plot_mass_magnitude('mag160w')

    # Finally, let's test the isochronePhot file generation
    metal_value = 0.
    metal_sign = 'm' if metal_value < 0 else 'p'
    assert os.path.exists(f'{iso_dir}/iso_{logAge:.2f}_{AKs:4.2f}_{str(distance).zfill(5)}_{metal_sign}{metal_value:.2f}.fits')

    # Check 1: If we try to remake the isochrone, does it read the file rather than
    # making a new one
    iso_new = syn.IsochronePhot(
        logAge,
        AKs,
        distance,
        evo_model=evo_model,
        atm_func=atm_func,
        red_law=redlaw,
        filters=filt_list,
        mass_sampling=mass_sampling,
        iso_dir=iso_dir
    )

    assert iso_new.recalc == False

    # Check 2: Confirm that adding a new column to an existing isochrone works properly.
    #    Does the new filter get added to the isochrone? And the old ones still there?
    #    Does the computed data for the new filter match the same result if you fully regenerate the isochrone?
    iso_new_addfilt = syn.IsochronePhot(
        logAge,
        AKs,
        distance,
        evo_model=evo_model,
        atm_func=atm_func,
        red_law=redlaw,
        filters=filt_list+['2mass,Ks'],
        mass_sampling=mass_sampling,
        iso_dir=iso_dir
    )

    assert iso_new_addfilt.recalc == False
    assert 'm_2mass_Ks' in iso_new_addfilt.points.colnames
    assert 'm_nirc2_J' in iso_new_addfilt.points.colnames

    iso_new_3filt = syn.IsochronePhot(
        logAge,
        AKs,
        distance,
        evo_model=evo_model,
        atm_func=atm_func,
        red_law=redlaw,
        filters=filt_list+['2mass,Ks'],
        mass_sampling=mass_sampling,
        recomp=True,
        iso_dir=iso_dir
    )
    np.testing.assert_almost_equal(iso_new_addfilt.points['m_2mass_Ks'], iso_new_3filt.points['m_2mass_Ks'])
    assert iso_new_3filt.recalc==True

    # Check 3: If we change evo model, atmo model, or redlaw,
    # does IsochronePhot regenerate the isochrone and overwrite the existing one?
    evo2 = evolution.MergedBaraffePisaEkstromParsec()
    mass_sampling=20

    iso_new = syn.IsochronePhot(
        logAge,
        AKs,
        distance,
        evo_model=evo2,
        atm_func=atm_func,
        red_law=redlaw,
        filters=filt_list,
        mass_sampling=mass_sampling,
        iso_dir=iso_dir
    )

    assert iso_new.recalc == True

    redlaw2 = reddening.RedLawHosek18b()
    iso_new = syn.IsochronePhot(
        logAge,
        AKs,
        distance,
        evo_model=evo2,
        atm_func=atm_func,
        red_law=redlaw2,
        filters=filt_list,
        mass_sampling=mass_sampling,
        iso_dir=iso_dir
    )

    assert iso_new.recalc == True

    atm2 = atmospheres.get_castelli_atmosphere
    iso_new = syn.IsochronePhot(
        logAge,
        AKs,
        distance,
        evo_model=evo2,
        atm_func=atm2,
        red_law=redlaw2,
        filters=filt_list,
        mass_sampling=mass_sampling,
        iso_dir=iso_dir
    )

    assert iso_new.recalc == True

    return

def test_ResolvedCluster():
    # Define cluster parameters
    logAge = 6.7
    AKs = 2.4
    distance = 4000
    cluster_mass = 10**5.
    mass_sampling=5
    iso_dir = f'{spisea_path}/tests/isochrones'

    # Test filters
    filt_list = ['nirc2,J', 'nirc2,Kp']

    startTime = time.time()

    evo = evolution.MergedBaraffePisaEkstromParsec()
    atm_func = atmospheres.get_merged_atmosphere

    red_law = reddening.RedLawNishiyama09()

    iso = syn.IsochronePhot(
        logAge,
        AKs,
        distance,
        evo_model=evo,
        atm_func=atm_func,
        red_law=red_law,
        filters=filt_list,
        mass_sampling=mass_sampling,
        iso_dir=iso_dir,
        recomp=True,
        verbose=False
    )

    print('Constructed isochrone: %d seconds' % (time.time() - startTime))

    # Now to create the cluster.
    imf_mass_limits = np.array([0.01, 0.05, 0.22, 0.55, 8, 120])
    imf_powers = np.array([-0.6, -0.25, -1.3, -2.3, -2.35])

    ##########
    # Start without multiplicity
    ##########
    my_imf1 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=None)
    print('Constructed IMF: %d seconds' % (time.time() - startTime))

    cluster1 = syn.ResolvedCluster(iso, my_imf1, cluster_mass)
    clust1 = cluster1.star_systems
    print('Constructed cluster: %d seconds' % (time.time() - startTime))

    # Check that stars are returned
    assert len(clust1) > 0

    # Check that the total mass in stars is less than requested (no compact objects).
    cluster_mass_out = clust1['systemMass'].sum()
    assert cluster_mass_out < cluster_mass

    plt.figure(3)
    plt.clf()
    plt.plot(clust1['m_nirc2_J'] - clust1['m_nirc2_Kp'], clust1['m_nirc2_J'], 'r.')
    plt.plot(iso.points['m_nirc2_J'] - iso.points['m_nirc2_Kp'], iso.points['m_nirc2_J'], 'c.')
    plt.gca().invert_yaxis()

    # *** Visual Inspections: ***
    #  - check that points (red) fall between isochrone points (blue)

    ##########
    # Test with multiplicity
    ##########
    multi = multiplicity.MultiplicityUnresolved()
    my_imf2 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=multi)
    print('Constructed IMF with multiples: %d seconds' % (time.time() - startTime))

    cluster2 = syn.ResolvedCluster(iso, my_imf2, cluster_mass)
    clust2 = cluster2.star_systems
    print('Constructed cluster with multiples: %d seconds' % (time.time() - startTime))

    assert len(clust2) > 0
    assert len(cluster2.companions) > 0
    assert np.sum(clust2['N_companions']) == len(cluster2.companions)

    ##########
    # Plots
    ##########
    # Plot an IR CMD and compare cluster members to isochrone.
    plt.figure(1)
    plt.clf()
    plt.plot(clust1['m_nirc2_J'] - clust1['m_nirc2_Kp'], clust1['m_nirc2_J'], 'r.')
    plt.plot(clust2['m_nirc2_J'] - clust2['m_nirc2_Kp'], clust2['m_nirc2_J'], 'b.')
    plt.plot(iso.points['m_nirc2_J'] - iso.points['m_nirc2_Kp'], iso.points['m_nirc2_J'], 'c-')
    plt.gca().invert_yaxis()
    plt.xlabel('J - Kp (mag)')
    plt.ylabel('J (mag)')

    # Plot a mass-magnitude relationship.
    plt.figure(2)
    plt.clf()
    plt.semilogx(clust1['mass'], clust1['m_nirc2_J'], 'r.')
    plt.semilogx(clust2['mass'], clust2['m_nirc2_J'], 'r.')
    plt.gca().invert_yaxis()
    plt.xlabel('Mass (Msun)')
    plt.ylabel('J (mag)')

    # # Plot the spectrum of the most massive star
    # idx = cluster.mass.argmax()
    # plt.clf()
    # plt.plot(cluster.stars[idx].wave, cluster.stars[idx].flux, 'k.')

    # # Plot an integrated spectrum of the whole cluster.
    # wave, flux = cluster.get_integrated_spectrum()
    # plt.clf()
    # plt.plot(wave, flux, 'k.')

    return

def test_ResolvedClusterDiffRedden():
    logAge = 6.7
    AKs = 2.4
    distance = 4000
    cluster_mass = 10**5.
    deltaAKs = 0.05
    mass_sampling=5
    iso_dir = f'{spisea_path}/tests/isochrones'

    # Test filters
    filt_list = ['nirc2,J', 'nirc2,Kp']

    startTime = time.time()

    evo = evolution.MergedBaraffePisaEkstromParsec()
    atm_func = atmospheres.get_merged_atmosphere

    red_law = reddening.RedLawNishiyama09()

    iso = syn.IsochronePhot(
        logAge,
        AKs,
        distance,
        evo_model=evo,
        atm_func=atm_func,
        red_law=red_law,
        filters=filt_list,
        mass_sampling=mass_sampling,
        iso_dir=iso_dir
    )

    print('Constructed isochrone: %d seconds' % (time.time() - startTime))

    imf_mass_limits = np.array([0.07, 0.5, 1, np.inf])
    imf_powers = np.array([-1.3, -2.3, -2.3])

    ##########
    # Start without multiplicity
    ##########
    my_imf1 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=None)
    print('Constructed IMF: %d seconds' % (time.time() - startTime))

    cluster1 = syn.ResolvedClusterDiffRedden(iso, my_imf1, cluster_mass, deltaAKs)
    clust1 = cluster1.star_systems
    print('Constructed cluster: %d seconds' % (time.time() - startTime))

    assert len(clust1) > 0

    plt.figure(3)
    plt.clf()
    plt.plot(clust1['m_nirc2_J'] - clust1['m_nirc2_Kp'], clust1['m_nirc2_J'], 'r.')
    plt.plot(iso.points['m_nirc2_J'] - iso.points['m_nirc2_Kp'], iso.points['m_nirc2_J'], 'c.')
    plt.gca().invert_yaxis()

    # *** Visual Inspections: ***
    #  - check that points (red) fall between isochrone points (blue)

    ##########
    # Test with multiplicity
    ##########
    multi = multiplicity.MultiplicityUnresolved()
    my_imf2 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=multi)
    print('Constructed IMF with multiples: %d seconds' % (time.time() - startTime))

    cluster2 = syn.ResolvedClusterDiffRedden(iso, my_imf2, cluster_mass, deltaAKs)
    clust2 = cluster2.star_systems
    print('Constructed cluster with multiples: %d seconds' % (time.time() - startTime))

    assert len(clust2) > 0
    assert len(cluster2.companions) > 0
    assert np.sum(clust2['N_companions']) == len(cluster2.companions)

    ##########
    # Plots
    ##########
    # Plot an IR CMD and compare cluster members to isochrone.
    plt.figure(1)
    plt.clf()
    plt.plot(clust1['m_nirc2_J'] - clust1['m_nirc2_Kp'], clust1['m_nirc2_J'], 'r.')
    plt.plot(clust2['m_nirc2_J'] - clust2['m_nirc2_Kp'], clust2['m_nirc2_J'], 'b.')
    plt.plot(iso.points['m_nirc2_J'] - iso.points['m_nirc2_Kp'], iso.points['m_nirc2_J'], 'c-')
    plt.gca().invert_yaxis()
    plt.xlabel('J - Kp (mag)')
    plt.ylabel('J (mag')

    # Plot a mass-magnitude relationship.
    plt.figure(2)
    plt.clf()
    plt.semilogx(clust1['mass'], clust1['m_nirc2_J'], 'r.')
    plt.semilogx(clust2['mass'], clust2['m_nirc2_J'], 'r.')
    plt.gca().invert_yaxis()
    plt.xlabel('Mass (Msun)')
    plt.ylabel('J (mag)')

    return

def test_UnresolvedCluster():
    log_age = 6.7
    AKs = 0.0
    distance = 4000
    metallicity=0
    cluster_mass = 10**4.

    startTime = time.time()
    multi = multiplicity.MultiplicityUnresolved()
    imf_in = imf.Salpeter_Kirkpatrick_2024(multiplicity=multi)
    evo = evolution.MergedBaraffePisaEkstromParsec()
    atm_func = atmospheres.get_merged_atmosphere
    iso = syn.Isochrone(log_age, AKs, distance, metallicity=metallicity,
                            evo_model=evo, atm_func=atm_func, mass_sampling=10)
    print('Made Isochrone: %d seconds' % (time.time() - startTime))

    cluster = syn.UnresolvedCluster(iso, imf_in, cluster_mass)
    print('Constructed unresolved cluster: %d seconds' % (time.time() - startTime))

    # Plot an integrated spectrum of the whole cluster.
    wave = cluster.wave_trim
    flux = cluster.spec_trim
    plt.clf()
    plt.plot(wave, flux, 'k.')

    return

def test_ifmr_multiplicity():
    # Define cluster parameters
    logAge = 9.7
    AKs = 0.0
    distance = 1000
    cluster_mass = 1e6
    mass_sampling = 5
    iso_dir = f'{spisea_path}/tests/isochrones'

    # Test all filters
    filt_list = ['nirc2,Kp', 'nirc2,H', 'nirc2,J']

    startTime = time.time()
    
    evo = evolution.MergedPhillipsBaraffePisaEkstromParsec()
    atm_func = atmospheres.get_merged_atmosphere
    ifmr_obj = ifmr.IFMR_Raithel18()

    red_law = reddening.RedLawNishiyama09()

    iso = syn.IsochronePhot(
        logAge,
        AKs,
        distance,
        evo_model=evo,
        atm_func=atm_func,
        red_law=red_law,
        filters=filt_list,
        mass_sampling=mass_sampling,
        iso_dir=iso_dir,
        recomp=True
    )

    print('Constructed isochrone: %d seconds' % (time.time() - startTime))

    # Now to create the cluster.
    imf_mass_limits = np.array([0.01, 0.07, 0.5, 1, np.inf])
    imf_powers = np.array([-0.3, -1.3, -2.3, -2.3])

    ##########
    # Start without multiplicity and IFMR
    ##########
    my_imf1 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=None)
    print('Constructed IMF: %d seconds' % (time.time() - startTime))

    cluster1 = syn.ResolvedCluster(iso, my_imf1, cluster_mass, ifmr=ifmr_obj)
    clust1 = cluster1.star_systems
    print('Constructed cluster: %d seconds' % (time.time() - startTime))


    ##########
    # Test with multiplicity and IFMR
    ##########
    multi = multiplicity.MultiplicityUnresolved()
    my_imf2 = imf.Salpeter_Kirkpatrick_2024(multiplicity=multi)
    print('Constructed IMF with multiples: %d seconds' % (time.time() - startTime))

    cluster2 = syn.ResolvedCluster(iso, my_imf2, cluster_mass, ifmr=ifmr_obj)
    clust2 = cluster2.star_systems
    comps2 = cluster2.companions
    print('Constructed cluster with multiples: %d seconds' % (time.time() - startTime))

    ##########
    # Tests
    ##########

    # Check that we have black holes, neutron stars, and white dwarfs in both.
    assert len(np.where(clust1['phase'] == 101)) > 0   # WD
    assert len(np.where(clust2['phase'] == 101)) > 0
    assert len(np.where(clust1['phase'] == 102)) > 0   # NS
    assert len(np.where(clust2['phase'] == 102)) > 0
    assert len(np.where(clust1['phase'] == 103)) > 0   # BH
    assert len(np.where(clust2['phase'] == 103)) > 0
    assert len(np.where(clust1['phase'] == 90)) > 0   # BD
    assert len(np.where(clust2['phase'] == 90)) > 0

    # Now check that we have companions that are WDs, NSs, BHs, and BDs
    assert len(np.where(comps2['phase'] == 101)) > 0
    assert len(np.where(comps2['phase'] == 102)) > 0
    assert len(np.where(comps2['phase'] == 103)) > 0
    assert len(np.where(comps2['phase'] == 90)) > 0

    # Make sure no funky phase designations (due to interpolation effects)
    # slipped through

    bd_masses = np.where((clust1['mass'] >= 0.01) & (clust1['mass'] <= 0.08))
    bd_mask = (clust1['mass'] >= 0.01) & (clust1['mass'] <= 0.08)
    print(clust1[['mass', 'phase']][bd_masses])
    idx = np.where( (clust1['phase'] > 5) & (clust1['phase'] < 90) & (clust1['phase'] != 9) )
    assert len(idx[0]) == 0

    """
    07/2024: Added more testing criteria for brown dwarf stars to ensure they are labeled appropriately for masses from 0.01 - 0.08 M_sun.
    """
    # For cluster objects
    MIN_MASS = 0.1
    BD_MIN_MASS = 0.01
    BD_MAX_MASS = 0.08

    wd_idx = np.where(clust1['phase'] == 101)
    ns_idx = np.where(clust1['phase'] == 102)
    bh_idx = np.where(clust1['phase'] == 103)
    bd_idx = np.where(clust1['phase'] == 90)

    print(clust1[bd_idx])
    assert np.all(clust1['mass'][wd_idx] > MIN_MASS)
    assert np.all(clust1['mass'][ns_idx] > MIN_MASS)
    assert np.all(clust1['mass'][bh_idx] > MIN_MASS)
    assert np.all(clust1['mass'][bd_idx] < MIN_MASS)

    # For companion objects
    comp_wd_idx = np.where(comps2['phase'] == 101)
    comp_ns_idx = np.where(comps2['phase'] == 102)
    comp_bh_idx = np.where(comps2['phase'] == 103)
    comp_bd_idx = np.where(comps2['phase'] == 90)
    assert np.all(comps2['mass'][comp_wd_idx] > MIN_MASS)
    assert np.all(comps2['mass'][comp_ns_idx] > MIN_MASS)
    assert np.all(comps2['mass'][comp_bh_idx] > MIN_MASS)
    assert np.all(comps2['mass'][comp_bd_idx] < MIN_MASS)

    # Ensure brown dwarfs are within 0.01 and 0.08 solar masses
    assert np.all((clust1['mass'][bd_idx] >= BD_MIN_MASS) &
                  (clust1['mass'][bd_idx] <= BD_MAX_MASS))
    assert np.all((comps2['mass'][comp_bd_idx] >= BD_MIN_MASS) &
                  (comps2['mass'][comp_bd_idx] <= BD_MAX_MASS))

    # Ensure no other objects are in the brown dwarf mass range
    non_bd_idx = np.where((clust1['phase'] != 90) &
                          (clust1['mass'] >= BD_MIN_MASS) &
                          (clust1['mass'] <= BD_MAX_MASS))
    assert len(non_bd_idx[0]) == 0  # asserting no non-brown dwarf objects in BD mass range

    comp_non_bd_idx = np.where((comps2['phase'] != 90) &
                               (comps2['mass'] >= BD_MIN_MASS) &
                               (comps2['mass'] <= BD_MAX_MASS))
    print(comps2[comp_non_bd_idx]['phase', 'mass'])
    assert len(comp_non_bd_idx[0]) == 0  # asserting no non-brown dwarf companions in BD mass range

    # Ensure BD temperature assignment is working correctly
    assert np.all(clust1['Teff'][bd_idx] != np.nan)
    assert np.all(comps2['Teff'][comp_bd_idx] != np.nan)

    return

def test_metallicity():
    """
    Test isochrone generation at different metallicities
    """
    # Define isochrone parameters
    logAge = np.log10(5*10**6.)
    AKs = 0.8
    dist = 4000
    evo_model = evolution.MISTv1()
    atm_func = atmospheres.get_phoenixv16_atmosphere
    red_law = reddening.RedLawHosek18b()
    filt_list = ['wfc3,ir,f127m', 'wfc3,ir,f139m', 'wfc3,ir,f153m']
    iso_dir = f'{spisea_path}/tests/isochrones'

    # Start with a solar metallicity isochrone
    metallicity= 0.
    metal_sign = 'm' if metallicity < 0 else 'p'

    # Make Isochrone object, with high mass_sampling to decrease compute time
    my_iso = syn.IsochronePhot(
        logAge,
        AKs,
        dist,
        metallicity=metallicity,
        evo_model=evo_model,
        atm_func=atm_func,
        red_law=red_law,
        filters=filt_list,
        mass_sampling=10,
        iso_dir=iso_dir
    )

    # Test isochrone properties
    assert my_iso.points.meta['METAL_IN'] == 0.0
    assert os.path.exists(f'{iso_dir}/iso_6.70_0.80_04000_p0.00.fits')

    # Now for non-solar metallicity
    metallicity= -1.5

    # Make Isochrone object, with high mass_sampling to decrease compute time
    my_iso = syn.IsochronePhot(
        logAge,
        AKs,
        dist,
        metallicity=metallicity,
        evo_model=evo_model,
        atm_func=atm_func,
        red_law=red_law,
        filters=filt_list,
        mass_sampling=10,
        iso_dir=iso_dir
    )

    # MIST model sub-directory names changed in SPISEA v2.1.4 update;
    # changing what "metal_act" value was. version 1 of MIST grid
    # names mistakenly used 0.015 for Zsolar, which was fixed for
    # version 2 onward. Note that this only effects this variable,
    # and not the grid output itself.
    # So, here the expected value for metal_act changes depending
    # on the grid value
    grid_version = evolution.get_installed_grid_num('{0}/evolution/'.format(os.environ['SPISEA_MODELS']))

    if grid_version == 1:
        metal_act = np.log10(0.00047 / 0.0142) # For Mist isochrones
    else:
        metal_act = np.log10(0.00045 / 0.0142) # For Mist isochrones

    # Test isochrone properties
    assert my_iso.points.meta['METAL_IN'] == -1.5
    assert np.isclose(my_iso.points.meta['METAL_ACT'], metal_act)
    assert os.path.exists(f'{iso_dir}/iso_6.70_0.80_04000_m1.50.fits')

    return

def test_cluster_mass():
    # Define cluster parameters
    logAge = 6.7
    AKs = 2.4
    distance = 4000
    cluster_mass = 10**5.
    mass_sampling = 5
    iso_dir = f'{spisea_path}/tests/isochrones'

    # Test filters
    filt_list = ['nirc2,J', 'nirc2,Kp']

    startTime = time.time()

    # Define evolution/atmosphere models and extinction law
    evo = evolution.MISTv1()
    atm_func = atmospheres.get_merged_atmosphere
    red_law = reddening.RedLawHosek18b()

    iso = syn.IsochronePhot(
        logAge,
        AKs,
        distance,
        evo_model=evo,
        atm_func=atm_func,
        red_law=red_law,
        filters=filt_list,
        mass_sampling=mass_sampling,
        iso_dir=iso_dir
    )

    print('Constructed isochrone: %d seconds' % (time.time() - startTime))

    # Now to create the cluster.
    imf_mass_limits = np.array([0.2, 0.5, 1, 120.0])
    imf_powers = np.array([-1.3, -2.3, -2.3])

    # IFMR
    my_ifmr = ifmr.IFMR_Raithel18()


    ##########
    # Start without multiplicity
    ##########
    my_imf1 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=None)
    print('Constructed IMF: %d seconds' % (time.time() - startTime))

    cluster1 = syn.ResolvedCluster(iso, my_imf1, cluster_mass, ifmr=my_ifmr)
    clust1 = cluster1.star_systems
    print('Constructed cluster: %d seconds' % (time.time() - startTime))

    # Check that the total mass is within tolerance of input mass
    cluster_mass_out = clust1['systemMass'].sum()
    assert np.abs(cluster_mass_out - cluster_mass) < 200.0   # within 200 Msun of desired mass.
    print('Cluster Mass: IN = ', cluster_mass, " OUT = ", cluster_mass_out)

    cluster3 = syn.ResolvedCluster(iso, my_imf1, cluster_mass, ifmr=my_ifmr, keep_low_mass_stars=True)
    clust3 = cluster3.star_systems
    print('Constructed cluster w/ keep_low_mass_stars: %d seconds' % (time.time() - startTime))

    # Check that the total mass is within tolerance of input mass
    cluster_mass_out = clust3['systemMass'].sum()
    assert np.abs(cluster_mass_out - cluster_mass) < 200.0   # within 200 Msun of desired mass.
    print('Cluster Mass: IN = ', cluster_mass, " OUT = ", cluster_mass_out)

    ##########
    # Test with multiplicity
    ##########
    multi = multiplicity.MultiplicityUnresolved()
    my_imf2 = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=multi)
    print('Constructed IMF with multiples: %d seconds' % (time.time() - startTime))

    cluster2 = syn.ResolvedCluster(iso, my_imf2, cluster_mass, ifmr=my_ifmr)
    clust2 = cluster2.star_systems
    print('Constructed cluster with multiples: %d seconds' % (time.time() - startTime))

    # Check that the total mass is within tolerance of input mass
    cluster_mass_out = clust2['systemMass'].sum()
    assert np.abs(cluster_mass_out - cluster_mass) < 200.0   # within 200 Msun of desired mass.
    print('Cluster Mass: IN = ', cluster_mass, " OUT = ", cluster_mass_out)

    return

def test_keep_low_mass_stars():
    """
    Test "keep_low_mass_stars = True" functionality introduced in v2.2
    """
    # Define cluster parameters, pulling on an isochrone generated in an earlier test (since
    # we don't care about isochrone generation here
    logAge = 6.7
    AKs = 2.4
    distance = 4000
    cluster_mass = 10**5.
    mass_sampling = 5
    iso_dir = f'{spisea_path}/tests/isochrones'

    # Test filters
    filt_list = ['nirc2,J', 'nirc2,Kp']

    # Define evolution/atmosphere models and extinction law
    evo = evolution.MISTv1()
    atm_func = atmospheres.get_merged_atmosphere
    red_law = reddening.RedLawHosek18b()

    iso = syn.IsochronePhot(
        logAge,
        AKs,
        distance,
        evo_model=evo,
        atm_func=atm_func,
        red_law=red_law,
        filters=filt_list,
        mass_sampling=mass_sampling,
        iso_dir=iso_dir,
        recomp=True
    )

    # Get the minimum mass in the isochrones. This should be the lowest
    # mass psosbile when keep_low_mass_stars == False.
    # Make sure this min mass is low enough for a reasonalbe test
    min_mass_iso = np.min(iso.points['mass'])
    assert min_mass_iso >= 0.05

    # Define IMF + IFMR. Make sure IMF goes to really low masses,
    # below the 0.08 Msun limit of the MIST isochrones
    imf_min = 0.01
    imf_mass_limits = np.array([imf_min, 0.5, 1, 120.0])
    imf_powers = np.array([-1.3, -2.3, -2.3])
    my_ifmr = ifmr.IFMR_Raithel18()
    my_imf = imf.IMF_broken_powerlaw(imf_mass_limits, imf_powers,
                                      multiplicity=multiplicity.MultiplicityUnresolved())

    # Define 2 clusters: one without keep_low_mass_stars and the other
    # with keep_low_mass_stars
    clust_remove = syn.ResolvedCluster(iso, my_imf, cluster_mass, ifmr=my_ifmr, keep_low_mass_stars=False)
    clust_keep = syn.ResolvedCluster(iso, my_imf, cluster_mass, ifmr=my_ifmr, keep_low_mass_stars=True)

    # Check the star_systems tables: clust_remove loweset mass should match
    # min_mass_iso, while clust_keep should stretch to IMF limit (1%)
    assert np.isclose(np.min(clust_keep.star_systems['mass']), imf_min, rtol=0.01, atol=10**-8)
    assert (np.isclose(np.min(clust_remove.star_systems['mass']), imf_min, rtol=0.01, atol=10**-8) == False)

    # NOTE: companion table always had the low-mass stars, this functionality is only on the star_systems table

    # Finally, make sure these low-mass stars have properly assigned phase = 98
    idx = np.where(clust_keep.star_systems['mass'] < min_mass_iso)
    assert ( (len(idx[0]) > 0) & (np.all(clust_keep.star_systems['phase'][idx] == 98)) )

    return


def test_compact_object_companions():

    # Define cluster parameters
    logAge = 6.7
    AKs = 2.4
    distance = 4000
    cluster_mass = 10**4.
    mass_sampling=5
    iso_dir = f'{spisea_path}/tests/isochrones'

    # Test filters
    filt_list = ['nirc2,J', 'nirc2,Kp']

    startTime = time.time()

    evo = evolution.MergedBaraffePisaEkstromParsec()
    atm_func = atmospheres.get_merged_atmosphere

    red_law = reddening.RedLawNishiyama09()

    iso = syn.IsochronePhot(
        logAge,
        AKs,
        distance,
        evo_model=evo,
        atm_func=atm_func,
        red_law=red_law,
        filters=filt_list,
        mass_sampling=mass_sampling,
        iso_dir=iso_dir
    )

    print('Constructed isochrone: %d seconds' % (time.time() - startTime))

    clust_multiplicity = multiplicity.MultiplicityResolvedDK()

    massLimits = np.array([0.2, 0.5, 1, 120]) # mass segments
    powers = np.array([-1.3, -2.3, -2.3]) # power-law exponents
    clust_imf_Mult = imf.IMF_broken_powerlaw(massLimits, powers, multiplicity=clust_multiplicity)
    #clust_imf_Mult = imf.Kroupa_2001(multiplicity=clust_multiplicity)
    clust_Mult = syn.ResolvedCluster(iso, clust_imf_Mult, cluster_mass, ifmr=ifmr.IFMR_Raithel18())

    # Makes sure compact object companions not including MIST WDs
    #(i.e. those with no luminosity) are being given phases
    nan_lum_companions = clust_Mult.companions[np.isnan(clust_Mult.companions['L'])]

    assert (len(nan_lum_companions) == 0) | (all(np.isnan(nan_lum_companions['phase'])) == False)

def _require_cosmic():
    """
    Skip the calling test if the optional `cosmic` package is not installed,
    emitting a warning so the skip is visible (rather than silent).
    """
    if importlib.util.find_spec('cosmic') is None:
        msg = ('COSMIC integration test skipped: the optional `cosmic` package '
               'is not installed. Run in an environment with COSMIC (e.g. '
               '`astro_cosmic`) to exercise these tests.')
        warnings.warn(msg)
        pytest.skip(msg)

    return

def test_COSMIC_evolve():
    """
    Test the COSMIC external evolution model's evolve() method directly on a
    small, hand-built set of star systems and companions. Uses a young, low-mass
    population so no compact remnants/disruptions occur, keeping the run fast
    and avoiding the merger/disruption branches.

    Skipped (with a warning) if the optional `cosmic` package is not installed.
    """
    _require_cosmic()

    from astropy.table import Table

    # Build a minimal star_systems table: 2 binaries + 1 single
    star_systems = Table()
    star_systems['mass'] = np.array([1.0, 0.9, 0.5])
    star_systems['isMultiple'] = np.array([True, True, False])
    star_systems['N_companions'] = np.array([1, 1, 0])
    star_systems['systemMass'] = np.array([1.5, 1.3, 0.5])

    # Companions for the first two systems
    companions = Table()
    companions['system_idx'] = np.array([0, 1])
    companions['mass'] = np.array([0.5, 0.4])
    companions['log_a'] = np.array([1.0, 1.5])   # log10(AU)
    companions['e'] = np.array([0.1, 0.2])
    companions['i'] = np.array([30.0, 60.0])      # degrees
    companions['Omega'] = np.array([0.0, 0.0])
    companions['omega'] = np.array([0.0, 0.0])

    evo = evolution.COSMIC(keep_disrupted_companions=False, keep_COSMIC_tables=True)
    ss, comp = evo.evolve(star_systems, companions, logAge=8.0, metallicity=0.0)

    # Check that evolve() populated the expected output columns
    expected_cols = ['mass_current', 'Teff', 'L', 'logg', 'phase',
                     'kick', 'kick_x', 'kick_y', 'kick_z']
    for col in expected_cols:
        assert col in ss.colnames, 'star_systems missing column {0}'.format(col)
        assert col in comp.colnames, 'companions missing column {0}'.format(col)

    # Core bookkeeping invariant: companion count must match companion table length
    assert np.sum(ss['N_companions']) == len(comp)

    # Scalar kick must be the magnitude of the kick vector components
    np.testing.assert_allclose(
        ss['kick'], np.sqrt(ss['kick_x']**2 + ss['kick_y']**2 + ss['kick_z']**2))
    np.testing.assert_allclose(
        comp['kick'], np.sqrt(comp['kick_x']**2 + comp['kick_y']**2 + comp['kick_z']**2))

    # Young, low-mass stars should remain stellar (no compact remnants here).
    # Compact-object phase codes are 101 (WD), 102 (NS), 103 (BH).
    assert np.all(ss['phase'] < 100)
    assert np.all(comp['phase'] < 100)

    # keep_COSMIC_tables=True should store the raw COSMIC output tables
    for attr in ['bpp', 'bcm', 'initC', 'kick_info']:
        assert hasattr(evo, attr), 'COSMIC model missing table {0}'.format(attr)

    return

def test_COSMIC_ResolvedCluster():
    """
    Test the full COSMIC cluster pipeline, mirroring the Cluster_w_COSMIC
    tutorial: build an IsochronePhotExternalEvolution, then a ResolvedCluster
    with binary multiplicity (CSF_max=1, companion_max=True), and check the
    output tables.

    Skipped (with a warning) if the optional `cosmic` package is not installed.
    """
    _require_cosmic()

    # Cluster/isochrone parameters (kept small/cheap)
    logAge = 9.0
    AKs = 0.0
    distance = 4000
    metallicity = 0.0
    cluster_mass = 10**3.
    mass_sampling = 10
    atm_grid_dir = f'{spisea_path}/tests/atm_cosmic'

    filt_list = ['ubv,V']

    # External evolution model (keep tables so we can verify them)
    evo = evolution.COSMIC(keep_COSMIC_tables=True)
    atm_func = atmospheres.get_merged_atmosphere_w_bb_supplement
    red_law = reddening.RedLawCardelli(3.1)

    iso = syn.IsochronePhotExternalEvolution(
        logAge,
        AKs,
        distance,
        metallicity=metallicity,
        evo_model=evo,
        atm_func=atm_func,
        red_law=red_law,
        filters=filt_list,
        atm_grid_dir=atm_grid_dir,
        mass_sampling=mass_sampling,
        recomp=False
    )

    # COSMIC only supports binaries: resolved multiplicity, no higher-order systems
    clust_multiplicity = multiplicity.MultiplicityResolvedDK(CSF_max=1, companion_max=True)
    my_imf = imf.Kroupa_2001(multiplicity=clust_multiplicity)

    cluster = syn.ResolvedCluster(iso, my_imf, cluster_mass)
    star_systems = cluster.star_systems
    companions = cluster.companions

    # Basic sanity: stars and companions were produced
    assert len(star_systems) > 0
    assert len(companions) > 0

    # Core bookkeeping invariant
    assert np.sum(star_systems['N_companions']) == len(companions)

    # Kick columns present and self-consistent on both tables
    for tab in (star_systems, companions):
        for col in ['kick', 'kick_x', 'kick_y', 'kick_z']:
            assert col in tab.colnames
        np.testing.assert_allclose(
            tab['kick'], np.sqrt(tab['kick_x']**2 + tab['kick_y']**2 + tab['kick_z']**2))

    # Synthetic photometry column should exist
    assert 'm_ubv_V' in star_systems.colnames

    # All phase codes should be in the allowed set:
    # 0-9 stellar phases, 101 (WD), 102 (NS), 103 (BH)
    allowed_phases = set(range(0, 10)) | {101, 102, 103}
    ss_phases = set(np.unique(star_systems['phase']).astype(int).tolist())
    comp_phases = set(np.unique(companions['phase']).astype(int).tolist())
    assert ss_phases.issubset(allowed_phases), 'Unexpected star phases: {0}'.format(ss_phases - allowed_phases)
    assert comp_phases.issubset(allowed_phases), 'Unexpected companion phases: {0}'.format(comp_phases - allowed_phases)

    # keep_COSMIC_tables=True should expose the raw COSMIC tables on the evo model
    for attr in ['bpp', 'bcm', 'initC', 'kick_info']:
        assert hasattr(iso.evo_model, attr), 'COSMIC model missing table {0}'.format(attr)

    return

#=================================#
# Additional timing functions
#=================================#

def time_test_cluster():
    logAge = 6.7
    AKs = 2.7
    distance = 4000
    cluster_mass = 10**4
    iso_dir = f'{spisea_path}/tests/isochrones'

    startTime = time.time()

    evo = evolution.MergedBaraffePisaEkstromParsec()
    atm_func = atmospheres.get_merged_atmosphere
    red_law = reddening.RedLawNishiyama09()
    filt_list = ['nirc2,J', 'nirc2,Kp']

    iso = syn.IsochronePhot(
        logAge,
        AKs,
        distance,
        evo_model=evo,
        atm_func=atm_func,
        red_law=red_law,
        filters=filt_list,
        iso_dir=iso_dir
    )
    print('Constructed isochrone: %d seconds' % (time.time() - startTime))

    imf_limits = np.array([0.07, 0.5, 150])
    imf_powers = np.array([-1.3, -2.35])
    multi = multiplicity.MultiplicityUnresolved()
    my_imf = imf.IMF_broken_powerlaw(imf_limits, imf_powers, multiplicity=multi)
    print('Constructed IMF with multiples: %d seconds' % (time.time() - startTime))

    cluster = syn.ResolvedCluster(iso, my_imf, cluster_mass)
    print('Constructed cluster: %d seconds' % (time.time() - startTime))

    return

def model_young_cluster_object(resolved=False):
    log_age = 6.5
    AKs = 0.1
    distance = 8000.0
    cluster_mass = 10000.

    multi = multiplicity.MultiplicityUnresolved()
    imf_in = imf.Kroupa_2001(multiplicity=multi)
    evo = evolution.MergedBaraffePisaEkstromParsec()
    atm_func = atmospheres.get_merged_atmosphere
    iso = syn.Isochrone(
        log_age,
        AKs,
        distance,
        evo_model=evo,
        atm_func=atm_func,
        mass_sampling=10
    )

    if resolved:
        cluster = syn.ResolvedCluster(iso, imf_in, cluster_mass)
    else:
        cluster = syn.UnresolvedCluster(iso, imf_in, cluster_mass, wave_range=[19000,24000])

    # Plot the spectrum of the most massive star
    idx = cluster.mass_all.argmax()
    print('Most massive star is {0:f} M_sun.'.format(cluster.mass_all[idx]))
    #bigstar = cluster.spec_list_trim[idx]
    plt.figure(1)
    plt.clf()
    s = cluster.spec_list_trim[idx]
    plt.plot(s.waveset, s(s.waveset), 'k.')

    # Plot an integrated spectrum of the whole cluster.
    wave, flux = cluster.spec_list_trim[idx].waveset, cluster.spec_trim
    plt.figure(2)
    plt.clf()
    plt.plot(wave, flux, 'k.')

    return

def time_test_mass_match():
    log_age = 6.7
    AKs = 2.7
    distance = 4000
    cluster_mass = 5e3
    iso_dir = f'{spisea_path}/tests/isochrones'

    imf_in = imf.Kroupa_2001(multiplicity=None)

    start_time = time.time()
    iso = syn.IsochronePhot(log_age, AKs, distance, iso_dir=iso_dir)
    iso_masses = iso.points['mass']
    print('Generated iso masses in {0:.0f} s'.format(time.time() - start_time))

    start_time = time.time()
    star_masses, isMulti, compMass, sysMass = imf_in.generate_cluster(cluster_mass)
    print('Generated cluster masses in {0:.0f} s'.format(time.time() - start_time))

    def match_model_masses1(isoMasses, starMasses):
        indices = np.empty(len(starMasses), dtype=int)

        for ii in range(len(starMasses)):
            theMass = starMasses[ii]

            dm = np.abs(isoMasses - theMass)
            mdx = dm.argmin()

            # Model mass has to be within 10% of the desired mass
            if (dm[mdx] / theMass) > 0.1:
                indices[ii] = -1
            else:
                indices[ii] = mdx

        return indices


    def match_model_masses2(isoMasses, starMasses):
        isoMasses_tmp = isoMasses.reshape((len(isoMasses), 1))
        kdt = KDTree(isoMasses_tmp)

        starMasses_tmp = starMasses.reshape((len(starMasses), 1))
        q_results = kdt.query(starMasses_tmp, k=1)
        indices = q_results[1]

        dm_frac = np.abs(starMasses - isoMasses[indices]) / starMasses

        idx = np.where(dm_frac > 0.1)[0]
        indices[idx] = -1

        return indices

    print('Test #1 START')
    start_time = time.time()
    idx1 = match_model_masses1(iso_masses, star_masses)
    stop_time = time.time()
    print('Test #1 STOPPED after {0:.0f} seconds'.format(stop_time - start_time))

    print('Test #2 START')
    start_time = time.time()
    idx2 = match_model_masses2(iso_masses, star_masses)
    stop_time = time.time()
    print('Test #2 STOPPED after {0:.0f} seconds'.format(stop_time - start_time))

    return

def generate_Spera15_IFMR():
    """
    Make a set of objects based on the Spera15 IFMR for the purposes of testing
    Expect 28 total objects
    8 invalids, 8 WDs, 2 NSs, and 10 BHs
    """
    Spera = ifmr.IFMR_Spera15()

    def FeH_from_Z(Z):
        return np.log10(Z/0.019)

    metal = np.array([2.0e-4, 1.0e-3, 2.0e-3, 2.0e-2]) #ensure that all Spera metallicity regimes are represented

    FeH = FeH_from_Z(metal) #generate death mass takes metallicty as [Fe/H]
    #want to get a good range of masses for Spera, should expect 8 invalids, 8 WDs, 3 NSs, and 9 BHs
    ZAMS = np.array([-0.2*np.ones(len(FeH)), 0.2*np.ones(len(FeH)), 4.0*np.ones(len(FeH)), 9.2*np.ones(len(FeH)),
                    15.0*np.ones(len(FeH)), 30.0*np.ones(len(FeH)), 150.0*np.ones(len(FeH))])


    output_array = np.concatenate((Spera.generate_death_mass(ZAMS[0], FeH), Spera.generate_death_mass(ZAMS[1], FeH),
                                  Spera.generate_death_mass(ZAMS[2], FeH), Spera.generate_death_mass(ZAMS[3], FeH),
                                  Spera.generate_death_mass(ZAMS[4], FeH), Spera.generate_death_mass(ZAMS[5], FeH),
                                  Spera.generate_death_mass(ZAMS[6], FeH)), axis=1)

    #count up number of objects formed and ensure it matahces the number of stars input
    bad_idx = np.where(output_array[1] == -1)
    WD_idx = np.where(output_array[1] == 101)
    NS_idx = np.where(output_array[1] == 102)
    BH_idx = np.where(output_array[1] == 103)

    rem_mass = output_array[0]

    return bad_idx[0], WD_idx[0], NS_idx[0], BH_idx[0], rem_mass

def test_Spera15_IFMR1():
    """
    Check to make sure the total number of objects input matches the number of objects output (28)
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Spera15_IFMR()

    total = len(WD_idx) + len(NS_idx) + len(BH_idx) + len(bad_idx)
    assert total == 28 , "The # of input objects does not match the number of output objects for the Spera15 IFMR"

    return

def test_Spera15_IFMR2():
    """
    Check that all negative remnant masses have type code -1
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Spera15_IFMR()

    #check to make sure no unhandled negative remnant masses (ie without type code -1 assigned)
    neg_mass_idx = np.where(rem_mass < 0)
    assert set(bad_idx) == set(neg_mass_idx[0]), "There are unhandled negative remnant masses for the Spera15 IFMR"

    return

def test_Spera15_IFMR_3():
    """
    Check that there are no left over zeroes in the remnant mass array
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Spera15_IFMR()

    assert np.all(rem_mass != 0) , "There are left over zeros in the remnant mass array for the Spera15 IFMR"

    return

def test_Spera15_IFMR_4():
    """
    Check that the correct number of invalid objects were generated (8)
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Spera15_IFMR()

    assert len(bad_idx) == 8 , "There are not the right number of invalid objects for the Spera15 IFMR"

    return

def test_Spera15_IFMR_5():
    """
    Check that the correct number of WDs were generated (8)
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Spera15_IFMR()

    assert len(WD_idx) == 8 , "There are not the right number of WDs for the Spera15 IFMR"

    return

def test_Spera15_IFMR_6():
    """
    Check that the correct number of NSs were generated (8)
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Spera15_IFMR()

    assert len(NS_idx) == 2 , "There are not the right number of NSs for the Spera15 IFMR"

    return

def test_Spera15_IFMR_7():
    """
    Check that the correct number of invalid objects were generated (8)
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Spera15_IFMR()

    assert len(BH_idx) == 10 , "There are not the right number of BHs for the Spera15 IFMR"

    return

def generate_Raithel18_IFMR():
    """
    Make a set of objects using the Raithel18 IFMR for the purposes of testing
    Will make a total of 16 objects
    Will have 3 invalid objects and 2 WDs
    """

    Raithel = ifmr.IFMR_Raithel18()
    ZAMS = np.array([-0.2, 0.2, 1.0, 7.0, 10.0, 14.0, 16.0, 18.0, 18.6, 22.0, 26.0, 28.0, 50.0, 61.0, 119.0, 121.0])
    #3 invalid indices, 2 WDs, cannot make statements about #of BHs and NSs because the Raithel IFMR has some randomness

    output_array = Raithel.generate_death_mass(ZAMS)

    #count up number of objects formed and ensure it matahces the number of stars input
    bad_idx = np.where(output_array[1] == -1)
    WD_idx = np.where(output_array[1] == 101)
    NS_idx = np.where(output_array[1] == 102)
    BH_idx = np.where(output_array[1] == 103)

    rem_mass = output_array[0]

    return bad_idx[0], WD_idx[0], NS_idx[0], BH_idx[0], rem_mass

def test_Raithel18_IFMR_1():
    """
    Check to make sure that the number of input objects matches the number of output objects
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Raithel18_IFMR()

    total = len(WD_idx) + len(NS_idx) + len(BH_idx) + len(bad_idx)
    assert total == 16 , "The # of input objects does not match the number of output objects for the Raithel18 IFMR"

    return

def test_Raithel18_IFMR_2():
    """
    Check that all negative remnant masses have the type code -1
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Raithel18_IFMR()

    #check to make sure no unhandled negative remnant masses (ie without type code -1 assigned)
    neg_mass_idx = np.where(rem_mass < 0)
    assert set(bad_idx) == set(neg_mass_idx[0]) , "There are unhandled negative remnant masses for the Raithel18 IFMR"

    return

def test_Raithel18_IFMR_3():
    """
    Check that there are no left over zeros in the remnant mass array
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Raithel18_IFMR()

    #check to make sure there are no left over zeroes in the remnant mass array
    assert np.all(rem_mass != 0), "There are left over zeros in the remnant mass array for the Raithel18 IFMR"

    return

def test_Raithel18_IFMR_4():
    """
    Check that the right number of invalid objects are returned (3)
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Raithel18_IFMR()

    assert len(bad_idx) == 3 , "There are not the right number of invalid objects for the Raithel18 IFMR"

    return

def test_Raithel18_IFMR_5():
    """
    Check that the right number of WDs are returned (2)
    """
    bad_idx, WD_idx, NS_idx, BH_idx, rem_mass = generate_Raithel18_IFMR()

    assert len(WD_idx) == 2 , "There are not the right number of WDs for the Raithel18 IFMR"

    return

def test_ResolvedCluster_random_state():
    """
    Test that the random state is properly set in ResolvedCluster, such that two clusters with the same seed have the same stars.
    """
    log_age = 6.7
    AKs = 2.7
    distance = 4000
    cluster_mass = 10**4.
    iso_dir = f'{spisea_path}/tests/isochrones'

    evo = evolution.MergedBaraffePisaEkstromParsec()
    atm_func = atmospheres.get_merged_atmosphere
    red_law = reddening.RedLawNishiyama09()
    filt_list = ['nirc2,J', 'nirc2,Kp']

    iso = syn.IsochronePhot(
        log_age,
        AKs,
        distance,
        evo_model=evo,
        atm_func=atm_func,
        red_law=red_law,
        filters=filt_list,
        mass_sampling=10,
        iso_dir=iso_dir,
        recomp=True,
        verbose=False
    )

    imf_limits = np.array([0.07, 0.5, 150])
    imf_powers = np.array([-1.3, -2.35])
    imf_multi = multiplicity.MultiplicityUnresolved()
    imf_test = imf.IMF_broken_powerlaw(imf_limits, imf_powers, multiplicity=imf_multi)

    # Test that the same random seed produces the same cluster
    imf_test.rng = np.random.default_rng(seed=42)
    result1 = imf_test.generate_cluster(cluster_mass)
    imf_test.rng = np.random.default_rng(seed=42)
    result2 = imf_test.generate_cluster(cluster_mass)
    np.testing.assert_equal(result1, result2)

    # Test that two clusters generated with the same seed have the same star systems and companions
    cluster1 = syn.ResolvedCluster(iso, imf_test, cluster_mass, seed=42)
    cluster2 = syn.ResolvedCluster(iso, imf_test, cluster_mass, seed=42)
    np.testing.assert_array_equal(cluster1.star_systems, cluster2.star_systems)

    old_star_systems = Table.read(f'{spisea_path}/tests/test_data/star_systems.fits')
    old_companion = Table.read(f'{spisea_path}/tests/test_data/companions.fits')

    #pdb.set_trace()
    for key in old_star_systems.colnames:
        #np.testing.assert_array_equal(cluster1.star_systems[key], old_star_systems[key])
        np.testing.assert_allclose(cluster1.star_systems[key], old_star_systems[key], rtol=5e-2, atol=5e-2)

    for key in old_companion.colnames:
        # Require values are consistent within reasonable bounds
        # np.testing.assert_array_equal(cluster1.companions[key], old_companion[key])
        np.testing.assert_allclose(cluster1.companions[key], old_companion[key], rtol=5e-2, atol=5e-2)

    return


def test_ResolvedCluster_no_companions():
    """
    Test case where no companions get generated to
    make sure we don't get any errors. This relies on using
    a specific seed that results in no companions being generated.
    """
    
    # Define cluster parameters
    logAge = 6.7
    AKs = 2.4
    distance = 4000
    cluster_mass = 10**2
    iso_dir = f'{spisea_path}/tests/isochrones'

    # Test filters
    filt_list = ['nirc2,J', 'nirc2,Kp']
    evo = evolution.MergedBaraffePisaEkstromParsec()
    atm_func = atmospheres.get_merged_atmosphere
    red_law = reddening.RedLawNishiyama09()

    iso = syn.IsochronePhot(logAge, AKs, distance, metallicity=0,
                            evo_model=evo, atm_func=atm_func,
                            red_law=red_law, filters=filt_list,
                                iso_dir=iso_dir)
                                
    imf_multi = multiplicity.MultiplicityResolvedDK()
    my_imf = imf.Kroupa_2001(multiplicity=imf_multi)
    
    cluster = syn.ResolvedCluster(iso, my_imf, cluster_mass,
                    keep_low_mass_stars=True, seed=1074)
                    
    assert(~np.any(cluster.star_systems['isMultiple']))


@pytest.mark.parametrize(
    "filt_str, expected_ab_minus_vega",
    [
        # Regression baselines (synphot 1.x + stsynphot): m_AB - m_Vega for Vega,
        # using SPISEA get_filter_info bandpasses and SourceSpectrum.from_vega().
        ("ubv,V", 0.016),
        ("ubv,B", -0.12),
        ("2mass,J", 0.913),
        ("wfc3,ir,f125w", 0.923)
    ],
)
def test_calc_ab_vega_filter_conversion_known(filt_str, expected_ab_minus_vega):
    """Vega AB minus Vega magnitude offset matches tabulated values for a few filters."""
    with patch("builtins.print"):
        out = syn.calc_ab_vega_filter_conversion(filt_str)
    mag = float(out.value) if hasattr(out, "value") else float(out)
    np.testing.assert_allclose(mag, expected_ab_minus_vega, rtol=1e-2, atol=1e-2)

    return

@pytest.mark.parametrize(
    "filt_str, expected_st_minus_vega",
    [
        # Regression baselines (synphot 1.x + stsynphot): m_ST - m_Vega for Vega,
        # using SPISEA get_filter_info bandpasses and SourceSpectrum.from_vega().
        ("ubv,B", -0.61),
        ("ubv,V", 0.019),
        ("2mass,J", 2.686),
        ("wfc3,ir,f125w", 2.703)
    ],
)
def test_calc_st_vega_filter_conversion_known(filt_str, expected_st_minus_vega):
    """Vega AB minus Vega magnitude offset matches tabulated values for a few filters."""
    with patch("builtins.print"):
        out = syn.calc_st_vega_filter_conversion(filt_str)
    mag = float(out.value) if hasattr(out, "value") else float(out)
    np.testing.assert_allclose(mag, expected_st_minus_vega, rtol=1e-2, atol=1e-2)

    return

def test_mag_in_filter():
    """Test that the mag_in_filter function works correctly."""
    old_stars = Table.read(f'{spisea_path}/tests/test_data/star_systems.fits')

    idx = np.where((old_stars['isMultiple'] == False) & (old_stars['isWR'] == False))[0]
    ii = idx[2]
    print(f'ii: {ii}')
    star = old_stars[ii]

    atm_func = atmospheres.get_merged_atmosphere
    red_law = reddening.RedLawNishiyama09()
    filt_list = ['nirc2,J', 'nirc2,Kp']
    wave_range = [3000, 52000]
    distance = 4000
    AKs = 2.7
    L = (star['L'] * units.Watt).to('erg/s')
    Teff = star['Teff'] * units.K
    R = np.sqrt(L / (4.0 * np.pi * c.sigma_sb * Teff**4)).to('pc')
    logg = star['logg']

    # Convert to expected units
    L = L.value
    Teff = Teff.value
    R = R.value

    spec = atm_func(temperature=Teff, gravity=logg)
    spec = syn.trim_spectrum(spec, wave_range[0], wave_range[1])
    spec *= (R / distance)**2  # in erg s^-1 cm^-2 A^-1
    spec *= red_law.extinction_curve(AKs, spec.waveset)

    filt_J = syn.get_filter_info(filt_list[0])
    filt_Kp = syn.get_filter_info(filt_list[1])

    mag_J = syn.mag_in_filter(spec, filt_J)
    mag_Kp = syn.mag_in_filter(spec, filt_Kp)

    print(mag_J, star['m_nirc2_J'])
    print(mag_Kp, star['m_nirc2_Kp'])

    assert np.isclose(mag_J, star['m_nirc2_J'], rtol=5e-2, atol=5e-2)
    assert np.isclose(mag_Kp, star['m_nirc2_Kp'], rtol=5e-2, atol=5e-2)

    return
