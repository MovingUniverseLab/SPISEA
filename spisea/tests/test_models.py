# Test functions for the different stellar evolution and atmosphere models
from spisea import evolution, atmospheres, synthetic
import numpy as np
import pdb

def test_evo_model_grid_num():
    """
    Make sure evolution models have both evo_grid_num
    and evo_grid_min (e.g., make sure these functions
    are working). Try it on one evolution model here;
    we'll test on all evo models in another function.
    """
    # Make MIST evolution model, check evo grid variables
    evo = evolution.MISTv1()
    assert isinstance(evo.evo_grid_min, float)

    return

def test_evolution_models():
    """
    Test to make sure the different evolution models work
    """
    # Age ranges to test
    age_young_arr = [6.7, 7.9]
    age_all_arr = [6.7, 8.0, 9.7]
    age_all_MIST_arr = [5.2, 6.7, 9.7, 10.13]
    bd_test = [6.0, 6.5, 7.4, 8.4, 10.0]

    # Metallicity ranges to test (if applicable)
    metal_range = [-2.5, -1.5, 0, 0.25, 0.4]
    metal_solar = [0]
    metal_Marley = [-0.5, 0.0, 0.5]

    # Array of evolution models to test
    evo_models = [evolution.MISTv1(version=1.2), evolution.MergedBaraffePisaEkstromParsec(), 
                      evolution.Parsec(), evolution.Baraffe15(), evolution.Ekstrom12(), evolution.Pisa(),
                      evolution.Phillips2020(), evolution.Marley2021(),
                      evolution.MergedPhillipsBaraffePisaEkstromParsec()]


    # Array of age_ranges for the specific evolution models to test
    age_vals = [age_all_MIST_arr, age_all_arr, age_all_arr, age_young_arr, age_young_arr, age_young_arr, age_all_arr, age_all_arr, bd_test]

    # Array of metallicities for the specific evolution models to test
    metal_vals = [metal_range, metal_solar, metal_solar, metal_solar, metal_solar, metal_solar, metal_solar, metal_Marley, metal_solar]

    assert len(evo_models) == len(age_vals) == len(metal_vals)

    # Loop through models, testing if them work
    for ii in range(len(evo_models)):
        evo = evo_models[ii]

        # Loop through metallicities
        for jj in metal_vals[ii]:
            # Loop through ages
            for kk in age_vals[ii]:
                try:
                    test = evo.isochrone(age=10**kk, metallicity=jj)

                    # Make sure the actual isochrone metallicity taken is
                    # indeed the closest to the desired metallicity (e.g., closest point
                    # in evo grid)
                    z_ratio = np.log10(np.array(evo.z_list) / evo.z_solar)
                    closest_idx = np.where( abs(z_ratio - jj) == min(abs(z_ratio - jj)) )[0][0]
                    expected_val = z_ratio[closest_idx]

                    assert np.isclose(test.meta['metallicity_act'], expected_val, atol=0.01)

                except:
                    raise Exception('EVO TEST FAILED: {0}, age = {1}, metal = {2}'.format(evo, kk, jj))

        print('Done {0}'.format(evo))

    return

def test_synthpop_MIST_extension():
    """
    Testing the synthpop MIST extension to consistently lower masses
    """
    evo1_grid = evolution.MISTv1(version=1.2, synthpop_extension=False)
    evo2_grid = evolution.MISTv1(version=1.2, synthpop_extension=True)

    # Extract same isochrone from these two models
    logAge = 10.12
    evo1 = evo1_grid.isochrone(10**logAge, metallicity=0)
    evo2 = evo2_grid.isochrone(10**logAge, metallicity=0)

    # I expect evo2 extends to low masses below evo1
    assert len(evo2) > len(evo1)
    assert np.min(evo2['mass']) < np.min(evo1['mass'])

    return

def test_COSMIC_init():
    """
    Test the COSMIC external evolution model constructor: default flags,
    default BSEDict, and that user-supplied options are stored.
    """
    # Default construction
    evo = evolution.COSMIC()
    assert evo.external_evol is True
    assert evo.z_solar == 0.02
    assert evo.model_version_name == 'COSMIC'
    assert evo.keep_disrupted_companions is True
    assert evo.keep_COSMIC_tables is False

    # Default BSEDict should be a populated dictionary of BSE parameters
    assert isinstance(evo.BSEDict, dict)
    assert len(evo.BSEDict) > 0

    # User-supplied options should be stored
    custom_dict = {'windflag': 3, 'neta': 0.5}
    evo2 = evolution.COSMIC(BSEDict=custom_dict, keep_disrupted_companions=False,
                            keep_COSMIC_tables=True)
    assert evo2.BSEDict == custom_dict
    assert evo2.keep_disrupted_companions is False
    assert evo2.keep_COSMIC_tables is True

    return

def test_COSMIC_calc_logg():
    """
    Test the COSMIC.calc_logg helper. For the Sun (M=1 Msun, R=1 Rsun)
    the surface gravity should be logg ~ 4.438 (cgs).
    """
    evo = evolution.COSMIC()

    # Scalar solar value
    assert np.isclose(evo.calc_logg(1.0, 1.0), 4.438, atol=0.01)

    # Array input should be handled element-wise
    masses = np.array([1.0, 2.0])
    radii = np.array([1.0, 2.0])
    logg = evo.calc_logg(masses, radii)
    assert np.isclose(logg[0], 4.438, atol=0.01)
    # logg scales as log10(M/R^2); doubling both M and R lowers logg by log10(2)
    assert np.isclose(logg[0] - logg[1], np.log10(2.0), atol=0.01)

    return

def test_COSMIC_get_kick_differential():
    """
    Test the COSMIC.get_kick_differential helper. The transformation is a
    pure rotation (Rz(theta) * Rx(phi)) of the kick vector, so it must
    preserve the vector magnitude. A zero kick must map to a zero kick.
    """
    evo = evolution.COSMIC()

    # Zero kick in -> zero kick out
    zeros = np.zeros(3)
    phase = np.array([0.3, 1.1, 2.0])
    incl = np.array([0.5, 1.5, 2.5])
    kd_zero = evo.get_kick_differential(zeros, zeros, zeros, phase=phase, inclination=incl)
    assert np.allclose(kd_zero.d_x.value, 0.0)
    assert np.allclose(kd_zero.d_y.value, 0.0)
    assert np.allclose(kd_zero.d_z.value, 0.0)

    # Magnitude is preserved under the rotation
    vx = np.array([10.0, -5.0, 3.0])
    vy = np.array([2.0, 7.0, -1.0])
    vz = np.array([-4.0, 1.0, 8.0])
    kd = evo.get_kick_differential(vx, vy, vz, phase=phase, inclination=incl)

    mag_in = np.sqrt(vx**2 + vy**2 + vz**2)
    mag_out = np.sqrt(kd.d_x.value**2 + kd.d_y.value**2 + kd.d_z.value**2)
    np.testing.assert_allclose(mag_out, mag_in, rtol=1e-10)

    return

def test_atmosphere_models():
    """
    Test the rebinned atmosphere models used for synthetic photometry
    """
    # Array of atmospheres
    atm_arr = [
        atmospheres.get_merged_atmosphere,
        atmospheres.get_castelli_atmosphere,
        atmospheres.get_phoenixv16_atmosphere,
        atmospheres.get_BTSettl_2015_atmosphere,
        atmospheres.get_BTSettl_atmosphere,
        atmospheres.get_kurucz_atmosphere,
        atmospheres.get_phoenix_atmosphere,
        atmospheres.get_wdKoester_atmosphere,
        atmospheres.get_Phillips2020_atmosphere,
        atmospheres.get_Meisner2023_atmosphere
    ]

    # Array of metallicities
    metals_range = [-2.0, 0, 0.15]
    bd_metals_range = [-1.0, -0.5, 0, 0.3]
    metals_solar = [0]
    metals_arr = [metals_solar, metals_range, metals_range, metals_solar, metals_range, metals_range, metals_range,
                  metals_solar, metals_solar, bd_metals_range]

    assert len(atm_arr) == len(metals_arr)

    # Loop through models, testing if them work
    for ii in range(len(atm_arr)):
        atm_func = atm_arr[ii]

        # Loop through metallicities
        for jj in metals_arr[ii]:
            try:
                test = atm_func(metallicity=jj)
            except:
                raise Exception('ATM TEST FAILED: {0}, metal = {1}'.format(atm_func, jj))

        print('Done {0}'.format(atm_func))

    # Test get_merged_atmospheres at different temps
    temp_range = [250, 1000, 2000, 3500, 4000, 5250, 6000, 12000]
    atm_func = atmospheres.get_merged_atmosphere
    for ii in bd_metals_range:
        for jj in temp_range:
            try:
                test = atm_func(metallicity=ii, temperature=jj, verbose=True)
            except:
                raise Exception('ATM TEST FAILED: {0}, metal = {1}, temp = {2}'.format(atm_func, ii, jj))


    print('get_merged_atmosphere: all temps/metallicities passed')

    # Test get_bb_atmosphere at different temps
    # This func only requests temp
    temp_range = [1000, 2000, 3500, 4000, 5250, 6000, 12000]
    atm_func = atmospheres.get_bb_atmosphere
    for jj in temp_range:
        try:
            test = atm_func(temperature=jj, verbose=True)
        except:
            raise Exception('ATM TEST FAILED: {0}, temp = {2}'.format(atm_func, jj))

    print('get_bb_atmosphere: all temps passed')

    # Test get_bd_atmosphere at different temps
    # This func only requests temp
    temp_range = [250, 400, 500, 750, 950, 1200]
    atm_func = atmospheres.get_bd_atmosphere
    for jj in temp_range:
        try:
            test = atm_func(temperature=jj, verbose=True)
        except:
            raise Exception('ATM TEST FAILED: {0}, temp = {1}'.format(atm_func, jj))

    print('get_bd_atmosphere: all temps passed')

    return

def test_filters():
    """
    Test to make sure all of the filters work as expected
    """
    # Define vega spectrum
    vega = synthetic.Vega()

    # Filter list to test
    filt_list = ['wfc3,ir,f127m','acs,wfc1,f814w',
                     '2mass,J', '2mass,H','2mass,Ks',
                     'ctio_osiris,K', 'ctio_osiris,H',
                     'ubv,U', 'ubv,B', 'ubv,V', 'ubv,R',
                     'ubv,I', 'jg,J', 'jg,H', 'jg,K',
                     'decam,Y', 'decam,i', 'decam,z',
                     'decam,u', 'decam,g', 'decam,r',
                     'gaia,dr1,G', 'gaia,dr1,Gbp', 'gaia,dr1,Grp',
                     'gaia,dr2,G', 'gaia,dr2,Gbp', 'gaia,dr2,Grp',
                     'gaia,dr2_rev,G', 'gaia,dr2_rev,Gbp', 'gaia,dr2_rev,Grp',
                     'gaia,edr3,G', 'gaia,edr3,Gbp', 'gaia,edr3,Grp',
                     'jwst,F070W', 'jwst,F090W', 'jwst,F115W', 'jwst,F140M',
                     'jwst,F150W', 'jwst,F150W2', 'jwst,F162M', 'jwst,F164N',
                     'jwst,F182M', 'jwst,F187N', 'jwst,F200W', 'jwst,F212N',
                     'jwst,F210M','jwst,F250M', 'jwst,F277W', 'jwst,F300M',
                     'jwst,F322W2', 'jwst,F323N', 'jwst,F335M', 'jwst,F356W',
                     'jwst,F360M', 'jwst,F405N', 'jwst,F410M', 'jwst,F430M',
                     'jwst,F444W', 'jwst,F460M', 'jwst,F466N', 'jwst,F470N',
                     'jwst,F480M', 'naco,J', 'naco,H', 'naco,Ks',
                     'naco,IB_2.00','naco,IB_2.03', 'naco,IB_2.06', 'naco,IB_2.24',
                     'naco,IB_2.27', 'naco,IB_2.30', 'naco,IB_2.33', 'naco,IB_2.36',
                     'nirc1,K', 'nirc1,H', 'nirc2,J', 'nirc2,H',
                     'nirc2,Kp', 'nirc2,K', 'nirc2,Lp', 'nirc2,Hcont',
                     'nirc2,FeII', 'nirc2,Brgamma', 'ps1,z',
                     'ps1,g', 'ps1,r','ps1,i', 'ps1,y', 'ps1,w',
                     'ukirt,Z','ukirt,Y','ukirt,J', 'ukirt,H', 'ukirt,K',
                     'vista,Y', 'vista,Z', 'vista,J',
                     'vista,H',  'vista,Ks', 'ztf,g', 'ztf,r', 'ztf,i',
                     'hawki,J', 'hawki,H', 'hawki,Ks', 'roman,wfi,f062',
                     'roman,wfi,f087', 'roman,wfi,f106', 'roman,wfi,f129',
                     'roman,wfi,f158', 'roman,wfi,f146', 'roman,wfi,f213',
                     'roman,wfi,f184', 'rubin,g', 'rubin,i', 'rubin,r',
                     'rubin,u', 'rubin,z', 'rubin,y',
                     'euclid,VIS', 'euclid,Y', 'euclid,J', 'euclid,H',
                     'nsfcam,L', 'tess,tess',
                     'washington,C', 'washington,M', 'washington,T1', 'washington,T2',
                     'hipparcos,Hp', 'tycho,B', 'tycho,V',
                     'kepler,Kp', 'ogle,Rw', 'ogle,I', 'ogle,V',
                     'subaru,hsc,g','subaru,hsc,r','subaru,hsc,i','subaru,hsc,z','subaru,hsc,Y',
                     'subaru,hsc,nb387', 'subaru,hsc,nb468', 'subaru,hsc,nb515', 'subaru,hsc,nb527',
                     'subaru,hsc,nb656', 'subaru,hsc,nb718', 'subaru,hsc,nb816', 'subaru,hsc,nb921',
                     'subaru,hsc,nb926', 'subaru,hsc,nb973',
                     'bessell,U', 'bessell,B', 'bessell,V', 'bessell,R', 'bessell,I']

    # Loop through filters to test that they work: get_filter_info
    for ii in filt_list:
        filt = synthetic.get_filter_info(ii, rebin=True, vega=vega)

    print('get_filter_info pass')

    # Loop through filters to test that they work: get_obs_str
    for ii in filt_list:
        # Test going from col_name to obs_str
        col_name = synthetic.get_filter_col_name(ii)
        obs_str = synthetic.get_obs_str('m_{0}'.format(col_name))
        # Does the obs_str work?
        filt_info = synthetic.get_filter_info(obs_str)

    print('get_obs_str pass')
    print('Filters done')

    return
