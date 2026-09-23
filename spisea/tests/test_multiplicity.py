import numpy as np
import time
import os
import spisea
from spisea.imf import imf, multiplicity

def test_create_MultiplicityUnresolved():
    """
    Tests creating and accessing a MultiplicityUnresolved object.
    """
    # All default parameters -- check their values
    mu1 = multiplicity.MultiplicityUnresolved()
    assert mu1.MF_amp == 0.44
    assert mu1.MF_pow == 0.51
    assert mu1.CSF_amp == 0.50
    assert mu1.CSF_pow == 0.45
    assert mu1.CSF_max == 3
    assert mu1.q_pow == -0.4
    assert mu1.q_min == 0.01

    # Test setting different parameters
    mu2 = multiplicity.MultiplicityUnresolved(MF_amp=0.4,
                                             MF_power=0.4,
                                             CSF_amp=0.4,
                                             CSF_power=0.4,
                                             CSF_max=4,
                                             q_power=0.4,
                                             q_min=0.04)
    assert mu2.MF_amp == 0.4
    assert mu2.MF_pow == 0.4
    assert mu2.CSF_amp == 0.4
    assert mu2.CSF_pow == 0.4
    assert mu2.CSF_max == 4
    assert mu2.q_pow == 0.4
    assert mu2.q_min == 0.04


def test_multiplicity_fraction():
    """
    Test creating a MultiplicityUnresolved object and getting
    the multiplicity fraction out.
    """
    # First set of multiplicity parameters
    mu1 = multiplicity.MultiplicityUnresolved()

    mf1_1 = mu1.multiplicity_fraction(1.0)
    np.testing.assert_almost_equal(mf1_1, 0.44, decimal=2)

    mf1_2 = mu1.multiplicity_fraction(10.0)
    np.testing.assert_almost_equal(mf1_2, 1.0, decimal=2)

    mf1_3 = mu1.multiplicity_fraction(0.1)
    np.testing.assert_almost_equal(mf1_3, 0.136, decimal=2)

    # Second set of multiplicity parameters
    mu2 = multiplicity.MultiplicityUnresolved(MF_amp=0.4, MF_power=0.4,
                                             CSF_amp=0.4, CSF_power=0.4, CSF_max=4,
                                             q_power=0.4, q_min=0.04)

    mf2_1 = mu2.multiplicity_fraction(1.0)
    np.testing.assert_almost_equal(mf2_1, 0.4, decimal=2)

    mf2_2 = mu2.multiplicity_fraction(10.0)
    np.testing.assert_almost_equal(mf2_2, 1.0, decimal=2)

    mf2_3 = mu2.multiplicity_fraction(0.1)
    np.testing.assert_almost_equal(mf2_3, 0.159, decimal=2)

    # Test brown dwarf mass fractions
    mf_bd1 = mu1.multiplicity_fraction(0.07)  # near upper BD limit
    mf_bd2 = mu1.multiplicity_fraction(0.04)  # mid BD
    mf_bd3 = mu1.multiplicity_fraction(0.01)  # lower BD limit
    assert np.isclose(mf_bd1, 0.16, atol=0.01)
    assert np.isclose(mf_bd2, 0.08, atol=0.01)
    assert np.isclose(mf_bd3, 0.0, atol=1e-6)


def test_multiplicity_fraction_array():
    """
    Test multiplicity_fraction() on the MultiplicityUnresolved object
    where the inputs and outputs are in array form.
    """
    # First set of multiplicity parameters
    mu1 = multiplicity.MultiplicityUnresolved()

    mass_array = np.array([1.0, 10.0, 0.1, 0.07, 0.04, 0.01])
    mf_array = mu1.multiplicity_fraction(mass_array)

    # Stellar regime checks
    np.testing.assert_almost_equal(mf_array[0], 0.44, decimal=2)
    np.testing.assert_almost_equal(mf_array[1], 1.0, decimal=2)
    np.testing.assert_almost_equal(mf_array[2], 0.136, decimal=2)

    # BD regime checks
    # interpolation between values implies lower masses --> lower mf
    assert mf_array[3] < mf_array[2]
    assert mf_array[4] <= mf_array[3]
    assert mf_array[5] <= mf_array[4]

    # Ensure mf stars within reasonable bound (upper limit is 0.2)
    assert np.all(mf_array[3:] >= 0.0)
    assert np.all(mf_array[3:] <= 0.2)

    
def test_companion_star_fraction():
    """
    Test the companion_star fraction on the MultiplicityUnresolved object.
    """
    # First set of multiplicity parameters
    mu1 = multiplicity.MultiplicityUnresolved()

    csf1_1 = mu1.companion_star_fraction(1.0)
    np.testing.assert_almost_equal(csf1_1, 0.5, decimal=2)

    csf1_2 = mu1.companion_star_fraction(70.0)
    np.testing.assert_almost_equal(csf1_2, 3.0, decimal=2)

    csf1_3 = mu1.companion_star_fraction(0.1)
    np.testing.assert_almost_equal(csf1_3, 0.177, decimal=2)

    # Second set of multiplicity parameters
    mu2 = multiplicity.MultiplicityUnresolved(MF_amp=0.4, MF_power=0.4,
                                              CSF_amp=0.4, CSF_power=0.4, CSF_max=2,
                                              q_power=0.4, q_min=0.04)

    # csf2_1 = mu1.companion_star_fraction(1.0)
    # np.testing.assert_almost_equal(csf2_1, 0.4, decimal=2)

    # csf2_2 = mu1.companion_star_fraction(70.0)
    # np.testing.assert_almost_equal(csf2_2, 2.0, decimal=2)

    # csf2_3 = mu1.companion_star_fraction(0.1)
    # np.testing.assert_almost_equal(csf2_3, 0.159, decimal=2)

    # BD CSF follows the stellar power law (no CSF=MF mass cut).
    csf_bd1 = mu1.companion_star_fraction(0.07)
    csf_bd2 = mu1.companion_star_fraction(0.04)
    csf_bd3 = mu1.companion_star_fraction(0.01)
    np.testing.assert_almost_equal(csf_bd1, 0.50 * 0.07 ** 0.45, decimal=4)
    np.testing.assert_almost_equal(csf_bd2, 0.50 * 0.04 ** 0.45, decimal=4)
    np.testing.assert_almost_equal(csf_bd3, 0.50 * 0.01 ** 0.45, decimal=4)


def test_resolvedmult():
    """
    Test creating a MultiplicityResolvedDK object
    and that the parameters it's populated with are correct.
    Updated to test for specific brown dwarf characteristics.
    """
    from spisea import synthetic, evolution, atmospheres, reddening, ifmr
    # Fetch isochrone
    logAge = 6.70 # Age in log(years)
    AKs = 1.0 # Ks filter extinction in mags
    dist = 4000 # distance in parsecs
    metallicity = 0 # metallicity in [M/H]
    atm_func = atmospheres.get_merged_atmosphere
    evo_merged = evolution.MergedPhillipsBaraffePisaEkstromParsec()
    redlaw = reddening.RedLawCardelli(3.1) # Rv = 3.1
    filt_list = ['nirc2,J', 'nirc2,Kp']

    startTime = time.time()

    iso_merged = synthetic.IsochronePhot(logAge, AKs, dist, metallicity=metallicity,
                                 evo_model=evo_merged, atm_func=atm_func,
                                 filters=filt_list, red_law=redlaw,
                                 mass_sampling=3, iso_dir=f'{spisea.__path__[0]}/tests/isochrones')
    print('Constructed isochrone: %d seconds' % (time.time() - startTime))

    # Now we can make the cluster.
    clust_mtot = 10**4.
    clust_multiplicity = multiplicity.MultiplicityResolvedDK()

    # Multiplicity is defined in the IMF object
    clust_imf_Mult = imf.Salpeter_Kirkpatrick_2024(multiplicity=clust_multiplicity)
    
    # Make clusters
    clust_Mult = synthetic.ResolvedCluster(iso_merged, clust_imf_Mult, clust_mtot)

    clust_Mult_ss = clust_Mult.star_systems

    print('Constructed cluster: %d seconds' % (time.time() - startTime))

    #check if columns were created
    assert 'log_a' in clust_Mult.companions.colnames
    assert 'e' in clust_Mult.companions.colnames
    assert 'i' in clust_Mult.companions.colnames
    assert 'Omega' in clust_Mult.companions.colnames
    assert 'omega' in clust_Mult.companions.colnames

    #check values are in correct range
    assert all(10**i<= 2000 and 10**i>= 0 for i in clust_Mult.companions['log_a']) #max separation is 2000 AU
    assert all(i<= 1 and i>= 0 for i in clust_Mult.companions['e'])
    assert all(i<= 180 and i>= 0 for i in clust_Mult.companions['i'])
    assert all(i<= 360 and i>= 0 for i in clust_Mult.companions['omega'])
    assert all(i<= 360 and i>= 0 for i in clust_Mult.companions['Omega'])

    #checks sign for inclination is being randomly genarated
    assert any(i > 90 for i in clust_Mult.companions['i']) and any(i < 90 for i in clust_Mult.companions['i'])

    #checks eccentricity follows f(e) = 2e pdf
    n, bins = np.histogram(clust_Mult.companions['e'], density = True)
    bin_centers = 0.5*(bins[1:] + bins[:-1])
    assert all(np.abs(i) < 0.3 for i in 2*bin_centers - n)

    #checks shape of inclination histogram is sin(i)
    n, bins = np.histogram(clust_Mult.companions['i'])
    bin_centers = 0.5*(bins[1:] + bins[:-1])
    assert all(np.abs(i) < 0.15 for i in n/max(n) - np.sin(np.pi*bin_centers/180))

    #checks for brown dwarf specific features
    bd_idx = np.where(clust_Mult.star_systems['mass'] < 0.08)[0]

    #check there is only one possible companion per BD
    assert all(clust_Mult.star_systems['N_companions'][bd_idx] <= 1), \
    "Brown dwarf primaries have >1 companion."

    comp_rows = []
    start = 0
    for ii, N in enumerate(clust_Mult.star_systems['N_companions']):
        if ii in bd_idx and N > 0:
            comp_rows.extend(range(start, start+N))
        start += N

    bd_companions = clust_Mult.companions[comp_rows]

    if len(bd_companions) > 30:  # only test if enough BD binaries
        mean_log_a = np.mean(bd_companions['log_a'])
        std_log_a = np.std(bd_companions['log_a'])

        bd_masses = clust_Mult.star_systems['mass'][bd_idx]
        expected_sigma = np.mean(
            np.interp(
                np.log10(bd_masses),
                [np.log10(0.01), np.log10(0.08)],
                [0.25, 0.5]
            )
        )

        #expect lognormal centered near log10(2.9 AU), width ~0.21
        assert abs(mean_log_a - np.log10(2.9)) < 0.25, \
            f"BD mean log(a) off: {mean_log_a:.2f}"

        assert abs(std_log_a - expected_sigma) < 0.15, \
            f"BD sigma log(a) off: {std_log_a:.2f}"

    return


# ---------------------------------------------------------------------------
# Offner et al. 2023 (Table 1) multiplicity
# ---------------------------------------------------------------------------

# Published Table 1 MF (%) converted to fraction, CF, and 1-sigma MF error.
# Masses are geometric means of the tabulated M1 intervals.
_OFFNER_TABLE1 = [
    # name, M_lo, M_hi, MF, MF_err, CF
    ('Fontanive+2018', 0.019, 0.058, 0.08, 0.06, 0.08),
    ('Burgasser 2007', 0.05, 0.08, 0.15, 0.04, 0.16),
    ('Close+2003', 0.080, 0.095, 0.19, 0.07, 0.19),
    ('Allen+2007', 0.06, 0.15, 0.20, 0.04, 0.20),
    ('Winters+2019 late-M', 0.075, 0.15, 0.19, 0.03, 0.21),
    ('Winters+2019 mid-M', 0.15, 0.30, 0.23, 0.02, 0.27),
    ('Winters+2019 early-M', 0.3, 0.6, 0.30, 0.02, 0.38),
    ('Raghavan+2010', 0.75, 1.25, 0.46, 0.03, 0.60),
    ('Tokovinin 2014b', 0.85, 1.5, 0.47, 0.03, 0.62),
    ('Moe & Kratter 2021', 1.6, 2.4, 0.68, 0.07, 0.99),
    ('Moe & Di Stefano 2017 3-5', 3.0, 5.0, 0.81, 0.06, 1.28),
    ('Moe & Di Stefano 2017 5-8', 5.0, 8.0, 0.89, 0.05, 1.55),
    ('Moe & Di Stefano 2017 8-17', 8.0, 17.0, 0.93, 0.04, 1.80),
    ('Sana et al. 17-50', 17.0, 50.0, 0.96, 0.04, 2.10),
]


def _geom_mean_mass(row):
    """
    Geometric-mean primary mass of a Table 1 M1 interval.

    Parameters
    ----------
    row : sequence
        Table 1 row whose indices 1 and 2 are the low and high
        edges of the primary-mass interval, in solar masses (Msun).

    Returns
    -------
    mass : float
        Geometric mean ``sqrt(M_lo * M_hi)`` in solar masses (Msun).
    """
    return np.sqrt(row[1] * row[2])


def test_piecewise_powerlaw_api():
    """Custom piecewise MF/CSF is vectorized and clips MF to [0, 1]."""
    mass_limits = np.array([0.1, 1.0, 10.0])
    # First segment: MF = 0.4 * M^0 → 0.4; second: 0.4 * M^1 so MF(10)=4 → clip 1
    mp = multiplicity.MultiplicityPiecewisePowerLaw(
        mass_limits,
        MF_amps=[0.4, 0.4], MF_powers=[0.0, 1.0],
        CSF_amps=[0.4, 0.5], CSF_powers=[0.0, 0.5])
    assert mp.multiplicity_fraction(0.2) == 0.4
    assert mp.multiplicity_fraction(1.0) == 0.4
    np.testing.assert_almost_equal(mp.multiplicity_fraction(10.0), 1.0)
    masses = np.array([0.2, 1.0, 10.0])
    mf = mp.multiplicity_fraction(masses)
    np.testing.assert_allclose(mf, [mp.multiplicity_fraction(m) for m in masses])


def test_logistic_api():
    """Custom logistic MF/CSF clips, vectorizes, and keeps CSF >= MF."""
    ml = multiplicity.MultiplicityLogistic(
        MF_A=0.1, MF_B=1.5, MF_M0=1.0, MF_k=2.0,
        CSF_A=0.2, CSF_B=4.0, CSF_M0=2.0, CSF_k=1.0,
        CSF_max=2.0)
    # Low-mass asymptote A for a very low-mass primary (not a missing mass)
    np.testing.assert_almost_equal(ml.multiplicity_fraction(1e-8), 0.1, decimal=4)
    # High-mass MF saturates at B then clips to 1
    np.testing.assert_almost_equal(ml.multiplicity_fraction(1e6), 1.0)
    # High-mass CSF clips to CSF_max
    np.testing.assert_almost_equal(ml.companion_star_fraction(1e6), 2.0)
    # No CSF=MF mass cut: low-mass CSF can exceed MF.
    assert ml.companion_star_fraction(0.05) > ml.multiplicity_fraction(0.05)
    masses = np.array([0.05, 1.0, 100.0])
    mf = ml.multiplicity_fraction(masses)
    np.testing.assert_allclose(mf, [ml.multiplicity_fraction(m) for m in masses])
    csf = ml.companion_star_fraction(masses)
    np.testing.assert_allclose(
        csf, [ml.companion_star_fraction(m) for m in masses])
    assert np.all(csf >= mf - 1e-12)


def test_offner2023_logistic_coefficients():
    """Offner stores the equal-weight logistic-in-log-mass coefficients."""
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    assert isinstance(multi, multiplicity.MultiplicityLogistic)
    assert not isinstance(multi, multiplicity.MultiplicityPiecewisePowerLaw)
    np.testing.assert_allclose(multi.MF_A, 0.14)
    np.testing.assert_allclose(multi.MF_B, 0.99)
    np.testing.assert_allclose(multi.MF_M0, 1.41)
    np.testing.assert_allclose(multi.MF_k, 1.25)
    np.testing.assert_allclose(multi.CSF_A, 0.12)
    np.testing.assert_allclose(multi.CSF_B, 2.35)
    np.testing.assert_allclose(multi.CSF_M0, 3.57)
    np.testing.assert_allclose(multi.CSF_k, 0.96)
    np.testing.assert_allclose(multi.q_A, 6.6)
    np.testing.assert_allclose(multi.q_B, -1.77)
    np.testing.assert_allclose(multi.q_M0, 0.0651)
    np.testing.assert_allclose(multi.q_k, 0.629)
    assert not hasattr(multi, 'a_mup')
    assert not hasattr(multi, 'sig_A')
    assert not hasattr(multi, 'log_a_mean')
    assert not hasattr(multi, 'a_mean')
    assert not hasattr(multi, 'sigma_log_a')
    import inspect
    unres_sig = inspect.signature(
        multiplicity.MultiplicityUnresolvedOffner2023.__init__)
    for name in ('a_mup', 'a_mp', 'a_alphaL', 'a_alphaR', 'a_s', 'a_min',
                 'sig_A', 'sig_B', 'sig_M0', 'sig_k'):
        assert name not in unres_sig.parameters
    assert not hasattr(multiplicity, 'FONTANIVE2018_BD_Q_POWER')
    assert not any(n.startswith('OFFNER2023_') for n in dir(multiplicity))
    resolved = multiplicity.MultiplicityResolvedOffner2023()
    np.testing.assert_allclose(resolved.MF_A, 0.14)
    np.testing.assert_allclose(resolved.a_mup, 44.46)
    np.testing.assert_allclose(resolved.a_mp, 0.819)
    np.testing.assert_allclose(resolved.a_alphaL, 1.005)
    np.testing.assert_allclose(resolved.a_alphaR, -0.308)
    np.testing.assert_allclose(resolved.a_s, 0.10)
    np.testing.assert_allclose(resolved.a_min, 0.1)
    np.testing.assert_allclose(resolved.sig_A, 0.7)
    np.testing.assert_allclose(resolved.sig_B, 1.5)
    np.testing.assert_allclose(resolved.sig_M0, 0.354)
    np.testing.assert_allclose(resolved.sig_k, 6.05)
    res_sig = inspect.signature(
        multiplicity.MultiplicityResolvedOffner2023.__init__)
    assert res_sig.parameters['MF_A'].default == 0.14
    assert res_sig.parameters['sig_k'].default == 6.05


def test_offner2023_mf_smooth():
    """MF is continuous and nearly C1 around 0.08 and 1.5 Msun."""
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    eps = 1e-8
    for m in (0.08, 1.5):
        mf_left = multi.multiplicity_fraction(m - eps)
        mf_right = multi.multiplicity_fraction(m + eps)
        mf_at = multi.multiplicity_fraction(m)
        np.testing.assert_allclose(mf_left, mf_right, atol=1e-6, rtol=0)
        np.testing.assert_allclose(mf_at, mf_right, atol=1e-6, rtol=0)
        d_left = (mf_at - mf_left) / eps
        d_right = (mf_right - mf_at) / eps
        np.testing.assert_allclose(d_left, d_right, atol=1e-3, rtol=0)
    for m in (0.04, 0.3, 1.0, 10.0):
        expected = multiplicity._logistic_in_logm(
            m, multi.MF_A, multi.MF_B, multi.MF_M0, multi.MF_k,
            clip_min=0.0, clip_max=1.0)
        np.testing.assert_allclose(multi.multiplicity_fraction(m), expected)


def test_offner2023_table1_mf():
    """
    Logistic MF matches Offner et al. 2023 Table 1 at geom-mean M1.

    Fontanive (8±6%) sits ~0.07 below the curve (~15%); other rows,
    including A/B stars, stay close.
    """
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    for row in _OFFNER_TABLE1:
        name, mlo, mhi, mf_tab, mf_err, cf_tab = row
        m = _geom_mean_mass(row)
        mf = multi.multiplicity_fraction(m)
        tol = max(0.08, 2.0 * mf_err)
        assert abs(mf - mf_tab) <= tol, \
            '{0}: MF({1:.3f})={2:.3f} vs Table 1 {3:.2f} ± {4:.2f}'.format(
                name, m, mf, mf_tab, mf_err)
        assert 0.0 <= mf <= 1.0


def test_offner2023_table1_csf():
    """CSF tracks Table 1 CF; CSF >= MF at all masses (no BD CSF=MF cut)."""
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    for row in _OFFNER_TABLE1:
        name, mlo, mhi, mf_tab, mf_err, cf_tab = row
        m = _geom_mean_mass(row)
        csf = multi.companion_star_fraction(m)
        mf = multi.multiplicity_fraction(m)
        if mlo >= 1.6:
            # Logistic CF tracks A/B; Moe & Kratter residual ~0.1 is ok
            tol = max(0.12, 0.12 * cf_tab)
        else:
            tol = max(0.08, 0.25 * cf_tab)
        assert abs(csf - cf_tab) <= tol, \
            '{0}: CSF({1:.3f})={2:.3f} vs Table 1 CF {3:.2f}'.format(
                name, m, csf, cf_tab)
        assert csf >= mf - 1e-12
        assert csf <= multi.CSF_max + 1e-12


def test_offner2023_array_vs_scalar():
    """Array and scalar MF/CSF evaluations agree."""
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    masses = np.array([_geom_mean_mass(row) for row in _OFFNER_TABLE1])
    mf_arr = multi.multiplicity_fraction(masses)
    csf_arr = multi.companion_star_fraction(masses)
    for i, m in enumerate(masses):
        np.testing.assert_allclose(mf_arr[i], multi.multiplicity_fraction(float(m)))
        np.testing.assert_allclose(csf_arr[i], multi.companion_star_fraction(float(m)))


def test_higher_order_multiples_at_all_masses():
    """
    BD and stellar primaries can have n_comp > 1 when CSF/MF allows.

    ``draw_n_companions`` is for systems already drawn as multiple,
    so n_comp >= 1. Single brown dwarfs are allowed via the MF
    coin-flip (``random_is_multiple``), not this helper.
    """
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    assert multi.multiplicity_fraction(0.04) < 1.0
    rng = np.random.default_rng(123)
    for m in (0.04, 1.0):
        masses = np.full(4000, m)
        mf = multi.multiplicity_fraction(masses)
        csf = np.full(len(masses), 2.5)
        n_comp = multi.draw_n_companions(masses, csf, mf, rng)
        assert np.any(n_comp > 1)
        assert np.all(n_comp >= 1)


def test_companion_max_caps_all_masses():
    """companion_max=True clips n_comp at CSF_max for BD and stellar."""
    multi = multiplicity.MultiplicityUnresolvedOffner2023(
        companion_max=True, CSF_max=1)
    rng = np.random.default_rng(0)
    for m in (0.04, 1.0):
        masses = np.full(2000, m)
        mf = multi.multiplicity_fraction(masses)
        csf = np.full(len(masses), 3.0)
        n_comp = multi.draw_n_companions(masses, csf, mf, rng)
        assert np.all(n_comp <= 1)
        assert np.all(n_comp >= 1)


def test_offner2023_q_more_equal_mass_for_bds():
    """BD mass ratios are more equal-mass (higher mean q) than solar-type."""
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    rng = np.random.default_rng(7)
    n = 20000
    q_bd = multi.draw_q(np.full(n, 0.04), rng=rng)
    q_sun = multi.draw_q(np.full(n, 1.0), rng=rng)
    assert np.mean(q_bd) > np.mean(q_sun) + 0.1
    # Err-wt logistic undershoots Fontanive 4.8 (~3.3 at 0.033 Msun)
    assert multi.q_power_at_mass(0.033) > 2.5
    assert multi.q_power_at_mass(1.0) < 0.5


def test_offner2023_bd_separations_peak_few_au():
    """BD lognormal separations peak at a few AU (μ(0.033)≈2.1 au)."""
    multi = multiplicity.MultiplicityResolvedOffner2023()
    rng = np.random.default_rng(0)
    log_a = multi.log_semimajoraxis(np.full(5000, 0.04), rng=rng)
    med_a = 10 ** np.median(log_a)
    assert 1.5 < med_a < 8.0, 'BD median a = {0:.2f} AU'.format(med_a)
    # Solar-type should be much wider (smooth-broken μ ~ 44 au)
    log_a_s = multi.log_semimajoraxis(np.full(5000, 1.0), rng=rng)
    med_a_s = 10 ** np.median(log_a_s)
    assert med_a_s > 10.0
    assert med_a_s > med_a


def test_lu2013_defaults_unchanged():
    """SPISEA v2.5 MultiplicityUnresolved defaults and stellar MF unchanged."""
    mu = multiplicity.MultiplicityUnresolved()
    assert mu.MF_amp == 0.44
    assert mu.MF_pow == 0.51
    assert mu.CSF_amp == 0.50
    assert mu.CSF_pow == 0.45
    np.testing.assert_almost_equal(mu.multiplicity_fraction(1.0), 0.44, decimal=2)
    np.testing.assert_almost_equal(mu.multiplicity_fraction(10.0), 1.0, decimal=2)
    np.testing.assert_almost_equal(mu.multiplicity_fraction(0.1), 0.136, decimal=2)
    assert mu.bd_q_power == 6.1
    assert np.isclose(mu.q_power_at_mass(0.04), 6.1)
    assert np.isclose(mu.q_power_at_mass(1.0), -0.4)
    # Scalar BD overrides (SPISEA v2.5 / Fontanive path)
    assert np.isclose(mu.multiplicity_fraction(0.07), 0.16, atol=0.01)
    assert np.isclose(mu.multiplicity_fraction(0.04), 0.08, atol=0.01)
    assert np.isclose(mu.multiplicity_fraction(0.01), 0.0, atol=1e-6)


def test_resolveddk_log_a_mean_and_sigma():
    """
    MultiplicityResolvedDK.a_mean / log_a_mean / sigma_log_a match the
    Duchêne & Kraus + Fontanive BD interp + sigmoid that
    log_semimajoraxis used to inline, and the draw uses those methods.
    """
    import astropy.modeling
    import inspect

    dk = multiplicity.MultiplicityResolvedDK()
    masses = np.array([0.01, 0.04, 0.08, 0.3, 1.0, 5.0, 50.0])

    def expected_log_a_mean(mass):
        mass = np.atleast_1d(np.asarray(mass, dtype=float))
        logm = np.log10(mass)
        a_mean_func = astropy.modeling.powerlaws.BrokenPowerLaw1D(
            amplitude=dk.a_amp, x_break=dk.a_break,
            alpha_1=dk.a_slope1, alpha_2=dk.a_slope2)
        log_a_mean_star = np.log10(a_mean_func(mass))
        log_a_mean_bd = np.interp(
            logm,
            [np.log10(0.01), np.log10(0.08)],
            [np.log10(2.5), np.log10(8.0)])
        w = 1.0 / (1.0 + np.exp(-(logm - np.log10(0.08)) / 0.15))
        return (1 - w) * log_a_mean_bd + w * log_a_mean_star

    def expected_sigma(mass):
        mass = np.atleast_1d(np.asarray(mass, dtype=float))
        logm = np.log10(mass)
        log_a_std_func = astropy.modeling.models.Linear1D(
            slope=dk.a_std_slope, intercept=dk.a_std_intercept)
        log_a_std_star = log_a_std_func(logm)
        log_a_std_star[mass >= 2.9] = log_a_std_func(np.log10(2.9))
        log_a_std_star = np.clip(log_a_std_star, 0.1, None)
        log_a_std_bd = np.interp(
            logm,
            [np.log10(0.01), np.log10(0.08)],
            [0.25, 0.5])
        w = 1.0 / (1.0 + np.exp(-(logm - np.log10(0.08)) / 0.15))
        return (1 - w) * log_a_std_bd + w * log_a_std_star

    np.testing.assert_allclose(dk.log_a_mean(masses), expected_log_a_mean(masses))
    np.testing.assert_allclose(dk.sigma_log_a(masses), expected_sigma(masses))
    np.testing.assert_allclose(dk.a_mean(masses), 10.0 ** dk.log_a_mean(masses))

    for m in masses:
        log_a = dk.log_a_mean(float(m))
        sig = dk.sigma_log_a(float(m))
        a_mean = dk.a_mean(float(m))
        assert isinstance(log_a, float)
        assert isinstance(sig, float)
        assert isinstance(a_mean, float)
        np.testing.assert_allclose(log_a, expected_log_a_mean(m)[0])
        np.testing.assert_allclose(sig, expected_sigma(m)[0])
        np.testing.assert_allclose(a_mean, 10.0 ** log_a)

    src = inspect.getsource(dk.log_semimajoraxis)
    assert 'self.log_a_mean' in src
    assert 'self.sigma_log_a' in src
    assert 'BrokenPowerLaw1D' not in src

    log_a_draw = dk.log_semimajoraxis(
        np.full(4000, 1.0), rng=np.random.default_rng(0))
    assert np.all(log_a_draw >= np.log10(0.01) - 1e-12)
    assert np.all(log_a_draw <= np.log10(2000) + 1e-12)
    np.testing.assert_allclose(
        np.median(log_a_draw), dk.log_a_mean(1.0), atol=0.15)


def test_offner_generate_cluster_companions():
    """IMF cluster generation with Offner multiplicity produces companions."""
    imf_multi = multiplicity.MultiplicityUnresolvedOffner2023()
    mass_limits = np.array([0.01, 0.08, 0.5, 120.0])
    powers = np.array([-0.3, -1.3, -2.3])
    my_imf = imf.IMF_broken_powerlaw(mass_limits, powers, imf_multi)
    my_imf.rng = np.random.default_rng(42)
    mass, is_multi, comp_mass, sys_mass = my_imf.generate_cluster(500.0)
    assert np.any(is_multi)
    assert np.abs(500.0 - sys_mass.sum()) < 500.0 * 0.05


def test_offner_bd_q_shallower_than_v25():
    """Offner BD γ and mean q are shallower than SPISEA v2.5 Fontanive 6.1."""
    offner = multiplicity.MultiplicityUnresolvedOffner2023()
    lu = multiplicity.MultiplicityUnresolved()
    q_off = offner.q_power_at_mass(0.04)
    q_lu = lu.q_power_at_mass(0.04)
    assert 2.0 <= q_off <= 5.5
    assert np.isclose(q_lu, 6.1)
    assert q_off != q_lu

    rng = np.random.default_rng(1)
    masses = np.full(3000, 0.04)
    is_mult = np.ones(len(masses), dtype=bool)
    mf = offner.multiplicity_fraction(masses)
    csf = offner.companion_star_fraction(masses)
    comp, _, _ = offner.draw_companion_masses(
        masses, is_mult, csf, mf, rng, mass_min=0.01)
    q = comp.compressed() / 0.04
    q_lu_draw = lu.draw_q(
        np.full(len(q), 0.04), rng=np.random.default_rng(1))
    # Offner BD gamma is shallower than Fontanive 6.1, so mean q is lower.
    assert np.mean(q) < np.mean(q_lu_draw)


def test_offner2023_mf_is_vectorized():
    """Offner MF is vectorized; SPISEA v2.5 scalar BD bins are not used here."""
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    masses = np.array([0.03, 0.06, 0.10, 1.0])
    mf = multi.multiplicity_fraction(masses)
    for i, m in enumerate(masses):
        np.testing.assert_allclose(mf[i], multi.multiplicity_fraction(float(m)))
    # Array path is not the SPISEA v2.5 stellar power law 0.44 * M**0.51
    lu_pl = 0.44 * masses ** 0.51
    assert not np.allclose(mf, np.clip(lu_pl, 0, 1), atol=0.02)


def _same_seed_rngs(seed=11):
    return np.random.default_rng(seed), np.random.default_rng(seed)


def test_draw_q_same_seed():
    """Same rng seed reproduces draw_q."""
    masses = np.array([0.04, 0.3, 1.0, 2.0] * 8)
    for cls in (multiplicity.MultiplicityUnresolved,
                multiplicity.MultiplicityUnresolvedOffner2023):
        multi = cls()
        rng1, rng2 = _same_seed_rngs()
        q1 = multi.draw_q(masses, rng=rng1, n_comp=2)
        q2 = multi.draw_q(masses, rng=rng2, n_comp=2)
        np.testing.assert_array_equal(q1, q2)


def test_unresolved_draw_q_keeps_bd_stellar_split():
    """v2.5 draw_q draws stars then BDs (two RNG calls)."""
    multi = multiplicity.MultiplicityUnresolved()
    masses = np.array([1.0, 0.04, 2.0, 0.03])
    rng = np.random.default_rng(5)
    q = multi.draw_q(masses, rng=rng, n_comp=2)
    rng2 = np.random.default_rng(5)
    star = masses > multiplicity.H_BURNING_MASS
    q_exp = np.empty((len(masses), 2))
    q_exp[star] = multiplicity._q_from_powerlaw(
        rng2.random((int(star.sum()), 2)),
        multi.q_power_at_mass(masses[star]),
        multi.q_min)
    q_exp[~star] = multiplicity._q_from_powerlaw(
        rng2.random((int((~star).sum()), 2)),
        multi.q_power_at_mass(masses[~star]),
        multi.q_min)
    np.testing.assert_array_equal(q, q_exp)


def test_offner_draw_q_is_single_shot():
    """Offner draw_q draws all primaries in one RNG call."""
    multi = multiplicity.MultiplicityUnresolvedOffner2023()
    masses = np.array([1.0, 0.04, 2.0, 0.03])
    rng = np.random.default_rng(5)
    q = multi.draw_q(masses, rng=rng, n_comp=2)
    rng2 = np.random.default_rng(5)
    q_exp = multiplicity._q_from_powerlaw(
        rng2.random((len(masses), 2)),
        multi.q_power_at_mass(masses),
        multi.q_min)
    np.testing.assert_array_equal(q, q_exp)


def test_random_q_is_deprecated():
    """v2.5 random_q warns, requires mass, and delegates to draw_q."""
    import inspect
    import warnings
    import pytest

    multi = multiplicity.MultiplicityUnresolved()

    with pytest.raises(TypeError):
        multi.random_q()

    masses = np.array([1.0, 0.04])
    rng1 = np.random.default_rng(5)
    rng2 = np.random.default_rng(5)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', DeprecationWarning)
        q = multi.random_q(masses, rng=rng1, n_comp=2)
    assert any(issubclass(w.category, DeprecationWarning) for w in caught)
    assert any('draw_q' in str(w.message) for w in caught)
    q_exp = multi.draw_q(masses, rng=rng2, n_comp=2)
    np.testing.assert_array_equal(q, q_exp)

    off_src = inspect.signature(
        multiplicity.MultiplicityUnresolvedOffner2023.__init__)
    assert 'q_power' not in off_src.parameters
    res_src = inspect.signature(
        multiplicity.MultiplicityResolvedOffner2023.__init__)
    assert 'q_power' not in res_src.parameters
    for cls in (multiplicity.MultiplicityUnresolvedOffner2023,
                multiplicity.MultiplicityResolvedOffner2023):
        with pytest.raises(TypeError, match='draw_q'):
            cls().random_q(0.5)

