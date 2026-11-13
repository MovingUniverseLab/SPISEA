import numpy as np
import time
import pdb
import cProfile
from spisea.imf import imf, multiplicity

def test_generate_cluster():
    # Make multiplicity object
    imf_multi = multiplicity.MultiplicityUnresolved()

    # Make IMF object; we'll use a broken power law with the parameters from Salpeter_Kirkpatrick+24 (updated with brown dwarf addition
    massLimits = np.array([0.01, 0.05, 0.22, 0.55, 8, 120]) # Define boundaries of each mass segement
    powers = np.array([-0.6, -0.25, -1.3, -2.3, -2.35]) # Power law slope associated with each mass segment
    my_imf = imf.IMF_broken_powerlaw(massLimits, powers, imf_multi)

    # Define total cluster mass
    M_cl = 10**5.

    mass, isMulti, compMass, sysMass = my_imf.generate_cluster(M_cl)

    # Make sure that the total mass is always within the expected
    # range of the requested mass (2%).
    assert np.abs(M_cl - sysMass.sum()) < M_cl*0.02

    # Check that enough companions were generated.
    # Should be greater than 25% of the stars with companions.
    mdx = np.where(isMulti)[0]
    for mm in mdx:
        assert len(compMass[mm]) > 0

    return


def test_generate_cluster_offner2023():
    """
    generate_cluster with Offner et al. 2023 multiplicity: vectorized MF,
    and Offner q (not Fontanive gamma=6.1) is used for BD companion
    masses.
    """
    imf_multi = multiplicity.MultiplicityUnresolvedOffner2023()
    massLimits = np.array([0.01, 0.05, 0.22, 0.55, 8, 120])
    powers = np.array([-0.6, -0.25, -1.3, -2.3, -2.35])
    my_imf = imf.IMF_broken_powerlaw(massLimits, powers, imf_multi)
    my_imf.rng = np.random.default_rng(7)

    M_cl = 2e3
    mass, isMulti, compMass, sysMass = my_imf.generate_cluster(M_cl)

    assert np.abs(M_cl - sysMass.sum()) < M_cl * 0.05
    bd = mass <= 0.08
    assert np.any(isMulti)

    # BD companions should be more equal-mass than SPISEA v2.5 stellar q_power=-0.4
    bd_mult = bd & isMulti
    q_bd = []
    for i in np.where(bd_mult)[0]:
        comps = compMass[i].compressed()
        if len(comps):
            q_bd.extend(list(comps / mass[i]))
    assert len(q_bd) >= 5
    assert np.mean(q_bd) > 0.5

    # Seeded golden values for rng seed 7 (same cluster as above).
    assert len(q_bd) == 85
    np.testing.assert_allclose(np.mean(q_bd), 0.8046575091308631)

    return


def test_prim_power():
    #mass_limits = np.array([0.1, 1.0, 100.0])
    #powers = np.array([-2.0, -1.8])
    mass_limits = np.array([1.0, 100.0])
    powers = np.array([-2.0])

    imf_tmp = imf.IMF_broken_powerlaw(mass_limits, powers)

    Mcl = 1e5
    (mass, is_multi, c_mass, s_mass) = imf_tmp.generate_cluster(Mcl)

    print('mass shape = ', mass.shape)
    print('mass array:')
    print(mass)
    print('is_multi array:')
    print(is_multi)
    print('c_mass array:')
    print(c_mass)
    print('s_mass array:')
    print(s_mass)
    print('np.isfinite(mass):')
    print(np.isfinite(mass))
    print('Mass Sum: ', mass.sum(), '  (should be {0})'.format(Mcl))

    return

def test_xi():
    #from pstats import SortKey
    pr = cProfile.Profile()

    mass_limits = np.array([0.1, 1.0, 10.0, 100.0])
    powers = np.array([-0.3, -1.5, -2.3])
    imf_tmp = imf.IMF_broken_powerlaw(mass_limits, powers)

    ##########
    #
    # Test validity of returned values.
    #
    ##########
    N_size = 10
    m = np.linspace(0.2, 20, N_size)
    val_good = np.array([1.6206566 , 0.26895718, 0.10135922, 0.05639448, 0.03703704,
                         0.0243668 , 0.01613091, 0.01137139, 0.00839526, 0.00642142])

    val_test = np.zeros(len(m), dtype=float)
    for ii in range(len(m)):
        val_test[ii] = imf_tmp.xi(m[ii])

    # print(repr(val_test))
    np.testing.assert_almost_equal(val_test, val_good)

    ##########
    #
    # Performance testing
    #
    ##########
    t1 = time.time()
    # pr.enable()

    # Run a time test
    N_size = int(1e4)
    m = np.random.uniform(1.1, 99.0, size=N_size)
    foo1 = imf_tmp.xi(m)

    # pr.disable()
    t2 = time.time()

    print('test_xi() runtime = {0:.3f} s for {1:d} masses'.format(t2 - t1, N_size))

    # s = io.StringIO()
    # sortby = SortKey.CUMULATIVE
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # ps.print_stats()
    # print(s.getvalue())

    return

def test_xi2():
    """
    Test that xi() produces the correct probability for
    a given slope.
    """
    #from pstats import SortKey
    pr = cProfile.Profile()

    # Negative slope
    mass_limits = np.array([1.0, 10.0])
    powers = np.array([-2.0])
    imf_tmp = imf.IMF_broken_powerlaw(mass_limits, powers)
    imf_tmp.normalize(1e4)

    tmp = np.random.rand(100)
    masses = imf_tmp.dice_star_cl(tmp)

    pdf = imf_tmp.xi(masses)

    sdx = masses.argsort()

    pdf_sorted = pdf[sdx]

    print('Testing negative slope')
    assert pdf_sorted[0] > pdf_sorted[-1]
    assert pdf_sorted[0] > pdf_sorted[1]
    assert pdf.max() == pdf_sorted[0]

    # Positive slope
    mass_limits2 = np.array([1.0, 10.0])
    powers2 = np.array([2.0])
    imf_tmp2 = imf.IMF_broken_powerlaw(mass_limits2, powers2)
    imf_tmp2.normalize(1e4)

    tmp2 = np.random.rand(100)
    masses2 = imf_tmp2.dice_star_cl(tmp2)

    pdf2 = imf_tmp2.xi(masses2)

    sdx2 = masses2.argsort()

    pdf_sorted2 = pdf2[sdx2]

    print('Testing positive slope')
    assert pdf_sorted2[0] < pdf_sorted2[-1]
    assert pdf_sorted2[0] < pdf_sorted2[1]
    assert pdf2.min() == pdf_sorted2[0]


    import pylab as plt
    plt.clf()
    plt.loglog(masses[sdx], pdf_sorted, 'k.')
    plt.axis('equal')

    return




def test_mxi():
    #from pstats import SortKey
    pr = cProfile.Profile()

    mass_limits = np.array([0.1, 1.0, 10.0, 100.0])
    powers = np.array([-0.3, -1.5, -2.3])
    imf_tmp = imf.IMF_broken_powerlaw(mass_limits, powers)

    ##########
    #
    # Test validity of returned values.
    #
    ##########
    N_size = 10
    m = np.linspace(0.2, 20, N_size)
    val_good = np.array([0.32413132, 0.64549722, 0.4662524 , 0.38348249, 0.33333333,
                         0.27290819, 0.21615424, 0.17739363, 0.14943557, 0.12842838])
    val_test = np.zeros(len(m), dtype=float)
    for ii in range(len(m)):
        val_test[ii] = imf_tmp.m_xi(m[ii])

    # print(repr(val_test))
    np.testing.assert_almost_equal(val_test, val_good)
    assert val_test.shape == m.shape

    ##########
    #
    # Performance testing
    #
    ##########
    t1 = time.time()
    # pr.enable()

    # Run a time test
    N_size = int(1e4)
    m = np.random.uniform(1.1, 99.0, size=N_size)
    foo1 = imf_tmp.m_xi(m)
    assert foo1.shape == m.shape

    # pr.disable()
    t2 = time.time()

    print('test_mxi() runtime = {0:.3f} s for {1:d} masses'.format(t2 - t1, N_size))

    # s = io.StringIO()
    # sortby = SortKey.CUMULATIVE
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # ps.print_stats()
    # print(s.getvalue())

    return

def test_theta_closed():
    #from pstats import SortKey

    mass_limits = np.array([0.1, 1.0, 10.0, 100.0])
    powers = np.array([-0.3, -1.5, -2.3])
    imf_tmp = imf.IMF_broken_powerlaw(mass_limits, powers)

    ##########
    #
    # Test validity of returned values.
    #
    ##########
    N_size = 10
    m = np.linspace(0.2, 20, N_size)
    val_good = np.array([[1., 0., 0.],
       [1., 1., 0.],
       [1., 1., 0.],
       [1., 1., 0.],
       [1., 1., 0.],
       [1., 1., 1.],
       [1., 1., 1.],
       [1., 1., 1.],
       [1., 1., 1.],
       [1., 1., 1.]])

    val_test = np.zeros((len(m), len(powers)), dtype=float)
    for ii in range(len(m)):
        val_test[ii] = imf.theta_closed(m[ii] - imf_tmp._m_limits_low)

    np.testing.assert_equal(val_test, val_good)



    ##########
    #
    # Speed tests and performance profiling.
    #
    ##########
    N_size = 10000
    m = np.linspace(1.1, 99, N_size)

    tmp = np.zeros((len(m), len(powers)), dtype=float)

    t1 = time.time()
    # pr = cProfile.Profile()
    # pr.enable()

    for ii in range(len(m)):
        tmp[ii] = imf.theta_closed(m[ii] - imf_tmp._m_limits_low)

    # pr.disable()
    t2 = time.time()

    print('Runtime = {0:.3f} s for {1:d} masses'.format(t2 - t1, N_size))

    # s = io.StringIO()
    # sortby = SortKey.CUMULATIVE
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # ps.print_stats()
    # print(s.getvalue())


    return

