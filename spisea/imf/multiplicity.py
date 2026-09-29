import numpy as np
import scipy.integrate
import astropy.modeling
from astropy.table import Table
from scipy.stats import truncnorm
import os
import pdb

defaultMF_amp = 0.44
defaultMF_power = 0.51
defaultCSF_amp = 0.50
defaultCSF_power = 0.45
defaultCSF_max = 3
defaultq_power = -0.4
defaultq_min = 0.01
default_aMean = 100.0 # log (AU)
default_aSigma = 0.1  # log (AU)

class MultiplicityUnresolved(object):
    """
    The properties of stellar companions (see notes below). 
    The default parameters are as described in 
    `Lu et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013ApJ...764..155L/abstract>`_.
    These parameters are most appropriate for stellar populations
    with ages <10 Myr.

    Notes
    -----
    The number of stellar companions, their masses, and separations
    are described by the following functions:

    **Multiplicity Fraction** -- the number of stellar systems that host 
    multiple stars. In other words, the number of primary stars with
    companions. The multiplicity fraction (MF) is typically described
    as::
                            B + T + Q + ...
                MF =     ---------------------
                          S + B + T + Q + ...

    where S = single, B = binary, T = triple, Q = quadruple, etc.
    The MF also changes with mass and this dependency can be 
    described as a power-law::
            
                MF(mass) = MF_amp * (mass ** MF_power)

    However, in the brown dwarf mass regime, it is currently recognized
    that only binaries are possible, and the MF decreases dissimilarly
    to higher masses (> 0.08 solar masses). The values for this range
    are given by Aberasturi et al. (2014) and Fontanive et al. (2023).

    **Companion Star Fraction** -- the expected number of companions in
    a multiple system. The companion star fraction (CSF) also 
    changes with mass and this dependency can be described as
    a power-law::
                
                CSF(mass) = CSF_amp * (mass ** CSF_power)

    The companion star fraction is clipped to some maximum
    value, CSF_max. The actual number of companions is drawn 
    from a Poisson distribution with an expectation value of CSF.

    In the brown dwarf regime we impose an assumption that only
    binary systems are possible due to current literature trends.

    **Mass Ratio (Q)** -- The ratio between the companion star 
    mass and primary star mass, Q = (m_comp / m_prim ) has
    a probability density function described by a powerlaw::

                P(Q) = Q ** q_power  for q_min <= Q <= 1

    Current observations show no significant mass dependence.
        
    Parameters
    ----------
    MF_amp : float, optional
        The amplitude of the power-law describing the Multiplicity 
        Fraction as a function of stellar mass. 

    MF_power : float, optional
        The power of the power-law describing the Multiplicity
        Fraction as a function of stellar mass.

    CSF_amp : float, optional
        The amplitude of the power-law describing the companion star 
        fraction as a function of stellar mass.

    CSF_power : float, optional
        The power of the power-law describing the companion star 
        fraction as a function of stellar mass.

    CSF_max : float, optional
        The maximum allowed companion star fraction, which is the
        expectation value for the number of companion stars. Given
        a CSF_max = 3, some systems will still have more than 3 
        companions.

    q_power : float, optional
        The power of the power-law describing the probability
        density function for the mass ratio for stellar primaries.
        
    q_power_bd : float, optional
        The power of the power-law describing the probability
        density function for the mass ratio for brown dwarf primaries.

    q_min : float, optional
        The minimum allowed Q value for the probability
        density function of the mass ratio.
    
    companion_max : bool, optional
        Sets CSF_max as the max number of companions.
        Default False.
    """
    def __init__(self, 
                 MF_amp=0.44, MF_power=0.51,
                 CSF_amp=0.50, CSF_power=0.45, CSF_max=3,
                 q_power=-0.4, q_min=0.01, q_power_bd=6.1,
                 companion_max=False):
         
        self.is_resolved = False
        self.MF_amp = MF_amp
        self.MF_pow = MF_power
        self.CSF_amp = CSF_amp
        self.CSF_pow = CSF_power
        self.CSF_max = CSF_max
        self.q_pow = q_power
        self.q_min = q_min
        self.q_pow_bd = q_power_bd
        self.companion_max = companion_max

    def multiplicity_fraction(self, mass):
        """
        Given a star's mass, determine the probability that the star is in a
        multiple system (multiplicity fraction = MF).

        Modified to allow binary fraction to decrease in brown dwarf regime.
        Supported by Aberasturi et al. (2014) and Fontanive et al. (2018/2023).

        Parameters
        ----------
        mass : float or numpy array
            Mass of primary star.

        Returns
        -------
        mf : float or numpy array
            Multiplicity Fraction, the fraction of stars at this mass
            that will have one or more companions.
        """
        mass = np.atleast_1d(mass)
        mf = self.MF_amp * (mass ** self.MF_pow)
        mf = np.minimum(mf, 1.0)

        # Brown dwarf overrides (Aberasturi+14, Fontanive+18/23)
        bd1 = (mass <= 0.08) & (mass > 0.06)
        bd2 = (mass <= 0.06) & (mass > 0.02)
        bd3 = (mass < 0.02)

        mf[bd1] = 0.16
        mf[bd2] = 0.08
        mf[bd3] = 0.0

        return float(mf[0]) if np.isscalar(mass) else mf

    def companion_star_fraction(self, mass):
        """
        Given a star's mass, determine the average number of
        companion stars (companion star fraction = CSF). For
        brown dwarfs we impose a hard limit of one companion.

        Parameters
        ----------
        mass : float or numpy array
            Mass of primary star.

        Returns
        -------
        csf : float or numpy array
            Companion Star Fraction, the expected number of companions
            for a star at this mass.
        """
        mass = np.atleast_1d(mass)
        csf = self.CSF_amp * (mass ** self.CSF_pow)
        csf = np.minimum(csf, self.CSF_max)

        # Enforce single-companion limit in brown dwarf regime
        bd = mass <= 0.08
        if np.any(bd):
            csf[bd] = self.multiplicity_fraction(mass[bd])

        return float(csf[0]) if np.isscalar(mass) else csf

    def random_q(self, x, primary_mass):
        """
        Generative function for companion mass ratio (inverse CDF).

            `q = m_companion / m_primary`
            `P(q) = q ** beta` for q_min <= q <= 1

        Parameters
        ----------
        x : float or array_like
            Random number(s) uniform on [0, 1].
        primary_mass : float or array_like
            Primary mass(es) in solar masses.

        Returns
        -------
        q : float or numpy.ndarray
            Companion mass ratio(s).
        """
        x_arr = np.asarray(x)
        m1_arr = np.asarray(primary_mass)

        # Align primary_mass dimensions with x for multi-companion matrices (e.g. N_stars x N_comps)
        if x_arr.ndim > m1_arr.ndim:
            m1_arr = np.expand_dims(m1_arr, axis=-1)

        # Exponent b = 1 + beta
        b_stellar = 1.0 + self.q_pow
        b_bd = 1.0 + self.q_pow_bd

        # Assign exponent based on primary mass threshold (0.08 M_sun)
        b = np.where(m1_arr <= 0.08, b_bd, b_stellar)

        # Inverse CDF calculation
        q = (x_arr * (1.0 - self.q_min ** b) + self.q_min ** b) ** (1.0 / b)

        return float(q) if (np.isscalar(x) and np.isscalar(primary_mass)) else q

    def get_companions(self, mass1, rng=np.random.default_rng()):
        """
        Generate companion masses and orbital parameters for multiplicity models.

        Returns a six-tuple of masked arrays:
        (compMasses, compLoga, compEcc, compI, compOmega, compomega)

        For unresolved models, the orbital parameter arrays are fully masked,
        while resolved models populate them according to the subclass physics.
        """
        mass1 = np.atleast_1d(mass1)
        n_primaries = len(mass1)

        mf = self.multiplicity_fraction(mass1)
        csf = self.companion_star_fraction(mass1)
        is_multi = self.random_is_multiple(rng.random(n_primaries), mf)

        system_idx = np.where(is_multi)[0]
        if len(system_idx) == 0:
            empty_ma = np.ma.masked_all((n_primaries, 1))
            return (
                empty_ma,
                np.ma.masked_all((n_primaries, 1)),
                np.ma.masked_all((n_primaries, 1)),
                np.ma.masked_all((n_primaries, 1)),
                np.ma.masked_all((n_primaries, 1)),
                np.ma.masked_all((n_primaries, 1)),
            )

        primary = mass1[system_idx]
        comp_nums = 1 + rng.poisson((csf[system_idx] / mf[system_idx]) - 1)

        if self.companion_max:
            comp_nums = np.minimum(comp_nums, self.CSF_max)

        max_comp = np.max(comp_nums)
        compMasses = np.zeros((n_primaries, max_comp))

        for c_num in np.unique(comp_nums):
            group_mask = comp_nums == c_num
            group_p_idx = system_idx[group_mask]
            group_primaries = primary[group_mask]

            rand_u = rng.random((len(group_p_idx), c_num))
            q_vals = self.random_q(rand_u, group_primaries)
            compMasses[group_p_idx, :c_num] = q_vals * group_primaries[:, None]

        mask = compMasses == 0
        compMasses_ma = np.ma.MaskedArray(compMasses, mask=mask)
        compLoga_ma = np.ma.masked_all((n_primaries, max_comp))
        compEcc_ma = np.ma.masked_all((n_primaries, max_comp))
        compI_ma = np.ma.masked_all((n_primaries, max_comp))
        compOmega_ma = np.ma.masked_all((n_primaries, max_comp))
        compomega_ma = np.ma.masked_all((n_primaries, max_comp))

        return compMasses_ma, compLoga_ma, compEcc_ma, compI_ma, compOmega_ma, compomega_ma

    def random_is_multiple(self, x, MF):
        """
        Helper function: determine if star is in multiple system.
        """
        return x < MF

    def random_companion_count(self, x, CSF, MF):
        """
        Helper function: calculate number of companions.
        """
        if MF <= 0:
            return 0

        n_comp = 1 + np.random.poisson((CSF / MF) - 1)
        
        if self.companion_max:
            if n_comp > self.CSF_max:
                n_comp = self.CSF_max
            
        return n_comp


class MultiplicityResolvedDK(MultiplicityUnresolved):
    """
    Sub-class of MultiplicityUnresolved that adds semimajor axis and eccentricity information 
    for multiple objects from distributions described in Duchene and Kraus 2013.

    For brown dwarf regime, mean separation and std are given by Fontanive et al. (2018).
    
    Parameters
    --------------
    a_amp: float, optional
        Amplitude of the broken power law describing the log_semimajoraxis

    a_break: float, optional
        Break location on the x-axis of the broken power law describing the log_semimajoraxis

    a_slope1: float, optional
        Slope of the left side of the broken power law describing the log_semimajoraxis

    a_slope2: float, optional
        Slope of the right side of the broken power law describing the log_semimajoraxis

    a_std_slope: float, optional
        Slope of the line that fit sigma_log_semimajoraxis vs log_mass

    a_std_intercept: float, optional
        Intercept of the line that fit sigma_log_semimajoraxis vs log_mass
    """
    def __init__(self, a_amp=379.79953034, a_break=4.90441533, a_slope1=-1.80171539, 
                 a_slope2=4.23325571, a_std_slope=1.19713084, a_std_intercept=1.28974264, **kwargs):
        super(MultiplicityResolvedDK, self).__init__(**kwargs)
        self.is_resolved = True
        self.a_amp = a_amp
        self.a_break = a_break
        self.a_slope1 = a_slope1
        self.a_slope2 = a_slope2
        self.a_std_slope = a_std_slope
        self.a_std_intercept = a_std_intercept
    
    def log_semimajoraxis(self, mass):
        """
        Generate the semimajor axis for a given mass. The mean and standard deviation of a given mass are determined 
        by fitting the data from fitting the semimajor axis data as a function of mass in table 1 of Duchene and Kraus 2013.
        Then a random semimajor axis is drawn from a log normal distribution with that mean and standard deviation.

        The brown dwarf range is covered by mass-dependent scaling of both the characteristic separation and dispersion
        matching trends described in Fontanive et al. (2018).
        
        Parameters
        ----------
        mass : array-like
            Mass array of primary star

        Returns
        -------
        log_semimajoraxis : array-like
            Log of the semimajor axis/separation between the stars in units of AU
        """
        is_scalar = np.isscalar(mass) 
        mass = np.atleast_1d(mass)
        logm = np.log10(mass)

        # Stellar mean and std (Duchêne & Kraus 2013)
        a_mean_func = astropy.modeling.powerlaws.BrokenPowerLaw1D(
            amplitude=self.a_amp, x_break=self.a_break, alpha_1=self.a_slope1, alpha_2=self.a_slope2
        )
        log_a_mean_star = np.log10(a_mean_func(mass))
        log_a_std_func = astropy.modeling.models.Linear1D(slope=self.a_std_slope, intercept=self.a_std_intercept)
        log_a_std_star = log_a_std_func(logm)
        log_a_std_star[mass >= 2.9] = log_a_std_func(np.log10(2.9))
        log_a_std_star = np.maximum(log_a_std_star, 0.1)

        # BD mean and std (Fontanive+18): interpolated over substellar range
        log_a_mean_bd = np.interp(
            logm,
            [np.log10(0.01), np.log10(0.08)],
            [np.log10(2.5), np.log10(8.0)]
        )
        log_a_std_bd = np.interp(
            logm,
            [np.log10(0.01), np.log10(0.08)],
            [0.25, 0.5]
        )

        # Sigmoid blend: smoothly transitions from BD to stellar regime at 0.08 M_sun
        w = 1.0 / (1.0 + np.exp(-(logm - np.log10(0.08)) / 0.15))
        log_a_mean = (1 - w) * log_a_mean_bd + w * log_a_mean_star
        log_a_std = (1 - w) * log_a_std_bd + w * log_a_std_star

        # Truncated normal distribution between log10(0.01) AU and log10(2000) AU
        log_a_lower = np.log10(0.01)
        log_a_upper = np.log10(2000)

        # Convert bounds to standard normal space
        a_lower_std = (log_a_lower - log_a_mean) / log_a_std
        a_upper_std = (log_a_upper - log_a_mean) / log_a_std

        log_semimajoraxis = truncnorm.rvs(a_lower_std, a_upper_std, loc=log_a_mean, scale=log_a_std)            
        return float(log_semimajoraxis[0]) if is_scalar else log_semimajoraxis

    def random_e(self, x):
        """
        Generate random eccentricity from the inverse of the CDF where the PDF is f(e) = 2e from Duchene and Kraus 2013
        
        Parameters
        ----------
        x : float or array_like
            Random number between 0 and 1.

        Returns
        -------
        e : float or array_like
            companion mass ratio(s)
        """
        return np.sqrt(x)

    def get_companions(self, mass1, rng=np.random.default_rng()):
        """
        Generate companion masses, semimajor axes (AU), and eccentricities
        as parallel 2D MaskedArrays. Supports higher-order multiple systems 
        (triples, quadruples, etc.) based on CSF and MF.

        Parameters
        ----------
        mass1 : array-like or float
            Primary masses in solar masses.

        Returns
        -------
        compMasses : np.ma.MaskedArray
            2D masked array of companion masses (shape: N x max_comp).
        compLoga : np.ma.MaskedArray
            2D masked array of semimajor axes in AU (shape: N x max_comp).
        compEcc : np.ma.MaskedArray
            2D masked array of orbital eccentricities (shape: N x max_comp).
        compI : np.ma.MaskedArray
            2D masked array of inclinations in degrees (shape: N x max_comp).
        compOmega : np.ma.MaskedArray
            2D masked array of longitudes of ascending node in degrees (shape: N x max_comp).
        compomega : np.ma.MaskedArray
            2D masked array of arguments of periastron in degrees (shape: N x max_comp).
        """
        mass1 = np.atleast_1d(mass1)
        n_primaries = len(mass1)

        mf = self.multiplicity_fraction(mass1)
        csf = self.companion_star_fraction(mass1)
        is_multi = self.random_is_multiple(rng.random(n_primaries), mf)
        
        system_idx = np.where(is_multi)[0]
        if len(system_idx) == 0:
            empty_ma = np.ma.masked_all((n_primaries, 1))
            return empty_ma, empty_ma, empty_ma, empty_ma, empty_ma, empty_ma

        primary = mass1[system_idx]
        
        # Calculate number of companions per multiple system using Poisson draw
        comp_nums = 1 + rng.poisson((csf[system_idx] / mf[system_idx]) - 1)
        if self.companion_max:
            comp_nums = np.minimum(comp_nums, self.CSF_max)

        max_comp = np.max(comp_nums)
        compMasses = np.zeros((n_primaries, max_comp))
        compLoga = np.zeros((n_primaries, max_comp))
        compEcc = np.zeros((n_primaries, max_comp))
        compI = np.zeros((n_primaries, max_comp))
        compOmega = np.zeros((n_primaries, max_comp))
        compomega = np.zeros((n_primaries, max_comp))

        for c_num in np.unique(comp_nums):
            group_mask = comp_nums == c_num
            group_p_idx = system_idx[group_mask]
            group_primaries = primary[group_mask]
            n_group = len(group_p_idx)

            # Sample mass ratios
            rand_q = rng.random((n_group, c_num))
            q_vals = self.random_q(rand_q, group_primaries)
            
            # Sample semimajor axes for each companion
            prim_repeated = np.repeat(group_primaries, c_num)
            log_a_flat = self.log_semimajoraxis(prim_repeated)
            log_a_vals = log_a_flat.reshape(n_group, c_num)
            
            # Sample eccentricities
            ecc_vals = self.random_e(rng.random((n_group, c_num)))
            i_vals, Omega_vals, omega_vals = purely_random_keplerian_parameters(rng=rng, n=c_num * n_group)
            i_vals = np.asarray(i_vals).reshape(n_group, c_num)
            Omega_vals = np.asarray(Omega_vals).reshape(n_group, c_num)
            omega_vals = np.asarray(omega_vals).reshape(n_group, c_num)

            compMasses[group_p_idx, :c_num] = q_vals * group_primaries[:, None]
            compLoga[group_p_idx, :c_num] = log_a_vals
            compEcc[group_p_idx, :c_num] = ecc_vals
            compI[group_p_idx, :c_num] = i_vals
            compOmega[group_p_idx, :c_num] = Omega_vals
            compomega[group_p_idx, :c_num] = omega_vals

        mask = compMasses == 0
        return (
            np.ma.MaskedArray(compMasses, mask=mask),
            np.ma.MaskedArray(compLoga, mask=mask),
            np.ma.MaskedArray(compEcc, mask=mask),
            np.ma.MaskedArray(compI, mask=mask),
            np.ma.MaskedArray(compOmega, mask=mask),
            np.ma.MaskedArray(compomega, mask=mask),
        )

    
def purely_random_keplerian_parameters(rng=None, n=1):
    """
    Generate random inclination and angles for a binary or multiple system.

    Parameters
    ----------
    rng : np.random.Generator, optional
        Random number generator used for the draws.
    n : int, optional
        Number of random orbital parameter sets to generate.

    Returns
    -------
    inclination : float or numpy.ndarray
        Angle of inclination in degrees.

    Omega : float or numpy.ndarray
        Longitude of ascending node in degrees.

    omega : float or numpy.ndarray
        Argument of periastron in degrees.
    """
    if rng is None:
        rng = np.random.default_rng()

    x = rng.random(n)
    y = rng.random(n)
    z = rng.random(n)

    sign = rng.choice([-1, 1], size=len(x))
    inclination = np.arccos(sign * x) * 180.0 / np.pi
    Omega = 360.0 * y
    omega = 360.0 * z

    # if n == 1:
    #     return float(inclination[0]), float(Omega[0]), float(omega[0])
    return inclination, Omega, omega


class Multiplicity_MoeDiStefano(MultiplicityUnresolved):
    """
    Multiplicity model based on Max Moe's IDL function and the Python
    implementation from COSMIC, adapted for SPISEA.
    """
    def __init__(self, regenerate_grid=False, **kwargs):
        super(Multiplicity_MoeDiStefano, self).__init__(**kwargs)
        self.is_resolved = True

        self.M2min = 0.01

        script_dir = os.path.dirname(os.path.abspath(__file__))
        grid_path = os.path.join(script_dir, 'moe_destefano_grid.fits')

        # Build FITS lookup table automatically if it doesn't exist yet
        if (not os.path.exists(grid_path)) or regenerate_grid:
            generate_moe_destefano_grid(grid_path)

        # Load grid directly from the script's directory
        grid = Table.read(grid_path)
        self.M1v = grid['M1v'][0]
        self.logPv = grid['logPv'][0]
        self.qv = grid['qv'][0]
        self.ev = grid['ev'][0]
        self.cumqdist = grid['cumqdist'][0]
        self.cumedist = grid['cumedist'][0]
        self.cumPbindist = grid['cumPbindist'][0]

        self.mf_m_grid = np.concatenate([[0.0, 0.02,0.04, 0.07],self.M1v])
        self.mf_f_grid = np.concatenate([[0.0, 0.0, 0.08, 0.16],self.cumPbindist[-1,:]])
        self.mf_interp = scipy.interpolate.interp1d(self.mf_m_grid, self.mf_f_grid, bounds_error=False, fill_value=(0.0, 1.0))

    def multiplicity_fraction(self, mass):
        """Return the multiplicity fraction for input primary mass(es)."""
        is_scalar = np.isscalar(mass)
        mass = np.atleast_1d(mass)
        mfs = self.mf_interp(mass)  
        return float(mfs[0]) if is_scalar else mfs

    def companion_star_fraction(self, mass):
        """Companion star fraction equals multiplicity fraction in this model."""
        return self.multiplicity_fraction(mass)

    def get_companions(self, mass1, rng=np.random.default_rng()):
        """
        Generate companion masses, orbital periods (days), and eccentricities
        as parallel 2D MaskedArrays using table lookup from Moe & Di Stefano (2017).

        Parameters
        ----------
        mass1 : array-like or float
            Primary masses in solar masses.

        Returns
        -------
        compMasses : np.ma.MaskedArray
            2D masked array of companion masses (shape: N x 1).
        compLoga : np.ma.MaskedArray
            2D masked array of orbital periods in days (shape: N x 1).
        compEcc : np.ma.MaskedArray
            2D masked array of orbital eccentricities (shape: N x 1).
        compI : np.ma.MaskedArray
            2D masked array of orbital inclinations in degrees (shape: N x 1).
        compOmega : np.ma.MaskedArray
            2D masked array of longitudes of ascending node in degrees (shape: N x 1).
        compomega : np.ma.MaskedArray
            2D masked array of arguments of periastron in degrees (shape: N x 1).
        """
        mass1 = np.atleast_1d(mass1)
        n_primaries = len(mass1)

        mf = self.multiplicity_fraction(mass1)
        is_bin = self.random_is_multiple(rng.random(n_primaries), mf)
        system_idx = np.where(is_bin)[0]
        n_binaries = len(system_idx)

        compMasses = np.zeros((n_primaries, 1))
        compLoga = np.zeros((n_primaries, 1))
        compEcc = np.zeros((n_primaries, 1))
        compI = np.zeros((n_primaries, 1))
        compOmega = np.zeros((n_primaries, 1))
        compomega = np.zeros((n_primaries, 1))

        if n_binaries > 0:
            for k in range(n_binaries):
                idx = system_idx[k]
                myM1 = mass1[idx]
                if myM1 >= 0.08:
                    indM1 = np.argmin(np.abs(myM1 - self.M1v))

                    mycumPbindist_flat = self.cumPbindist[:, indM1].flatten()
                    if myM1 <= 0.8:
                        mycumPbindist_flat = mycumPbindist_flat * np.interp(
                                np.log10(myM1), np.log10([0.08, 0.8]), [0.0, 1.0])

                    mybinfrac = np.max(mycumPbindist_flat)
                    myrand = rng.random() * mybinfrac

                    mylogP = np.interp(myrand, mycumPbindist_flat, self.logPv)
                    indlogP = np.argmin(np.abs(mylogP - self.logPv))

                    mye = np.interp(rng.random(), self.cumedist[:, indlogP, indM1].flatten(), self.ev)

                    mycumqdist = self.cumqdist[:, indlogP, indM1].flatten()
                    if myM1 < self.M2min*10:
                        q_min = self.M2min / myM1 # TODO: I CHANGED THIS, ARE WE OK WITH IT
                        cum_qmin = np.interp(q_min, self.qv, mycumqdist)
                        mycumqdist = mycumqdist - cum_qmin
                        max_q = np.max(mycumqdist)
                        if max_q > 0:
                            mycumqdist = mycumqdist / max_q
                        indq = np.where(self.qv <= q_min)
                        mycumqdist[indq] = 0.0

                    myq = np.interp(rng.random(), mycumqdist, self.qv)
                    p_days = 10.0 ** mylogP
                    # Convert periods to log_a (period in days, a in AU)
                    log_a = np.log10((p_days**2 * myM1*(1+myq) * 7.496e-6)**(1.0/3)) # 

                    compMasses[idx, 0] = myq * myM1
                    compLoga[idx, 0] = log_a
                    compEcc[idx, 0] = mye

        # Handle BDs (Fontanive+18)
        # Find binary BDs
        bd_comp_idxs = np.where((mass1<0.08) & is_bin)[0]
        logm = np.log10(mass1[bd_comp_idxs])
        # MASS FIRST
        b = 1.0 + self.q_pow_bd
        # Inverse CDF calculation
        q_bds = (rng.random(len(bd_comp_idxs)) * (1.0 - self.q_min ** b) + self.q_min ** b) ** (1.0 / b)
        # SEMI-MAJOR AXIS NEXT
        # Calculate mean and standard deviation semi-major axes
        log_a_mean = np.interp(
            logm,
            [np.log10(0.01), np.log10(0.08)],
            [np.log10(2.5), np.log10(8.0)]
        )
        log_a_std = np.interp(
            logm,
            [np.log10(0.01), np.log10(0.08)],
            [0.25, 0.5]
        )
        # Truncated normal distribution between log10(0.01) AU and log10(2000) AU
        log_a_lower = np.log10(0.01)
        log_a_upper = np.log10(2000)
        # Convert bounds to standard normal space
        a_lower_std = (log_a_lower - log_a_mean) / log_a_std
        a_upper_std = (log_a_upper - log_a_mean) / log_a_std
        # Draw log_a
        log_a_bds = truncnorm.rvs(a_lower_std, a_upper_std, loc=log_a_mean, scale=log_a_std)
        # LAST: ECCENTRICITY
        ecc_bds = np.sqrt(rng.random(len(bd_comp_idxs)))
        # SAVE PROPERTIES TO ARRAYS
        compMasses[bd_comp_idxs,0] = mass1[bd_comp_idxs]*q_bds
        compLoga[bd_comp_idxs,0] = log_a_bds
        compEcc[bd_comp_idxs,0] = ecc_bds
        
        # Generate all orbital angles at once for all binaries
        all_orbital_angles = purely_random_keplerian_parameters(rng=rng, n=n_binaries)
        compI[system_idx, 0] = all_orbital_angles[0]
        compOmega[system_idx, 0] = all_orbital_angles[1]
        compomega[system_idx, 0] = all_orbital_angles[2]

        #pdb.set_trace()

        mask = compMasses == 0
        return (
            np.ma.MaskedArray(compMasses, mask=mask),
            np.ma.MaskedArray(compLoga, mask=mask),
            np.ma.MaskedArray(compEcc, mask=mask),
            np.ma.MaskedArray(compI, mask=mask),
            np.ma.MaskedArray(compOmega, mask=mask),
            np.ma.MaskedArray(compomega, mask=mask),
        )


def generate_moe_destefano_grid(grid_path):
    """Pre-computes lookup table grids and saves them to moe_destefano_grid.fits in this file's folder."""
    numM1 = 101
    bwlogP = 0.05
    numq = 91
    nume = 100

    porb_lo, porb_hi = 0.15, 8.0
    M1v = np.logspace(np.log10(0.8), np.log10(40.0), numM1)
    logPv = np.arange(porb_lo, porb_hi + bwlogP, bwlogP)
    qv = np.linspace(0.1, 1.0, numq)
    ev = np.linspace(0.0, 0.99, nume) + 0.0001

    numlogP = len(logPv)
    cumqdist = np.zeros([numq, numlogP, numM1])
    cumedist = np.zeros([nume, numlogP, numM1])
    cumPbindist = np.zeros([numlogP, numM1])

    flogP_sq = np.zeros([numlogP, numM1])
    probbin = np.zeros([numlogP, numM1])

    alpha = 0.018
    DlogP = 0.7

    H = np.zeros(numq)
    H[qv >= 0.95] = 1.0
    H = H / idl_tabulate(qv, H)

    indlq = np.where(qv >= 0.3)
    indsq = np.where(qv < 0.3)
    indq0p3 = np.min(indlq)

    for i in range(numM1):
        myM1 = M1v[i]

        FtwinlogPle1 = 0.3 - 0.15 * np.log10(myM1)
        logPtwin = 1.5 if myM1 >= 6.5 else 8.0 - myM1

        flogPle1 = 0.020 + 0.04 * np.log10(myM1) + 0.07 * (np.log10(myM1)) ** 2.0
        flogPeq2p7 = 0.039 + 0.07 * np.log10(myM1) + 0.01 * (np.log10(myM1)) ** 2.0
        flogPeq5p5 = 0.078 - 0.05 * np.log10(myM1) + 0.04 * (np.log10(myM1)) ** 2.0

        for j in range(numlogP):
            mylogP = logPv[j]

            if mylogP <= 1.0:
                Ftwin = FtwinlogPle1
            elif mylogP >= logPtwin:
                Ftwin = 0.0
            else:
                Ftwin = FtwinlogPle1 * (1.0 - (mylogP - 1.0) / (logPtwin - 1.0))

            gl_1p2 = -0.5 if mylogP <= 5.0 else -0.5 - 0.3 * (mylogP - 5.0)

            if mylogP <= 1.0:
                gl_3p5 = -0.5
            elif mylogP <= 4.5:
                gl_3p5 = -0.5 - 0.2 * (mylogP - 1.0)
            elif mylogP <= 6.5:
                gl_3p5 = -1.2 - 0.4 * (mylogP - 4.5)
            else:
                gl_3p5 = -2.0

            if mylogP <= 1.0:
                gl_6 = -0.5
            elif mylogP <= 2.0:
                gl_6 = -0.5 - 0.9 * (mylogP - 1.0)
            elif mylogP <= 4.0:
                gl_6 = -1.4 - 0.3 * (mylogP - 2.0)
            else:
                gl_6 = -2.0

            if myM1 <= 1.2:
                gl = gl_1p2
            elif myM1 <= 3.5:
                gl = np.interp(np.log10(myM1), np.log10([1.2, 3.5]), [gl_1p2, gl_3p5])
            elif myM1 <= 6.0:
                gl = np.interp(np.log10(myM1), np.log10([3.5, 6.0]), [gl_3p5, gl_6])
            else:
                gl = gl_6

            gs_1p2 = 0.3
            if mylogP <= 2.5:
                gs_3p5 = 0.2
            elif mylogP <= 5.5:
                gs_3p5 = 0.2 - 0.3 * (mylogP - 2.5)
            else:
                gs_3p5 = -0.7 - 0.2 * (mylogP - 5.5)

            if mylogP <= 1.0:
                gs_6 = 0.1
            elif mylogP <= 3.0:
                gs_6 = 0.1 - 0.15 * (mylogP - 1.0)
            elif mylogP <= 5.6:
                gs_6 = -0.2 - 0.50 * (mylogP - 3.0)
            else:
                gs_6 = -1.5

            if myM1 <= 1.2:
                gs = gs_1p2
            elif myM1 <= 3.5:
                gs = np.interp(np.log10(myM1), np.log10([1.2, 3.5]), [gs_1p2, gs_3p5])
            elif myM1 <= 6.0:
                gs = np.interp(np.log10(myM1), np.log10([3.5, 6.0]), [gs_3p5, gs_6])
            else:
                gs = gs_6

            fq = qv ** gl
            fq = fq / idl_tabulate(qv[indlq], fq[indlq])
            fq = fq * (1.0 - Ftwin) + H * Ftwin
            fq[indsq] = fq[indq0p3] * (qv[indsq] / 0.3) ** gs
            cumfq = np.cumsum(fq) - fq[0]
            cumqdist[:, j, i] = cumfq / np.max(cumfq)

            q_factor = idl_tabulate(qv, fq)

            if mylogP >= 0.7:
                eta_3 = 0.6 - 0.7 / (mylogP - 0.5)
                eta_7 = 0.9 - 0.2 / (mylogP - 0.5)
            else:
                eta_3, eta_7 = -2.9, -0.1

            if myM1 <= 3.0:
                eta = eta_3
            elif myM1 <= 7.0:
                eta = np.interp(np.log10(myM1), np.log10([3.0, 7.0]), [eta_3, eta_7])
            else:
                eta = eta_7

            if 10 ** mylogP <= 2.0:
                fe = ev ** (-3.2)
            else:
                fe = ev ** eta
                e_max = 1.0 - (10 ** mylogP / 2.0) ** (-2.0 / 3.0)
                fe[ev >= e_max] = 0.0
                ind_trans = np.where((ev >= 0.8 * e_max) & (ev <= 1.0 * e_max))
                if len(ind_trans[0]) > 0:
                    ind_cont = np.min(ind_trans) - 1
                    fe[ind_trans] = np.interp(
                        ev[ind_trans], [0.8 * e_max, 1.0 * e_max], [fe[ind_cont], 0.0]
                    )

            cumfe = np.cumsum(fe) - fe[0]
            cumedist[:, j, i] = cumfe / np.max(cumfe)

            if mylogP <= 1.0:
                flogP = flogPle1
            elif mylogP <= 2.7 - DlogP:
                flogP = flogPle1 + (mylogP - 1.0) / (1.7 - DlogP) * (
                    flogPeq2p7 - flogPle1 - alpha * DlogP
                )
            elif mylogP <= 2.7 + DlogP:
                flogP = flogPeq2p7 + alpha * (mylogP - 2.7)
            elif mylogP <= 5.5:
                flogP = flogPeq2p7 + alpha * DlogP + (mylogP - 2.7 - DlogP) / (2.8 - DlogP) * (
                    flogPeq5p5 - flogPeq2p7 - alpha * DlogP
                )
            else:
                flogP = flogPeq5p5 * np.exp(-0.3 * (mylogP - 5.5))

            flogP_sq[j, i] = flogP * q_factor

            probbin[j, i] = (
                max(0.0, 1.0 - 0.11 * (mylogP - 1.5) ** 1.43 * (myM1 / 10.0) ** 0.56)
                if mylogP > 1.5
                else 1.0
            )

        mycumPbindist = np.cumsum(flogP_sq[:, i] * probbin[:, i]) - flogP_sq[0, i] * probbin[0, i]
        cumPbindist[:, i] = (
            mycumPbindist
            / np.max(mycumPbindist)
            * idl_tabulate(logPv, flogP_sq[:, i] * probbin[:, i])
        )

    grid_table = Table({
        'M1v': [M1v],
        'logPv': [logPv],
        'qv': [qv],
        'ev': [ev],
        'cumqdist': [cumqdist],
        'cumedist': [cumedist],
        'cumPbindist': [cumPbindist]
    })

    grid_table.write(grid_path, format='fits', overwrite=True)


def idl_tabulate(x, f, p=5):
    """Replicates IDL int_tabulated function via a p-point Newton-Cotes integration."""
    def newton_cotes(x, f):
        if x.shape[0] < 2:
            return 0
        rn = (x.shape[0] - 1) * (x - x[0]) / (x[-1] - x[0])
        weights = scipy.integrate.newton_cotes(rn)[0]
        return (x[-1] - x[0]) / (x.shape[0] - 1) * np.dot(weights, f)

    ret = 0
    for idx in range(0, x.shape[0], p - 1):
        ret += newton_cotes(x[idx: idx + p], f[idx: idx + p])
    return ret
