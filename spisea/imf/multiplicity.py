import warnings

import numpy as np
import astropy.modeling
from scipy.stats import truncnorm

# Hydrogen-burning limit used for BD-primary (binaries-only) logic.
# Offner et al. 2023 use M_comp > 0.075 Msun as the MS companion cut;
# SPISEA keeps 0.08 Msun for consistency with existing BD handling.
H_BURNING_MASS = 0.08


class _ResolvedOrbitalMixin(object):
    """Eccentricity and Keplerian angles shared by resolved multiplicity classes."""

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
            eccentricity
        """
        e = np.sqrt(x)
        
        return e
    
    def random_keplarian_parameters(self, x, y, z, rng=None):
        """
        Generate random inclination and angles of a binary system.

        Inclination uses the inverse CDF of isotropic orientations
        (``i = arccos(s * x)`` with a random sign ``s``). ``Omega``
        and ``omega`` are uniform on [0, 360) deg from ``y`` and ``z``.

        Parameters
        ----------
        x : float or array_like
            Random number between 0 and 1, used for inclination.
        y : float or array_like
            Random number between 0 and 1, used for Omega.
        z : float or array_like
            Random number between 0 and 1, used for omega.
        rng : numpy.random.Generator, optional
            Random number generator used for the inclination sign.
            Default is a new ``numpy.random.default_rng()``.

        Returns
        -------
        inclination : float or array_like
            Inclination angle in degrees.
        Omega : float or array_like
            Longitude of the ascending node in degrees.
        omega : float or array_like
            Argument of periastron in degrees.
        """
        if rng is None:
            rng = np.random.default_rng()
        sign = rng.choice([-1, 1], size=len(x))
        x = sign * x
        inclination = np.arccos(x) * 180 / np.pi

        Omega = 360 * y
        omega = 360 * z

        return inclination, Omega, omega


class MultiplicityUnresolved(object):
    """
    SPISEA v2.5 default unresolved multiplicity (companions, no orbits).

    The default parameters are as described in
    `Lu et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013ApJ...764..155L/abstract>`_.
    These parameters are most appropriate for stellar populations
    with ages < 10 Myr. This is the unresolved default for backwards
    compatibility. For a scientifically preferred model that includes
    brown dwarfs, see :class:`MultiplicityUnresolvedOffner2023`
    (opt-in; not the default).

    Notes
    -----
    The number of stellar companions and their masses are described
    by the following functions.

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

    Defaults are MF_amp = 0.44, MF_power = 0.51 (Lu et al. 2013).
    MF is clipped to [0, 1].

    In the brown-dwarf regime, only binaries are expected and the MF
    does not follow the stellar power law. For **scalar** mass,
    :meth:`multiplicity_fraction` applies a staircase from
    Aberasturi et al. (2014) and Fontanive et al. (2018)::

        M < 0.02 Msun          MF = 0
        0.02 < M <= 0.06 Msun  MF = 0.08
        0.06 < M <= 0.08 Msun  MF = 0.16

    Array masses use the stellar power law only over the whole mass range, 
    including brown dwarf primaries. .
    
    **Companion Star Fraction** -- the expected number of companions in
    a multiple system. The companion star fraction (CSF) also
    changes with mass and this dependency can be described as
    a power-law::

                CSF(mass) = CSF_amp * (mass ** CSF_power)

    Defaults are CSF_amp = 0.50, CSF_power = 0.45. CSF is clipped
    to CSF_max (default 3). The actual number of companions is drawn
    from a Poisson distribution with an expectation value of CSF.
    Higher-order multiples (triples+) are allowed at all masses,
    including brown dwarfs. If ``companion_max`` is True, the
    companion count is capped at CSF_max at all masses.

    **Mass Ratio (Q)** -- The ratio between the companion star
    mass and primary star mass, Q = (m_comp / m_prim) has
    a probability density function described by a power law::

                P(Q) = Q ** q_power  for q_min <= Q <= 1

    Stellar primaries use q_power = −0.4 (Lu et al. 2013).
    Brown-dwarf primaries (M <= 0.08 Msun) use γ = 6.1 from
    Fontanive et al. (2018) via :meth:`q_power_at_mass` /
    :meth:`draw_q`. That 0.08 Msun switch is a mass-ratio index
    only, not a companion-count policy. :meth:`random_q` is
    deprecated; use :meth:`draw_q`.

    Parameters
    ----------
    MF_amp : float, optional
        Amplitude of the MF power law, dimensionless
        (units of MF / Msun**MF_power). Default 0.44.
    MF_power : float, optional
        Power-law index of MF(M), dimensionless. Default 0.51.
    CSF_amp : float, optional
        Amplitude of the CSF power law, dimensionless mean
        companion count / Msun**CSF_power. Default 0.50.
    CSF_power : float, optional
        Power-law index of CSF(M), dimensionless. Default 0.45.
    CSF_max : float, optional
        Maximum companion star fraction, dimensionless mean
        companion count (not bounded by 1). Default 3.
    q_power : float, optional
        Stellar mass-ratio power-law index γ, dimensionless.
        Default −0.4 (Lu et al. 2013).
    q_min : float, optional
        Minimum mass ratio m_comp/m_prim, dimensionless.
        Default 0.01.
    companion_max : bool, optional
        If True, cap companion counts at CSF_max at all masses.
        Default False.
    bd_q_power : float, optional
        Mass-ratio power-law index γ for brown-dwarf primaries
        (M <= 0.08 Msun), dimensionless. Default 6.1
        (Fontanive et al. 2018).
    """
    def __init__(self,
                 MF_amp=0.44, MF_power=0.51,
                 CSF_amp=0.50, CSF_power=0.45, CSF_max=3,
                 q_power=-0.4, q_min=0.01, companion_max=False,
                 bd_q_power=6.1):

        self.MF_amp = MF_amp
        self.MF_pow = MF_power
        self.CSF_amp = CSF_amp
        self.CSF_pow = CSF_power
        self.CSF_max = CSF_max
        self.q_pow = q_power
        self.q_min = q_min
        self.companion_max = companion_max
        self.bd_q_power = bd_q_power

    def multiplicity_fraction(self, mass):
        """
        Given a star's mass, determine the probability that the star is in a
        multiple system (multiplicity fraction = MF).

        Arrays use ``MF = MF_amp * M**MF_power`` clipped to 1.
        Scalars also apply the Aberasturi/Fontanive brown-dwarf
        staircase (0 / 8% / 16% below 0.02 / 0.06 / 0.08 Msun).

        Parameters
        ----------
        mass : float or array_like
            Primary mass n solar masses (Msun).

        Returns
        -------
        mf : float or ndarray
            Multiplicity fraction, dimensionless, in [0, 1].
            The fraction of stars at this mass that will have one or
            more companions. Python float if ``mass`` is scalar,
            ndarray otherwise.
        """
        # Multiplicity Fraction
        mf = self.MF_amp * mass ** self.MF_pow

        if np.isscalar(mf):
            if mf > 1:
                mf = 1
            # physically override mf for brown dwarfs
            if (mass <= 0.08) & (mass > 0.06):
                mf = 0.16
            if (mass <= 0.06) & (mass > 0.02):
                mf = 0.08
            if (mass < 0.02):
                mf = 0
        else:
            mf[mf > 1] = 1

        return mf

    def companion_star_fraction(self, mass):
        """
        Given a star's mass, determine the average number of
        companion stars (companion star fraction = CSF).

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        csf : float or ndarray
            Companion star fraction, the expected number of companions
            for a star at this mass. Dimensionless mean companion count
            (not bounded by 1). Python float if ``mass`` is scalar,
            ndarray otherwise.
        """
        # Companion Star Fraction
        csf = self.CSF_amp * mass ** self.CSF_pow

        if np.isscalar(csf):
            if csf > self.CSF_max:
                csf = self.CSF_max
        else:
            csf[csf > self.CSF_max] = self.CSF_max

        return csf

    def q_power_at_mass(self, mass):
        """
        Mass-ratio power-law index, P(q) ∝ q ** q_power.

        Lu et al. (2013) use a single ``q_power`` for stellar primaries.
        Brown-dwarf primaries (M <= 0.08 Msun) use ``bd_q_power``
        (default 6.1, Fontanive et al. 2018), matching the
        companion-mass draw previously special-cased in
        ``imf.calc_multi``.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        q_power : float or ndarray
            Mass-ratio power-law index γ, dimensionless.
            Python float if ``mass`` is scalar, ndarray otherwise.
        """
        mass_arr = np.atleast_1d(np.asarray(mass, dtype=float))
        q_pow = np.full(mass_arr.shape, self.q_pow, dtype=float)
        q_pow[mass_arr <= H_BURNING_MASS] = self.bd_q_power
        if np.isscalar(mass):
            return float(q_pow[0])
        return q_pow

    def random_q(self, mass, rng=None, n_comp=1):
        """
        Deprecated wrapper for companion mass ratios.

        Use :meth:`draw_q`. ``mass`` is required; γ comes from
        :meth:`q_power_at_mass`. Random draws use ``rng``.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).
        rng : numpy.random.Generator, optional
            Random generator. Default ``numpy.random.default_rng()``.
        n_comp : int, optional
            Companions per primary. Default 1.

        Returns
        -------
        q : float or ndarray
            Mass ratios in [q_min, 1]. Same shape convention as
            :meth:`draw_q`.
        """
        warnings.warn(
            "random_q is deprecated; use draw_q(mass, rng=...).",
            DeprecationWarning,
            stacklevel=2)

        return self.draw_q(mass, rng=rng, n_comp=n_comp)

    def draw_q(self, mass, rng=None, n_comp=1):
        """
        Draw companion mass ratios q = m_comp / m_prim.

        Uses ``q_power_at_mass(mass)`` for γ. ``mass`` is required.
        Stellar and brown-dwarf primaries are drawn in two separate
        RNG calls (star_mask then bd_mask, split at
        ``H_BURNING_MASS``) so the SPISEA v2.5 ``imf.calc_multi``
        random sequence is unchanged.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).
        rng : numpy.random.Generator, optional
            Random generator. Default ``numpy.random.default_rng()``.
        n_comp : int, optional
            Companions per primary. Default 1.

        Returns
        -------
        q : float or ndarray
            Mass ratios in [q_min, 1]. Scalar if ``mass`` is scalar
            and ``n_comp`` is 1; shape ``(n_mass,)`` or
            ``(n_mass, n_comp)`` otherwise.
        """
        if rng is None:
            rng = np.random.default_rng()
        n_comp = int(n_comp)
        return_scalar = np.isscalar(mass) and n_comp == 1
        mass_arr = np.atleast_1d(np.asarray(mass, dtype=float))
        q_values = np.empty((len(mass_arr), n_comp), dtype=float)
        bd_mask = mass_arr <= H_BURNING_MASS
        star_mask = ~bd_mask

        if np.any(star_mask):
            q_pow = self.q_power_at_mass(mass_arr[star_mask])
            x = rng.random((int(star_mask.sum()), n_comp))
            q_values[star_mask] = _q_from_powerlaw(x, q_pow, self.q_min)

        if np.any(bd_mask):
            q_pow = self.q_power_at_mass(mass_arr[bd_mask])
            x = rng.random((int(bd_mask.sum()), n_comp))
            q_values[bd_mask] = _q_from_powerlaw(x, q_pow, self.q_min)

        return _finalize_q_draw(q_values, return_scalar, n_comp)

    def random_is_multiple(self, x, MF):
        """
        Helper function: determine if star is in multiple system.
        """
        return x < MF

    def random_companion_count(self, x, CSF, MF, mass=None, rng=None):
        """
        Number of companions for primaries already identified as multiple.

        The count is drawn from a Poisson with expectation CSF/MF - 1,
        then 1 is added so every multiple has at least one companion.
        ``x`` is unused and kept for API compatibility.

        Parameters
        ----------
        x : float or array_like
            Unused (historical signature). Dimensionless uniform
            draw in [0, 1] if provided.
        CSF : float or array_like
            Companion star fraction, dimensionless mean companion
            count (not bounded by 1).
        MF : float or array_like
            Multiplicity fraction, dimensionless, in [0, 1].
        mass : float or array_like, optional
            Primary mass in solar masses (Msun). 
        rng : numpy.random.Generator, optional
            Random generator. If omitted, uses ``numpy.random`` (the
            historical scalar helper).

        Returns
        -------
        n_comp : int or ndarray of int
            Number of companions, integer count. Python int if ``CSF``
            and ``MF`` are scalar, ndarray otherwise.
        """
        return_scalar = np.isscalar(CSF) and np.isscalar(MF)
        if return_scalar and rng is None:
            if MF <= 0:
                return 0
            n_comp = 1 + np.random.poisson((CSF / MF) - 1)
            if self.companion_max and n_comp > self.CSF_max:
                n_comp = self.CSF_max
            return int(n_comp)

        CSF = np.atleast_1d(np.asarray(CSF, dtype=float))
        MF = np.atleast_1d(np.asarray(MF, dtype=float))
        if rng is None:
            n_comp = 1 + np.random.poisson((CSF / MF) - 1)
        else:
            n_comp = 1 + rng.poisson((CSF / MF) - 1)

        if self.companion_max:
            n_comp = np.minimum(n_comp, self.CSF_max)

        if return_scalar:
            return int(n_comp[0])

        return n_comp

    def draw_n_companions(self, mass, CSF, MF, rng):
        """
        Vectorized companion counts for primaries that are already
        identified as multiple. Delegates to
        :meth:`random_companion_count`. If ``companion_max`` is True,
        counts are capped at CSF_max at all masses.

        Parameters
        ----------
        mass : array_like
            Primary masses of systems already identified as multiple.
            Must be positive, in solar masses (Msun).
        CSF : array_like
            Companion star fraction at each primary, dimensionless
            mean companion count (not bounded by 1).
        MF : array_like
            Multiplicity fraction at each primary, dimensionless,
            in [0, 1].
        rng : numpy.random.Generator
            Random generator used for the Poisson companion-count draw.

        Returns
        -------
        n_comp : ndarray of int
            Number of companions per primary, integer count, shape
            matching ``mass``.
        """
        n_comp = self.random_companion_count(None, CSF, MF, mass=mass, rng=rng)

        return np.atleast_1d(n_comp)

    def draw_companion_masses(self, primary_masses, is_multiple, CSF, MF,
                              rng, mass_min):
        """
        Assign companion masses for a set of primaries.

        This is the multiplicity-object entry point used by
        ``IMF.calc_multi``. Companion-mass draws use :meth:`draw_q`
        and live here rather than in ``imf.py``.

        Parameters
        ----------
        primary_masses : array_like
            Primary masses must be positive, in solar masses (Msun).
        is_multiple : array_like of bool
            True for primaries drawn as multiple systems.
        CSF : array_like
            Companion star fraction at each primary, dimensionless
            mean companion count (not bounded by 1).
        MF : array_like
            Multiplicity fraction at each primary, dimensionless,
            in [0, 1].
        rng : numpy.random.Generator
            Random generator.
        mass_min : float
            Minimum companion mass in solar masses (Msun); lighter
            companions are masked.

        Returns
        -------
        comp_masses : numpy.ma.MaskedArray
            Companion masses in solar masses (Msun), shape
            (n_primaries, max_n_comp).
        system_masses : ndarray
            Primary plus unmasked companion mass, in solar masses
            (Msun).
        is_multiple : ndarray of bool
            Updated multiplicity flags after masking sub-minimum
            companions.
        """
        primary_masses = np.asarray(primary_masses, dtype=float)
        is_multiple = np.asarray(is_multiple, dtype=bool)
        CSF = np.asarray(CSF, dtype=float)
        MF = np.asarray(MF, dtype=float)

        system_masses = primary_masses.copy()
        multiple_idx = np.where(is_multiple)[0]
        primary = primary_masses[multiple_idx]
        n_comp = self.draw_n_companions(primary, CSF[multiple_idx], MF[multiple_idx], rng)

        if len(multiple_idx) == 0:
            comp_masses = np.zeros((len(primary_masses), 1), dtype=float)
            comp_masses = np.ma.MaskedArray(comp_masses, mask=comp_masses < mass_min)

            return comp_masses, system_masses, is_multiple

        n_unique = np.unique(n_comp)
        n_indices = [np.where(n_comp == i)[0] for i in n_unique]
        comp_masses = np.zeros((len(primary_masses), int(np.max(n_unique))))

        for n_c, idx in zip(n_unique, n_indices):
            prim_subset = primary[idx]
            q_values = np.asarray(
                self.draw_q(prim_subset, rng=rng, n_comp=int(n_c)),
                dtype=float)
            if q_values.ndim == 1:
                q_values = q_values[:, np.newaxis]
            m_comp = q_values * prim_subset[:, np.newaxis]
            comp_masses[multiple_idx[idx], :int(n_c)] = m_comp

        comp_masses = np.ma.MaskedArray(comp_masses, mask=comp_masses < mass_min)
        system_masses[multiple_idx] += comp_masses[multiple_idx].sum(axis=1)
        is_multiple = np.any(~comp_masses.mask, axis=1)

        return comp_masses, system_masses, is_multiple
    

class MultiplicityResolvedDK(MultiplicityUnresolved, _ResolvedOrbitalMixin):
    """
    SPISEA v2.5 default resolved multiplicity.

    Same MF/CSF/q as :class:`MultiplicityUnresolved` (Lu et al. 2013
    stellar power laws, scalar BD staircase, Fontanive γ = 6.1), plus
    semimajor axis and eccentricity from Duchêne & Kraus (2013).
    This is the resolved default for backwards compatibility. It is
    **not** meant for the brown-dwarf regime the way
    :class:`MultiplicityResolvedOffner2023` is.

    Notes
    -----
    Characteristic μ(a) is a broken power law in a vs M
    (``astropy.modeling.powerlaws.BrokenPowerLaw1D`` with the fitted
    ``a_amp``, ``a_break``, ``a_slope1``, ``a_slope2``).
    σ(log10 a) is a linear fit in log M (``a_std_slope``,
    ``a_std_intercept``) that saturates above 2.9 Msun and is
    clipped to ≥ 0.1 dex. The dip near 0.08 Msun is the BD/stellar
    blend: Fontanive et al. (2018) interpolation of the BD mean
    (2.5–8 AU) and width (0.25–0.5 dex) over 0.01–0.08 Msun,
    combined with the stellar law by a sigmoid at 0.08 Msun.

    :meth:`log_a_mean`, :meth:`a_mean`, and :meth:`sigma_log_a`
    return that same characteristic mean and width.
    :meth:`log_semimajoraxis` draws from a truncated lognormal
    with ``loc = log_a_mean(mass)`` and
    ``scale = sigma_log_a(mass)``, truncated to 0.01–2000 AU.
    Eccentricity follows f(e) = 2e; inclination and angles are
    random (Duchêne & Kraus 2013).

    Parameters
    ----------
    a_amp : float, optional
        Amplitude of the broken power law for μ(a), in AU.
        Default 379.79953034.
    a_break : float, optional
        Break mass of the broken power law, in solar masses (Msun).
        Default 4.90441533.
    a_slope1 : float, optional
        Power-law index below ``a_break``, dimensionless.
        Default −1.80171539.
    a_slope2 : float, optional
        Power-law index above ``a_break``, dimensionless.
        Default 4.23325571.
    a_std_slope : float, optional
        Slope of σ(log10 a) vs log10(M / Msun), in dex per dex.
        Default 1.19713084.
    a_std_intercept : float, optional
        Intercept of σ(log10 a) vs log10(M / Msun), in dex.
        Default 1.28974264.
    **kwargs
        Passed to :class:`MultiplicityUnresolved` (MF/CSF/q).
    """
    def __init__(self, a_amp = 379.79953034, a_break = 4.90441533, a_slope1 = -1.80171539, 
                 a_slope2 = 4.23325571, a_std_slope = 1.19713084, a_std_intercept = 1.28974264, **kwargs):
        super(MultiplicityResolvedDK, self).__init__(**kwargs)
        self.a_amp = a_amp
        self.a_break = a_break
        self.a_slope1 = a_slope1
        self.a_slope2 = a_slope2
        self.a_std_slope = a_std_slope
        self.a_std_intercept = a_std_intercept

        return

    def log_a_mean(self, mass):
        """
        Characteristic log10(a/AU) used as the loc of
        :meth:`log_semimajoraxis`.

        Duchêne & Kraus (2013) broken power law for stellar primaries,
        Fontanive et al. (2018) interpolation for brown dwarfs, and a
        sigmoid blend at 0.08 Msun.

        Parameters
        ----------
        mass : float or array_like
            Primary mass must be positive, in solar masses (Msun).

        Returns
        -------
        log_a_mean : float or ndarray
            Characteristic log10(a / 1 AU) in dex (not ln, not AU).
            Python float if ``mass`` is scalar, ndarray otherwise.
        """
        return_scalar = np.isscalar(mass)
        mass = np.atleast_1d(np.asarray(mass, dtype=float))
        logm = np.log10(mass)

        a_mean_func = astropy.modeling.powerlaws.BrokenPowerLaw1D(
            amplitude=self.a_amp, x_break=self.a_break,
            alpha_1=self.a_slope1, alpha_2=self.a_slope2)
        log_a_mean_star = np.log10(a_mean_func(mass))
        log_a_mean_bd = np.interp(
            logm,
            [np.log10(0.01), np.log10(0.08)],
            [np.log10(2.5), np.log10(8.0)]
        )
        w = 1.0 / (1.0 + np.exp(-(logm - np.log10(0.08)) / 0.15))
        log_a_mean = (1 - w) * log_a_mean_bd + w * log_a_mean_star
        if return_scalar:
            return float(log_a_mean[0])
        return log_a_mean

    def a_mean(self, mass):
        """
        Characteristic μ(a) in AU, ``10 ** log_a_mean(mass)``.

        Parameters
        ----------
        mass : float or array_like
            Primary mass must be positive, in solar masses (Msun).

        Returns
        -------
        a_mean : float or ndarray
            Characteristic separation μ(a) in AU. Python float if
            ``mass`` is scalar, ndarray otherwise.
        """
        log_a = self.log_a_mean(mass)
        if np.isscalar(log_a):
            return 10.0 ** log_a
        return 10.0 ** np.asarray(log_a, dtype=float)

    def sigma_log_a(self, mass):
        """
        σ(log10 a) used as the scale of :meth:`log_semimajoraxis`.

        Linear fit in log-mass for stellar primaries (saturates above
        2.9 Msun, clipped to ≥ 0.1), Fontanive et al. (2018)
        interpolation for brown dwarfs, and a sigmoid blend at
        0.08 Msun.

        Parameters
        ----------
        mass : float or array_like
            Primary mass must be positive, in solar masses (Msun).

        Returns
        -------
        sigma_log_a : float or ndarray
            Standard deviation of log10(a / 1 AU), in dex.
            Python float if ``mass`` is scalar, ndarray otherwise.
        """
        return_scalar = np.isscalar(mass)
        mass = np.atleast_1d(np.asarray(mass, dtype=float))
        logm = np.log10(mass)

        log_a_std_func = astropy.modeling.models.Linear1D(
            slope=self.a_std_slope, intercept=self.a_std_intercept)
        log_a_std_star = log_a_std_func(logm)
        log_a_std_star[mass >= 2.9] = log_a_std_func(np.log10(2.9))
        log_a_std_star = np.clip(log_a_std_star, 0.1, None)
        log_a_std_bd = np.interp(
            logm,
            [np.log10(0.01), np.log10(0.08)],
            [0.25, 0.5]
        )
        w = 1.0 / (1.0 + np.exp(-(logm - np.log10(0.08)) / 0.15))
        log_a_std = (1 - w) * log_a_std_bd + w * log_a_std_star
        if return_scalar:
            return float(log_a_std[0])
        return log_a_std
    
    def log_semimajoraxis(self, mass, rng=None):
        """
        Draw log10(a/AU) from a mass-dependent truncated lognormal.

        The mean and standard deviation at a given mass come from
        fitting Table 1 of Duchêne and Kraus 2013. The brown dwarf
        range uses mass-dependent scaling of both the characteristic
        separation and dispersion (Fontanive et al. 2018).

        Draws use ``loc = log_a_mean(mass)`` and
        ``scale = sigma_log_a(mass)``.

        Parameters
        ----------
        mass : array-like
            Primary mass must be positive, in solar masses (Msun).
        rng : numpy.random.Generator, optional
            Random number generator passed to ``truncnorm.rvs``.
            Default is a new ``numpy.random.default_rng()``.

        Returns
        -------
        log_semimajoraxis : array-like
            Drawn log10(a / 1 AU) in dex, truncated so a is between
            0.01 AU and 2000 AU.
        """
        if rng is None:
            rng = np.random.default_rng()
        mass = np.atleast_1d(mass)
        log_a_mean = np.atleast_1d(np.asarray(self.log_a_mean(mass),
                                             dtype=float))
        log_a_std = np.atleast_1d(np.asarray(self.sigma_log_a(mass),
                                            dtype=float))

        # Trunc normal between log10(0.01) AU and log10(2000) AU
        log_a_lower = np.log10(0.01)
        log_a_upper = np.log10(2000)
        a_lower_std = (log_a_lower - log_a_mean) / log_a_std
        a_upper_std = (log_a_upper - log_a_mean) / log_a_std

        log_semimajoraxis = truncnorm.rvs(
            a_lower_std, a_upper_std, loc=log_a_mean, scale=log_a_std,
            random_state=rng)
        return log_semimajoraxis


class MultiplicityPiecewisePowerLaw(MultiplicityUnresolved):
    """
    Generic helper: multiplicity as a piecewise power law in primary mass.

    Use this when MF and CSF should change slope at user-specified
    mass edges (survey-specific segments), rather than a single
    power law or a logistic. Offner et al. 2023 does **not** use
    this class; it uses :class:`MultiplicityLogistic`.

    Notes
    -----
    On each mass segment i, with edges
    ``mass_limits[i] <= M < mass_limits[i+1]``::

        MF(M)  = MF_amp[i]  * M ** MF_power[i]
        CSF(M) = CSF_amp[i] * M ** CSF_power[i]

    MF is clipped to [0, 1]. CSF is clipped to [0, CSF_max] and
    raised to at least MF so the Poisson companion-count draw is
    well defined. Higher-order multiples are allowed at all masses.
    If ``companion_max`` is True, counts are capped at CSF_max at
    all masses. Mass-ratio draws use a single ``q_power`` (same as
    :class:`MultiplicityUnresolved`) unless a subclass overrides
    :meth:`q_power_at_mass`. There is no scalar-only brown-dwarf
    staircase; each segment is evaluated for both scalar and array
    mass.

    Parameters
    ----------
    mass_limits : array_like
        Segment edges in solar masses (Msun), length N+1, strictly
        increasing.
    MF_amps : array_like
        Length-N amplitudes for the multiplicity fraction,
        dimensionless (units of MF / Msun**MF_power).
    MF_powers : array_like
        Length-N power-law indices for the multiplicity fraction,
        dimensionless.
    CSF_amps : array_like
        Length-N amplitudes for the companion star fraction,
        dimensionless (mean companion count / Msun**CSF_power).
    CSF_powers : array_like
        Length-N power-law indices for the companion star fraction,
        dimensionless.
    CSF_max : float, optional
        Maximum companion star fraction, dimensionless mean companion
        count (not bounded by 1). Default 3.
    q_power : float, optional
        Mass-ratio power-law index, dimensionless. Default -0.4.
    q_min : float, optional
        Minimum mass ratio m_comp/m_prim, dimensionless. Default 0.01.
    companion_max : bool, optional
        If True, cap companion counts at CSF_max at all masses.
        Default False.
    """
    def __init__(self, mass_limits, MF_amps, MF_powers, CSF_amps, CSF_powers,
                 CSF_max=3, q_power=-0.4, q_min=0.01, companion_max=False):
        mass_limits = np.asarray(mass_limits, dtype=float)
        MF_amps = np.asarray(MF_amps, dtype=float)
        MF_powers = np.asarray(MF_powers, dtype=float)
        CSF_amps = np.asarray(CSF_amps, dtype=float)
        CSF_powers = np.asarray(CSF_powers, dtype=float)

        nseg = len(MF_amps)

        if len(mass_limits) != nseg + 1:
            raise ValueError('len(mass_limits) must be len(MF_amps) + 1')

        if not (len(MF_powers) == len(CSF_amps) == len(CSF_powers) == nseg):
            raise ValueError('MF/CSF amplitude and power arrays must have equal length')

        if np.any(np.diff(mass_limits) <= 0):
            raise ValueError('mass_limits must be strictly increasing')

        super(MultiplicityPiecewisePowerLaw, self).__init__(
            MF_amp=MF_amps[-1], MF_power=MF_powers[-1],
            CSF_amp=CSF_amps[-1], CSF_power=CSF_powers[-1],
            CSF_max=CSF_max, q_power=q_power, q_min=q_min,
            companion_max=companion_max)

        self.mass_limits = mass_limits
        self.MF_amps = MF_amps
        self.MF_powers = MF_powers
        self.CSF_amps = CSF_amps
        self.CSF_powers = CSF_powers

        return

    def multiplicity_fraction(self, mass):
        """
        Multiplicity fraction as a piecewise power law in primary mass.
        Clipped to [0, 1].

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        mf : float or ndarray
            Multiplicity fraction, dimensionless, in [0, 1].
            Python float if ``mass`` is scalar, ndarray otherwise.
        """
        mf = _piecewise_powerlaw(mass, self.mass_limits, self.MF_amps, self.MF_powers, clip_min=0.0, clip_max=1.0)

        return mf

    def companion_star_fraction(self, mass):
        """
        Companion star fraction as a piecewise power law in primary mass.

        Clipped to [0, CSF_max] and raised to at least MF.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        csf : float or ndarray
            Companion star fraction, dimensionless mean companion
            count (not bounded by 1). Python float if ``mass`` is
            scalar, ndarray otherwise.
        """
        return_scalar = np.isscalar(mass)
        mass_arr = np.atleast_1d(np.asarray(mass, dtype=float))
        csf = _piecewise_powerlaw(
            mass_arr, self.mass_limits, self.CSF_amps, self.CSF_powers,
            clip_min=0.0, clip_max=self.CSF_max)
        mf = _piecewise_powerlaw(
            mass_arr, self.mass_limits, self.MF_amps, self.MF_powers,
            clip_min=0.0, clip_max=1.0)
        csf = np.maximum(csf, mf)
        if return_scalar:
            return float(csf[0])
        return csf


class MultiplicityLogistic(MultiplicityUnresolved):
    """
    Generic helper: multiplicity as a logistic in log primary mass.

    Use this for a C-infinity smooth MF/CSF that saturates at low
    and high mass. :class:`MultiplicityUnresolvedOffner2023` is this
    class with coefficients fitted to Offner et al. (2023) Table 1.

    Notes
    -----
    The same 4-parameter logistic is used for MF and CSF with
    independent coefficients::

        f(M) = A + (B - A) / (1 + (M / M0)**(-k))

    As M → 0, f → A; as M → ∞, f → B. The curve is not a piecewise
    interpolation of survey knots. MF is clipped to [0, 1]. CSF is
    clipped to [0, CSF_max] and raised to at least MF. Higher-order
    multiples are allowed at all masses. If ``companion_max`` is
    True, counts are capped at CSF_max at all masses. Mass-ratio
    draws use a single ``q_power`` unless a subclass overrides
    :meth:`q_power_at_mass`.

    Parameters
    ----------
    MF_A, MF_B : float
        Low-mass and high-mass asymptotes of the multiplicity-fraction
        logistic, dimensionless (MF).
    MF_M0 : float
        Characteristic mass of the MF logistic, in solar masses (Msun).
    MF_k : float
        MF logistic slope, dimensionless.
    CSF_A, CSF_B : float
        Low-mass and high-mass asymptotes of the companion-star-fraction
        logistic, dimensionless mean companion count (not bounded by 1).
    CSF_M0 : float
        Characteristic mass of the CSF logistic, in solar masses (Msun).
    CSF_k : float
        CSF logistic slope, dimensionless.
    CSF_max : float, optional
        Maximum companion star fraction, dimensionless mean companion
        count (not bounded by 1). Default 3.
    q_power : float, optional
        Mass-ratio power-law index, dimensionless. Default -0.4.
    q_min : float, optional
        Minimum mass ratio m_comp/m_prim, dimensionless. Default 0.01.
    companion_max : bool, optional
        If True, cap companion counts at CSF_max at all masses.
        Default False.
    """
    def __init__(self, MF_A, MF_B, MF_M0, MF_k,
                 CSF_A, CSF_B, CSF_M0, CSF_k,
                 CSF_max=3, q_power=-0.4, q_min=0.01, companion_max=False):
        super(MultiplicityLogistic, self).__init__(
            MF_amp=1.0, MF_power=0.0, CSF_amp=1.0, CSF_power=0.0,
            CSF_max=CSF_max, q_power=q_power, q_min=q_min,
            companion_max=companion_max)
        self.MF_A = float(MF_A)
        self.MF_B = float(MF_B)
        self.MF_M0 = float(MF_M0)
        self.MF_k = float(MF_k)
        self.CSF_A = float(CSF_A)
        self.CSF_B = float(CSF_B)
        self.CSF_M0 = float(CSF_M0)
        self.CSF_k = float(CSF_k)

        return

    def multiplicity_fraction(self, mass):
        """
        Multiplicity fraction as a logistic in log primary mass.
        Clipped to [0, 1].

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        mf : float or ndarray
            Multiplicity fraction, dimensionless, in [0, 1].
            Python float if ``mass`` is scalar, ndarray otherwise.
        """
        mf = _logistic_in_logm(mass, self.MF_A, self.MF_B, self.MF_M0, self.MF_k,
                               clip_min=0.0, clip_max=1.0)

        return mf

    def companion_star_fraction(self, mass):
        """
        Companion star fraction as a logistic in log primary mass.

        Clipped to [0, CSF_max] and raised to at least MF.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        csf : float or ndarray
            Companion star fraction, dimensionless mean companion
            count (not bounded by 1). Python float if ``mass`` is
            scalar, ndarray otherwise.
        """
        return_scalar = np.isscalar(mass)
        mass_arr = np.atleast_1d(np.asarray(mass, dtype=float))

        # Calculate the multiplicity fraction
        mf = _logistic_in_logm(mass_arr, self.MF_A, self.MF_B, self.MF_M0, self.MF_k,
                               clip_min=0.0, clip_max=1.0)

        # Calculate the companion star fraction
        csf = _logistic_in_logm(mass_arr, self.CSF_A, self.CSF_B, self.CSF_M0, self.CSF_k,
                                clip_min=0.0, clip_max=self.CSF_max)
        
        # Ensure the companion star fraction is at least the multiplicity fraction
        csf = np.maximum(csf, mf)

        if return_scalar:
            return float(csf[0])

        return csf


class MultiplicityUnresolvedOffner2023(MultiplicityLogistic):
    """
    Opt-in unresolved multiplicity derived from data in Offner et al. 2023 Table 1,
    including brown dwarfs.

    Scientifically preferred over the MultiplicityUnresolved default, but **not**
    the default (backwards compatibility). Companions only; for orbits
    use :class:`MultiplicityResolvedOffner2023`.

    Citation: Offner, S. S. R., Moe, M., Kratter, K. M., Sadavoy, S. I.,
    Jensen, E. L. N., & Tobin, J. J. 2023, in Protostars and Planets VII,
    ASP Conf. Ser. 534, 275 (`arXiv:2203.10066
    <https://arxiv.org/abs/2203.10066>`_; ADS
    `2023ASPC..534..275O
    <https://ui.adsabs.harvard.edu/abs/2023ASPC..534..275O>`_).
    Table 1 data: Zenodo `10.5281/zenodo.6628915
    <https://doi.org/10.5281/zenodo.6628915>`_.

    Notes
    -----
    MF and CSF are a **4-parameter logistic in log-mass**, fitted with
    equal weight to the geom-mean MF/CF columns of Table 1::

        f(M) = A + (B - A) / (1 + (M / M0)**(-k))

    with (A, B, M0, k) = (0.14, 0.99, 1.41, 1.25) for MF and
    (0.12, 2.35, 3.57, 0.96) for CSF/CF. The curve is C-infinity
    smooth (not a broken power law).
    The curve saturates at B ~ 1 for MF so A/B stars stay near the
    Raghavan/MDS/Sana points, and has a low-mass floor of A ~ 0.14.
    The MF is clipped to [0, 1] and the CSF is clipped to
    [0, CSF_max]. Higher-order multiples
    are allowed at all masses, including brown dwarfs. If
    ``companion_max`` is True, the number of companions is
    capped at CSF_max at all masses.

    SPISEA draws companions down to ``q_min`` (default 0.01), 
    so brown-dwarf secondaries of stellar primaries are generated. 

    **Mass-ratio index.** :meth:`q_power_at_mass` is an error-weighted
    logistic in log-mass fitted to Table 1 γ_trunc (1–100 au)::

        γ(M) = A + (B - A) / (1 + (M / M0)**(-k))

    with (A, B, M0, k) = (6.6, −1.77, 0.0651, 0.629). The q API
    is :meth:`draw_q` / :meth:`q_power_at_mass`; ``random_q`` is
    not supported. BD companions are still more equal-mass than
    solar-type companions. The err-weighted fit undershoots
    Fontanive 4.8 ± 2.2 (~3.3 at 0.033 Msun); that is the fit,
    not a bug.

    This class is companions-only (MF/CSF/q). Orbital methods
    (:meth:`~MultiplicityResolvedOffner2023.log_a_mean`,
    :meth:`~MultiplicityResolvedOffner2023.a_mean`,
    :meth:`~MultiplicityResolvedOffner2023.sigma_log_a`) live on
    :class:`MultiplicityResolvedOffner2023`.

    Parameters
    ----------
    MF_A, MF_B : float, optional
        Low- and high-mass MF logistic asymptotes, dimensionless.
        Defaults 0.14, 0.99 (Offner et al. 2023 Table 1,
        equal-weight geom-mean MF fit).
    MF_M0 : float, optional
        Characteristic mass of the MF logistic, in solar masses
        (Msun). Default 1.41.
    MF_k : float, optional
        MF logistic slope, dimensionless. Default 1.25.
    CSF_A, CSF_B : float, optional
        Low- and high-mass CSF logistic asymptotes, dimensionless
        mean companion count. Defaults 0.12, 2.35 (Table 1
        equal-weight geom-mean CF fit).
    CSF_M0 : float, optional
        Characteristic mass of the CSF logistic, in solar masses
        (Msun). Default 3.57.
    CSF_k : float, optional
        CSF logistic slope, dimensionless. Default 0.96.
    q_A, q_B : float, optional
        Low- and high-mass γ logistic asymptotes, dimensionless.
        Defaults 6.6, −1.77 (Table 1 γ_trunc, error-weighted).
    q_M0 : float, optional
        Characteristic mass of the γ logistic, in solar masses
        (Msun). Default 0.0651.
    q_k : float, optional
        γ logistic slope, dimensionless. Default 0.629.
    CSF_max : float, optional
        Maximum companion star fraction, dimensionless mean companion
        count (not bounded by 1). Default 3.
    q_min : float, optional
        Minimum mass ratio m_comp/m_prim, dimensionless, in [q_min, 1].
        Default 0.01.
    companion_max : bool, optional
        If True, cap companion counts at CSF_max at all masses.
        Default False.
    """
    def __init__(self, MF_A=0.14, MF_B=0.99, MF_M0=1.41, MF_k=1.25,
                 CSF_A=0.12, CSF_B=2.35, CSF_M0=3.57, CSF_k=0.96,
                 q_A=6.6, q_B=-1.77, q_M0=0.0651, q_k=0.629,
                 CSF_max=3, q_min=0.01, companion_max=False):

        super(MultiplicityUnresolvedOffner2023, self).__init__(
            MF_A=MF_A, MF_B=MF_B, MF_M0=MF_M0, MF_k=MF_k,
            CSF_A=CSF_A, CSF_B=CSF_B, CSF_M0=CSF_M0, CSF_k=CSF_k,
            CSF_max=CSF_max, q_min=q_min,
            companion_max=companion_max)
        self.q_A = float(q_A)
        self.q_B = float(q_B)
        self.q_M0 = float(q_M0)
        self.q_k = float(q_k)

        return

    def random_q(self, mass, rng=None, n_comp=1):
        """
        Not supported on Offner classes.

        Use :meth:`draw_q` with a primary mass.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun). Unused; this method
            always raises.
        rng : numpy.random.Generator, optional
            Unused; this method always raises.
        n_comp : int, optional
            Unused; this method always raises.

        Raises
        ------
        TypeError
            Always. Offner q draws require ``draw_q(mass, rng=...)``.
        """
        raise TypeError(
            "Offner multiplicity does not support random_q; "
            "use draw_q(mass, rng=...)")

    def q_power_at_mass(self, mass):
        """
        Mass-ratio power-law index γ(M), P(q) ∝ q^γ on [q_min, 1].

        Error-weighted logistic in log-mass fitted to Table 1
        γ_trunc (1–100 au). Not an interpolation of the Table 1
        knots::

            γ(M) = A + (B - A) / (1 + (M / M0)**(-k))

        with (A, B, M0, k) = (6.6, −1.77, 0.0651, 0.629).
        Undershoots Fontanive 4.8 ± 2.2 (~3.3 at 0.033 Msun).

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        gamma : float or ndarray
            Mass-ratio power-law index γ, dimensionless.
            Python float if ``mass`` is scalar, ndarray otherwise.
        """
        gamma = _logistic_in_logm(
            mass, self.q_A, self.q_B, self.q_M0, self.q_k)
        
        return gamma

    def draw_q(self, mass, rng=None, n_comp=1):
        """
        Draw companion mass ratios q = m_comp / m_prim.

        Uses ``q_power_at_mass(mass)`` for the Table 1 γ logistic.
        ``mass`` is required. All primaries are drawn in one RNG
        call; γ(M) is already mass-dependent.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).
        rng : numpy.random.Generator, optional
            Random generator. Default ``numpy.random.default_rng()``.
        n_comp : int, optional
            Companions per primary. Default 1.

        Returns
        -------
        q : float or ndarray
            Mass ratios in [q_min, 1]. Scalar if ``mass`` is scalar
            and ``n_comp`` is 1; shape ``(n_mass,)`` or
            ``(n_mass, n_comp)`` otherwise.
        """
        if rng is None:
            rng = np.random.default_rng()
        n_comp = int(n_comp)
        return_scalar = np.isscalar(mass) and n_comp == 1
        mass_arr = np.atleast_1d(np.asarray(mass, dtype=float))
        q_pow = self.q_power_at_mass(mass_arr)
        x = rng.random((len(mass_arr), n_comp))
        q_values = _q_from_powerlaw(x, q_pow, self.q_min)
        return _finalize_q_draw(q_values, return_scalar, n_comp)


class MultiplicityResolvedOffner2023(MultiplicityUnresolvedOffner2023,
                                     _ResolvedOrbitalMixin):
    """
    Opt-in resolved Offner et al. 2023 multiplicity.

    Same MF/CSF/q as :class:`MultiplicityUnresolvedOffner2023`
    (Table 1 logistic MF/CSF, error-weighted γ logistic, Table 1
    companion-cut caveat). Adds mass-dependent separations.
    Higher-order multiples are allowed at all masses.
    Scientifically preferred over :class:`MultiplicityResolvedDK`
    in the brown-dwarf regime, but **not** the default.

    Notes
    -----
    :meth:`log_a_mean` / :meth:`a_mean` are a smooth broken power
    law in log10 a vs log10 M (FGK-pulled, s = 0.1 dex),
    C-infinity via a stable logcosh (not ``log(cosh x)``)::

        v = log10(M / Mp),   yp = log10(μp)
        log10 a = yp + 0.5*(αL+αR)*v + 0.5*(αR-αL)*s * logcosh(v/s)

    with μp = 44.46 AU, Mp = 0.819 Msun, αL = 1.005, αR = −0.308,
    s = 0.10. Linear-space a is clipped above 0.1 AU.

    :meth:`sigma_log_a` is a 2-parameter logistic pinned at
    0.7 / 1.5::

        σ(M) = 0.7 + 0.8 / (1 + (M / 0.354)**(-6.05))

    :meth:`log_semimajoraxis` draws log10(a/AU) from a truncated
    lognormal with those loc / scale values. Brown-dwarf binaries
    peak near a few AU (μ(0.033) ≈ 2 AU). Truncation is
    0.01–2000 AU, same limits as :class:`MultiplicityResolvedDK`.

    Eccentricity and Keplerian angles still follow Duchêne & Kraus
    (2013): f(e) = 2e, random inclination and angles. Same mixin as
    :class:`MultiplicityResolvedDK`.

    Parameters
    ----------
    MF_A, MF_B : float, optional
        Low- and high-mass MF logistic asymptotes, dimensionless.
        Defaults 0.14, 0.99 (Offner et al. 2023 Table 1,
        equal-weight geom-mean MF fit).
    MF_M0 : float, optional
        Characteristic mass of the MF logistic, in solar masses
        (Msun). Default 1.41.
    MF_k : float, optional
        MF logistic slope, dimensionless. Default 1.25.
    CSF_A, CSF_B : float, optional
        Low- and high-mass CSF logistic asymptotes, dimensionless
        mean companion count. Defaults 0.12, 2.35 (Table 1
        equal-weight geom-mean CF fit).
    CSF_M0 : float, optional
        Characteristic mass of the CSF logistic, in solar masses
        (Msun). Default 3.57.
    CSF_k : float, optional
        CSF logistic slope, dimensionless. Default 0.96.
    q_A, q_B : float, optional
        Low- and high-mass γ logistic asymptotes, dimensionless.
        Defaults 6.6, −1.77 (Table 1 γ_trunc, error-weighted).
    q_M0 : float, optional
        Characteristic mass of the γ logistic, in solar masses
        (Msun). Default 0.0651.
    q_k : float, optional
        γ logistic slope, dimensionless. Default 0.629.
    a_mup : float, optional
        Smooth-broken-PL pivot μ(a), in AU. Default 44.46.
    a_mp : float, optional
        Smooth-broken-PL pivot mass, in solar masses (Msun).
        Default 0.819.
    a_alphaL, a_alphaR : float, optional
        Smooth-broken-PL slopes in log10 a vs log10 M,
        dimensionless. Defaults 1.005, −0.308.
    a_s : float, optional
        Smooth-broken-PL smoothing scale, in dex. Default 0.10.
    a_min : float, optional
        Minimum characteristic a, in AU. Default 0.1.
    sig_A, sig_B : float, optional
        Low- and high-mass σ(log10 a) logistic values, in dex.
        Defaults 0.7, 1.5 (Table 2 pins).
    sig_M0 : float, optional
        Characteristic mass of the σ logistic, in solar masses
        (Msun). Default 0.354.
    sig_k : float, optional
        σ logistic slope, dimensionless. Default 6.05.
    CSF_max : float, optional
        Maximum companion star fraction, dimensionless mean companion
        count (not bounded by 1). Default 3.
    q_min : float, optional
        Minimum mass ratio m_comp/m_prim, dimensionless. Default 0.01.
    companion_max : bool, optional
        If True, cap companion counts at CSF_max at all masses.
        Default False.
    sep_sig : array_like, optional
        Table 2 σ(log10 a) knots, in dex. Default (0.7, 1.3, 1.5).
    sep_sig_mass : array_like, optional
        Table 2 primary-mass knots for the σ comparison plot, in
        solar masses (Msun). Defaults are geom-mean M1 of the
        late-M, early-M, and FGK bins.
    """
    def __init__(self, MF_A=0.14, MF_B=0.99, MF_M0=1.41, MF_k=1.25,
                 CSF_A=0.12, CSF_B=2.35, CSF_M0=3.57, CSF_k=0.96,
                 q_A=6.6, q_B=-1.77, q_M0=0.0651, q_k=0.629,
                 a_mup=44.46, a_mp=0.819, a_alphaL=1.005,
                 a_alphaR=-0.308, a_s=0.10, a_min=0.1,
                 sig_A=0.7, sig_B=1.5, sig_M0=0.354, sig_k=6.05,
                 CSF_max=3, q_min=0.01, companion_max=False,
                 sep_sig=(0.7, 1.3, 1.5),
                 sep_sig_mass=(np.sqrt(0.075 * 0.15),
                               np.sqrt(0.3 * 0.6),
                               np.sqrt(0.75 * 1.25))):
        super(MultiplicityResolvedOffner2023, self).__init__(
            MF_A=MF_A, MF_B=MF_B, MF_M0=MF_M0, MF_k=MF_k,
            CSF_A=CSF_A, CSF_B=CSF_B, CSF_M0=CSF_M0, CSF_k=CSF_k,
            q_A=q_A, q_B=q_B, q_M0=q_M0, q_k=q_k,
            CSF_max=CSF_max, q_min=q_min,
            companion_max=companion_max)
        self.a_mup = float(a_mup)
        self.a_mp = float(a_mp)
        self.a_alphaL = float(a_alphaL)
        self.a_alphaR = float(a_alphaR)
        self.a_s = float(a_s)
        self.a_min = float(a_min)
        self.sig_A = float(sig_A)
        self.sig_B = float(sig_B)
        self.sig_M0 = float(sig_M0)
        self.sig_k = float(sig_k)
        # Table 2 σ knots for the comparison plot; draws use sigma_log_a.
        self.sep_sig_mass = np.array(sep_sig_mass, dtype=float)
        self.sep_sig = np.array(sep_sig, dtype=float)

        return

    def log_a_mean(self, mass):
        """
        Characteristic log10(a/AU) from the smooth broken power law.

        FGK-pulled, s = 0.1 dex, C-infinity (stable logcosh; not
        ``log(cosh x)`` and not a hard ``where`` break)::

            v = log10(M / Mp),   yp = log10(μp)
            log10 a = yp + 0.5*(αL+αR)*v + 0.5*(αR-αL)*s * logcosh(v/s)

        with μp = 44.46 AU, Mp = 0.819 Msun, αL = 1.005,
        αR = −0.308, s = 0.10. Linear-space a is clipped to 0.1 AU.
        Uses ``logcosh x = |x| + log(1 + e**(-2|x|)) - log 2``.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        log_a_mean : float or ndarray
            Characteristic log10(a / 1 AU) in dex (not ln, not AU).
            Python float if ``mass`` is scalar, ndarray otherwise.
        """
        log_a_mean = _smooth_broken_loglog(
            mass, self.a_mup, self.a_mp,
            self.a_alphaL, self.a_alphaR, self.a_s,
            a_min=self.a_min)

        return log_a_mean

    def a_mean(self, mass):
        """
        Characteristic μ(a) in AU, ``10 ** log_a_mean(mass)``.

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        a_mean : float or ndarray
            Characteristic separation μ(a) in AU. Python float if
            ``mass`` is scalar, ndarray otherwise.
        """
        log_a = self.log_a_mean(mass)

        if np.isscalar(log_a):
            a_mean = 10.0 ** log_a
            return a_mean

        a_mean = 10.0 ** np.asarray(log_a, dtype=float)

        return a_mean

    def sigma_log_a(self, mass):
        """
        σ(log10 a) from a 2-parameter logistic in log-mass.

        Floors/ceilings pinned at 0.7 / 1.5; clipped to ≥ 0.1::

            σ(M) = 0.7 + 0.8 / (1 + (M / 0.354)**(-6.05))

        i.e. (A, B, M0, k) = (0.7, 1.5, 0.354, 6.05).

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).

        Returns
        -------
        sigma_log_a : float or ndarray
            Standard deviation of log10(a / 1 AU), in dex.
            Python float if ``mass`` is scalar, ndarray otherwise.
        """
        sigma_log_a = _logistic_in_logm(
            mass, self.sig_A, self.sig_B, self.sig_M0, self.sig_k,
            clip_min=0.1)

        return sigma_log_a

    def log_semimajoraxis(self, mass, rng=None):
        """
        Draw log10(a/AU) from a mass-dependent truncated lognormal.

        Uses ``loc = log_a_mean(mass)`` and
        ``scale = sigma_log_a(mass)``. Truncated so a is between
        0.01 AU and 2000 AU (same limits as
        :class:`MultiplicityResolvedDK`).

        Parameters
        ----------
        mass : float or array_like
            Primary mass in solar masses (Msun).
        rng : numpy.random.Generator, optional
            Random number generator passed to ``truncnorm.rvs``.
            Default is a new ``numpy.random.default_rng()``.

        Returns
        -------
        log_semimajoraxis : ndarray
            Drawn log10(a / 1 AU) in dex (not ln, not AU), truncated
            so a is between 0.01 AU and 2000 AU.
        """
        if rng is None:
            rng = np.random.default_rng()
        mass = np.atleast_1d(np.asarray(mass, dtype=float))
        log_a_mean = np.atleast_1d(np.asarray(self.log_a_mean(mass),
                                             dtype=float))
        log_a_std = np.atleast_1d(np.asarray(self.sigma_log_a(mass),
                                            dtype=float))

        log_a_lower = np.log10(0.01)
        log_a_upper = np.log10(2000)
        a_lower_std = (log_a_lower - log_a_mean) / log_a_std
        a_upper_std = (log_a_upper - log_a_mean) / log_a_std

        # Draw the log10(a / 1 AU) from a truncated normal
        log_a = truncnorm.rvs(a_lower_std, a_upper_std,
                             loc=log_a_mean, scale=log_a_std,
                             random_state=rng)

        return log_a


def _piecewise_powerlaw(mass, mass_limits, amps, powers, clip_min=None,
                        clip_max=None):
    """
    Evaluate y = A_i * M**alpha_i on mass segments.

    Segment i applies for mass_limits[i] <= M < mass_limits[i+1].
    The first segment also covers M below the lowest limit; the last
    segment is closed on the right.

    Parameters
    ----------
    mass : float or array_like
        Primary mass in solar masses (Msun).
    mass_limits : array_like
        Segment edges in solar masses (Msun), length N+1, strictly
        increasing.
    amps : array_like
        Length-N amplitudes A_i, in units of y / Msun**alpha_i.
    powers : array_like
        Length-N power-law indices alpha_i, dimensionless.
    clip_min, clip_max : float or None, optional
        Optional lower/upper clips on y, in the same units as y.
        ``None`` means no clip on that side.

    Returns
    -------
    y : float or ndarray
        Piecewise power-law value, in the same units as
        ``amps * mass**powers``. Python float if ``mass`` is scalar,
        ndarray otherwise.
    """
    return_scalar = np.isscalar(mass)
    mass_arr = np.atleast_1d(np.asarray(mass, dtype=float))
    out = np.empty(mass_arr.shape, dtype=float)
    nseg = len(amps)

    # Evaluate the piecewise power-law for each segment
    for i in range(nseg):
        lo = mass_limits[i]
        hi = mass_limits[i + 1]

        # Determine the mask for the current segment
        if i == 0:
            mask = mass_arr < hi
        elif i == nseg - 1:
            mask = mass_arr >= lo
        else:
            mask = (mass_arr >= lo) & (mass_arr < hi)
        out[mask] = amps[i] * np.power(mass_arr[mask], powers[i])

    # Apply the clips
    if clip_min is not None:
        out = np.maximum(out, clip_min)

    if clip_max is not None:
        out = np.minimum(out, clip_max)

    # Return the result
    if return_scalar:
        return float(out[0])

    return out


def _logistic_in_logm(mass, A, B, M0, k, clip_min=None, clip_max=None):
    """
    Evaluate y = A + (B - A) / (1 + (M / M0)**(-k)).

    This is a logistic in log-mass: as M -> 0+, y -> A; as M -> inf,
    y -> B.

    Parameters
    ----------
    mass : float or array_like
        Primary mass in solar masses (Msun).
    A, B : float
        Low-mass and high-mass asymptotes, in the same units as y.
        Dimensionless for MF and γ; mean companion count for CSF;
        dex for σ(log10 a).
    M0 : float
        Characteristic mass in solar masses (Msun).
    k : float
        Logistic slope, dimensionless.
    clip_min, clip_max : float or None, optional
        Optional lower/upper clips on y, in the same units as y.
        ``None`` means no clip on that side.

    Returns
    -------
    y : float or ndarray
        Logistic value in the same units as ``A`` and ``B``.
        Python float if ``mass`` is scalar, ndarray otherwise.
    """
    return_scalar = np.isscalar(mass)
    mass_arr = np.atleast_1d(np.asarray(mass, dtype=float))

    # Evaluate the logistic in log-mass
    out = A + (B - A) / (1.0 + np.power(mass_arr / M0, -k))

    # Apply the clips
    if clip_min is not None:
        out = np.maximum(out, clip_min)
    if clip_max is not None:
        out = np.minimum(out, clip_max)

    # Return the result
    if return_scalar:
        return float(out[0])

    return out


def _logcosh(x):
    """
    Numerically stable log(cosh(x)).

    Uses |x| + log1p(exp(-2|x|)) - log(2) rather than np.log(np.cosh(x)),
    which overflows for |x| ≳ 700.

    Parameters
    ----------
    x : float or array_like
        Argument of cosh, dimensionless (for the smooth broken power
        law this is v/s, with v and s in dex).

    Returns
    -------
    logcosh_x : float or ndarray
        log(cosh(x)), dimensionless. Same shape as ``x``.
    """
    ax = np.abs(np.asarray(x, dtype=float))

    # Calculate the log(cosh(x))
    logcosh_x = ax + np.log1p(np.exp(-2.0 * ax)) - np.log(2.0)

    # Return the result
    return logcosh_x


def _smooth_broken_loglog(mass, mup, Mp, alpha_L, alpha_R, s, a_min=0.1):
    """
    Smooth broken power law in log10(a) vs log10(M).

        v  = log10(M / Mp)
        yp = log10(mup)
        log10(a) = yp + 0.5*(αL+αR)*v + 0.5*(αR-αL)*s * logcosh(v/s)

    ``s`` is the smoothing scale in dex (C-infinity; logcosh). The
    linear-space value is clipped to ``a_min``.

    Parameters
    ----------
    mass : float or array_like
        Primary mass in solar masses (Msun).
    mup : float
        Characteristic separation at the break mass, in AU.
    Mp : float
        Break mass in solar masses (Msun).
    alpha_L, alpha_R : float
        Power-law indices below and above ``Mp``, dimensionless.
    s : float
        Smoothing scale in dex of log10(M / 1 Msun).
    a_min : float, optional
        Minimum linear-space separation in AU. Default 0.1 AU.

    Returns
    -------
    log_a : float or ndarray
        log10(a / 1 AU) in dex (not ln, not AU). Python float if
        ``mass`` is scalar, ndarray otherwise.
    """
    return_scalar = np.isscalar(mass)
    mass_arr = np.atleast_1d(np.asarray(mass, dtype=float))
    v = np.log10(mass_arr / float(Mp))
    yp = np.log10(float(mup))

    # Calculate the log10(a / 1 AU) using the smooth broken power law
    log_a = (yp
             + 0.5 * (alpha_L + alpha_R) * v
             + 0.5 * (alpha_R - alpha_L) * s * _logcosh(v / s))
    a = np.maximum(10.0 ** log_a, float(a_min))
    log_a = np.log10(a)

    # Return the result
    if return_scalar:
        return float(log_a[0])

    return log_a


def _q_from_powerlaw(x, q_pow, q_min):
    """
    Inverse CDF of P(q) ∝ q**q_pow for q_min <= q <= 1.

    ``q_pow`` may be a scalar or an array broadcastable to ``x``.
    The q_pow = -1 (b = 0) limit is q = q_min**(1 - x).

    Parameters
    ----------
    x : float or array_like
        Uniform random draw, dimensionless, in [0, 1].
    q_pow : float or array_like
        Mass-ratio power-law index γ, dimensionless. Broadcastable
        to ``x``.
    q_min : float
        Minimum mass ratio m_comp/m_prim, dimensionless, in (0, 1].

    Returns
    -------
    q : ndarray
        Companion mass ratio m_comp/m_prim, dimensionless, in
        [q_min, 1]. Same shape as the broadcast of ``x`` and ``q_pow``.
    """
    x = np.asarray(x, dtype=float)
    q_pow = np.asarray(q_pow, dtype=float)

    # Broadcast the arrays if necessary
    if x.ndim > q_pow.ndim:
        q_pow = q_pow.reshape(q_pow.shape + (1,) * (x.ndim - q_pow.ndim))

    b = 1.0 + q_pow
    b, x = np.broadcast_arrays(b, x)

    # Create an empty array to store the result
    out = np.empty(x.shape, dtype=float)

    # Determine the mask for values near zero
    near_zero = np.abs(b) < 1e-12

    # Determine the mask for values far from zero
    ok = ~near_zero

    # Calculate the mass ratio for values near zero
    if np.any(near_zero):
        out[near_zero] = q_min ** (1.0 - x[near_zero])
    if np.any(ok):
        out[ok] = (x[ok] * (1.0 - q_min ** b[ok]) + q_min ** b[ok]) ** (1.0 / b[ok])

    # Return the result
    return out


def _finalize_q_draw(q_values, return_scalar, n_comp):
    """
    Shape a (n_mass, n_comp) q array for :meth:`draw_q`.

    Parameters
    ----------
    q_values : ndarray
        Mass ratios, dimensionless, shape (n_mass, n_comp).
    return_scalar : bool
        If True, return a Python float.
    n_comp : int
        Companions per primary, integer count.

    Returns
    -------
    q : float or ndarray
        Scalar, shape (n_mass,), or shape (n_mass, n_comp).
    """
    q_values = np.asarray(q_values, dtype=float)
    if return_scalar:
        return float(q_values[0, 0])
    if n_comp == 1:
        return q_values[:, 0]
    return q_values

