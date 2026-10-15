##################################################
# J. Lu
#
# Original code was taken from libimf package written by Jan Pflamm-Altenburg
# and has been modified only marginally. The libimf code was licensed under
# a GNU General Public License.
#
# When I use this code, I should cite Pflamm-Altenburg & Kroupa 2006
#
# Unfortunately, the code was almost completely un-commented, so all
# comments are mine. I have also added substantially to the code to
# make more convenient and general purpose functions.
#
##################################################

import numpy as np
import time
import pdb
import logging

log = logging.getLogger('imf')

class IMF(object):
    """
    The IMF base class. The mass sampling and multiplicity
    implementation is here.


    Parameters
    ----------
    massLimits : 2 element array; optional
        Define the minimum and maximum stellar masses in the IMF, in
        solar masses. First element is taken as the min, second element
        the max (e.g. `massLimits` = [min_mass, max_mass]).

    multiplicity : Multiplicity object or None
        If None, no multiplicity is assumed. Otherwise, use
        multiplicity object to create multiple star systems.

    seed : int, optional
        Seed for the random generator numpy.random.default_rng(seed).
        All random functions in the class will use this generator, by default None.
        Behavior:
        ::

            imf = IMF(..., seed=42)
            result1 = imf.generate_cluster()
            result2 = imf.generate_cluster()
            imf = IMF(..., seed=42)
            result3 = imf.generate_cluster()
            result4 = imf.generate_cluster()

        result1==result3, result2==result4, but result1≠result2, result3≠result4.
        This is the same behavior as
        ::

            rng = np.random.default_rng(seed=42)
            result1 = rng.random(1)
            result2 = rng.random(1)
            rng = np.random.default_rng(seed=42)
            result3 = rng.random(1)
            result4 = rng.random(1)

        If identical output is desired over each run, the random state can be reset before running the function, e.g.
        ::

            imf.rng = np.random.default_rng(seed=42)
            result1 = imf.generate_cluster()
            imf.rng = np.random.default_rng(seed=42)
            result2 = imf.generate_cluster()

        In this case, result1==result2

    Notes
    -----
    Code author: J. Lu.

    Original code was taken from libimf package written by Jan Pflamm-Altenburg
    (`Pflamm-Altenburg & Kroupa 2006 <https://ui.adsabs.harvard.edu/abs/2006MNRAS.373..295P/abstract>`_)
    and has been modified only marginally, though more convinient and general purpose
    functions have been added. The libimf code was licensed under
    a GNU General Public License.

    """
    def __init__(self, massLimits=np.array([0.01,150]), multiplicity=None, seed=None):
        self._multi_props = multiplicity
        self._mass_limits = np.atleast_1d(massLimits)
        self.rng = np.random.default_rng(seed)

        if multiplicity:
            self.make_multiples = True
        else:
            self.make_multiples = False

        return

    def generate_cluster(self, totalMass):
        """
        Generate a cluster of stellar systems with the specified IMF.

        Randomly sample from an IMF with specified mass
        limits until the desired total mass is reached. The maximum
        stellar mass is not allowed to exceed the total cluster mass.
        The simulated total mass will not be exactly equivalent to the
        desired total mass; but we will take one star above or below
        (whichever brings us closer to the desired total) the desired
        total cluster mass break point.

        Primary stars are sampled from the IMF, companions are generated
        based on the multiplicity properties provided.

        Parameters
        ----------
        totalMass : float
            The total mass of the cluster (including companions) in solar masses.

        Returns
        -------
        masses : numpy float array
            Array of primary star masses.

        isMultiple : numpy boolean array
            Array of booleans with True for each primary star that is in a multiple
            system and False for each single star.

        companionMasses : numpy masked array
            Masked array of companion masses. Each row corresponds to a primary star, and each column corresponds to a companion. The mask is True for entries that are not valid companions (e.g. for single stars or for companions that are below the minimum mass limit).

        systemMasses : numpy float array
            Array of total system masses (primary + companions) for each primary star.

        """
        initial_mass_limit = self._mass_limits[-1]

        if (self._mass_limits[-1] > totalMass):
            log.info('sample_imf: Setting maximum allowed mass to %d' %
                      (totalMass))
            self._mass_limits[-1] = totalMass

        # Estimate the mean number of stars expected.
        self.normalize(totalMass)
        mean_number = self.int_xi(self._mass_limits[0], self._mass_limits[-1])
        newStarCount = np.round(mean_number)
        if self._multi_props == None:
            newStarCount *= 1.1

        # Generate output arrays.
        masses = np.array([], dtype=float)
        isMultiple = np.array([], dtype=bool)
        # compMasses = {} # Hashmap for index -> compMasses for faster lookup
        compMasses = []
        systemMasses = np.array([], dtype=float)

        # Loop through and add stars to the cluster until we get to
        # the desired total cluster mass.
        totalMassTally = 0
        loopCnt = 0

        # start_while = time.time()
        while totalMassTally < totalMass:
            # Generate a random number array.
            uniX = self.rng.random(newStarCount.astype(int))
            # Convert into the IMF from the inverted CDF
            newMasses = self.dice_star_cl(uniX)

            # Testing for Nans produced in masses
            test = np.isnan(newMasses)
            if np.sum(test) > 0:
                raise ValueError('Nan detected in cluster mass')

            # Dealing with multiplicity
            if self._multi_props:
                # newCompMasses = np.empty((len(newMasses),), dtype=object)
                # newCompMasses.fill([])
                # Determine the multiplicity of every star
                MF = self._multi_props.multiplicity_fraction(newMasses)
                CSF = self._multi_props.companion_star_fraction(newMasses)

                newIsMultiple = self.rng.random(newStarCount.astype(int)) < MF

                # Function to calculate multiple systems more efficiently
                # start_calc = time.time()
                newCompMasses, newSystemMasses, newIsMultiple = self.calc_multi(newMasses, newIsMultiple, CSF, MF)
                # end_calc = time.time()
                # print('Time taken for calc_multi: ', end_calc - start_calc)
                newTotalMassTally = newSystemMasses.sum()
                isMultiple = np.append(isMultiple, newIsMultiple)
                systemMasses = np.append(systemMasses, newSystemMasses)
                compMasses.append(newCompMasses)

            else:
                newTotalMassTally = newMasses.sum()
            # end_while = time.time()
            # print('Time taken for while loop: ', end_while - start_while)
            # Append to our primary masses array
            masses = np.append(masses, newMasses)

            if (loopCnt >= 0):
                log.info('sample_imf: Loop %d added %.2e Msun to previous total of %.2e Msun' %
                         (loopCnt, newTotalMassTally, totalMassTally))

            totalMassTally += newTotalMassTally
            newStarCount = mean_number * 0.1  # increase by 20% each pass
            loopCnt += 1

        # Make a running sum of the system masses
        if self._multi_props:
            # Concatenate the companion masses
            if len(compMasses) > 1:
                max_cols = max(compMass.shape[1] for compMass in compMasses)

                # Pad each array to have the same number of columns
                padded_arrays = [
                    np.ma.masked_all((compMass.shape[0], max_cols)) for compMass in compMasses
                ]

                for i, compMass in enumerate(compMasses):
                    padded_arrays[i][:, :compMass.shape[1]] = compMass

                # Vertically stack the padded arrays
                compMasses = np.ma.vstack(padded_arrays)

            else:
                compMasses = compMasses[0]

            # Make a running sum of the system masses
            massCumSum = systemMasses.cumsum()
        else:
            massCumSum = masses.cumsum()

        # Find the index where we are closest to the desired
        # total mass.
        idx = np.abs(massCumSum - totalMass).argmin()

        masses = masses[:idx+1]

        if self._multi_props:
            systemMasses = systemMasses[:idx+1]
            isMultiple = isMultiple[:idx+1]
            compMasses = compMasses[:idx+1]
        else:
            isMultiple = np.zeros(len(masses), dtype=bool)
            systemMasses = masses

        self._mass_limits[-1] = initial_mass_limit

        return (masses, isMultiple, compMasses, systemMasses)

    def calc_multi(self, newMasses, newIsMultiple, CSF, MF):
        """
        Helper function to calculate multiples more efficiently.
        Companion counts and companion-mass draws (including brown-dwarf
        q distributions and binaries-only BD systems) are delegated to
        the multiplicity object.
        """
        return self._multi_props.draw_companion_masses(
            newMasses, newIsMultiple, CSF, MF,
            rng=self.rng, mass_min=self._mass_limits[0])


class IMF_broken_powerlaw(IMF):
    """
    Initialize a multi-part power-law with N parts. Each part of the
    power-law is described with a probability density function:

        P(m) ∝ m ** power[n]

    for mass_limits[n] < m <= mass_limits[n+1].

    Parameters
    ----------
    mass_limits : numpy array
        Array of length (N + 1) with lower and upper limits of
        the power-law segments.

    powers : numpy array
        Array of length N that contains the powers for each
        power-law segment.

    multiplicity : Multiplicity object or None
        If None, no multiplicity is assumed. Otherwise, use
        multiplicity object to create multiple star systems.
    """
    def __init__(self, mass_limits, powers, multiplicity=None, seed=None):
        super().__init__(massLimits=mass_limits, multiplicity=multiplicity, seed=seed)
        powers = np.atleast_1d(powers)
        if len(mass_limits) != len(powers) + 1:
            msg = 'Incorrect specification of multi-part powerlaw.\n'
            msg += '    len(massLimts) != len(powers)+1\n'
            msg += '    len(massLimits) = \n' + str(len(mass_limits))
            msg += '    len(powers) = \n' + str(len(powers))

            raise RuntimeError(msg)
        mass_limits = np.atleast_1d(mass_limits)
        self._m_limits_low = mass_limits[0:-1]
        self._m_limits_high = mass_limits[1:]
        self._powers = np.atleast_1d(powers)

        # Calculate the coeffs to make the function continuous
        nterms = len(self._powers)
        coeffs = np.ones(nterms, dtype=float)

        # First term is just 1.0
        # Subsequent terms are products of previous terms and then some.
        for i in range(1, nterms):
            y = self._m_limits_low[i] ** self._powers[i-1]
            z = self._m_limits_low[i] ** self._powers[i]

            coeffs[i] *= coeffs[i-1] * y / z

        self.nterms = nterms
        self.coeffs = coeffs
        self.k = 1

    def xi(self, m):
        """
        Probability density describing the IMF.

        Input:
        m - mass of a star

        Output:
        xi - probability of measuring that mass.
        """
        returnFloat = isinstance(m, (float, np.floating))

        m = np.atleast_1d(m)

        # Temporary arrays
        y = np.zeros(len(m), dtype=float)
        z = np.ones(len(m), dtype=float)

        # Loop through the different segments of the power law.
        for i in range(self.nterms): # For i = 0 --> n, where n is the number of intervals
            aux = m - self._m_limits_low[i] #---Should this be i - 1?

            # Only continue for those entries that are in later segments
            idx = np.where(aux >= 0)[0]

            # Maybe we are all done?
            if len(idx) == 0:
                break

            m_tmp = m[idx]
            aux_tmp = aux[idx]

            y_i = gamma_closed(m_tmp, self._m_limits_low[i], self._m_limits_high[i])
            y_i *= self.coeffs[i] * m_tmp**self._powers[i]

            # Save results into the y array
            y[idx] += y_i

            z *= delta(m - self._m_limits_high[i])

        xi = self.k * z * y

        if returnFloat:
            return xi[0]
        else:
            return xi

    def m_xi(self, m):
        """
        Mass-weighted probability m*xi
        """
        returnFloat = isinstance(m, (float, np.floating))
        m = np.atleast_1d(m)
        mxi = np.zeros(len(m), dtype=float)

        # Temporary arrays
        y = np.zeros(len(m), dtype=float)
        z = np.ones(len(m), dtype=float)

        # Loop through the different segments of the power law.
        for i in range(self.nterms): # For i = 0 --> n, where n is the number of intervals
            aux = m - self._m_limits_low[i] #---Should this be i - 1?

            # Only continue for those entries that are in later segments
            idx = np.where(aux >= 0)[0]

            # Maybe we are all done?
            if len(idx) == 0:
                break

            m_tmp = m[idx]
            aux_tmp = aux[idx]

            y_i = gamma_closed(m_tmp, self._m_limits_low[i], self._m_limits_high[i])
            y_i *= self.coeffs[i] * m_tmp**(self._powers[i] + 1)

            # Save results into the y array
            y[idx] += y_i

            z *= delta(m - self._m_limits_high[i])

        mxi = self.k * z * y

        if returnFloat:
            return mxi[0]
        else:
            return mxi


    def getProbabilityBetween(self, massLo, massHi):
        """Return the integrated probability between some low and high
        mass value.
        """
        return self.int_xi(massLo, massHi)

    def int_xi(self, massLo, massHi):
        """Return the integrated probability between some low and high
        mass value.
        """
        return self.prim_xi(massHi) - self.prim_xi(massLo)

    def getMassBetween(self, massLo, massHi):
        """Return the integrated mass between some low and high
        mass value.
        """
        return self.int_mxi(massLo, massHi)

    def int_mxi(self, massLo, massHi):
        """Return the integrated total mass between some low and high stellar
        mass value. Be sure to normalize the IMF instance beforehand.
        """
        return self.prim_mxi(massHi) - self.prim_mxi(massLo)

    def prim_xi(self, a):
        """
        Helper function
        """
        returnFloat = isinstance(a, (float, np.floating))

        a = np.atleast_1d(a)
        val = np.zeros(len(a), dtype=float)

        for i in range(len(val)):
            t1 = theta_open(a[i] - self._m_limits_high) * self.coeffs
            t2 = prim_power(self._m_limits_high, self._powers)
            t3 = prim_power(self._m_limits_low, self._powers)
            y1 = (t1 * (t2 - t3)).sum()

            t1 = gamma_closed(a[i], self._m_limits_low, self._m_limits_high)
            t1 *= self.coeffs
            t2 = prim_power(a[i], self._powers)
            t3 = prim_power(self._m_limits_low, self._powers)
            y2 = (t1 * (t2 - t3)).sum()

            val[i] = self.k * (y1 + y2)

        if returnFloat:
            return val[0]
        else:
            return val

    def prim_mxi(self, a):
        """
        Helper function
        """
        returnFloat = isinstance(a, (float, np.floating))

        a = np.atleast_1d(a)
        val = np.zeros(len(a), dtype=float)

        for i in range(len(val)):
            t1 = theta_open(a[i] - self._m_limits_high) * self.coeffs
            t2 = prim_power(self._m_limits_high, self._powers+1)
            t3 = prim_power(self._m_limits_low, self._powers+1)
            y1 = (t1 * (t2 - t3)).sum()

            t1 = gamma_closed(a[i], self._m_limits_low, self._m_limits_high)
            t1 *= self.coeffs
            t2 = prim_power(a[i], self._powers+1)
            t3 = prim_power(self._m_limits_low, self._powers+1)
            y2 = (t1 * (t2 - t3)).sum()

            val[i] = self.k * (y1 + y2)

        if returnFloat:
            return val[0]
        else:
            return val

    def normalize(self, Mcl, Mmin=None, Mmax=None):
        """
        Normalize the IMF to a total cluster mass within a specified
        minimum and maximum stellar mass range.
        """
        self.k = 1.0
        self.Mcl = Mcl

        if Mmax == None:
            Mmax = self._m_limits_high[-1]

        if Mmin == None:
            Mmin = self._m_limits_low[0]

        if Mmax > Mcl:
            Mmax = Mcl

        if Mmax > self._m_limits_high[-1]:
            Mmax = self._m_limits_high[-1]

        if Mmin < self._m_limits_low[0]:
            Mmin = self._m_limits_low[0]

        self.norm_Mmin = Mmin
        self.norm_Mmax = Mmax

        int_mxi = self.int_mxi(self.norm_Mmin, self.norm_Mmax)

        self.k = Mcl / int_mxi
        self.lamda = self.int_xi_cl(self._m_limits_low[0], self._mass_limits)

        return

    def norm_cl_wk04(self, Mcl, Mmax=None, Mmin=None):
        """
        Helper function
        """
        self.k = 1.0
        self.Mcl = Mcl

        if Mmax == None:
            Mmax = self._m_limits_high[-1]

        if Mmin == None:
            Mmin = self._m_limits_low[0]

        if Mmax > Mcl:
            Mmax = Mcl

        if Mmax > self._m_limits_high[-1]:
            Mmax = self._m_limits_high[-1]

        if Mmin < self._m_limits_low[0]:
            Mmin = self._m_limits_low[0]

        a = Mmin
        c = Mmax
        b = (c + a) / 2.0
        while (((c/b)-(a/b)) > 0.00001):
            mb = self.int_mxi(Mmin, b) / self.int_xi(b, Mmax)
            if mb < Mcl:
                a = b
            else:
                c = b
            b = (c + a) / 2.0

        Mmax = b
        self.norm_Mmin = Mmin
        self.norm_Mmax = Mmax

        self.k = Mcl / self.int_mxi(Mmin, Mmax)
        self.lamda = self.int_xi_cl(self._m_limits_low[0], self._mass_limits)

    def xi_cl(self, m):
        """
        Helper function
        """
        return theta_closed(self.norm_Mmax - m) * self.xi(m)

    def mxi_cl(self, m):
        """
        Helper function
        """
        return theta_closed(self.norm_Mmax - m) * self.m_xi(m)

    def int_xi_cl(self, left, right):
        """
        Helper function
        """
        t1 = self.prim_xi(right)
        t2 = theta_closed(right - self.norm_Mmax)
        t3 = self.int_xi(self.norm_Mmax, right)
        t4 = self.prim_xi(left)
        t5 = theta_closed(left - self.norm_Mmax)
        t6 = self.int_xi(self.norm_Mmax, left)

        return (t1 - t2*t3) - (t4 - t5*t6)

    def int_mxi_cl(self, left, right):
        t1 = self.prim_mxi(right)
        t2 = theta_closed(right - self.norm_Mmax)
        t3 = self.int_mxi(self.norm_Mmax, right)
        t4 = self.prim_mxi(left)
        t5 = theta_closed(left - self.norm_Mmax)
        t6 = self.int_mxi(self.norm_Mmax, left)

        return (t1 - t2*t3) - (t4 - t5*t6)

    def dice_star_cl(self, r):
        """
        Given a list of random numbers (r), return a list of masses
        selected from the IMF.
        """
        returnFloat = isinstance(r, (float, np.floating))
        r = np.atleast_1d(r)  # Make sure it is an array

        x = r * self.lamda[-1]
        y = np.zeros_like(r)
        z = np.ones_like(r)

        # Loop through the different parts of the power law.
        for i in range(self.nterms): #-----For i = 0 --> n, where n is the number of intervals
            aux = x - self.lamda[i] #---Should this be i - 1?

            # Only continue for those entries that are in later segments
            idx = aux >= 0

            # Maybe we are all done?
            if sum(idx) == 0:
                break

            x_tmp = x[idx]
            aux_tmp = aux[idx]

            t1 = aux_tmp / (self.coeffs[i] * self.k)
            t1 += prim_power(self._m_limits_low[i], self._powers[i])
            y_i = gamma_closed(x_tmp, self.lamda[i], self.lamda[i+1])
            y_i *= inv_prim_power(t1, self._powers[i])

            # Save results into the y array
            y[idx] += y_i

            z *= delta(x - self.lamda[i+1]) # new version
            #z *= delta(x - self.lamda[i])   # old version

        if returnFloat:
            return y[0] * z[0]
        else:
            return y * z

class IMFSalpeter1955(IMF_broken_powerlaw):
    """
    Define IMF from `Salpeter (1955) <https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/abstract>`_.
    Mass range is 0.4 M_sun - 10 M_sun.
    """
    def __init__(self, multiplicity=None):

        massLimits = np.array([0.40, 10.0])
        powers = np.array([-2.3])

        IMF_broken_powerlaw.__init__(self, massLimits, powers,
                                     multiplicity=multiplicity)


class Miller_Scalo_1979(IMF_broken_powerlaw):
    """
    Define IMF from `Miller & Scalo (1979) <https://ui.adsabs.harvard.edu/abs/1979ApJS...41..513M/abstract>`_.
    Mass range is 0.1 M_sun - inf M_sun.
    """
    def __init__(self, multiplicity=None):
        massLimits = np.array([0.1, 1, 10, np.inf])
        powers = np.array([-1.4, -2.5, -3.3])

        IMF_broken_powerlaw.__init__(self, massLimits, powers,
                                     multiplicity=multiplicity)

class Kennicutt_1983(IMF_broken_powerlaw):
    """
    Define IMF from `Kennicutt (1983) <https://ui.adsabs.harvard.edu/abs/1983ApJ...272...54K/abstract>`_.
    Mass range is 0.1 M_sun - inf M_sun.
    """
    def __init__(self, multiplicity=None):
        massLimits = np.array([0.1, 1, np.inf])
        powers = np.array([-1.4, -2.5])

        IMF_broken_powerlaw.__init__(self, massLimits, powers,
                                     multiplicity=multiplicity)

class Kroupa_2001(IMF_broken_powerlaw):
    """
    Define IMF from `Kroupa (2001) <https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/abstract>`_.
    Mass range is 0.01 M_sun - inf M_sun.
    """
    def __init__(self, multiplicity=None):
        massLimits = np.array([0.01, 0.08, 0.5, 1, np.inf])
        powers = np.array([-0.3, -1.3, -2.3, -2.3])

        IMF_broken_powerlaw.__init__(self, massLimits, powers,
                                     multiplicity=multiplicity)

class Weidner_Kroupa_2004(IMF_broken_powerlaw):
    """
    Define IMF from `Weidner & Kroupa (2004) <https://ui.adsabs.harvard.edu/abs/2004MNRAS.348..187W/abstract>`_.
    Mass range is 0.01 M_sun - inf M_sun.
    """
    def __init__(self, multiplicity=None):
        massLimits = np.array([0.01, 0.08, 0.5, 1, 120])
        powers = np.array([-0.3, -1.3, -2.3, -2.35])

        IMF_broken_powerlaw.__init__(self, massLimits, powers,
                                     multiplicity=multiplicity)

class Salpeter_Kirkpatrick_2024(IMF_broken_powerlaw):
    """
    Define combined IMF from Kirkpatrick (2024) and Salpeter (1955) to allow
    inclusion of the brown dwarf mass range.
    Mass range:
        * 0.01 M_sun - 8 M_sun: Kirkpatrick 2024
        <https://ui.adsabs.harvard.edu/abs/2024ApJS..271...55K/abstract>`_.
        * 8 M_sun - 120 M_sun: Salpeter 1955
        <https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/abstract>`_.
    """
    def __init__(self, multiplicity=None):
        massLimits = np.array([0.01, 0.05, 0.22, 0.55, 8, 120])
        powers = np.array([-0.6, -0.25, -1.3, -2.3, -2.35])

        IMF_broken_powerlaw.__init__(self, massLimits, powers,
                                     multiplicity=multiplicity)

##################################################
#
# Generic functions -- see if we can move these up.
#
##################################################
def prim_power(m, power):
    """
    Takes floats or arrays, but returns arrays.
    returns: m**(power + 1) / (power + 1) and handles the case when power = -1
    """
    returnFloat = isinstance(m, (float, np.floating)) and isinstance(power, (float, np.floating))

    m = np.atleast_1d(m)
    power = np.atleast_1d(power)

    if (len(m) == 1) and (len(power) > 1):
        m = np.repeat(m, len(power))
    if (len(power) == 1) and (len(m) > 1):
        power = np.repeat(power, len(m))

    z = 1.0 + power
    val = np.empty_like(m)
    valid_idx = power != -1
    val[valid_idx] = (m[valid_idx]**z[valid_idx]) / z[valid_idx]
    val[~valid_idx] = np.log(m[~valid_idx])

    if returnFloat:
        return val[0]
    else:
        return val

def inv_prim_power(x, power):
    """
    returns ((1+power) * x)**(1.0 / (1 + power)) and handles the case
    when power == -1.
    """
    returnFloat = isinstance(x, (float, np.floating)) and isinstance(power, (float, np.floating))

    x = np.atleast_1d(x)
    power = np.atleast_1d(power)

    if (len(x) == 1) and (len(power) > 1):
        x = np.repeat(x, len(power))
    if (len(power) == 1) and (len(x) > 1):
        power = np.repeat(power, len(x))

    if x.shape != power.shape:
        raise ValueError('spisea.imf.inv_prim_power: Dimension mismatch, x and power must have the same shape')

    z = 1.0 + power
    val = np.empty_like(x)
    valid_idx = power != -1
    val[valid_idx] = (z[valid_idx] * x[valid_idx])**(1.0 / z[valid_idx])
    val[~valid_idx] = np.exp(x[~valid_idx])

    if returnFloat:
        return val[0]
    else:
        return val


def log_normal(m, mean_logm, sigma_logm):
    returnFloat = isinstance(m, (float, np.floating)) and isinstance(mean_logm, (float, np.floating)) and isinstance(sigma_logm, (float, np.floating))

    m = np.atleast_1d(m)
    mean_logm = np.atleat_1d(mean_logm)
    sigma_logm = np.atleat_1d(sigma_logm)

    z = np.log10(m) - mean_logm
    val = np.exp(-z**2 / (2.0 * sigma_logm**2)) / m

    if returnFloat:
        return val[0]
    else:
        return val

def prim_log_normal(m, mean_logm, sigma_logm):
    returnFloat = isinstance(m, (float, np.floating)) and isinstance(mean_logm, (float, np.floating)) and isinstance(sigma_logm, (float, np.floating))

    m = np.atleast_1d(m)
    mean_logm = np.atleat_1d(mean_logm)
    sigma_logm = np.atleat_1d(sigma_logm)

    mu = (np.log10(m) - mean_logm) / (1.4142135623731 * sigma_logm)
    val = 2.88586244942136 * sigma_logm * error(mu)

    if returnFloat:
        return val[0]
    else:
        return val

def inv_prim_log_normal(m, mean_logm, sigma_logm):
    returnFloat = isinstance(m, (float, np.floating)) and isinstance(mean_logm, (float, np.floating)) and isinstance(sigma_logm, (float, np.floating))

    m = np.atleast_1d(m)
    mean_logm = np.atleat_1d(mean_logm)
    sigma_logm = np.atleat_1d(sigma_logm)

    mu = inv_error(0.346516861952484 * m / sigma_logm)
    val = 10.0**(1.4142135623731 * sigma_logm * mu + mean_logm)

    if returnFloat:
        return val[0]
    else:
        return val

def mlog_normal(m, mean_logm, sigma_logm):
    returnFloat = isinstance(m, (float, np.floating)) and isinstance(mean_logm, (float, np.floating)) and isinstance(sigma_logm, (float, np.floating))

    m = np.atleast_1d(m)
    mean_logm = np.atleat_1d(mean_logm)
    sigma_logm = np.atleat_1d(sigma_logm)

    z = np.log10(m) - mean_logm
    val = np.exp(-z**2 / (2.0 * sigma_logm**2))

    if returnFloat:
        return val[0]
    else:
        return val

def prim_mlog_normal(m, mean_logm, sigma_logm):
    returnFloat = isinstance(m, (float, np.floating)) and isinstance(mean_logm, (float, np.floating)) and isinstance(sigma_logm, (float, np.floating))

    m = np.atleast_1d(m)
    mean_logm = np.atleat_1d(mean_logm)
    sigma_logm = np.atleat_1d(sigma_logm)

    eta = np.log10(m) - mean_logm - (sigma_logm**2 * 2.30258509299405)
    eta /= 1.4142135623731 * sigma_logm

    t1 = (1.15129254649702 * sigma_logm**2) + mean_logm

    val = error(eta)
    val *= 2.88586244942136 * sigma_logm * np.exp(2.30258509299405 * t1)

    if returnFloat:
        return val[0]
    else:
        return val


def theta_closed(x):
    """
    Pass in float or array of floats (x) and return 0 for x<0
    and 1 for everything else.
    """
    isFloat = type(x) == float

    x = np.atleast_1d(x)
    val = (x >= 0).astype('float')

    if isFloat:
        return val[0]
    else:
        return val

def theta_open(x):
    """
    Pass in float or array of floats (x) and return 1 for x>0
    and 0 for everything else.
    """
    isFloat = type(x) == float

    x = np.atleast_1d(x)
    val = (x > 0).astype('float')

    if isFloat:
        return val[0]
    else:
        return val

def delta(x):
    """
    Pass in float or array of floats (x) and return 0.5 for x==0
    and 1.0 for everything else.
    """
    isFloat = type(x) == float

    x = np.atleast_1d(x)
    val = np.ones(len(x), dtype=float)
    val[x == 0] = 0.5

    if isFloat:
        return val[0]
    else:
        return val

def gamma_closed(m, left, right):
    """

    """
    return theta_closed(m - left) * theta_closed(right - m)


def error(x):
    x2 = x**2
    ax2 = 0.140012288686666 * x2

    val = np.sqrt(1.0 - np.exp(-x2*(1.27323954473516+ax2)/(1+ax2)))

    if x >=0:
        return val
    else:
        return -val

def inv_error(x):
    x2 = x**2
    lnx2 = np.log(1.0 - x2)
    aux = 4.54688497944829 + (lnx2 / 2.0)
    y = -aux + np.sqrt(aux**2 - (lnx2 / 0.140012288686666))

    val = np.sqrt(y)

    if x>=0:
        return y
    else:
        return -y
