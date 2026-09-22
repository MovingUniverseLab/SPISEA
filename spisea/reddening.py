#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Reddening laws.
"""
import pylab as py
import numpy as np
from scipy import interpolate
from astropy import units as u
from astropy.table import Table
from scipy.linalg import solve_banded
import numbers

from synphot import exceptions
from synphot import units as su
from synphot.reddening import ExtinctionCurve, ReddeningLaw, ExtinctionModel1D

import ast
import os
import pdb


class RedLawBase(ExtinctionCurve):
    """
    Base class for reddening or extinction laws in SPISEA.

    Subclass of synphot `ExtinctionCurve` at :math:`A_{Ks} = 1` mag. For arbitrary
    :math:`A_{Ks}` and wavelength sampling, use `extinction_at` or `ExtinctionAtAKs`.

    Parameters
    ----------
    waveset : numpy array
        Wavelength grid in Angstrom.
    A_lambda_over_AKs : numpy array
        Extinction in magnitudes at each wavelength.
    name : str
        Name of the reddening law.
    litref : str
        Literature reference for the reddening law.

    Attributes
    ----------
    wave : numpy array
        Wavelength grid in Angstrom.
    A_lambda_over_AKs : numpy array
        Extinction in magnitudes at each wavelength.
    obscuration : numpy array
        Same as ``A_lambda_over_AKs`` (:math:`A_\\lambda/A_{Ks}` on the native grid).
    _init_name : str
        Name of the reddening law.
    _litref : str
        Literature reference for the reddening law.
    """

    def __init__(self, waveset, A_lambda_over_AKs, name="", litref=""):
        """
        Parameters
        ----------
        waveset : numpy array
            Wavelength grid in Angstrom.
        A_lambda_over_AKs : numpy array
            Extinction in magnitudes at each wavelength.
        name : str
            Name of the reddening law.
        litref : str
            Literature reference for the reddening law.

        Returns
        -------
        self : `RedLawBase`
            The reddening law object.
        """

        # Handle input units
        if isinstance(waveset, u.Quantity):
            waveset = waveset.to(u.AA)
        else:
            waveset = waveset * u.AA

        # Convert to throughput for application to flux.
        thru = 10 ** (-0.4 * A_lambda_over_AKs) * su.THROUGHPUT

        # Define some meta data for the extinction curve
        meta = {}
        if name:
            meta["name"] = name
        if litref:
            meta["litref"] = litref

        # Initialize the extinction curve
        super().__init__(ExtinctionModel1D,
                         points=waveset,
                         lookup_table=thru,
                         meta=meta)

        # Save this in the obscuration attribute as well.
        self.wave = waveset.value  # Stored in Angstroms
        self.A_lambda_over_AKs = A_lambda_over_AKs
        self.obscuration = self.A_lambda_over_AKs
        self.name = name
        self.litref = litref

        return
    

    def extinction_curve(self, AKs, wavelengths=None):
        """
        Return the extinction curve for this law at a given :math:`A_{Ks}`
        and at the specified wavelength grid.

        Parameters
        ----------
        AKs : float
            Total extinction in AKs, in mags
        wavelengths : numpy or Quantity array
            Wavelength array in Angstrom or Quantity. If not provided, the native 
            wavelength array is used, which is in Angstrom..

        Returns
        -------
        ext_at_wave : `ExtinctionCurve`
            Extinction curve at the specified wavelengths with :math:`A_{Ks}` = `AKs`.

        Raises
        ------
        synphot.exceptions.SynphotError
            If ``AKs`` is not a real number or a magnitude Quantity.
        """
        if isinstance(AKs, u.Quantity) and AKs.unit.decompose() == u.mag:
            AKs = AKs.value
        elif not isinstance(AKs, numbers.Real):
            raise exceptions.SynphotError('AKs={0} is invalid.'.format(AKs))

        if not isinstance(AKs, u.Quantity):
            AKs *= u.mag

        if wavelengths is None:
            wavelengths = self.wave

        x = self._validate_wavelengths(wavelengths)

        #if (self.low_lim is not None) and (x.min() < self.low_lim):
        #    raise ValueError(f'{self.name}: wavelength values below valid range of {self.low_lim}')
        #if (self.high_lim is not None) and (x.max() > self.high_lim):
        #    raise ValueError(f'{self.name}: wavelength values above valid range of {self.high_lim}')

        self.meta['A_Ks'] =  AKs

        # Get the A_lambda extinction at the original wavelengths.
        A_lambda_orig = self.A_lambda_over_AKs * AKs

        # Interpolate over the curve to the new wavelengths.
        spline_interp = interpolate.splrep(self.wave, A_lambda_orig.value, k=3, s=0)
        A_lambda_new = interpolate.splev(x, spline_interp)

        # Throughputs but only at the input wavelengths. Now rescaled to AKs.
        thru = 10 ** (-0.4 * A_lambda_new) * su.THROUGHPUT

        return ExtinctionCurve(ExtinctionModel1D, 
                               points=x, 
                               lookup_table=thru,
                               meta=self.meta)


def get_red_law(str):
    """
    Given a reddening law name, return the reddening
    law object.

    Parameters:
    ----------
    str: str
        Reddening law name and additional params (comma-separated).
        Name must match
    """
    # Parse the string, extracting redlaw name and other params
    tmp = str.split(',')
    name = tmp[0]

    # How we split this up changes for the broken power law EL
    # versus the other ELs (since we have arrays for broken power law EL)
    if name == 'broken_pl':
        # resplit the str, given that arrays are in it
        tmp = str.split('[')
        tmp2 = str.split(']')
        param1 = ast.literal_eval('[{0}'.format(tmp[1][:-1])) # Funny formatting to remove ,
        param2 = ast.literal_eval('{0}]'.format(tmp2[1][1:])) # Funny formatting to remove ,
        param3 = float(str.split(',')[-1])

        params = (param1, param2, param3)
    else:
        params = ()
        if len(tmp) > 1:
            for ii in range(len(tmp) - 1):
                params = params + (float(tmp[ii+1]),)

    # Define dictionary connecting redlaw names to the redlaw classes
    name_dict = {'pl': RedLawPowerLaw,
                     'broken_pl': RedLawBrokenPowerLaw,
                     'RL85': RedLawRiekeLebofsky,
                     'C89': RedLawCardelli,
                     'I05': RedLawIndebetouw05,
                     'RZ07': RedLawRomanZuniga07,
                     'F09': RedLawFitzpatrick09,
                     'N09':RedLawNishiyama09,
                     'S10': RedLawSchoedel10,
                     'F11': RedLawFritz11,
                     'D16': RedLawDamineli16,
                     'DM16': RedLawDeMarchi16,
                     'S16': RedLawSchlafly16,
                     'H18b': RedLawHosek18b,
                     'NL18': RedLawNoguerasLara18,
                     'NL20': RedLawNoguerasLara20,
                     'SODC': RedLawSODC
                    }

    # Make reddening law object, including params if necessary.
    # This is not great coding, but I really strugged to generalize this...
    if len(params) == 0:
        red_law = name_dict[name]()
    elif len(params) == 1:
        red_law = name_dict[name](params[0])
    elif len(params) == 2:
        red_law = name_dict[name](params[0], params[1])
    elif len(params) == 3:
        red_law = name_dict[name](params[0], params[1], params[2])
    elif len(params) == 4:
        red_law = name_dict[name](params[0], params[1], params[2], params[3])
    else:
        mes = 'Redlaw contains more params than reddening.get_red_law currently supports'
        raise ValueError(mes)

    return red_law


class RedLawNishiyama09(RedLawBase):
    """
    The extinction law towards the Galactic Center
    from `Nishiyama et al. 2009
    <https://ui.adsabs.harvard.edu/abs/2009ApJ...696.1407N/abstract>`_,
    combined with the Av / AKs value from `Nishiyama et al. 2008
    <https://ui.adsabs.harvard.edu/abs/2008ApJ...680.1174N/abstract>`_.
    This law is defined between 0.551 - 8.0 microns.

    This law is constructed in 3 segments:

    * 0.551 -- 1.24 microns: a linear interpolation in log(1/lambda) vs log(A/AKs) space
      between the Av/AKs and AJ/AKs values
    * 1.25 -- 2.14 microns: a power law with index of 2.0
    * 2.14 -- 8.0 microns: a spline interpolation between the observed extinction values

    A_lambda / A_Ks = 1 when lambda = 2.14 microns

    This is the default extinction law.
    """
    def __init__(self):
        # Fetch the extinction curve, pre-interpolate across 3-8 microns
        wave = np.linspace(0.5, 8.0, 7501) * u.micron

        # This will eventually be scaled by AKs when you
        # call reddening(). Right now, calc for AKs=1
        wave_vals, Alambda_scaled = RedLawNishiyama09._derive_nishiyama09(wave)

        RedLawBase.__init__(
            self,
            wave_vals,
            Alambda_scaled,
            name='Nishiyama09',
            litref='Nishiyama+ 2009',
        )

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave_vals)
        self.high_lim = max(wave_vals)

        # Other info
        self.scale_lambda = 2.14
        self.name = 'N09'


    @staticmethod
    def _derive_nishiyama09(wavelength):
        """
        Calculate the N09 extinction law as defined in the paper:
        a A_lambda/AKs = power law of exponent -2.0 between JHK. Then
        use a *linear* interpolation in 1/lambda space to go from J to the V-band observation,
        in order to avoid imposing more structure. A cublic spline interpolation
        across the wavelength points is used longward of K-band

        Parameters
        ----------
        wavelength : float or array of Quantity
            in microns
        AKs : float
            in magnitudes
        """
        # Handle input units
        if isinstance(wavelength, u.Quantity):
            wavelength = wavelength.to(u.micron)
        else:
            wavelength = wavelength * u.micron

        # Check if wavelength is already a numpy array, if not, convert to numpy array
        if not isinstance(wavelength, np.ndarray):
            wavelength = np.atleast_1d(wavelength)

        #-----Define power law extinction law between JHK----#
        jhk_idx = np.where( (wavelength >= 1.25 * u.micron) & (wavelength <= 2.14 * u.micron) )

        alpha = 2.0
        wave_jhk = wavelength[jhk_idx]
        idx_scale = np.argmin(np.abs(wave_jhk - 2.14*u.micron))

        A_jhk = (wave_jhk / (1 * u.micron))**(-1.0*alpha)
        A_Ks_jhk = A_jhk / A_jhk[-1]

        #----Now do a linear interpolation (in log(1/lambda) vs log(A/AKs) space) between 1.25 microns and 0.551 microns---#
        jv_idx = np.where( (wavelength < 1.25 * u.micron) & (wavelength > 0.551 * u.micron) )
        Av = 16.13
        func = interpolate.interp1d(np.log10(np.array([1.0/1.25, 1.0/0.551])), np.log10(np.array([A_Ks_jhk[0], Av])),
                                        kind='linear')
        A_Ks_jv = func(np.log10(1.0 / wavelength[jv_idx].to('micron').value))

        # Convert back to linear space
        A_Ks_jv = 10**A_Ks_jv

        #---Do a spline interpolation for the rest of the (long-wavelength) law---#
        # We do this since no other function form is given
        long_idx = np.where(wavelength > 2.14 * u.micron)
        wave = np.array([0.551, 1.25, 1.63, 2.14, 3.545, 4.442, 5.675, 7.760]) * u.micron
        A_AKs = np.array([16.13, 3.02, 1.73, 1.00, 0.500, 0.390, 0.360, 0.430])
        interp_idx = np.where(wave > 2.1 * u.micron)

        spline_interp = interpolate.splrep(wave[interp_idx].value, A_AKs[interp_idx], k=3, s=0)
        A_AKs_long = interpolate.splev(wavelength[long_idx].value, spline_interp)

        # Stitch together sections for the final law
        wave_vals = np.concatenate((wavelength[jv_idx[0]], wavelength[jhk_idx[0]]))
        A_AKs_vjhk = np.concatenate((A_Ks_jv, A_Ks_jhk))

        # Now add the long-wavelength law
        wave_vals = np.concatenate((wave_vals, wavelength[long_idx[0]]))
        A_AKs_final = np.concatenate((A_AKs_vjhk, A_AKs_long))

        return wave_vals, A_AKs_final

    def Nishiyama09(self, wavelength, AKs):
        """
        Return the extinction at a given wavelength assuming the
        extinction law and an overall `AKs` value.

        Parameters
        ----------
        wavelength : float or array of astropy Quantity
            Wavelength to return extinction for, in microns
        AKs : float
            Total extinction in AKs, in mags
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron

        if not isinstance(AKs, u.Quantity):
            AKs *= u.mag

        # If input entry is a single float, turn it into an array
        wavelength = np.atleast_1d(wavelength)

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < self.low_lim) | (max(wavelength) > self.high_lim)):
            raise ValueError(f'{self.name}: wavelength values beyond interpolation range')

        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * u.AA
        law = self.A_lambda_over_AKs

        # Find the value of the law at the closest points to wavelength
        wave_idxs = np.argmin(np.abs(wave-wavelength[:,np.newaxis]),axis=1)

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = law[wave_idxs] * AKs

        return A_at_wave

    def plot_Nishiyama09(self):
        """
        Plot Nishiayama+09 law vs. the observed measurements
        from their Table 1.
        """
        # Extract the law
        wave = self.wave * u.AA
        law = self.obscuration

        # Measured extinction values from Nishiyama+09, Table 1
        wave_obs = np.array([0.551, 1.25, 1.63, 2.14, 3.545, 4.442, 5.675, 7.760]) * u.micron
        A_AKs = np.array([16.13, 3.02, 1.73, 1.00, 0.500, 0.390, 0.360, 0.430]) # unitless
        A_AKs_err = np.array([0.21, 0.04, 0.03, 0.0, 0.01, 0.01, 0.01, 0.01]) # unitless

        # Make plot
        py.figure(figsize=(10,10))
        py.plot(wave.to('micron'), law, 'r-', label='EL Function')
        py.errorbar(wave_obs, A_AKs, yerr=A_AKs_err, fmt='k.', ms=10,
                        label='Measured')
        py.xlabel('Wavelength (microns)')
        py.ylabel(r'Extinction (A$_{\lambda}$)')
        py.title('Nishiyama+09 EL')
        py.gca().set_xscale('log')
        py.gca().set_yscale('log')
        py.legend()
        py.savefig('nishiyama09_el.png')

        return

class RedLawCardelli(RedLawBase):
    r"""
    Defines the extinction law from
    `Cardelli et al. 1989 <https://ui.adsabs.harvard.edu/abs/1989ApJ...345..245C/abstract>`_.
    The law is defined from 0.3 - 3 microns, and in terms
    of :math:`A_{\lambda} / A_{Ks}`, where Ks is 2.174 microns.

    Parameters
    ----------
    Rv : float
        Ratio of absolute to selective extinction, :math:`A(V) / E(B-V)`.
        The standard value for the diffuse ISM is 3.1.
    """
    def __init__(self, Rv):
        # Fetch the extinction curve, pre-interpolate across 0.3-3 microns
        wave = np.linspace(0.3, 3.0, 2701) * u.micron

        # This will eventually be scaled by AKs when you
        # call reddening(). Produces A_lambda for AKs = 1, which will be
        # scaled later. Expects wavelength in microns
        Alambda_scaled = RedLawCardelli._derive_cardelli(wave, Rv)

        super().__init__(
            wave,
            Alambda_scaled,
            name='Cardelli89',
            litref='Cardelli+ 2009',
        )

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)

        # other info
        self.scale_lambda = 2.174 * u.micron
        self.name = 'C89,{0}'.format(Rv)

    @staticmethod
    def _derive_cardelli(wavelength, Rv):
        """
        Cardelli extinction law. This produces extinction values expected
        for AKs = 1.

        Parameters
        ----------
        wavelength : float or array of astropy Quantity
            Wavelength to return extinction for, in microns
        Rv : float
            Ratio of absolute to selective extinction, :math:`A(V) / E(B-V)`.
            The standard value for the diffuse ISM is 3.1.

        Returns
        -------
        output : float or array of astropy Quantity
            Extinction in magnitudes at each wavelength.
        """
        # Handle input units
        if isinstance(wavelength, u.Quantity):
            wavelength = wavelength.to(u.micron)
        else:
            wavelength = wavelength * u.micron

        wavelength = np.atleast_1d(wavelength)

        x = 1.0 * u.micron / wavelength

        # check for applicability
        if (np.min(x) < 0.3):
            print( 'wavelength is longer than applicable range for Cardelli law')
            return None

        if (np.max(x) > 8.0):
            print( 'wavelength is shorter than applicable range for Cardelli law')
            return None

        # Set up some arrays for coefficients that we will need
        a = np.zeros(len(x), dtype=float)
        b = np.zeros(len(x), dtype=float)

        y = x - 1.82

        # Calculate coefficients for long wavelengths (low wavenumber)
        # Wavenumger <= 1.1 (Eq. 2a, 2b)
        idx = np.where(x <= 1.1)[0]
        a[idx] =  0.574 * x[idx] ** 1.61
        b[idx] = -0.527 * x[idx] ** 1.61

        # Calculate coefficients for intermediate wavelengths
        # 1.1 < wavenumber <= 3.3 (Eq. 3a, 3b)
        idx = np.where((x > 1.1) & (x <= 3.3))[0]
        yy = y[idx]
        a[idx] = 1 + (0.17699 * yy) - (0.50447 * yy ** 2) - \
            (0.02427 * yy ** 3) + (0.72085 * yy ** 4) + \
            (0.01979 * yy ** 5) - (0.77530 * yy ** 6) + \
            (0.32999 * yy ** 7)
        b[idx] = (1.41338 * yy) + (2.28305 * yy ** 2) + \
            (1.07233 * yy ** 3) - (5.38434 * yy ** 4) - \
            (0.62251 * yy ** 5) + (5.30260 * yy ** 6) - \
            (2.09002 * yy ** 7)

        # Calculate the long wavelength
        # 3.3 < wavenumber < 5.9 (Eq. 4a, 4b)
        idx = np.where((x > 3.3) & (x < 5.9))[0]
        xx = x[idx]
        a[idx] = 1.752 - (0.316 * xx) - (0.104/((xx - 4.67) ** 2 + 0.341))
        b[idx] = -3.090 + (1.825 * xx) + (1.206/((xx - 4.62) ** 2 + 0.263))

        # Calculate the longest wavelength
        # 5.9 <= wavenumber (Eq. 4a, 4b)
        idx = np.where(x >= 5.9)[0]
        xx = x[idx]
        a[idx] = 1.752 - (0.316 * xx) - (0.104/((xx - 4.67) ** 2 + 0.341)) + \
            (-0.04473 * (xx - 5.9) ** 2) - (0.009779 * (xx - 5.9) ** 3)
        b[idx] = -3.090 + (1.825 * xx) + (1.206/((xx - 4.62) ** 2 + 0.263)) + \
            (0.2130 * (xx - 5.9) ** 2) + (0.1207 * (xx - 5.9) ** 3)

        # A(lam) / A(V), from Eq. 1
        extinction = a + b/Rv

        # Now, want to produce A_lambda / AKs, to match other laws
        k_ind = np.argmin(np.abs(x-0.46))
        Aks_Av = a[k_ind] + b[k_ind]/Rv # Aks / Av
        Av_Aks = 1.0 / Aks_Av # Av / Aks

        output = extinction * Av_Aks # (A(lamb) / Av) * (Av / Aks) = (A(lamb) / Aks)

        return output

    def Cardelli89(self, wavelength, AKs):
        """
        Return the extinction at a given wavelength assuming the
        extinction law and an overall `AKs` value.

        Parameters
        ----------
        wavelength : float or array of astropy Quantity
            Wavelength to return extinction for, in microns
        AKs : float
            Total extinction in AKs, in mags
        """
        # Handle input units
        if isinstance(wavelength, u.Quantity):
            wavelength = wavelength.to(u.micron)
        else:
            wavelength = wavelength * u.micron

        if not isinstance(AKs, u.Quantity):
            AKs *= u.mag

        # If input entry is a single float, turn it into an array
        wavelength = np.atleast_1d(wavelength)

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < (self.low_lim)) | (max(wavelength) > (self.high_lim))):
            raise ValueError(f'{self.name}: wavelength values beyond interpolation range')

        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * u.AA
        law = self.A_lambda_over_AKs

        # Find the value of the law at the closest points to wavelength
        wave_idxs = np.argmin(np.abs(wave-wavelength[:,np.newaxis]),axis=1)

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = law[wave_idxs] * AKs

        return A_at_wave

class RedLawSODC(RedLawBase):
    r"""
    Defines the SODC extinction law from SynthPop, described by
    `Klüter & Huston et al. (2025) <https://ui.adsabs.harvard.edu/abs/2025AJ....169..317K/abstract>`_.
    It is based on the `Cardelli et al. (1989) <https://ui.adsabs.harvard.edu/abs/1989ApJ...345..245C/abstract>`_
    formulation with an updated optical end from `O'Donnell et al (1994) 
    <https://ui.adsabs.harvard.edu/abs/1994ApJ...422..158O/abstract>`_ and infrared side adjusted to match
    `Surot et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020A%26A...644A.140S/abstract>`_.
    The law is defined from 0.25 - 3.5 microns, and in terms
    of :math:`A_{\lambda} / A_{Ks}`, where Ks is 2.174 microns.

    Parameters
    ----------
    Rv : float
        Ratio of absolute to selective extinction, :math:`A(V) / E(B-V)`.
        The standard value for the diffuse ISM is 3.1. Toward the Galactic
        bulge, 2.5 is more typical.
    """
    def __init__(self, Rv):
        # Fetch the extinction curve, pre-interpolate across 0.3-3 microns
        wave = np.linspace(0.25, 3.5, 3251) * u.micron

        # This will eventually be scaled by AKs when you
        # call reddening(). Produces A_lambda for AKs = 1, which will be
        # scaled later. Expects wavelength in microns
        Alambda_scaled = RedLawSODC._derive_sodc(wave, Rv)

        super().__init__(
            wave,
            Alambda_scaled,
            name='SODC',
            litref='Klüter & Huston + 2025',
        )

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)

        # other info
        self.scale_lambda = 0.549 * u.micron
        self.name = 'SODC,{0}'.format(Rv)

    @staticmethod
    def _derive_sodc(wavelength, Rv):
        """
        SODC extinction law. This produces extinction values expected
        for AKs = 1
        """
        # Handle input units
        if isinstance(wavelength, u.Quantity):
            wavelength = wavelength.to(u.micron)
        else:
            wavelength = wavelength * u.micron

        wavelength = np.atleast_1d(wavelength)

        x = 1.0 * u.micron / wavelength

        # check for applicability
        if (np.min(wavelength.value) < 0.25):
            print( 'wavelength is shorter than applicable range for SODC law')
            return None

        if (np.max(wavelength.value) > 3.5):
            print( 'wavelength is longer than applicable range for SODC law')
            return None

        # Set up some arrays for coefficients that we will need
        a = np.zeros(len(x), dtype=float)
        b = np.zeros(len(x), dtype=float)

        y = x - 1.82

        # Calculate coefficients for long wavelengths (low wavenumber)
        # Wavenumger <= 1.1 
        idx = np.where(x <= 1.1)[0]
        a[idx] =  0.53974  * x[idx] ** 2.255
        b[idx] = -0.495567 * x[idx] ** 2.255

        # Calculate coefficients for short wavelengths
        # 1.1 < wavenumber
        idx = np.where((x > 1.1))[0]
        yy = y[idx]
        a[idx] = 1 + (0.104 * yy) - (0.609 * yy ** 2) + \
            (0.701 * yy ** 3) + (1.137* yy ** 4) - \
            (1.718 * yy ** 5) - (0.827 * yy ** 6) + \
            (1.647 * yy ** 7) - (0.505 * yy ** 8)
        b[idx] = (1.952 * yy) + (2.908 * yy ** 2) - \
            (3.989 * yy ** 3) - (7.985 * yy ** 4) + \
            (11.102 * yy ** 5) + (5.491 * yy ** 6) - \
            (10.805 * yy ** 7) + (3.347 * yy ** 8)

        # A(lam) / A(V), from Eq. 1
        extinction = a + b/Rv

        # Now, want to produce A_lambda / AKs, to match other laws
        k_ind = np.argmin(np.abs(x-0.46))
        Aks_Av = a[k_ind] + b[k_ind]/Rv # Aks / Av
        Av_Aks = 1.0 / Aks_Av # Av / Aks

        output = extinction * Av_Aks # (A(lamb) / Av) * (Av / Aks) = (A(lamb) / Aks)

        return output

    def SODC(self, wavelength, AKs):
        """
        Return the extinction at a given wavelength assuming the
        extinction law and an overall `AKs` value.

        Parameters
        ----------
        wavelength : float or array of astropy Quantity
            Wavelength to return extinction for, in microns
        AKs : float
            Total extinction in AKs, in mags
        """
        # Handle input units
        if isinstance(wavelength, u.Quantity):
            wavelength = wavelength.to(u.micron)
        else:
            wavelength = wavelength * u.micron

        if not isinstance(AKs, u.Quantity):
            AKs *= u.mag

        # If input entry is a single float, turn it into an array
        wavelength = np.atleast_1d(wavelength)

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < (self.low_lim)) | (max(wavelength) > (self.high_lim))):
            raise ValueError(f'{self.name}: wavelength values beyond interpolation range')

        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * u.AA
        law = self.A_lambda_over_AKs

        # Find the value of the law at the closest points to wavelength
        wave_idxs = np.argmin(np.abs(wave-wavelength[:,np.newaxis]),axis=1)

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = law[wave_idxs] * AKs

        return A_at_wave

class RedLawRomanZuniga07(RedLawBase):
    """
    Defines extinction law from `Roman-Zuniga et al. 2007
    <https://ui.adsabs.harvard.edu/abs/2007ApJ...664..357R/abstract>`_
    for the dense cloud core Barnard 59. The law is a cubic spline fit
    to the values of A_lambda / A_Ks derived using the color-color
    diagrams slopes in their Table 1. It is defined between 1.0 - 8.0
    microns.

    A_lambda / A_Ks = 1 when lambda = 2.164 microns
    """
    def __init__(self):
        # Fetch the extinction curve, pre-interpolate across 1-8 microns
        wave = np.linspace(1.0, 8.0, 7001) * u.micron
        wave = wave.to(u.AA)

        # This will eventually be scaled by AKs when you
        # call reddening(). Right now, calc for AKs=1
        A_lambda_over_AKs = RedLawRomanZuniga07._derive_romanzuniga07(wave)

        super().__init__(wave, A_lambda_over_AKs , name='RomanZuniga07', litref='Roman-Zuniga+ 2007')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)
        self.name = 'RZ07'

    @staticmethod
    def _derive_romanzuniga07(wavelength):
        """
        Measurements taken from C-C column in Table 1 of RZ07
        """
        # Convert to numpy array
        wavelength = np.atleast_1d(wavelength)

        # Handle input units
        if isinstance(wavelength, u.Quantity):
            wavelength = wavelength.to(u.micron)
        else:
            wavelength = wavelength * u.micron

        # Define the filters and their corresponding wavelengths and A/AKs values
        filters = ['J', 'H', 'Ks', '[3.6]', '[4.5]', '[5.8]', '[8.0]']
        wave =      np.array([1.240, 1.664, 2.164, 3.545, 4.442, 5.675, 7.760]) * u.micron
        A_AKs =     np.array([2.299, 1.550, 1.000, 0.618, 0.525, 0.462, 0.455]) # unitless
        A_AKs_err = np.array([0.530, 0.080, 0.000, 0.077, 0.063, 0.055, 0.059]) # unitless

        # Interpolate over the curve
        spline_interp = interpolate.splrep(wave.value, A_AKs, k=3, s=0)
        A_AKs_at_wave = interpolate.splev(wavelength.value, spline_interp)

        return A_AKs_at_wave

    def RomanZuniga07(self, wavelength, AKs):
        """
        Return the extinction at a given wavelength assuming the
        extinction law and an overall `AKs` value.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to return extinction for, in microns
        AKs : float
            Total extinction in AKs, in mags
        """
        # Handle input units
        if isinstance(wavelength, u.Quantity):
            wavelength = wavelength.to(u.micron)
        else:
            wavelength = wavelength * u.micron

        if not isinstance(AKs, u.Quantity):
            AKs *= u.mag

        # If input entry is a single float, turn it into an array
        wavelength = np.atleast_1d(wavelength)

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < self.low_lim) | (max(wavelength) > self.high_lim)):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))

        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * u.AA
        law = self.obscuration

        # Find the value of the law at the closest points to wavelength
        wave_idxs = np.argmin(np.abs(wave-wavelength[:,np.newaxis]),axis=1)

        # Now multiply by AKs (since law assumes AKs = 1)
        A_lambda = law[wave_idxs] * AKs

        return A_lambda

    def plot_RomanZuniga07(self):
        """
        Plot law against measured values from Roman-Zuniga+07
        """
        wave = self.wave * u.AA
        law = self.obscuration

        # Get the observed values
        wave_obs = np.array([1.240, 1.664, 2.164, 3.545, 4.442, 5.675, 7.760]) * u.micron
        A_AKs = np.array([2.299, 1.550, 1.000, 0.618, 0.525, 0.462, 0.455])
        A_AKs_err = np.array([0.530, 0.080, 0.000, 0.077, 0.063, 0.055, 0.059])

        # Make plot
        py.figure(figsize=(10,10))
        py.plot(wave.to('micron'), law, 'r-', label='EL Function')
        py.errorbar(wave_obs, A_AKs, yerr=A_AKs_err, fmt='k.', ms=10,
                        label='Measured')
        py.xlabel('Wavelength (microns)')
        py.ylabel(r'Extinction (A$_{\lambda}$)')
        py.title('Roman-Zuniga+07 EL')
        py.gca().set_xscale('log')
        py.gca().set_yscale('log')
        py.legend()
        py.savefig('romanzuniga07_el.png')
        return

class RedLawRiekeLebofsky(RedLawBase):
    """
    Defines the extinction law from `Rieke & Lebofsky 1985
    <https://ui.adsabs.harvard.edu/abs/1985ApJ...288..618R/abstract>`_
    for the Galactic Center. The law is defined between 1.0 - 5.0 microns.

    A_lambda / A_Ks = 1 when lambda = 2.2 microns
    """
    def __init__(self):
        # Define the wavelength range of the extinction law
        wave = np.linspace(1.0, 5.0, 4001) * u.micron

        # This will eventually be scaled by AKs when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawRiekeLebofsky._derive_RiekeLebofsky(wave)

        super().__init__(waveset=wave, 
                         A_lambda_over_AKs=Alambda_scaled, 
                         name='RiekeLebofsky', 
                         litref='Rieke+Lebovsky 1985')

        self.low_lim = min(wave)
        self.high_lim = max(wave)

        return

    @staticmethod
    def _derive_RiekeLebofsky(wavelength):
        """
        Calculate the resulting extinction for an array of wavelengths,
        using the values calculated in RL85. Extinction is originally
        in A_lambda / A_v, but is re-normalized to A_lambda / A_Ks

        Data pulled from Rieke+Lebofsky 1985, Table 3. Note that their
        Table contains values from 0.3 - 13 microns, but only 1 - 5 microns
        is measured directly.

        Wavelengths for filters I - M are from Rieke+89, table 4.
        """
        # Arrays with the values from the paper
        filters = ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'L', 'M',
                   '[8.0]', '[8.5]', '[9.0]', '[9.5]', '[10.0]', '[10.5]',
                   '[11.0]', '[11.5]', '[12.0]', '[12.5]', '[13.0]']
        wave = np.array([0.365, 0.445, 0.551, 0.658, 0.9, 1.25, 1.60, 2.2,
                         3.50, 4.8, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0,
                        11.5, 12.0, 12.5, 13.0]) * u.micron
        A_Av = np.array([1.531, 1.324, 1.00, 0.748, 0.482, 0.282, 0.175, 0.112,
                         0.058, 0.023, 0.02, 0.043, 0.074, 0.087, 0.083,
                         0.074, 0.060, 0.047, 0.037, 0.030, 0.027])

        # Want to change this from A/Av to A/AK
        k_ind = np.where(np.array(filters) == 'K')
        Ak_Av = A_Av[k_ind]
        Av_Ak = 1.0 / Ak_Av

        A_Ak = A_Av * Av_Ak

        # Only interpolate between I and M
        idx1 = np.where(np.array(filters) == 'I')[0][0]
        idx2 = np.where(np.array(filters) == 'M')[0][0]

        wave_interp = wave[idx1:idx2+1]
        A_Ak_interp = A_Ak[idx1:idx2+1]
        assert len(wave_interp) == 6

        # Interpolate over the curve over desired wavelength range
        spline_interp = interpolate.splrep(wave_interp.value, A_Ak_interp, k=3, s=0)
        A_Ak_at_wave = interpolate.splev(wavelength.value, spline_interp)

        return A_Ak_at_wave

    def RiekeLebofsky85(self, wavelength, AKs):
        """
        Return the extinction at a given wavelength assuming the
        extinction law and an overall `AKs` value.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to return extinction for, in microns
        AKs : float
            Total extinction in AKs, in mags
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron

        if not isinstance(AKs, u.Quantity):
            AKs *= u.mag

        # If input entry is a single float, turn it into an array
        wavelength = np.atleast_1d(wavelength)

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < self.low_lim) | (max(wavelength) > self.high_lim)):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))

        # Extract wave and A/AKs from law, turning wave into micron units.
        # Start with the input law wave (AA) and A/AKs.
        wave = self.wave * u.AA
        law = self.A_lambda_over_AKs

        # Find the value of the law at the closest points to wavelength
        wave_idxs = np.argmin(np.abs(wave-wavelength[:,np.newaxis]),axis=1)

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = law[wave_idxs] * AKs

        return A_at_wave

    def plot_RiekeLebofsky85(self):
        """
        Plot law against measured values from RL+85, from their
        Table 3
        """
        wave = self.wave * u.AA # Angstroms
        law = self.obscuration

        # Get the observed values from their Table 3. Note only JHKLM is
        # measured directly by RL85, other values come from elsewhere.
        # Wavelengths are from Rieke+89, Table 4
        filters = ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'L', 'M',
                   '[8.0]', '[8.5]', '[9.0]', '[9.5]', '[10.0]', '[10.5]',
                   '[11.0]', '[11.5]', '[12.0]', '[12.5]', '[13.0]']
        wave_obs = np.array([0.365, 0.445, 0.551, 0.658, 0.9, 1.25, 1.60, 2.2,
                         3.50, 4.8, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0,
                        11.5, 12.0, 12.5, 13.0]) * u.micron
        A_Av = np.array([1.531, 1.324, 1.00, 0.748, 0.482, 0.282, 0.175, 0.112,
                         0.058, 0.023, 0.02, 0.043, 0.074, 0.087, 0.083,
                         0.074, 0.060, 0.047, 0.037, 0.030, 0.027])

        # Want to change this from A/Av to A/AK
        k_ind = np.where(np.array(filters) == 'K')
        Ak_Av = A_Av[k_ind]
        Av_Ak = 1.0 / Ak_Av

        A_Ak = A_Av * Av_Ak

        # RL law is only defined between I and M. Restrict to this range
        idx1 = np.where(np.array(filters) == 'I')[0][0]
        idx2 = np.where(np.array(filters) == 'M')[0][0]

        wave_obs_f = wave_obs[idx1:idx2+1]
        A_Ak_f = A_Ak[idx1:idx2+1]

        # Make plot
        py.figure(figsize=(10,10))
        py.plot(wave.to('micron'), law, 'r-', label='EL Function')
        py.plot(wave_obs_f.to('micron'), A_Ak_f, 'k.', ms=10,
                        label='Measured')
        py.xlabel('Wavelength (microns)')
        py.ylabel(r'Extinction (A$_{\lambda}$)')
        py.title('Rieke+Lebofsky+85 EL')
        py.gca().set_xscale('log')
        py.gca().set_yscale('log')
        py.legend()
        py.savefig('rl85_el.png')

        return


class RedLawDamineli16(RedLawBase):
    """
    Defines the extinction law of `Damineli et al. 2016
    <https://ui.adsabs.harvard.edu/abs/2016MNRAS.463.2653D/abstract>`_,
    derived for the Wd1 cluster. The law is derived between
    0.4 - 4.8 microns.

    Law is scaled such that A_lambda / A_Ks = 1 at 2.156 microns
    """
    def __init__(self):
        # Fetch the extinction curve, pre-interpolate across 0.4-4.8 microns
        wave = np.linspace(0.4, 4.8, 4401) * u.micron

        # This will eventually be scaled by AKs when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawDamineli16._derive_Damineli16(wave)

        super().__init__(waveset=wave, 
                         A_lambda_over_AKs=Alambda_scaled, 
                         name='Damineli16', 
                         litref='Damineli+ 2016')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)

        self.scale_lambda = 2.159 * u.micron
        self.name = 'D16'

        return

    @staticmethod
    def _derive_Damineli16(wavelength):
        """
        Calculate the Damineli+16 extinction law using their equation 19

        Parameters
        ----------
        wavelength : float
            in microns
        AKs : float
            in magnitudes
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        # From their eq 19
        x = np.log10(2.159 * u.micron / wavelength)
        log_A_AKs = -0.015 + 2.33*x + 0.522*x**2. - 3.001*x**3. + 2.034*x**4.

        # Now to convert this back to linear space
        A_AKs_at_wave = 10**log_A_AKs

        return A_AKs_at_wave

    def Damineli16(self, wavelength, AKs):
        """
        Return the extinction at a given wavelength assuming the
        extinction law and an overall `AKs` value.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to return extinction for, in microns
        AKs : float
            Total extinction in AKs, in mags
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        # If input entry is a single float, turn it into an array
        wavelength = np.atleast_1d(wavelength)

        # Return error if any wavelength is beyond interpolation range of
        # extinction law2
        if ((min(wavelength) < self.low_lim) | (max(wavelength) > self.high_lim)):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))

        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * u.AA
        law = self.obscuration

        # Find the value of the law at the closest points to wavelength
        wave_idxs = np.argmin(np.abs(wave-wavelength[:,np.newaxis]),axis=1)

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = law[wave_idxs] * AKs

        return A_at_wave

    def plot_Damineli16(self):
        """
        Plot law against measured values from Damineli+16
        from their Table 1
        """
        wave = self.wave * u.AA # Angstroms
        law = self.obscuration

        # Get the observed values from their Table 1 (Wd1 + RC)
        wave_obs = np.array([0.442, 0.537, 0.664, 0.805, 0.878, 1.021,
                                 1.244, 1.651, 2.159, 3.295, 4.4809]) * u.micron
        A_AKs = np.array([21.43, 14.95, 11.25, 8.72, 7.23, 5.10, 3.23, 1.77,
                              1.0, 0.39, 0.26])

        # Make plot
        py.figure(figsize=(10,10))
        py.plot(wave.to('micron'), law, 'r-', label='EL Function')
        py.errorbar(wave_obs.to('micron'), A_AKs, fmt='k.', ms=10,
                        label='Measured')
        py.xlabel('Wavelength (microns)')
        py.ylabel(r'Extinction (A$_{\lambda}$)')
        py.title('Damineli+16 EL')
        py.gca().set_xscale('log')
        py.gca().set_yscale('log')
        py.legend()
        py.savefig('damineli16_el.png')

        return

class RedLawDeMarchi16(RedLawBase):
    """
    Defines extinction law from `De Marchi et al. 2016
    <https://ui.adsabs.harvard.edu/abs/2016MNRAS.455.4373D/abstract>`_
    derived for 30 Dorodus. The law is defined between 0.3 - 8.0 microns.
    """
    def __init__(self):
        # Fetch the extinction curve, pre-interpolate across 1-8 microns
        wave = np.linspace(0.3, 8.0, 7701) * u.micron

        # This will eventually be scaled by AK when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawDeMarchi16._derive_DeMarchi16(wave)

        super().__init__(waveset=wave, 
                         A_lambda_over_AKs=Alambda_scaled, 
                         name='DeMarchi16', 
                         litref='DeMarchi+ 2016')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)
        self.name = 'DM16'

    @staticmethod
    def _derive_DeMarchi16(wavelength):
        """
        Calculate the resulting extinction for an array of wavelengths.
        The extinction is normalized with A_Ks.

        Data pulled from DeMarchi+16, Table 3

        Note: Authors measure R_VI (V) = 3.09 +/- 0.15,
        so we use this to calculate the absolute extinctions in all
        of the other bands. This corresponds to A_I/ A_V = 0.676

        Note that they extrapolate their curve to get to K-band

        Parameters
        ----------
        wavelength : float
            in microns
        AKs : float
            in magnitudes
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        wavelength = np.atleast_1d(wavelength)

        AI_AV = 0.676

        # Extracting the values from the paper
        filters = ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']
        wave = np.array([0.365, 0.445, 0.551, 0.658, 0.806, 1.22, 1.63, 2.19]) * u.micron
        R_VI = np.array([4.41, 3.78, 3.09, 2.58, 2.09, 1.26, 0.84, 0.52])
        R_VI_err = np.array([0.18, 0.15, 0.15, 0.13, 0.17, 0.18, 0.12, 0.08])

        # We'll calculate A_AKs from R_VI
        A_Av = R_VI * (1. - AI_AV)
        AK_Av = A_Av[-1]
        A_AK = A_Av / AK_Av

        # Interpolate over the curve
        spline_interp = interpolate.splrep(wave.value, A_AK, k=3, s=0)
        A_AK_at_wave = interpolate.splev(wavelength.value, spline_interp)

        return A_AK_at_wave

    def DeMarchi16(self, wavelength, AK):
        """
        Return the extinction at a given wavelength assuming the
        extinction law and an overall `AKs` value.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to return extinction for, in microns
        AKs : float
            Total extinction in AKs, in mags
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        wavelength = np.atleast_1d(wavelength)

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < self.low_lim) | (max(wavelength) > self.high_lim)):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))

        # Extract wave and A/AKs from law, turning wave into micron units
        law = self.obscuration
        wave = self.wave * u.micron

        # Find the value of the law at the closest points to wavelength
        wave_idxs = np.argmin(np.abs(wave-wavelength[:,np.newaxis]),axis=1)

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = law[wave_idxs] * AKs

        return A_at_wave

class RedLawFitzpatrick09(RedLawBase):
    """
    Defines the extinction law from
    `Fitzpatrick et al. 2009 <https://ui.adsabs.harvard.edu/abs/2009ApJ...699.1209F/abstract>`_.
    The law is defined between 0.5 -- 3 microns.

    The extinction law is as defined in their equation 5, and has two
    free parameters: :math:`\alpha` and R(V). Averaged over 14 sight-lines,
    the authors generally find either :math:`alpha` ~ 2.5, R(V) ~ 3, or
    :math:`alpha` ~ 1.8, R(V) ~ 5 (their Figure 6).

    A_lambda / A_K = 1 at lambda = 2.18

    Parameters
    ----------
    alpha : float
         alpha parameter for extinction law.

    RV : float
        R(V) parameter for extinction law.
    """
    def __init__(self, alpha, RV):
        # Fetch the extinction curve, pre-interpolate across 1-8 microns
        wave = np.linspace(0.5, 3.0, 2501) * u.micron

        # This will eventually be scaled by AK when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawFitzpatrick09._derive_Fitzpatrick09(wave, alpha, RV)

        super().__init__(waveset=wave, 
                         A_lambda_over_AKs=Alambda_scaled, 
                         name='Fitzpatrick09', 
                         litref='Fitzpatrick+ 2009')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)

        self.scale_lambda = 2.18 * u.micron
        self.name = 'F09,{0},{1}'.format(alpha, RV)

    @staticmethod
    def _derive_Fitzpatrick09(wavelength, alpha, RV):
        """
        Calculate the resulting extinction for an array of wavelengths.
        The extinction is normalized with A_Ks.

        Data pulled from Fitzpactrick09, equation 5

        Parameters
        ----------
        wavelength : float
            in microns

        alpha: float
            Free parameter alpha

        RV: float
            Free parameter RV
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        wavelength = np.atleast_1d(wavelength)

        alpha = float(alpha)
        RV = float(RV)

        # First we'll calculate k(lambda - V) = E(lambda - V) / E(B - V),
        # directly from equation 5
        k = (0.349 + 2.087*RV) * (1.0 / (1.0 + (wavelength.value / 0.507)**alpha)) - RV

        # We'll calculate Alam/Av from K + Rv
        Alam_Av = (k / RV) + 1.

        # Finally, to get A_lambda/Aks we need to divide Alam_Av by AKs_Av.
        # We'll assume a wavelength of 2.18 for Ks, since it is the wavelength
        # they report for K-band
        idx = np.argmin(np.abs(wavelength.value - 2.18))

        A_AKs_at_wave = Alam_Av / Alam_Av[idx]

        return A_AKs_at_wave

    def Fitzpatrick09(self, wavelength, AKs):
        """
        Return the extinction at a given wavelength assuming the
        extinction law and an overall `AKs` value.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to return extinction for, in microns
        AKs : float
            Total extinction in AKs, in mags
        """
        # If input entry is a single float, turn it into an array
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        wavelength = np.atleast_1d(wavelength)

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < self.low_lim) | (max(wavelength) > self.high_lim)):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))

        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * u.AA
        law = self.A_lambda_over_AKs

        # Find the value of the law at the closest points to wavelength
        wave_idxs = np.argmin(np.abs(wave-wavelength[:,np.newaxis]),axis=1)

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = law[wave_idxs] * AKs

        return A_at_wave

class RedLawSchlafly16(RedLawBase):
    """
    Defines the extinction law from `Schlafly et al. 2016
    <https://ui.adsabs.harvard.edu/abs/2016ApJ...821...78S/abstract>`_.
    The law is defined between 0.5 - 4.8 microns.

    The law is scaled such that A_lambda / A_Ks = 1 when lambda = 2.151

    Parameters
    ----------
    AH_AKs : float
        Ratio of A_H / A_Ks, which sets the normalization of the law (see Schlafly+16)
    x : float
        Free parameter in extinction law (see Schlafly+16, Eqn 6)
    """
    def __init__(self, AH_AKs, x):
        # Fetch the extinction curve, pre-interpolate across 0.5-4.8 microns
        wave = np.linspace(0.5, 4.8, 4301) * u.micron

        # This will eventually be scaled by AK when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawSchlafly16._derive_Schlafly16(wave, AH_AKs, x)

        super().__init__(waveset=wave, 
                         A_lambda_over_AKs=Alambda_scaled, 
                         name='Schlafly16', 
                         litref='Schlafly+ 2016')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)
        self.name = 'S16,{0},{1}'.format(AH_AKs,x)

    @staticmethod
    def _derive_Schlafly16(wavelength, AH_AKs, x):
        """
        Calculate Schalfly+16 extinction law according to
        code provided in appendix of the paper. AH_AKs sets the
        gray component while x sets the shape of the law in an
        Rv-like way
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        wavelength = np.atleast_1d(wavelength)

        # Use the function from the Schlafly+16 appendix to get the extinciton law
        # for given AH_AKs and x value. This is given in terms of A_lambda / A(5420)
        law_func = RedLawSchlafly16._Schlafly_appendix(x, AH_AKs)

        # Evaluate function for desired wavelengths (in angstroms)
        law = law_func(wavelength.to(u.AA).value)

        idx = np.argmin(np.abs(wavelength.value - 2.151))
        law_out = law / law[idx]

        return law_out

    @staticmethod
    def _Schlafly_appendix(x, rhk):
        """
        Build the Schlafly+16 extinction curve as a callable spline.

        This is the appendix implementation from Schlafly et al. 2016, modified
        so that both the shape parameter ``x`` and the gray-component ratio
        ``rhk`` = A(H)/A(K) are caller-controlled. Anchor wavelengths and
        extinction vectors are the Schlafly+16 values (PS1 g through WISE W2).

        Parameters
        ----------
        x : float
            Shape parameter of the extinction curve (R(V)-like; Schlafly+16
            Eqn 6).
        rhk : float
            Near-IR gray-component ratio A(H)/A(K) used to pin the H and K
            anchors. Schlafly+16 used 1.55 (Indebetouw+2005).

        Returns
        -------
        cs : `CubicSpline`
            Callable extinction curve. Evaluate as ``cs(wave_aa)``, where
            ``wave_aa`` is wavelength in Angstrom, to get
            :math:`A(\\lambda)/A(5420\\,\\mathrm{\\AA})`.

        Notes
        -----
        The underlying photometry spans roughly 5000--45000 Angstrom. Values
        outside that range are extrapolated. The gray component is relatively
        uncertain, and its variation with ``x`` is largely made up.
        """
        # Schlafly+2016
        ra = np.array([ 0.65373283,  0.39063843,  0.20197893,  0.07871701, -0.00476316,
                   -0.14213929, -0.23660605, -0.28522577, -0.321301  , -0.33503192])
        dra = np.array([-0.54278669,  0.03404903,  0.36841725,  0.42265873,  0.38247769,
                     0.14148814, -0.04020524, -0.13457319, -0.26883343, -0.36269229])

        # "isoreddening wavelengths" for extinction curve, at E(g-r) = 0.65 reddening
        # T_eff = 4500, Fe/H = 0, log g = 2.5
        lam = np.array([  5032.36441067,   6280.53335141,   7571.85928312,   8690.89321059,
                      9635.52560909,  12377.04268274,  16381.78146718,  21510.20523237,
                     32949.54009328,  44809.4919175 ])


        anchors = ra + x*dra
        # fix gray component so that A(H)/A(K) = 1.55
        anchors += (-anchors[6] + rhk*anchors[7])/(1 - rhk)
        cs0 = CubicSpline(lam, anchors, yp='3d=0')
        # normalize at 5420 angstroms
        cs = CubicSpline(lam, anchors/cs0(5420.), yp='3d=0')
        return cs

    def Schlafly16(self, wavelength, AKs):
        """
        Return the extinction at a given wavelength assuming the
        extinction law and an overall `AKs` value.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to return extinction for, in microns
        AKs : float
            Total extinction in AKs, in mags
        """
        # If input entry is a single float, turn it into an array
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        wavelength = np.atleast_1d(wavelength)

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < self.low_lim) | (max(wavelength) > self.high_lim)):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))

        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * u.AA
        law = self.A_lambda_over_AKs

        # Find the value of the law at the closest points to wavelength
        wave_idxs = np.argmin(np.abs(wave-wavelength[:,np.newaxis]),axis=1)

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = law[wave_idxs] * AKs

        return A_at_wave

class RedLawIndebetouw05(RedLawBase):
    """
    Defines the extinction law from `Indebetouw et al. 2005
    <https://ui.adsabs.harvard.edu/abs/2005ApJ...619..931I/abstract>`_.
    The law is defined between 1.25 - 8 microns using Equation 4
    in their paper.

    Note that A_lambda / A_K = 1 when lambda = 2.164 microns.
    """
    def __init__(self):
        # Get A_lambda / A_K values from Indebetouw+05
        wave = np.linspace(1.25, 8.0, 6751) * u.micron

        Alambda_scaled = RedLawIndebetouw05._derive_Indebetouw05(wave)

        super().__init__(waveset=wave, 
                         A_lambda_over_AKs=Alambda_scaled, 
                         name='Indebetouw05', 
                         litref='Indebetouw+ 2005')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)

        # Other info
        self.scale_lambda = 2.164 * u.micron
        self.name = 'I05'

    @staticmethod
    def _derive_Indebetouw05(wave):
        """
        Calculate Indebetouw+05 extinction law using equation 4 of their
        paper.
        """
        log_Alambda_AK = 0.61 - 2.22*np.log10(wave.value) + 1.21 * np.log10(wave.value)**2.

        # Return output in terms of A_lambda / A_K
        Alambda_AK = 10**log_Alambda_AK

        return Alambda_AK

    def Indebetouw05(self, wavelength, AKs):
        """
        Return the extinction at a given wavelength assuming the
        extinction law and an overall extinction at AK (2.164 microns)

        Parameters
        ----------
        wavelength : float or array
            Wavelength to return extinction for, in microns
        AKs : float
            Total extinction in AKs, in mags
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        wavelength = np.atleast_1d(wavelength)

        if not isinstance(AKs, u.Quantity):
            AKs *= u.mag

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < self.low_lim) | (max(wavelength) > self.high_lim)):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))

        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * u.AA
        law = self.A_lambda_over_AKs

        # Find the value of the law at the closest points to wavelength
        wave_idxs = np.argmin(np.abs(wave-wavelength[:,np.newaxis]),axis=1)

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = law[wave_idxs] * AKs

        return A_at_wave

    def plot_Indebetouw05(self):
        """
        Plot Indebetouw+05 extinciton curve versus their
        actual measured values (their Table 1).
        This is similar to their Figure 6.

        Saves plot as indebetouw05_el.png in cwd
        """
        # Get the extinction law
        wave = self.wave * u.AA # Angstroms
        law = self.obscuration # A_lambda / AK


        # Their average measurements across sight lines
        # from Table 1
        wave_arr = [1.240, 1.664, 2.164, 3.545, 4.442, 5.675, 7.760] * u.micron
        law_obs_arr = [2.50, 1.55, 1.0, 0.56, 0.43, 0.43, 0.43]
        law_obs_err_arr = [0.15, 0.08, 0.0, 0.06, 0.08, 0.10, 0.10]

        # Make plot
        py.figure(figsize=(10,10))
        py.plot(wave.to('micron'), law, 'r-', label='EL Function')
        py.errorbar(wave_arr, law_obs_arr, yerr=law_obs_err_arr, fmt='k.', ms=10,
                        label='Measured')
        py.xlabel('Wavelength (microns)')
        py.ylabel(r'Extinction (A$_{\lambda}$)')
        py.title('Indebetouw+05 EL')
        py.gca().set_xscale('log')
        py.gca().set_yscale('log')
        py.legend()
        py.savefig('indebetouw05_el.png')

        return

class RedLawPowerLaw(RedLawBase):
    r"""
    Extinction object that is a power-law extinction law:
    :math:`A_{\lambda} \propto \lambda^{-\alpha}`.

    For example, to create an extinction law between
    0.8 and 3 microns where :math:`\alpha = 2.21`,
    where :math:`A_{\lambda} / A_{Ks} = 1` at 2.12 microns:

    >>> from spisea import reddening
    >>> red_law = reddening.RedLawPowerLaw(2.21, 2.12, wave_min=0.8, wave_max=3.0)

    Parameters
    ----------
    alpha : float
        Exponent of the extinction power-law.

    K_wave : float
        Extinction law is normalized such that AKs = 1 at `K_wave`.

    wave_min : float; optional
        Minimum wavelength of the extinction law, in microns.
        Default is 0.5 microns.

    wave_max : float; optional
        Maximum wavelength of the extinction law, in microns.
        Default is 5.0 microns
    """
    def __init__(self, alpha, K_wave, wave_min=0.5, wave_max=5.0):
        # Fetch the extinction curve, pre-interpolate across wave_min to wave_max
        wave = np.linspace(wave_min, wave_max, int(round((wave_max-wave_min)*1000+1))) * u.micron

        # This will eventually be scaled by AK when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawPowerLaw._derive_powerlaw(wave, alpha, K_wave * u.micron)

        super().__init__(waveset=wave, 
                         A_lambda_over_AKs=Alambda_scaled, 
                         name='Power law',
                         litref='Power law')
        

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)
        self.name = 'pl,{0},{1},{2},{3}'.format(alpha,K_wave,wave_min,wave_max)

    @staticmethod
    def _derive_powerlaw(wavelength, alpha, K_wave):
        """
        Calculate the resulting extinction for an array of wavelengths.
        The extinction is normalized with A_Ks.

        Parameters
        ----------
        wavelength : float
            in microns

        alpha: float
            -1.0 * (power law exponent)

        K_wave: float
            Desired K-band wavelength, in microns
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        if not isinstance(K_wave, u.Quantity):
            K_wave = K_wave * u.micron
        else:
            K_wave = K_wave.to(u.micron)

        # Create extinction law (Power Law relative to 1 micron)
        law = (wavelength / u.micron)**(-1.0 * alpha)

        # We'll identify K-band as 2.14 microns
        idx = np.argmin(np.abs(wavelength - K_wave))

        A_AKs_at_wave = law / law[idx]

        return A_AKs_at_wave

    def powerlaw(self, wavelength, AKs):
        """
        Return the extinction at a given wavelength assuming the
        extinction law and an overall `AKs` value.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to return extinction for, in microns
        AKs : float
            Total extinction in AKs, in mags
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        if not isinstance(AKs, u.Quantity):
            AKs *= u.mag

        # If input entry is a single float, turn it into an array
        wavelength = np.atleast_1d(wavelength)

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < self.low_lim) | (max(wavelength) > self.high_lim)):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))

        # Extract wave and A/AKs from law, turning wave into micron units
        law = self.A_lambda_over_AKs

        # Find the value of the law at the closest points to wavelength
        wave_idxs = np.argmin(np.abs(wave-wavelength[:,np.newaxis]),axis=1)

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = law[wave_idxs] * AKs

        return A_at_wave

class RedLawBrokenPowerLaw(RedLawBase):
    r"""
    Extinction object that is a broken power-law extinction law:
    :math:`A_{\lambda} \propto \lambda^{\alpha[n]}`

    for:
    :math: `\lambda_{limits}[n] < \lambda <= \lambda_{limits}[n+1]`

    Note: lambda_limits must be continuous in wavelength and K_wave must be
    within one of the section defined by the lambda_limits array.
    Extinction law is only defined over lambda_limits

    Units of lambda_limits array is microns.

    Parameters
    ----------
    lambda_limits : numpy array
        Array of length (N + 1) with lower and upper wavelength limits of
        the power-law segments. Units are microns.

    alpha_vals : numpy array
        Array of length N that contains the powers for each
        power-law segment.

    K_wave : float
        Extinction law is normalized such that AKs = 1 at `K_wave`.
    """
    def __init__(self, lambda_limits, alpha_vals, K_wave):
        # Handle input units
        lambda_limits = np.atleast_1d(lambda_limits)

        if not isinstance(lambda_limits, u.Quantity):
            lambda_limits = lambda_limits * u.micron
        else:
            lambda_limits = lambda_limits.to(u.micron)

        if not isinstance(K_wave, u.Quantity):
            K_wave = K_wave * u.micron
        else:
            K_wave = K_wave.to(u.micron)

        # Fetch the extinction curve, pre-interpolate across defined wavelength range
        wave_min,wave_max = np.min(lambda_limits)/u.micron, np.max(lambda_limits)/u.micron
        wave = np.linspace(wave_min, wave_max, int(round((wave_max-wave_min)*1000+1))) * u.micron

        # Deal with pesky floating point issues that can artificially push the upper
        # value of wave above the max(lambda_limit)
        wave[-1] = np.max(lambda_limits/u.micron) * u.micron

        # Assert that K_wave is within lambda_limits
        try:
            assert (K_wave >= np.min(lambda_limits)) & (K_wave <= np.max(lambda_limits))
        except:
            raise Exception('K_wave not within lambda_limits bounds')

        # This will eventually be scaled by AK when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawBrokenPowerLaw._derive_broken_powerlaw(wave, lambda_limits, alpha_vals, K_wave)

        super().__init__(waveset=wave, 
                         A_lambda_over_AKs=Alambda_scaled, 
                         name='Broken Power law',
                         litref='Broken Power law')

        self.low_lim = np.min(lambda_limits)
        self.high_lim = np.max(lambda_limits)
        self.name = 'broken_pl,{0},{1},{2}'.format(lambda_limits, alpha_vals, K_wave)
        
        return

    @staticmethod
    def _derive_broken_powerlaw(wave, lambda_limits, alpha_vals, K_wave):
        r"""
        Calculate the resulting extinction for an array of wavelengths.
        The extinction is normalized with A_Ks.

        Parameters
        ----------
        wavelength : float
            in microns

        alpha: float
            The power-law exponent where the law scales as A_lambda \propto lambda ** -alpha.

        K_wave: float
            Desired K-band wavelength, in microns
        """
        # Handle input units
        if not isinstance(wave, u.Quantity):
            wave = wave * u.micron
        else:
            wave = wave.to(u.micron)

        if not isinstance(K_wave, u.Quantity):
            K_wave = K_wave * u.micron
        else:
            K_wave = K_wave.to(u.micron)

        if not isinstance(lambda_limits, u.Quantity):
            lambda_limits = lambda_limits * u.micron
        else:
            lambda_limits = lambda_limits.to(u.micron)


        # Create extinction law in segments
        law = np.ones(len(wave)) * np.nan
        for ii in range(len(alpha_vals)):
            wave_min = lambda_limits[ii]
            wave_max = lambda_limits[ii+1]
            alpha = alpha_vals[ii]

            # Find elements of wavelength array in this segment
            idx = np.where( (wave >= wave_min) & (wave <= wave_max))

            # Calculate coefficient for this segment to ensure
            # law is continuous
            coeff = 1
            if ii > 0:
                for jj in range(ii):
                    wave_connect = lambda_limits[jj+1]
                    val = ((wave_connect / u.micron) ** (-1*alpha_vals[jj])) / ((wave_connect / u.micron) ** (-1*alpha_vals[jj+1]))

                    coeff *= val

            law[idx] = coeff * ((wave[idx] / u.micron)**(-1.0 * alpha))

        # Let's make sure we didn't miss updating any parts of the law
        assert np.sum(np.isnan(law)) == 0

        # We'll identify K-band as 2.14 microns
        idx = np.argmin(np.abs(wave - K_wave))
        A_AKs_at_wave = law / law[idx]

        return A_AKs_at_wave

    def broken_powerlaw(self, wavelength, AKs):
        """
        Return the extinction at a given wavelength assuming the
        extinction law and an overall `AKs` value.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to return extinction for, in microns
        AKs : float
            Total extinction in AKs, in mags
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        if not isinstance(AKs, u.Quantity):
            AKs *= u.mag

        # If input entry is a single float, turn it into an array
        wavelength = np.atleast_1d(wavelength)

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < self.low_lim) | (max(wavelength) > self.high_lim)):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))

        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * u.AA
        law = self.A_lambda_over_AKs

        # Find the value of the law at the closest points to wavelength
        wave_idxs = np.argmin(np.abs(wave-wavelength[:,np.newaxis]),axis=1)

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = law[wave_idxs] * AKs

        return A_at_wave

class RedLawFritz11(RedLawBase):
    """
    Defines extinction law from `Fritz et al. 2011
    <https://ui.adsabs.harvard.edu/abs/2011ApJ...737...73F/abstract>`_
    for the Galactic Center. The law is defined from 1.0 -- 26 microns.

    By default, law is scaled such that A_lambda / A_2.166 microns = 1.
    According to Fritz+11, A_2.166 microns = 2.62 +/- 0.11 mag is the total
    extinction towards the inner 14"x20" of the MW.

    Parameters:
    -----------
    scale_labmda: float
        Wavelength at which extinction law is scaled (A_lambda / A_scale_lambda = 1),
        in microns. Default is 2.166 microns

    """
    def __init__(self, scale_lambda=2.166):
        # Read in the interpolated extinction curve from Fritz+11,
        # based on their Table 8. Wavelengths in microns, extinction in mags
        wave, ext, ext_err = RedLawFritz11._read_Fritz11()

        # Convert wavelength to angstroms
        wave = u.Quantity(wave, unit=u.micron)
        wave = wave.to(u.AA)

        # Rescale extinction law such that A_lambda / A_2.166 microns = 1
        idx = np.argmin(np.abs(wave - (scale_lambda * u.micron)))
        ext_scale = ext / ext[idx]

        # Make custom reddening law
        super().__init__(waveset=wave, 
                         A_lambda_over_AKs=ext_scale.data, 
                         name='Fritz11', 
                         litref='Fritz+2011')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)

        # Other useful variables
        self.scale_lambda = scale_lambda
        self.name = 'F11'

        return

    @staticmethod
    def _read_Fritz11():
        """
        Return the interpolated extinction curve from Fritz+11,
        as defined in their Table 8.

        Output:
        ------
        wave: array
            Wavelength in microns

        ext: array
            Extinction in mags

        ext_err: array
            Extinction error, in mags
        """
        # Read in file with Table 8 info (published with Fritz+11 paper)
        inpath = os.path.dirname(os.path.abspath(__file__))
        infile = os.path.join(inpath, 'el_files', 'fritz11_EL_table8.fits')

        t = Table.read(infile, format='fits')
        wave = t['lambda']
        ext = t['A']
        ext_err = t['e_A']

        return wave, ext, ext_err

    @staticmethod
    def _read_Fritz11_obs():
        """
        Return the Fritz+11 observed values, from their Table 2

        Output:
        -------
        wave: array
            Wavelength in microns

        ext: array
            Extinction in mags

        ext_err: array
            Extinction error, in mags
        """
        # Read in file with Table 8 info (published with Fritz+11 paper)
        sep = '/'
        inpath = sep.join(__file__.split('/')[:-1])
        infile = '{0}/el_files/fritz11_EL_table2.txt'.format(inpath)

        t = Table.read(infile, format='ascii')
        wave = t['wave']
        ext = t['ext']
        ext_err = t['ext_err']

        return wave, ext, ext_err

    def plot_Fritz11(self):
        """
        Plot Fritz+11 interpolated extinciton curve (their Fig 8)
        versus their actual measured values (their Table 2).
        This is similar to their Figure 8.

        Saves plot as fritz11_el.png in cwd
        """
        # Read in the Fritz+11 table 8
        wave, ext, ext_err = RedLawFritz11._read_Fritz11()

        # Convert wavelength to angstroms
        wave = u.Quantity(wave, unit=u.micron)
        wave = wave.to(u.AA)

        # Read in Fritz+11 measurements (Table 2)
        wave_obs, ext_obs, ext_obs_err = RedLawFritz11._read_Fritz11_obs()

        # Convert wavelength to angstroms
        wave_obs = u.Quantity(wave_obs, unit=u.micron)
        wave_obs = wave_obs.to(u.AA)

        # Now plot the scaled extinction law, scaled to the Fritz+11
        # extinction at 2.166 microns. Remember that this produces
        # throughput = 10^-0.4*Alambda
        ext_scaled = self.extinction_curve(2.62, self.wave * u.AA)
        law = self.obscuration

        # Make plot
        py.figure(figsize=(10,10))
        py.plot(wave.to(u.micron).value, ext, 'g-', label='Interpolated EL')
        py.fill_between(wave.to(u.micron).value, ext+ext_err, ext-ext_err, color='blue',
                            alpha=0.3)
        py.errorbar(wave_obs.to(u.micron).value, ext_obs, yerr=ext_obs_err, fmt='k.', ms=10,
                        label='Measured')
        py.plot(self.wave/1e4, law, 'r-', label='Scaled EL')
        py.xlabel('Wavelength (microns)')
        py.ylabel(r'Extinction (A$_{\lambda}$)')
        py.title('Fritz+11 EL')
        py.gca().set_xscale('log')
        py.gca().set_yscale('log')
        py.legend()
        py.savefig('fritz11_el.png')

        return

    def Fritz11(self, wavelength, A_scale_lambda):
        """
        Return the extinction at a given wavelength assuming the
        extinction law and a total extinction at the scale_lambda
        (the wavelength where the extinction law = 1)

        Parameters
        ----------
        wavelength : float or array
            Wavelength to return extinction for, in microns
        A_scale_lambda : float
            Total extinction at scale_lambda, in mags
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        if not isinstance(A_scale_lambda, u.Quantity):
            A_scale_lambda = A_scale_lambda * u.mag
        else:
            A_scale_lambda = A_scale_lambda.to(u.mag)

        # If input entry is a single float, turn it into an array
        wavelength = np.atleast_1d(wavelength)

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < self.low_lim) | (max(wavelength) > self.high_lim)):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))

        # Extract wave and A/A_scale_lambda from law, turning wave into micron units
        wave = self.wave * u.AA
        law = self.obscuration

        # Find the value of the law at the closest points to wavelength
        wave_idxs = np.argmin(np.abs(wave-wavelength[:,np.newaxis]),axis=1)

        # Now multiply by A_scale_lambda (since law assumes A_scale_lambda = 1)
        A_at_wave = law[wave_idxs] * A_scale_lambda

        return A_at_wave


class RedLawHosek18b(RedLawBase):
    """
    Defines extinction law from `Hosek et al. 2019
    <https://ui.adsabs.harvard.edu/abs/2019ApJ...870...44H/abstract>`_
    for the Arches cluster and Wd1.
    The law is derived between 0.7 - 3.54 microns
    """
    def __init__(self):
        # Fetch the extinction curve, pre-interpolate across 3-8 microns
        wave = np.linspace(0.7, 3.545, 2846) * u.micron

        # This will eventually be scaled by AKs when you
        # call reddening(). Right now, calc for AKs=1
        Alambda_scaled = RedLawHosek18b._derive_Hosek18b(wave)

        super().__init__(waveset=wave, 
                         A_lambda_over_AKs=Alambda_scaled, 
                         name='Hosek+18b', 
                         litref='Hosek+ 2018b')

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = min(wave)
        self.high_lim = max(wave)
        self.name = 'H18b'

    @staticmethod
    def _derive_Hosek18b(wavelength):
        """
        Derive the Hosek+18 extinction law, using the data from Table 4.

        Calculate the resulting extinction for an array of wavelengths.
        The extinction is normalized with A_Ks.

        Data pulled from Hosek+18, Table 4

        Parameters
        ----------
        wavelength : float
            Wavelength range to define extinction law over, in microns
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        wavelength = np.atleast_1d(wavelength)

        # Extinction law definition
        wave = np.array([0.8059, 0.962, 1.25, 1.53, 2.14, 3.545]) * u.micron
        A_AKs = np.array([7.943, 5.715, 3.142, 2.04, 1.0, 0.50])

        # Following Hosek+18, Interpolate over the curve with cubic spline interpolation
        spline_interp = interpolate.splrep(wave.to(u.AA).value, A_AKs, k=3, s=0)
        A_AKs_at_wave = interpolate.splev(wavelength.to(u.AA).value, spline_interp)

        # This curve already assumes A_Ks = 1.0, so we can go straight to
        # output
        return A_AKs_at_wave

    def Hosek18b(self, wavelength, AKs):
        """
        Return the extinction at a given wavelength assuming the
        extinction law and an overall `AKs` value.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to return extinction for, in microns
        AKs : float
            Total extinction in AKs, in mags
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        if not isinstance(AKs, u.Quantity):
            AKs *= u.mag

        wavelength = np.atleast_1d(wavelength)

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < self.low_lim) | (max(wavelength) > self.high_lim)):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))

        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * u.AA
        law = self.obscuration

        # Find the value of the law at the closest points to wavelength
        wave_idxs = np.argmin(np.abs(wave-wavelength[:,np.newaxis]),axis=1)

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = law[wave_idxs] * AKs

        return A_at_wave

class RedLawSchoedel10(RedLawBrokenPowerLaw):
    """
    Defines extinction law from `Schoedel et al. 2010
    <https://ui.adsabs.harvard.edu/abs/2010A%26A...511A..18S/abstract>`_
    for the Galactic Center. It is defined between 1.5 - 3.8 microns.

    Power law indices:

    * 1.677 - 2.168 microns ---> alpha = 2.21 +/- 0.24
    * 2.168 - 3.636 microns ---> alpha = 1.34 +/- 0.29

    Wavelengths come from effective wavelengths of observations (some buffer
    is added to either side of these values).

    Reddening law is scaled such that A_lambda / A_Ks = 1 at
    lambda = 2.168 microns.
    """
    def __init__(self):
        lambda_limits = [1.5, 2.168, 3.8]
        alpha_vals = [2.21, 1.34]
        K_wave = 2.168
        super().__init__(lambda_limits, alpha_vals, K_wave)

        # Set the upper/lower wavelength limits of law (in microns)
        self.low_lim = np.min(lambda_limits)
        self.high_lim = np.max(lambda_limits)

        # Other useful variables
        self.scale_lambda = K_wave
        self.name = 'S10'

        return

    def Schoedel10(self, wavelength, AKs):
        """
        Return the extinction at a given wavelength assuming the
        extinction law and a total extinction at scale_lambda
        (the wavelength where the extinction law = 1)

        Parameters
        ----------
        wavelength : float or array
            Wavelength to return extinction for, in microns
        AKs : float
            Total extinction at scale_lambda, in mags
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        if not isinstance(AKs, u.Quantity):
            AKs *= u.mag

        wavelength = np.atleast_1d(wavelength)

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < self.low_lim) | (max(wavelength) > self.high_lim)):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))

        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * u.AA
        law = self.obscuration

        # Find the value of the law at the closest points to wavelength
        wave_idxs = np.argmin(np.abs(wave-wavelength[:,np.newaxis]),axis=1)

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = law[wave_idxs] * AKs

        return A_at_wave


class RedLawNoguerasLara18(RedLawPowerLaw):
    """
    Defines extinction law from `Nogueras-Lara et al. 2018
    <https://ui.adsabs.harvard.edu/abs/2018A%26A...610A..83N/abstract>`_
    for the Galactic Center.  It is defined between 1.0 - 3.0 microns.
    Measurements were made in JHK, with effective wavelengths
    of 1.2685, 1.6506, and 2.1629 microns, respectively.

    This extinction law is a single power law with exponent
    of alpha = 2.3 +/- 0.08.

    Reddening law is scaled such that A_lambda / A_Ks = 1 at
    lambda = 2.163 microns (the observed K-band)
    """
    def __init__(self):
        wave_min = 1.0
        wave_max = 3.0
        K_wave = 2.163
        RedLawPowerLaw.__init__(self, 2.30, K_wave, wave_min=wave_min, wave_max=wave_max)

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = wave_min*10**4
        self.high_lim = wave_max*10**4

        self.scale_lambda = 2.163
        self.name = 'NL18'

    def NoguerasLara18(self, wavelength, AKs):
        """
        Return the extinction at a given wavelength assuming the
        extinction law and an overall `AKs` value.

        Parameters
        ----------
        wavelength : float or array
            Wavelength to return extinction for, in microns
        AKs : float
            Total extinction in AKs, in mags
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        if not isinstance(AKs, u.Quantity):
            AKs *= u.mag

        wavelength = np.atleast_1d(wavelength)

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < self.low_lim*u.angstrom) | (max(wavelength) > self.high_lim*u.angstrom)):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))

        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * u.AA
        law = self.obscuration

        # Find the value of the law at the closest points to wavelength
        wave_idxs = np.argmin(np.abs(wave-wavelength[:,np.newaxis]),axis=1)

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = law[wave_idxs] * AKs

        return A_at_wave

class RedLawNoguerasLara20(RedLawBrokenPowerLaw):
    """
    Defines extinction law from `Nogueras-Lara et al. 2020
    <https://ui.adsabs.harvard.edu/abs/2020A%26A...641A.141N/abstract>`_
    for the Galactic Center. It is defined between 1.0 -- 3 microns.
    Measurements were made in JHK, with effective wavelengths
    of 1.2685, 1.6506, and 2.1629 microns, respectively

    Measured power law indices:

    * 1.2685 - 1.6505 microns ---> alpha = 2.44 +/- 0.05
    * 1.6505 - 2.1629 microns ---> alpha = 2.23 +/- 0.05

    Wavelengths come from effective wavelengths of observations (some buffer
    is added to either side of these values).

    Reddening law is scaled such that A_lambda / A_Ks = 1 at
    lambda = 2.163 microns (the observed K-band)
    """
    def __init__(self):
        lambda_limits = [1.0, 1.6505, 3.0]
        alpha_vals = [2.23, 2.44]
        K_wave = 2.163
        super().__init__(lambda_limits, alpha_vals, K_wave)

        # Set the upper/lower wavelength limits of law (in angstroms)
        self.low_lim = np.min(lambda_limits) * u.micron
        self.high_lim = np.max(lambda_limits) * u.micron

        # Other useful variables
        self.scale_lambda = K_wave
        self.name = 'NL20'

        return

    def NoguerasLara20(self, wavelength, AKs):
        """
        Return the extinction at a given wavelength assuming the
        extinction law and a total extinction at scale_lambda
        (the wavelength where the extinction law = 1)

        Parameters
        ----------
        wavelength : float or array
            Wavelength to return extinction for, in microns
        AKs : float
            Total extinction at scale_lambda, in mags
        """
        # Handle input units
        if not isinstance(wavelength, u.Quantity):
            wavelength = wavelength * u.micron
        else:
            wavelength = wavelength.to(u.micron)

        if not isinstance(AKs, u.Quantity):
            AKs *= u.mag

        wavelength = np.atleast_1d(wavelength)

        # Return error if any wavelength is beyond interpolation range of
        # extinction law
        if ((min(wavelength) < self.low_lim) | (max(wavelength) > self.high_lim)):
            return ValueError('{0}: wavelength values beyond interpolation range'.format(self))

        # Extract wave and A/AKs from law, turning wave into micron units
        wave = self.wave * u.AA
        law = self.obscuration

        # Find the value of the law at the closest points to wavelength
        wave_idxs = np.argmin(np.abs(wave-wavelength[:,np.newaxis]),axis=1)

        # Now multiply by AKs (since law assumes AKs = 1)
        A_at_wave = law[wave_idxs] * AKs

        return A_at_wave

#---------------------------#
# Cubic spline function from Schalfly+16 appendix
#---------------------------#
def splint(spl, x):
    npts = len(spl.x)
    lo = np.searchsorted(spl.x, x)-1
    lo = np.clip(lo, 0, npts-2)
    hi = lo + 1
    dx = spl.x[hi] - spl.x[lo]
    a = (spl.x[hi] - x)/dx
    b = (x-spl.x[lo])/dx
    y = (a*spl.y[lo]+b*spl.y[hi]+
         ((a**3-a)*spl.y2[lo]+(b**3-b)*spl.y2[hi])*dx**2./6.)
    return y

class CubicSpline:
    def __init__(self, x, y, yp=None):
        npts = len(x)
        mat = np.zeros((3, npts))
        # enforce continuity of 1st derivatives
        mat[1,1:-1] = (x[2:  ]-x[0:-2])/3.
        mat[2,0:-2] = (x[1:-1]-x[0:-2])/6.
        mat[0,2:  ] = (x[2:  ]-x[1:-1])/6.
        bb = np.zeros(npts)
        bb[1:-1] = ((y[2:  ]-y[1:-1])/(x[2:  ]-x[1:-1]) -
                    (y[1:-1]-y[0:-2])/(x[1:-1]-x[0:-2]))
        if yp is None: # natural cubic spline
            mat[1,0] = 1.
            mat[1,-1] = 1.
            bb[0] = 0.
            bb[-1] = 0.
        elif yp == '3d=0':
            mat[1, 0] = -1./(x[1]-x[0])
            mat[0, 1] =  1./(x[1]-x[0])
            mat[1,-1] =  1./(x[-2]-x[-1])
            mat[2,-2] = -1./(x[-2]-x[-1])
            bb[ 0] = 0.
            bb[-1] = 0.
        else:
            mat[1, 0] = -1./3.*(x[1]-x[0])
            mat[0, 1] = -1./6.*(x[1]-x[0])
            mat[2,-2] =  1./6.*(x[-1]-x[-2])
            mat[1,-1] =  1./3.*(x[-1]-x[-2])
            bb[ 0] = yp[0]-1.*(y[ 1]-y[ 0])/(x[ 1]-x[ 0])
            bb[-1] = yp[1]-1.*(y[-1]-y[-2])/(x[-1]-x[-2])
        y2 = solve_banded((1,1), mat, bb)
        self.x, self.y, self.y2 = (x, y, y2)
    def __call__(self, x):
        return splint(self, x)
