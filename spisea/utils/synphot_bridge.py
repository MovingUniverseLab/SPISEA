"""
Helpers for synphot / stsynphot (replacing pysynphot) in SPISEA.
"""
from __future__ import annotations

import numpy as np
from astropy import units as u

from astropy.modeling import CompoundModel

from synphot import Observation
from synphot.models import Empirical1D
from synphot.spectrum import SourceSpectrum, SpectralElement


def trim_spectrum(sp, wmin, wmax):
    """
    Trim *sp* to *[wmin, wmax]* wavelengths in Angstroms using an 
    empirical mask on ``waveset``.

    Parameters:
    ----------
    sp : SourceSpectrum
        The spectrum to trim.
    wmin : float
        The minimum wavelength to trim to.
    wmax : float
        The maximum wavelength to trim to.
        
    Returns:
    -------
    SourceSpectrum
        A new SourceSpectrum object with the trimmed spectrum.
    """
    w = sp.waveset
    y = sp(w)
    w_aa = sp.waveset.to(u.AA).value
    mask = (w_aa >= wmin) & (w_aa <= wmax)
    wsub = w[mask]
    ysub = y[mask]
    
    return SourceSpectrum(
        Empirical1D,
        points=wsub,
        lookup_table=ysub,
        meta=getattr(sp, "meta", None) or {}
    )


def tabulate_if_needed(sp):
    """pysynphot CompositeSourceSpectrum.tabulate → empirical SourceSpectrum."""
    if hasattr(sp, "model") and isinstance(sp.model, CompoundModel):
        w = sp.waveset
        y = sp(w)
        return SourceSpectrum(
            Empirical1D,
            points=w,
            lookup_table=y,
            meta=getattr(sp, "meta", None) or {},
        )
    return sp


def resample_source_to(sp, wave_target_aa):
    """Resample spectrum onto *wave_target_aa* (Angstrom, 1-D array)."""
    sp = tabulate_if_needed(sp)
    w = np.asarray(wave_target_aa, dtype=np.float64)
    fp = sp(sp.waveset)
    f = np.interp(
        w,
        sp.waveset.to(u.AA).value,
        np.asarray(fp.value, dtype=np.float64).ravel(),
        left=np.nan,
        right=np.nan,
    )
    yu = sp(sp.waveset)
    unit = yu.unit
    return SourceSpectrum(
        Empirical1D,
        points=w * u.AA,
        lookup_table=f * unit,
        meta=getattr(sp, "meta", None) or {},
    )


def bandpass_from_stsyn(bp_raw, name=None):
    """Return stsynphot ``band()`` result; optional ``meta['expr']`` label."""
    if name is not None:
        meta = getattr(bp_raw, "meta", None)
        if meta is None:
            bp_raw.meta = {}
            meta = bp_raw.meta
        meta["expr"] = name
    return bp_raw


def resample_bandpass(bp, new_wave_aa):
    """Resample a `SpectralElement` onto a new wavelength grid (Angstrom)."""
    return _resample_spectral_element(bp, new_wave_aa)


def bandpass_wave_aa(bp):
    """Wavelength sampling for *bp* (Angstrom, 1-D `~numpy.ndarray`)."""
    wq = bp.waveset
    return np.asarray(wq.to(u.AA).value, dtype=np.float64)


def bandpass_throughput_array(bp):
    """Throughput samples for *bp* (dimensionless, matching ``bandpass_wave_aa``)."""
    wq = bp.waveset
    y = bp(wq)
    return np.asarray(y.value, dtype=np.float64).ravel()


def _resample_spectral_element(bp, new_wave_aa):
    wnew = np.asarray(new_wave_aa, dtype=np.float64)
    wq = bp.waveset
    wold = np.asarray(wq.to(u.AA).value, dtype=np.float64)
    told = np.asarray(bp(wq).value, dtype=np.float64).ravel()
    tnew = np.interp(wnew, wold, told, left=0.0, right=0.0)
    from synphot import units as su

    return SpectralElement(
        Empirical1D,
        points=wnew * u.AA,
        lookup_table=tnew * su.THROUGHPUT,
        meta=getattr(bp, "meta", None) or {},
    )



def make_observation(spec, band, binset_aa, force="taper"):
    """synphot Observation with binset in Angstrom."""
    binset = np.asarray(binset_aa, dtype=np.float64)
    return Observation(spec, band, binset=binset * u.AA, force=force)


def observation_bin_edges(obs):
    """Approximate bin widths like np.diff(binwave) with synphot ``binset``."""
    bs = np.asarray(obs.binset.to(u.AA).value, dtype=np.float64)
    diff = np.diff(bs)
    diff = np.append(diff, diff[-1])
    return diff


def rebin_spec(waveset_in, spec, waveset_out):
    """
    Resample *spec* onto *waveset_out* in a manner that conserves the flux.

    Parameters
    ----------
    waveset_in : np.array or Quantity
        Input wavelength grid (e.g. ``sp.waveset``), ``astropy.units.Quantity``.
        If not units are included, then it is assumed to be in Angstroms.
    spec : SourceSpectrum
        ``synphot.spectrum.SourceSpectrum`` defined on *waveset_in*.
    waveset_out : np.array or Quantity
        Target bin centers (e.g. ``sp_atlas.waveset``), ``Quantity``.
        If not units are included, then it is assumed to be in Angstroms.

    Returns
    -------
    np.array
        Rebinned flux values.
    """
    from synphot import units as su

    if not isinstance(waveset_in, u.Quantity):
        waveset_in = waveset_in * u.AA
    if not isinstance(waveset_out, u.Quantity):
        waveset_out = waveset_out * u.AA

    n = len(waveset_in)
    filt = SpectralElement(
        Empirical1D,
        points=waveset_in,
        lookup_table=np.ones(n) * su.THROUGHPUT,
    )

    obs_f = Observation(spec, filt, binset=waveset_out, force="taper")

    return obs_f.binflux
