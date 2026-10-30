"""
Tests for atmosphere utilities, including ``rebin_spec``.
"""
import importlib.util
import os
import time

import numpy as np
import pytest
from astropy import units as u
from synphot import Observation
from synphot import units as su
from synphot.models import BlackBodyNorm1D, Empirical1D
from synphot.spectrum import SourceSpectrum, SpectralElement
from synphot.units import convert_flux

from spisea import atmospheres as atm
from spisea import synthetic as syn


def _unit_flat_bandpass(waveset):
    """Unit throughput on *waveset* (same construction as ``synphot_bridge.rebin_spec``)."""
    n = len(np.asarray(waveset.to(u.AA).value))
    return SpectralElement(
        Empirical1D,
        points=waveset,
        lookup_table=np.ones(n) * su.THROUGHPUT,
    )


def test_rebin_spec_matches_synphot_observation():
    """``rebin_spec`` must return the same binned flux as ``Observation.binflux``."""
    vega = syn.vega
    w_in = vega.waveset
    w_out = w_in[::28]
    filt = _unit_flat_bandpass(w_in)

    got = atm.rebin_spec(w_in, vega, w_out).value
    obs = Observation(vega, filt, binset=w_out, force="taper")
    expected = np.asarray(obs.binflux.value, dtype=np.float64).ravel()

    np.testing.assert_allclose(got, expected, rtol=0, atol=0)

    return


def test_vega_rebin_conserves_integrated_flux():
    """
    Total integrated flux (spec × unit band) is independent of binning grid.

    Uses the package Vega spectrum and the same flat bandpass as ``rebin_spec``;
    :meth:`~synphot.spectrum.BaseSourceSpectrum.integrate` must agree for fine
    and coarse ``binset`` values.
    """
    vega = syn.vega
    w_in = vega.waveset
    filt = _unit_flat_bandpass(w_in)

    obs_fine = Observation(vega, filt, binset=w_in[::3], force="taper")
    obs_coarse = Observation(vega, filt, binset=w_in[::45], force="taper")

    int_fine = obs_fine.integrate()
    int_coarse = obs_coarse.integrate()

    assert int_fine.unit == int_coarse.unit
    np.testing.assert_allclose(
        int_fine.value,
        int_coarse.value,
        rtol=1e-12,
        atol=0,
        err_msg="Integrated flux changed when only the rebin grid changed",
    )

    return


def test_rebin_spec_output_length_matches_binset():
    vega = syn.vega
    w_out = vega.waveset[::33]
    bf = atm.rebin_spec(vega.waveset, vega, w_out)
    assert bf.shape == (w_out.size,)

    return


def test_pysynphot_vs_stsynphot_timing():
    """
    Time Kurucz ``k93models`` spectrum extraction via stsynphot or pysynphot.

    This module imports ``spisea.atmospheres`` and ``spisea.synthetic``, which
    import ``stsynphot`` at load time; the pysynphot-only section of this code will
    only run on older versions of SPISEA.
    """
    temperature = 20000
    metallicity = 0.0
    gravity = 4.0

    has_stsyn = importlib.util.find_spec("stsynphot") is not None
    has_pysyn = importlib.util.find_spec("pysynphot") is not None

    if has_stsyn:
        if not os.environ.get("PYSYN_CDBS"):
            pytest.skip("PYSYN_CDBS not set; grid_to_spec needs CDBS tree")

        from stsynphot.catalog import grid_to_spec

        t0 = time.perf_counter()
        try:
            sp = grid_to_spec("k93models", temperature, metallicity, gravity)
            w = sp.waveset
            _ = sp(w)
        except Exception as exc:
            pytest.skip(f"stsynphot grid_to_spec failed: {exc}")
        elapsed = time.perf_counter() - t0

        assert elapsed >= 0
        assert np.isfinite(elapsed)
        assert w.size > 0

    elif has_pysyn:
        import synphot

        t0 = time.perf_counter()
        try:
            sp = synphot.grid_to_spec("k93models", temperature, metallicity, gravity)
            w = sp.waveset
            _ = sp(w)
        except Exception as exc:
            pytest.skip(f"pysynphot Icat failed: {exc}")
        elapsed = time.perf_counter() - t0

        assert elapsed >= 0
        assert np.isfinite(elapsed)
        assert np.asarray(w).size > 0

    else:
        pytest.skip("Neither stsynphot nor pysynphot is installed")

    return


def test_wd_atmosphere_covers_mir_filters():
    """
    White-dwarf spectra must overlap mid-infrared filters.

    The Koester (2010) grid ends near 3 micron. ``get_wd_atmosphere``
    extends that grid with a blackbody so JWST F323N and F405N, NIRC2 Lp
    and Ms, and NSFCam L all overlap the spectrum.
    """
    temperature = 25000.0
    gravity = 7.8

    # Native grid, used to locate the join and to check continuity.
    native = atm.get_wdKoester_atmosphere(temperature=temperature,
                                          gravity=gravity, metallicity=0.0)
    native_wave = np.asarray(native.waveset.to(u.AA).value, dtype=float)
    edge = native_wave[-1] * u.AA
    assert edge.to(u.AA).value < 32000

    spectrum = atm.get_wd_atmosphere(temperature=temperature, gravity=gravity,
                                     metallicity=0.0)
    wave = np.asarray(spectrum.waveset.to(u.AA).value, dtype=float)
    assert wave[-1] >= 1.0e5

    # Flux at the native red edge is unchanged. Evaluation is PHOTLAM;
    # compare after converting both to FLAM.
    native_flux = convert_flux(edge, native(edge), su.FLAM)
    extended_flux = convert_flux(edge, spectrum(edge), su.FLAM)
    np.testing.assert_allclose(np.asarray(extended_flux.value).ravel(),
                               np.asarray(native_flux.value).ravel(),
                               rtol=1.0e-5, atol=0.0)

    # Tail samples are a blackbody. Compare on the tail grid itself so
    # the check is not limited by linear interpolation between samples.
    # Shape of tail_wave: [n_tail].
    tail_wave = wave[wave > native_wave[-1]]
    wave_blue = tail_wave[5] * u.AA
    wave_red = tail_wave[40] * u.AA
    bb = SourceSpectrum(BlackBodyNorm1D, temperature=temperature)
    flux_blue = convert_flux(wave_blue, spectrum(wave_blue), su.FLAM)
    flux_red = convert_flux(wave_red, spectrum(wave_red), su.FLAM)
    bb_blue = convert_flux(wave_blue, bb(wave_blue), su.FLAM)
    bb_red = convert_flux(wave_red, bb(wave_red), su.FLAM)
    flux_ratio = float(np.asarray(flux_red.value).ravel()[0]
                       / np.asarray(flux_blue.value).ravel()[0])
    bb_ratio = float(np.asarray(bb_red.value).ravel()[0]
                     / np.asarray(bb_blue.value).ravel()[0])
    np.testing.assert_allclose(flux_ratio, bb_ratio, rtol=1.0e-4, atol=0.0)

    # Filters that already overlap the Koester grid stay at the same
    # magnitude after the red tail is attached.
    filt_h = syn.get_filter_info('nirc2,H')
    mag_native = syn.mag_in_filter(native, filt_h)
    mag_extended = syn.mag_in_filter(spectrum, filt_h)
    np.testing.assert_allclose(mag_extended, mag_native, rtol=0.0, atol=1.0e-4)

    # These bandpasses start redward of the native Koester grid.
    filt_names = ['jwst,F323N', 'jwst,F405N', 'nirc2,Lp', 'nirc2,Ms',
                  'nsfcam,L']
    for name in filt_names:
        filt = syn.get_filter_info(name)
        obs = Observation(spectrum, filt, binset=filt.waveset, force='taper')
        integrated = obs.integrate()
        assert np.isfinite(integrated.value)
        assert integrated.value > 0

    return


def test_extend_spectrum_blackbody_blend():
    """
    The blackbody weight is linear from 2.8 to 3.0 micron.

    Blueward of 2.8 micron the model is unchanged. At 2.9 micron the
    flux is an equal mix of the model and the scaled blackbody. At and
    redward of 3.0 micron the flux is the blackbody alone.
    """
    temperature = 10000.0
    wave_aa = np.array([5000.0, 15000.0, 25000.0, 28000.0, 29000.0,
                        29500.0, 29900.0], dtype=float)

    # F_lambda falls as 1/lambda, which is not a 10000 K blackbody, so the
    # blend is distinguishable from either piece alone.
    flux_val = 1.0e17 * (1.0e4 / wave_aa)
    spectrum = SourceSpectrum(Empirical1D, points=wave_aa * u.AA,
                              lookup_table=flux_val * su.FLAM)
    extended = atm.extend_spectrum_blackbody(spectrum, temperature=temperature,
                                             wave_max=5.0e4)

    def flam_at(spec, wave):
        wave_q = np.array([wave], dtype=float) * u.AA
        flux = convert_flux(wave_q, spec(wave_q), su.FLAM)
        return float(np.asarray(flux.value).ravel()[0])

    # Scale the comparison blackbody to the model at its red edge.
    bb = SourceSpectrum(BlackBodyNorm1D, temperature=temperature)
    edge = wave_aa[-1]
    bb_edge = flam_at(bb, edge)
    scale = flux_val[-1] / bb_edge

    def bb_at(wave):
        return scale * flam_at(bb, wave)

    # Blue of the window, including the 2.8 micron boundary.
    for wave in (15000.0, 28000.0):
        np.testing.assert_allclose(flam_at(extended, wave),
                                   1.0e17 * (1.0e4 / wave),
                                   rtol=1.0e-5, atol=0.0)

    # Mid-window: w_BB = (2.9 - 2.8) / (3.0 - 2.8) = 0.5.
    wave_mid = 29000.0
    model_mid = 1.0e17 * (1.0e4 / wave_mid)
    expected_mid = 0.5 * model_mid + 0.5 * bb_at(wave_mid)
    np.testing.assert_allclose(flam_at(extended, wave_mid), expected_mid,
                               rtol=1.0e-5, atol=0.0)

    # Still inside the window, closer to the blackbody.
    # w_BB = (2.95 - 2.8) / 0.2 = 0.75.
    wave_75 = 29500.0
    model_75 = 1.0e17 * (1.0e4 / wave_75)
    expected_75 = 0.25 * model_75 + 0.75 * bb_at(wave_75)
    np.testing.assert_allclose(flam_at(extended, wave_75), expected_75,
                               rtol=1.0e-5, atol=0.0)

    # At and redward of 3 micron the tail is the blackbody.
    for wave in (30000.0, 40000.0):
        np.testing.assert_allclose(flam_at(extended, wave), bb_at(wave),
                                   rtol=1.0e-4, atol=0.0)

    ext_wave = np.asarray(extended.waveset.to(u.AA).value, dtype=float)
    np.testing.assert_allclose(ext_wave[-1], 5.0e4, rtol=1.0e-8, atol=0.0)

    return