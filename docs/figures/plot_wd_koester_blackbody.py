#!/usr/bin/env python
"""
Generate white-dwarf Koester plus blackbody figures for the docs.

Writes ``docs/figures/wd_koester_bb_sed.png`` and
``docs/figures/wd_koester_bb_join.png``.

Run from the repository root::

    python docs/figures/plot_wd_koester_blackbody.py
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from astropy import units as u
from synphot.models import BlackBodyNorm1D
from synphot.spectrum import SourceSpectrum
from synphot.units import convert_flux
from synphot import units as su


# Allow running without installing the package.
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from spisea import atmospheres
from spisea import synthetic


# log g = 8 is a typical white-dwarf surface gravity on the Koester grid.
_LOGG = 8.0
_TEMPERATURES = [8000.0, 15000.0, 25000.0, 40000.0]

# Mid-IR filters that start redward of the Koester grid, plus two NIR
# filters that still overlap it. Spans are the wavelengths where the
# throughput exceeds 1% of its peak.
_FILTERS = ['nirc2,H', 'nirc2,Kp', 'jwst,F323N', 'jwst,F405N',
            'nirc2,Lp', 'nirc2,Ms']

_KOESTER_COLOR = '#0072B2'
_BLEND_COLOR = '#009E73'
_TAIL_COLOR = '#D55E00'
_BB_COLOR = '0.45'
_BLEND_LO_UM = 2.8
_BLEND_HI_UM = 3.0


def flam_on_waveset(spectrum):
    """
    Evaluate a spectrum in FLAM on its own wavelength grid.

    Parameters
    ----------
    spectrum : synphot.spectrum.SourceSpectrum
        Atmosphere or blackbody spectrum.

    Returns
    -------
    wave_um : numpy.ndarray
        Wavelength in microns. Shape ``[N]``.
    flux : numpy.ndarray
        F_lambda in FLAM. Shape ``[N]``.
    """
    wave_aa = np.asarray(spectrum.waveset.to(u.AA).value, dtype=float)
    flux = convert_flux(spectrum.waveset, spectrum(spectrum.waveset), su.FLAM)
    flux = np.asarray(flux.value, dtype=float).ravel()

    # Shape of wave_um: [N].
    wave_um = wave_aa * 1.0e-4

    return wave_um, flux


def scaled_blackbody(temperature, wave_um, wave_edge, flux_edge):
    """
    Blackbody scaled to a model flux at one wavelength.

    Parameters
    ----------
    temperature : float
        Blackbody temperature in Kelvin.
    wave_um : array_like
        Wavelengths in microns. Shape ``[N]``.
    wave_edge : float
        Wavelength, in microns, where the blackbody is matched to the model.
    flux_edge : float
        Model F_lambda at ``wave_edge``, in FLAM.

    Returns
    -------
    flux : numpy.ndarray
        Scaled blackbody F_lambda in FLAM. Shape ``[N]``.
    """
    wave_aa = np.asarray(wave_um, dtype=float) * 1.0e4
    bb = SourceSpectrum(BlackBodyNorm1D, temperature=temperature)

    # Shape of the edge sample: [1].
    edge = np.array([wave_edge * 1.0e4]) * u.AA
    bb_edge = convert_flux(edge, bb(edge), su.FLAM)
    scale = flux_edge / float(np.asarray(bb_edge.value).ravel()[0])

    wave = wave_aa * u.AA
    bb_flux = convert_flux(wave, bb(wave), su.FLAM)
    flux = np.asarray(bb_flux.value, dtype=float).ravel() * scale

    return flux


def load_stitched_spectrum(temperature, gravity=_LOGG):
    """
    Load a Koester atmosphere and the blackbody-blended spectrum.

    Parameters
    ----------
    temperature : float
        Effective temperature in Kelvin.
    gravity : float, optional
        log g in cgs. Default is 8.

    Returns
    -------
    pieces : dict
        Koester and returned spectra, plus the flux used to normalize
        plots at 3 micron.
    """
    native = atmospheres.get_wdKoester_atmosphere(metallicity=0.0,
                                                  temperature=temperature,
                                                  gravity=gravity)
    stitched = atmospheres.get_wd_atmosphere(metallicity=0.0,
                                             temperature=temperature,
                                             gravity=gravity,
                                             verbose=False)

    native_wave, native_flux = flam_on_waveset(native)
    ext_wave, ext_flux = flam_on_waveset(stitched)

    # Normalize at 3 micron, where the returned spectrum is all blackbody.
    at_3um = np.argmin(np.abs(ext_wave - _BLEND_HI_UM))
    norm = float(ext_flux[at_3um])

    return {
        'temperature': temperature,
        'native_wave': native_wave,
        'native_flux': native_flux,
        'ext_wave': ext_wave,
        'ext_flux': ext_flux,
        'norm': norm,
    }


def _segment(wave, flux, wave_lo, wave_hi):
    """
    Slice a spectrum to a wavelength interval, in microns.

    Parameters
    ----------
    wave : array_like
        Wavelength in microns. Shape ``[N]``.
    flux : array_like
        Flux. Shape ``[N]``.
    wave_lo : float
        Inclusive blue edge, in microns.
    wave_hi : float
        Inclusive red edge, in microns.

    Returns
    -------
    wave_seg, flux_seg : numpy.ndarray
        Samples inside the interval.
    """
    keep = (np.asarray(wave) >= wave_lo) & (np.asarray(wave) <= wave_hi)
    return np.asarray(wave)[keep], np.asarray(flux)[keep]


def filter_spans(names=_FILTERS):
    """
    Wavelength span of each filter, in microns.

    Parameters
    ----------
    names : list of str, optional
        SPISEA filter names.

    Returns
    -------
    spans : list of tuple
        ``(name, wave_lo, wave_hi)`` with wavelengths in microns, limited
        to samples whose throughput is at least 1% of the peak.
    """
    spans = []
    for name in names:
        filt = synthetic.get_filter_info(name, rebin=False)
        wave_um = np.asarray(filt.waveset.to(u.AA).value, dtype=float) * 1.0e-4
        thru = np.asarray(filt(filt.waveset).value, dtype=float).ravel()
        keep = thru > (0.01 * np.nanmax(thru))
        spans.append((name, float(wave_um[keep][0]), float(wave_um[keep][-1])))

    return spans


def _style_spectrum_axis(ax):
    """
    Inward ticks on a spectrum panel.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to style.

    Returns
    -------
    None
    """
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.grid(False)

    return None


def _plot_returned_segments(ax, case, wave_lo, wave_hi):
    """
    Plot the returned spectrum in Koester, blend, and blackbody colors.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to draw on.
    case : dict
        One entry from ``load_stitched_spectrum``.
    wave_lo : float
        Blue edge of the plotted window, in microns.
    wave_hi : float
        Red edge of the plotted window, in microns.

    Returns
    -------
    None
    """
    norm = case['norm']
    wave = case['ext_wave']
    flux = case['ext_flux'] / norm

    # Blue of 2.8 micron is the Koester model. The middle segment is the
    # weighted average. Red of 3 micron is the blackbody alone.
    segments = [
        (wave_lo, _BLEND_LO_UM, _KOESTER_COLOR, 0.9),
        (_BLEND_LO_UM, _BLEND_HI_UM, _BLEND_COLOR, 1.6),
        (_BLEND_HI_UM, wave_hi, _TAIL_COLOR, 1.6),
    ]
    for seg_lo, seg_hi, color, lw in segments:
        w_seg, f_seg = _segment(wave, flux, max(seg_lo, wave_lo),
                                min(seg_hi, wave_hi))
        if len(w_seg) == 0:
            continue
        ax.plot(w_seg, f_seg, color=color, lw=lw, zorder=4, solid_capstyle='butt')

    return None


def plot_sed(cases, path):
    """
    Full 0.1–10 micron SEDs for several white-dwarf temperatures.

    Parameters
    ----------
    cases : list of dict
        Stitched spectra from ``load_stitched_spectrum``.
    path : str
        Output PNG path.

    Returns
    -------
    None
    """
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.6), sharex=True)
    fig.suptitle(r'White-dwarf spectra at $\log g = 8$', fontsize=13)

    # Comparison blackbody on a smooth grid across the plotted window.
    wave_bb = np.logspace(np.log10(0.1), np.log10(10.0), 500)

    # Shared log limits from the returned spectra, ignoring a few wild lines.
    stacked = []
    for case in cases:
        keep = ((case['ext_wave'] >= 0.1) & (case['ext_wave'] <= 10.0)
                & (case['ext_flux'] > 0))
        stacked.append(case['ext_flux'][keep] / case['norm'])
    stacked = np.concatenate(stacked)
    ymin = 10.0 ** np.floor(np.log10(np.percentile(stacked, 0.5)))
    ymax = 10.0 ** np.ceil(np.log10(np.percentile(stacked, 99.5)))

    for ax, case in zip(axes.ravel(), cases):
        # Blackbody scaled to the Koester red edge, drawn at all wavelengths.
        edge = float(case['native_wave'][-1])
        edge_flux = float(case['native_flux'][-1])
        bb_flux = scaled_blackbody(case['temperature'], wave_bb, edge, edge_flux)

        ax.axvspan(_BLEND_LO_UM, _BLEND_HI_UM, color='#d9f2ea', zorder=0)
        ax.plot(wave_bb, bb_flux / case['norm'], color=_BB_COLOR, ls='--',
                lw=1.1, zorder=2)
        _plot_returned_segments(ax, case, 0.1, 10.0)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.1, 10.0)
        ax.set_ylim(ymin, ymax)
        ax.set_title(r'$T_{\mathrm{eff}} = %d$ K' % case['temperature'])
        _style_spectrum_axis(ax)

    for ax in axes[:, 0]:
        ax.set_ylabel(r'$F_\lambda \,/\, F_\lambda(3\,\mu\mathrm{m})$')
    for ax in axes[1, :]:
        ax.set_xlabel(r'Wavelength ($\mu$m)')

    legend_handles = [
        Line2D([0], [0], color=_KOESTER_COLOR, lw=1.4,
               label=r'Koester, $\lambda < 2.8\,\mu$m'),
        Line2D([0], [0], color=_BLEND_COLOR, lw=1.6,
               label=r'Blend, $2.8$–$3.0\,\mu$m'),
        Line2D([0], [0], color=_TAIL_COLOR, lw=1.6,
               label=r'Blackbody, $\lambda > 3\,\mu$m'),
        Line2D([0], [0], color=_BB_COLOR, ls='--', lw=1.1,
               label='Blackbody at all wavelengths'),
        Patch(facecolor='#d9f2ea', edgecolor='none',
              label=r'Blend window'),
    ]
    fig.legend(handles=legend_handles, loc='lower center', ncol=2,
               frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.08, 1, 0.96))
    fig.savefig(path, dpi=160, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    return None


def _percent_from_blackbody(wave, flux, wave_bb, bb_flux):
    """
    Percent difference of a spectrum from a blackbody sampled on ``wave_bb``.

    Parameters
    ----------
    wave : array_like
        Wavelength in microns. Shape ``[N]``.
    flux : array_like
        F_lambda on ``wave``. Shape ``[N]``.
    wave_bb : array_like
        Blackbody wavelength grid in microns. Shape ``[M]``.
    bb_flux : array_like
        Blackbody F_lambda on ``wave_bb``. Shape ``[M]``.

    Returns
    -------
    percent : numpy.ndarray
        ``100 * (F / F_BB - 1)``. Shape ``[N]``.
    """
    bb_at_wave = np.interp(np.asarray(wave), np.asarray(wave_bb),
                           np.asarray(bb_flux))
    percent = 100.0 * (np.asarray(flux) / bb_at_wave - 1.0)

    return percent


def plot_join(cases, spans, path):
    """
    Zoom on the blend, as a percent difference from the scaled blackbody.

    The Koester and blackbody continua agree to about a percent near 3
    micron, so the weight is invisible on a linear F_lambda plot. This
    figure shows that residual.

    Parameters
    ----------
    cases : list of dict
        Stitched spectra from ``load_stitched_spectrum``.
    spans : list of tuple
        Filter ``(name, wave_lo, wave_hi)`` in microns.
    path : str
        Output PNG path.

    Returns
    -------
    None
    """
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.6), sharex=True)
    fig.suptitle(r'Blend weight near 3 $\mu$m, $\log g = 8$', fontsize=13)

    short_label = {
        'jwst,F323N': 'F323N',
        'jwst,F405N': 'F405N',
        'nirc2,Lp': 'Lp',
        'nirc2,Ms': 'Ms',
    }
    x_lo, x_hi = 2.50, 4.20
    wave_bb = np.linspace(x_lo, x_hi, 600)

    for ax, case in zip(axes.ravel(), cases):
        edge = float(case['native_wave'][-1])
        edge_flux = float(case['native_flux'][-1])
        bb_flux = scaled_blackbody(case['temperature'], wave_bb, edge, edge_flux)

        native_pct = _percent_from_blackbody(case['native_wave'],
                                             case['native_flux'],
                                             wave_bb, bb_flux)
        ext_pct = _percent_from_blackbody(case['ext_wave'], case['ext_flux'],
                                          wave_bb, bb_flux)

        ax.axvspan(_BLEND_LO_UM, _BLEND_HI_UM, color='#d9f2ea', zorder=0)
        for name, filt_lo, filt_hi in spans:
            if name not in short_label:
                continue
            if (filt_hi < x_lo) or (filt_lo > x_hi):
                continue
            ax.axvspan(max(filt_lo, x_lo), min(filt_hi, x_hi), color='#f6d3b8',
                       alpha=0.7, lw=0, zorder=1)

        native_w, native_y = _segment(case['native_wave'], native_pct, x_lo, x_hi)
        ax.axhline(0.0, color=_BB_COLOR, ls='--', lw=1.2, zorder=2)
        ax.plot(native_w, native_y, color=_KOESTER_COLOR, lw=1.3, zorder=3)
        # Returned spectrum only from the start of the blend, so the Koester
        # curve stays visible blueward of 2.8 micron.
        for seg_lo, seg_hi, color, lw in (
                (_BLEND_LO_UM, _BLEND_HI_UM, _BLEND_COLOR, 2.0),
                (_BLEND_HI_UM, x_hi, _TAIL_COLOR, 1.8)):
            w_seg, y_seg = _segment(case['ext_wave'], ext_pct, seg_lo, seg_hi)
            ax.plot(w_seg, y_seg, color=color, lw=lw, zorder=4)

        in_window = (native_w >= x_lo) & (native_w <= _BLEND_HI_UM)
        y_span = native_y[in_window]
        y_max = float(np.nanmax(y_span)) if len(y_span) else 1.0
        pad = max(0.15 * y_max, 0.05)
        ax.set_ylim(-0.35 * pad, y_max + pad)
        ax.set_xlim(x_lo, x_hi)
        ax.set_title(r'$T_{\mathrm{eff}} = %d$ K' % case['temperature'])
        _style_spectrum_axis(ax)
        ax.text(0.5 * (_BLEND_LO_UM + _BLEND_HI_UM),
                ax.get_ylim()[1] - 0.06 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
                r'$w=0 \rightarrow 1$', ha='center', va='top', fontsize=8,
                color='#066b4f')

    # Filter names sit inside the first panel, above the zero residual.
    top_left = axes[0, 0]
    y0, y1 = top_left.get_ylim()
    for name, filt_lo, filt_hi in spans:
        if name not in short_label:
            continue
        if (filt_hi < x_lo) or (filt_lo > x_hi):
            continue
        center = 0.5 * (max(filt_lo, x_lo) + min(filt_hi, x_hi))
        top_left.text(center, y1 - 0.04 * (y1 - y0), short_label[name],
                      ha='center', va='top', fontsize=7, color='0.25')

    for ax in axes[:, 0]:
        ax.set_ylabel(r'$(F_\lambda / F_{\mathrm{BB}} - 1)$ (%)')
    for ax in axes[1, :]:
        ax.set_xlabel(r'Wavelength ($\mu$m)')

    legend_handles = [
        Line2D([0], [0], color=_KOESTER_COLOR, lw=1.4, label='Koester model'),
        Line2D([0], [0], color=_BLEND_COLOR, lw=2.0,
               label=r'Returned blend, $2.8$–$3.0\,\mu$m'),
        Line2D([0], [0], color=_TAIL_COLOR, lw=1.8,
               label=r'Returned blackbody, $\lambda > 3\,\mu$m'),
        Line2D([0], [0], color=_BB_COLOR, ls='--', lw=1.2,
               label='Blackbody scaled to the Koester edge'),
        Patch(facecolor='#f6d3b8', edgecolor='none', label='Mid-IR filters'),
    ]
    fig.legend(handles=legend_handles, loc='lower center', ncol=2,
               frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.05))
    fig.tight_layout(rect=(0, 0.10, 1, 0.96))
    fig.savefig(path, dpi=160, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    return None


def main():
    """
    Write the two documentation figures.

    Returns
    -------
    None
    """
    cases = []
    for temperature in _TEMPERATURES:
        case = load_stitched_spectrum(temperature)
        cases.append(case)
        print('T=%d K  koester %.3f–%.3f um  returned %.3f–%.3f um' % (
            temperature, case['native_wave'][0], case['native_wave'][-1],
            case['ext_wave'][0], case['ext_wave'][-1]))

    spans = filter_spans()
    for name, wave_lo, wave_hi in spans:
        print('filter %s  %.3f–%.3f um' % (name, wave_lo, wave_hi))

    out_dir = os.path.dirname(__file__)
    sed_path = os.path.join(out_dir, 'wd_koester_bb_sed.png')
    join_path = os.path.join(out_dir, 'wd_koester_bb_join.png')
    plot_sed(cases, sed_path)
    plot_join(cases, spans, join_path)
    print('Wrote', sed_path)
    print('Wrote', join_path)

    return None


if __name__ == '__main__':
    main()
