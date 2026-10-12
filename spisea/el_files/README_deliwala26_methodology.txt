README: Methodology for creation of Deliwala+26 extinction curve in SPISEA

Near-IR

- Started with Fritz+11 (hereafter F11) interpolated wavelength grid in SPISEA.
- A single power law was close, but gave F182M = 1.31 instead of 1.21.
- Used a smooth spline correction to match our short-wavelength ratios, reconnecting with F11 at 2.7 μm.

F323N and F405N

- Raised the existing F11 values inside the narrow bandpasses by about 61% and 2%, respectively.
- Smoothly joined these adjustments back to F11 over 0.05 μm on either side; retained F11 elsewhere.

Matching the Deliwala+26 extinction ratios in the JWST filters

- Used a SPISEA 8 Gyr RC spectrum with [Fe/H] = +0.3, reddened by F11 to 2.5 mag at 2.166 μm
- Added a small amount of extinction using the trial curve and calculated the magnitude changes through our NIRCam filters.
- Adjusted the curve until those changes, relative to F212N, reproduced all five measured differential ratios.
- Changing the reference extinction from 1.8–3.2 mag gives F115W ratios of 3.15–2.95 without refitting; F323N and F405N barely change.
