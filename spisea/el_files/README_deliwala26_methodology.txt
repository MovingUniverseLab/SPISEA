README: Methodology for creation of Deliwala+26 extinction curve in SPISEA

Near-IR
- Started with F11's interpolated wavelength grid in SPISEA.
- Used a smooth spline correction to match our short-wavelength ratios, reconnecting with F11 at 2.7 μm.

F323N and F405N
- Raised the existing F11 values inside the narrow bandpasses by about 61% and 2%, respectively.
- Smoothly joined these adjustments back to F11 over 0.05 μm on either side; retained F11 elsewhere.

Matching the Deliwala+26 extinction ratios in the JWST filters
- I used an RC spectrum from an 8 Gyr MIST population with [Fe/H] = +0.3, selecting the model nearest the median RC temperature.
- I calculated its unreddened magnitudes through the filters using SPISEA.
- I applied 2.5 mag of the trial curve at 2.166\mu m and calculated the magnitudes again.
- I subtracted the unreddened magnitudes to obtain the extinction in each filter, then divided by the F212N extinction.
- I adjusted the curve until these total ratios matched our five measured values: 3.06, 2.21, 1.21, 1.05, and 0.39.
