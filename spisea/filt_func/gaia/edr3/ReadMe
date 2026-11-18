J/A+A/649/A3     Gaia Early Data Release 3 photometric passbands (Riello+, 2021)
================================================================================
Gaia Early Data Release 3: Photometric content and validation.
    Riello M., De Angeli F., Evans D.W., Montegriffo P., Carrasco J.M,
    Busso G., Palaversa L., Burgess P., Diener C., Davidson M., Rowell N.,
    Fabricius C., Jordi C., Bellazzini M., Pancino E., Harrison D.L.,
    Cacciari C., van Leeuwen F., Hambly N.C., Hodgkin S.T., Osborne P.J.,
    Altavilla G., Barstow M.A., Brown A.G.A., Castellani M., Cowell S.,
    De Luise F., Gilmore G., Giuffrida G., Hidalgo S., Holland G., Marinoni S.,
    Pagani C., Piersimoni A.M., Pulone L., Ragaini S., Rainer M., Richards P.J.,
    Sanna N., Walton N.A., Weiler M., Yoldas A.
    <Astron. Astrophys. 649, A3 (2021)>
    =2021A&A...649A...3R        (SIMBAD/NED BibCode)
================================================================================
ADC_Keywords: Surveys ; Photometry, G band
Keywords: catalogues - surveys - instrumentation: photometers -
          techniques: photometric - Galaxy: general

Abstract:
    Gaia Early Data Release 3 (Gaia EDR3) contains astrometry and
    photometry results for about 1.8 billion sources based on
    observations collected by the European Space Agency Gaia satellite
    during the first 34 months of its operational phase.

    In this paper, we focus on the photometric content, describing the
    input data, the algorithms, the processing, and the validation of the
    results. Particular attention is given to the quality of the data and
    to a number of features that users may need to take into account to
    make the best use of the Gaia EDR3 catalogue.

    The processing broadly followed the same procedure as for Gaia DR2,
    but with significant improvements in several aspects of the blue and
    red photometer (BP and RP) preprocessing and in the photometric
    calibration process. In particular, the treatment of the BP and RP
    background has been updated to include a better estimation of the
    local background, and the detection of crowding effects has been used
    to exclude affected data from the calibrations. The photometric
    calibration models have also been updated to account for flux loss
    over the whole magnitude range. Significant improvements in the
    modelling and calibration of the Gaia point and line spread functions
    have also helped to reduce a number of instrumental effects that were
    still present in DR2.

    Gaia EDR3 contains 1.806 billion sources with G-band photometry and
    1.540 billion sources with GBP and GRP photometry. The median
    uncertainty in the G-band photometry, as measured from the standard
    deviation of the internally calibrated mean photometry for a given
    source, is 0.2mmag at magnitude G=10 to 14, 0.8mmag at G~17, and
    2.6mmag at G~19. The significant magnitude term found in the Gaia DR2
    photometry is no longer visible, and overall there are no trends
    larger than 1mmag/mag. Using one passband over the whole colour and
    magnitude range leaves no systematics above the 1% level in magnitude
    in any of the bands, and a larger systematic is present for a very
    small sample of bright and blue sources. A detailed description of the
    residual systematic effects is provided. Overall the quality of the
    calibrated mean photometry in Gaia EDR3 is superior with respect to
    DR2 for all bands.

Description:
    These tabular data describes the photometric system defined by the G,
    G_BP_ and G_RP_ Gaia bands for Gaia Early Data Release 3.

    The tables provide the full passband and the corresponding zero point
    for each of the photometric bands. The zero points are available both
    in the VEGAMAG and AB systems.

    The passband calibration is based on the modelling of a set of
    corrections applied to pre-launch knowledge of the instrument, to find
    the best match between observed and synthetic photometry for a set of
    calibrators. For EDR3 a large set of calibrators covering a wide range
    of spectra types was used for the passband calibration.

    For these sources reconstructed SEDs were obtained from externally
    calibrated Gaia BP/RP spectra.

File Summary:
--------------------------------------------------------------------------------
 FileName      Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe            80        .   This file
passband.dat     103      781   G, G_BP_ anf G_RP_ passbands used to generate
                                 the magnitudes and astrophysical parameters
                                 included in Gaia EDR3
zeropt.dat        97        2   G, G_BP_ and G_RP_ zero points used to generate
                                 the magnitudes and astrophysical parameters
                                 included in Gaia EDR3
--------------------------------------------------------------------------------

See also:
        I/350 : Gaia EDR3 (Gaia Collaboration, 2020)
 J/A+A/649/A6 : Gaia Catalogue of Nearby Stars - GCNS (Gaia collaboration, 2021)
 J/A+A/649/A7 : MC structure and properties (Gaia Collaboration+, 2021)

Byte-by-byte Description of file: passband.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  7  F7.2  nm      lambda    Wavelength
  10- 23  E14.9 mag     GPb       ?=99.99 G transmissivity curve at the
                                   corresponding wavelength (1)
  26- 39  E14.9 mag   e_GPb       ?=99.99 Uncertainty on the G transmissivity
                                   curve (1)
  42- 55  E14.9 mag     BPPb      ?=99.99 BP transmissivity curve at the
                                   corresponding wavelength (1)
  58- 71  E14.9 mag   e_BPPb      ?=99.99 Uncertainty on the BP transmissivity
                                   curve (1)
  74- 87  E14.9 mag     RPPb      ?=99.99 RP transmissivity curve at the
                                   corresponding wavelength (1)
  90-103  E14.9 mag   e_RPPb      ?=99.99 Uncertainty on the RP transmissivity
                                   curve (1)
--------------------------------------------------------------------------------
Note (1): In correspondence to wavelength values where the passband is not
   defined, both transmissivity and uncertainty are set to 99.99.
--------------------------------------------------------------------------------

Byte-by-byte Description of file: zeropt.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   2- 14 F13.10 ---     GZp       G band zero point
  18- 29 F12.10 ---   e_GZp       G band zero point uncertainty
  32- 44 F13.10 ---     BPZp      G_BP_ band zero point
  48- 59 F12.10 ---   e_BPZp      G_BP_ band zero point uncertainty
  62- 74 F13.10 ---     RPZp      G_RP_ band zero point
  78- 89 F12.10 ---   e_RPZp      G_RP_ band zero point uncertainty
  91- 97  A7    ---     System    [VEGAMAG, AB] Photometric system
--------------------------------------------------------------------------------

Acknowledgements:
    Marco Riello, mriello(at)ast.cam.ac.uk
    Francesca De Angeli, fda(at)ast.cam.ac.uk

================================================================================
(End)     Francesca De Angeli [IoA, UK], Patricia Vannier [CDS]      05-Jan-2021
