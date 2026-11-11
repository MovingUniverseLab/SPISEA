.. _api:

==============
API Reference
==============

Class and function docstrings are rendered on the topic pages.
This page links to those references.

Isochrones and clusters
-----------------------

* :class:`synthetic.Isochrone`
* :class:`synthetic.IsochronePhot`
* :class:`synthetic.IsochronePhotExternalEvolution`
* :func:`synthetic.calc_ab_vega_filter_conversion`
* :func:`synthetic.calc_st_vega_filter_conversion`
* :class:`synthetic.Cluster`
* :class:`synthetic.ResolvedCluster`
* :class:`synthetic.ResolvedClusterDiffRedden`
* :class:`synthetic.UnresolvedCluster`

Evolution models
----------------

* :class:`evolution.StellarEvolution`
* :class:`evolution.MISTv1`
* :class:`evolution.MergedBaraffePisaEkstromParsec`
* :class:`evolution.MergedPhillipsBaraffePisaEkstromParsec`
* :class:`evolution.Baraffe15`
* :class:`evolution.Ekstrom12`
* :class:`evolution.Parsec`
* :class:`evolution.Pisa`
* :class:`evolution.Phillips2020`
* :class:`evolution.COSMIC`

Atmosphere models
-----------------

* :func:`atmospheres.get_merged_atmosphere`
* :func:`atmospheres.get_wd_atmosphere`
* :func:`atmospheres.extend_spectrum_blackbody`
* :func:`atmospheres.get_bb_atmosphere`
* :func:`atmospheres.get_castelli_atmosphere`
* :func:`atmospheres.get_phoenixv16_atmosphere`
* :func:`atmospheres.get_BTSettl_2015_atmosphere`
* :func:`atmospheres.get_BTSettl_atmosphere`
* :func:`atmospheres.get_wdKoester_atmosphere`
* :func:`atmospheres.get_kurucz_atmosphere`
* :func:`atmospheres.get_phoenix_atmosphere`
* :func:`atmospheres.get_Phillips2020_atmosphere`
* :func:`atmospheres.get_Meisner2023_atmosphere`

Initial mass function
---------------------

* :class:`imf.imf.IMF`
* :class:`imf.imf.IMF_broken_powerlaw`
* :class:`imf.imf.IMFSalpeter1955`
* :class:`imf.imf.Miller_Scalo_1979`
* :class:`imf.imf.Kennicutt_1983`
* :class:`imf.imf.Kroupa_2001`
* :class:`imf.imf.Weidner_Kroupa_2004`
* :class:`imf.imf.Salpeter_Kirkpatrick_2024`

Multiplicity
------------

* :class:`imf.multiplicity.MultiplicityUnresolved`
* :class:`imf.multiplicity.MultiplicityPiecewisePowerLaw`
* :class:`imf.multiplicity.MultiplicityLogistic`
* :class:`imf.multiplicity.MultiplicityUnresolvedOffner2023`
* :class:`imf.multiplicity.MultiplicityResolvedDK`
* :class:`imf.multiplicity.MultiplicityResolvedOffner2023`

Initial-final mass relation
---------------------------

* :class:`ifmr.IFMR_Spera15`
* :class:`ifmr.IFMR_Raithel18`
* :class:`ifmr.IFMR_N20_Sukhbold`

Extinction laws
---------------

* :class:`reddening.RedLawPowerLaw`
* :class:`reddening.RedLawBrokenPowerLaw`
* :class:`reddening.RedLawRiekeLebofsky`
* :class:`reddening.RedLawCardelli`
* :class:`reddening.RedLawIndebetouw05`
* :class:`reddening.RedLawRomanZuniga07`
* :class:`reddening.RedLawFitzpatrick09`
* :class:`reddening.RedLawNishiyama09`
* :class:`reddening.RedLawSchoedel10`
* :class:`reddening.RedLawFritz11`
* :class:`reddening.RedLawDamineli16`
* :class:`reddening.RedLawDeMarchi16`
* :class:`reddening.RedLawSchlafly16`
* :class:`reddening.RedLawHosek18b`
* :class:`reddening.RedLawNoguerasLara18`
* :class:`reddening.RedLawNoguerasLara20`
* :class:`reddening.RedLawSODC`
