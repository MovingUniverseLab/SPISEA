---

title: 'SPISEA v3.1: Synthetic stellar populations with updated photometry, brown dwarfs, and binary evolution'
tags:
- Python
- astronomy
- synthetic photometry
- stellar populations
- star clusters
authors:
- name: Mace J. Huston
  orcid: 0000-0003-4591-3201
  affiliation: 1
- name: Jessica R. Lu
  orcid: 0000-0001-9611-0009
  affiliation: 1
- name: Matthew W. Hosek Jr.
  affiliation: 2
- name: Lingfeng Wei
  affiliation: 2
- name: Natasha S. Abrams
  orcid: 0000-0002-0287-3783
  affiliation: 1
- name: Riley Patlak
  affiliation: 1
- name: Caitlin Begbie
  affiliation: 1
- name: Deval Deliwala
  affiliation: 1
- name: J. Lewis-Merrill
  affiliation: 3
affiliations:
- index: 1
  name: Department of Astronomy, University of California, Berkeley, CA 94720, USA
- index: 2
  name: Department of Physics and Astronomy, University of California, Los Angeles, CA 90025, USA
- index: 3
  name: Department of Earth and Planetary Sciences, University of California, Los Angeles, CA 90025, USA
date: 14 November 2026
bibliography: references.bib

---

# Summary

SPISEA (Stellar Population Interface for Stellar Evolution and Atmospheres) builds single-age, single-metallicity stellar populations with user control of the initial mass function (IMF), multiplicity, extinction, filters, and stellar evolution and atmosphere grids. The software was introduced by [@Hosek2020] and extended into the substellar regime by [@Begbie2026] (submitted). Version 3.1 modernizes the synthetic-photometry stack, expands brown-dwarf and binary modeling, and documents preferred model combinations for new work.

Astronomers use SPISEA to generate isochrones and resolved or unresolved clusters for Galactic and extragalactic star-cluster studies, including crowded fields where multiplicity and differential extinction matter. The software returns stellar parameters and multi-band photometry (and optionally compact remnants) in a reproducible table-based workflow.

# Statement of need

Observers and theorists need flexible, open tools that couple evolution tracks, atmosphere spectra, and filter throughputs without rebuilding a photometry pipeline for each project. Legacy SPISEA depended on `pysynphot`, which is no longer maintained and is incompatible with current NumPy stacks. Users also need brown-dwarf physics across different metallicities [@Begbie2026], optional analytic IMFs beyond broken power laws, updated multiplicity statistics, and a path to evolve interacting binaries.

SPISEA 3.1 addresses these gaps for researchers who simulate star clusters, calibrate CMDs, or forward-model surveys. It remains a focused cluster/isochrone engine rather than a full galaxy SPS code: the target audience is Galactic and Local Group stellar-population work where control of IMF, multiples, and extinction laws is primary. Typical workflows build an isochrone once, cache it as a FITS table, then draw many clusters with different IMFs or multiplicity assumptions—so a maintainable photometry backend and explicit model-grid versioning matter as much as new physics modules.

# State of the field

Galaxy-scale SPS packages (for example FSPS, BPASS-oriented tools, and related frameworks) emphasize integrated light and large grids. SPISEA instead emphasizes *resolved* populations, arbitrary filters via STScI throughputs, and plug-in evolution/atmosphere objects [@Hosek2020; @Begbie2026]. Relative to those earlier SPISEA papers, version 3.1 does not replace galaxy-scale SPS codes; it restores a maintainable photometry backend (`synphot` / `stsynphot`; [@STScIDevelopmentTeam2018; @STScIDevelopmentTeam2020]) and adds science capabilities that those tools do not package together for cluster work (metallicity-dependent brown dwarfs, Offner-style multiplicity, optional `initial_mass_function` draws, and COSMIC binary evolution [@Breivik2020]).

We build on [@Hosek2020] and [@Begbie2026] rather than contributing only patches to `pysynphot` because the research need is an end-to-end, versioned cluster simulator with documented model grids—not a throughput library alone.

# Software design

SPISEA keeps a small object graph: evolution and atmosphere callables feed `IsochronePhot`; an IMF (with optional multiplicity) and optional IFMR feed `ResolvedCluster` / related classes. Version 3.1 preserves that API while swapping the spectrum engine to Astropy-affiliated `synphot`/`stsynphot`, tabulating atmospheres and bandpasses for faster photometry without changing magnitudes for matched inputs.

Trade-offs for 3.1:

- **Preferred vs default models.** Documentation recommends `MergedFlameSkimmerRedDiamondbackMIST` with `get_merged_atmosphere_v3` (and Hosek18b extinction) for new work, while `IsochronePhot` still defaults to the legacy MISTv1 / atmosphere v1 / Nishiyama09 combination so existing scripts and paper-reproduction notebooks keep working.
- **Atmosphere packaging.** SPISEA-specific CDBS catalogs (R=250) ship separately from STScI Reference Atlases; optional high-resolution archives are large and unused by default.
- **Binaries.** COSMIC evolves mass-exchanging binaries outside the static isochrone grid, so it is slower and limited to binaries (not triples+), with explicit guards in the cluster constructor.
- **Optional IMF package.** `IMF_package` lazily imports `initial_mass_function` so core installs stay light while Chabrier-like forms become available [@Richardson2026].
- **Multiplicity.** New work can opt into Offner et al. (2023) unresolved or resolved multiplicity models, including brown-dwarf companions [@Offner2023], while retaining the SPISEA v2.5 multiplicity classes for reproducibility.

# Research impact statement

SPISEA underpins peer-reviewed cluster and Galactic-center studies building on the original release [@Hosek2020] and the submitted brown-dwarf extension [@Begbie2026]. Version 3.1 ships rendered tutorials and paper-reproduction notebooks, a public model-file server with `grid_version` 3.1 archives, and Read the Docs builds for the `dev`/`stable` lines. Near-term significance is a supportable dependency stack (post-`pysynphot`) plus documented preferred models so new analyses and survey forward-models can adopt 3.1 grids without reverse-engineering notebook paths.

# AI usage disclosure

Generative AI assistants (including Cursor) were used to help draft documentation, release checklists, and portions of this paper text, and to assist with routine code edits and test scaffolding during the 3.1 development cycle. Scientific content, model choices, and software behavior were reviewed and verified by the authors against SPISEA tests, tutorials, and published model references before release.

# Acknowledgements

We thank the SPISEA contributor community and STScI for synphot/stsynphot. Funding acknowledgements will be finalized before JOSS submission.

# References
