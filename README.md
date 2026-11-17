<p align="center">
  <img width="380" height="336" src="./SPISEA_logo_final.png">
</p>

# SPISEA: Stellar Population Interface for Stellar Evolution and Atmospheres

SPISEA is an python package that generates single-age, single-metallicity
populations (i.e. star clusters). It gives the user control over many parameters:

* Cluster characteristics (age, metallicity, mass, distance)
* Total extinction, differential extinction, and extinction law
* Stellar evolution and atmosphere models
* Stellar multiplicity and Initial Mass Function
* Initial-Final Mass Relation
* Photometric filters

Here is a brief list of things that SPISEA can do:

* make a cluster isochrone in many filters using different stellar models
* make a star cluster at any age with an unusual IMF and unresolved multiplicity
* make a spectrum of a star cluster in integrated light

Install with pip:

```bash
pip install SPISEA
```

You still need the CDBS and evolution/atmosphere model archives; see
[Installation](https://spisea.readthedocs.io/en/latest/getting_started.html)
for the full setup. For an editable clone, use `pip install -e .`.

See [documentation](https://github.com/MovingUniverseLab/spisea/tree/main/docs) for
details on installing and running SPISEA. The tutorials, in order, are
[making a cluster](https://github.com/MovingUniverseLab/spisea/blob/main/docs/notebooks/01_make_cluster.ipynb),
[brown dwarfs](https://github.com/MovingUniverseLab/spisea/blob/main/docs/notebooks/02_brown_dwarfs.ipynb),
[resolved multiples](https://github.com/MovingUniverseLab/spisea/blob/main/docs/notebooks/03_resolved_multiples.ipynb),
and [COSMIC evolution](https://github.com/MovingUniverseLab/spisea/blob/main/docs/notebooks/04_cluster_with_cosmic.ipynb).
Paper figures are in
[docs/notebooks/papers](https://github.com/MovingUniverseLab/spisea/tree/main/docs/notebooks/papers).

If you use SPISEA in your research, please cite
[Hosek et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020AJ....160..143H).
If you use the brown-dwarf / substellar capabilities (SPISEA 2.5+),
also cite [Begbie et al. (2026)](https://arxiv.org/abs/2607.14292).

## Change Log
SPISEA is actively supported and is growing in functionality, and
subsequently has had several updates since its
initial release. See the
[main documentation page](https://github.com/MovingUniverseLab/spisea/blob/main/docs/index.rst)
for the change log describing the updates in each release.

## Contributions
We welcome and encourage contributions to SPISEA!
 For feature additions, we ask that users make their
own fork of the repository, make their changes, and then submit a pull
request to the "dev" branch.

All contributions will be acknowledged on the
[contributors page](https://github.com/MovingUniverseLab/spisea/blob/main/docs/contributors.rst) (with permission).
Contributors with features used in code releases will be co-authors in future SPISEA software papers.

## License
This project is Copyright (c) Matthew Hosek Jr., Jessica Lu, Casey
Lam, Abhimat Gautam, Kelly Lockhart, Dongwon Kim, and Siyao Jia and licensed under
the terms of the GNU GPL v3+ license. This package is based upon
the [Astropy package template](https://github.com/astropy/package-template)
which is licensed under the BSD 3-clause license. See the licenses folder for
more information. This program is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  


<sub>Logo by Natasha Abrams</sub>
