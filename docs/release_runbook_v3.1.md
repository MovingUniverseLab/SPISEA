# SPISEA v3.1 release runbook

Work top to bottom. Check boxes as you go.

**Current tip (when this was written):** `dev` @ `07a3ec5` (synced with `origin/dev`; release PR `dev`→`main` still to open — `gh` token invalid locally).

---

## Phase 0 — Finish pre-release polish

These were requested but may still be open on `dev`:

### 0.1 Docs: which atmosphere grids need which CDBS tarballs

- [x] Add a table on `docs/atmo_models.rst` (and a short pointer in `docs/getting_started.rst` / data-download page) mapping:
  - **STScI** `synphot1_throughput-master.tar` → filters (`comp/`, `mtab/`)
  - **STScI** `synphot5_pheonix-models.tar` → base Phoenix / ATLAS-style grids (`phoenix`, `ck04models`, `k93models`, …)
  - **SPISEA** `spisea_cdbs.tar.gz` (R=250, required) → FlameSkimmer, Meisner, Phillips, BTSettl rebinned, WD Koester, `merged_*` used by v1/v3
  - **SPISEA** `spisea_cdbs_highres.tar.gz` (optional) → high-R versions of those
  - Call out grids that are **not** in the standard tarballs (e.g. `cmfgen*`, `nextgen`, `AMESdusty` from issue #170) as unsupported unless the user adds catalogs themselves
- [x] Commit + push to `dev` (`5d68cb8`)

### 0.2 README: Begbie citation

- [x] Update `README.md` cite line to Hosek+2020 **and** Begbie+2026 (BD / substellar use)
- [x] Mirror the same wording on `docs/index.rst` if desired
- [x] Commit + push to `dev`

### 0.3 JOSS paper (short)

- [x] Finish `docs/joss_paper/paper.md` (Summary, Statement of need, brief design, impact, AI disclosure)
- [x] Keep it JOSS-short: ~750–1000 words of prose; cite Hosek2020, Begbie2026, synphot/stsynphot, COSMIC as needed
- [x] Fill `docs/joss_paper/references.bib` for any new cites
- [x] Commit + push to `dev` (submission can wait until after the tag)

### 0.4 Decide defaults (one-time decision)

Pick **A** or **B** before tagging:

| | Choice |
|---|---|
| **A (safer)** | Leave code defaults as `MISTv1` / `get_merged_atmosphere_v1` / `Nishiyama09`; docs already say “preferred for new work” is Flame+v3+Hosek18b |
| **B (aligned)** | Change `synthetic.py` defaults to `MergedFlameSkimmerRedDiamondbackMIST` / `get_merged_atmosphere_v3` / `Hosek18b` and update tests |

- [x] **Decision: A (safer).** Code defaults stay `MISTv1` / `get_merged_atmosphere_v1` / `Nishiyama09`. Preferred for new work remains Flame+v3+Hosek18b in the docs. Recorded in the release PR / notes.
- [x] Confirmed in `spisea/synthetic.py`: `default_evo_model = evolution.MISTv1()`, `default_atm_func = atm.get_merged_atmosphere_v1`, `default_red_law = reddening.RedLawNishiyama09()`.

**Note:** SPISEA has never shipped a conda package (`conda install spisea`). Install remains git + `PYTHONPATH` (and Docker). Do not block the release on conda-forge.

---

## Phase 1 — Verify data + software

### 1.1 Model tarballs on the Berkeley server

- [ ] Confirm `spisea_models.tar.gz` and `spisea_cdbs.tar.gz` on https://w.astro.berkeley.edu/~jlu/spisea/ are **grid_version 3.1**
- [ ] Spot-check: FlameSkimmer / merged `flame_red_diamondback_mist` present after unpack
- [ ] Confirm `docs/for_model_fileserver/index.html` matches sizes/dates if you re-upload

### 1.2 Clean install smoke test

- [x] Used existing `py314` env (Python 3.14; `spisea-synphot` not present on this machine) with synphot 1.7 / stsynphot git + COSMIC
- [x] Set `PYTHONPATH` to repo, `PYSYN_CDBS=/Users/jlu/work/models/cdbs/`, `SPISEA_MODELS=/Users/jlu/work/models/` (local grid_version 3.1)
- [x] Tutorial 01-style cluster (default MISTv1 / atm v1 / N09) + short BD path (`MergedFlameSkimmerRedDiamondbackMIST` + `get_merged_atmosphere_v3` + Hosek18b)
- [x] Full tests: `conda run -n py314 python -m pytest spisea/tests -q` — all passed (includes COSMIC)

```bash
conda run -n <env> python -m pytest spisea/tests -q
```

(include COSMIC only if `cosmic-popsynth` is installed)

### 1.3 Docs build

- [x] Latest RTD `dev` build green (`2e2d115` / #35899389 Built)
- [x] Spot-check https://spisea.readthedocs.io/en/dev/ — nav (Installation / Quick Start / API / Development), Hosek+Begbie cites, 3.1 changelog present

### 1.4 Changelog freeze

- [x] Set `docs/index.rst` date to `3.1 (2026-11-14)`
- [x] Skimmed 3.1 bullets; fixed typo `systynphot` → `stsynphot`
- [x] Commit + push with release-prep changes

---

## Phase 2 — Release PR (`dev` → `main`)

### 2.1 Open the PR

```bash
gh pr create --base main --head dev \
  --title "Release SPISEA v3.1.0" \
  --body "$(cat <<'EOF'
## Summary
- SPISEA 3.1.0: synphot/stsynphot, brown-dwarf grids, Offner multiplicity, IMF_package, atmosphere v1/v3, docs/notebooks.

## Defaults decision (runbook 0.4)
- **Option A (safer):** leave code defaults as MISTv1 / get_merged_atmosphere_v1 / Nishiyama09.
  Preferred for new work remains Flame+v3+Hosek18b in the docs.

## Test plan
- [x] Local pytest on py314 (+ COSMIC)
- [x] Tutorial smoke test with grid_version 3.1 models
- [x] RTD dev build green

## Release notes
See docs/index.rst Change Log — 3.1.
EOF
)"
```

- [ ] PR opened — blocked: `gh` keyring token for `jluastro` is invalid (`gh auth login` needed). Compare URL: https://github.com/MovingUniverseLab/SPISEA/compare/main...dev?expand=1
- [ ] Review / fix anything that comes up
- [ ] Merge with **merge commit** (repo convention)

### 2.2 Sync local `main`

```bash
git fetch origin
git checkout main
git pull --ff-only origin main
```

- [ ] Local `main` matches release tip

---

## Phase 3 — Tag and GitHub Release

### 3.1 Tag

Confirm `pyproject.toml` / `__version__` are `3.1.0`, then:

```bash
git checkout main
git pull --ff-only
git tag -a v3.1.0 -m "SPISEA v3.1.0"
git push origin v3.1.0
```

- [ ] Tag pushed

### 3.2 GitHub Release

```bash
gh release create v3.1.0 --title "SPISEA v3.1.0" --notes "$(cat <<'EOF'
## SPISEA 3.1.0

Highlights:
- pysynphot → synphot/stsynphot
- Brown-dwarf evolution/atmospheres (MergedFlameSkimmerRedDiamondbackMIST, get_merged_atmosphere_v3)
- Offner+2023 multiplicity; optional IMF_package (Chabrier etc.)
- Docs: preferred models, tutorials, paper notebooks

**Data:** Requires grid_version 3.1 model archives from https://w.astro.berkeley.edu/~jlu/spisea/

**Install note:** Until stsynphot 1.7 is on PyPI, install the git URL in docs/getting_started.rst.

Cite Hosek et al. (2020); also cite Begbie et al. (2026) for brown-dwarf work.
EOF
)"
```

- [ ] Release published
- [ ] Zenodo (if you archive from GitHub tags) — confirm new DOI or concept DOI update

---

## Phase 4 — Read the Docs + Docker

### 4.1 RTD

- [ ] Ensure version **stable** builds from tag `v3.1.0` (or `main` after merge)
- [ ] Ensure **latest** tracks `main` (so it is no longer July-era content)
- [ ] Set project **default version** to `stable` (or `latest` once `main` is 3.1) so https://spisea.readthedocs.io/ is not stuck on old docs
- [ ] Verify https://spisea.readthedocs.io/en/stable/ and the bare URL

### 4.2 Docker (optional but recommended)

- [ ] Rebuild image from current `Dockerfile` (Python 3.14 + synphot stack + note on model download)
- [ ] Push e.g. `amigahub/spisea:v3.1` (and update `docs/getting_started.rst` if Hub tag changes)
- [ ] Or document that Hub `v1` is obsolete until refreshed

---

## Phase 5 — Announce + housekeeping

- [ ] Post short note (Slack / collaborators): tag, docs URL, data URL, stsynphot caveat, citations
- [ ] Close or comment on release-related issues (#170 can stay open as “document unsupported grids” if the table is enough for now)
- [ ] Submit JOSS when `docs/joss_paper/paper.md` is ready (separate from the software tag)
- [ ] Leave non-blockers for later: #165 origin column, #164 COSMIC grid, #155 CI, #109, #50, feature requests

---

## Quick “done” criteria

| Check | Pass when |
|---|---|
| Code | `main` has 3.1; tag `v3.1.0` exists |
| Data | Berkeley 3.1 tarballs verified |
| Docs | RTD default shows 3.1; preferred models + CDBS mapping live |
| Install | Fresh env + models can run tutorial 01 |
| Cite | README has Hosek + Begbie |
| Paper | JOSS draft complete (submit when ready) |

---

## Suggested next step

Phase **0.1 → 0.2 → 0.3** on `dev`, then **1.1–1.4**, then open the `dev` → `main` PR (Phase 2).
