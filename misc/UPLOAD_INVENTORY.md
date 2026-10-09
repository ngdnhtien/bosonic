# Inventory of the three source archives and where each file went

Three archives were consolidated into this repository. They overlapped heavily:
`spider.zip` is a superset of `bqec/bqec_spider.zip` (adds the compiled PDF and build files),
and `bqec/CAN.zip`, `bqec/phase_measurement.zip`, `bqec/rotation-symmetric_bosonic_code.zip`
share the same template, bibliography and five reference images (three different notes).
20 byte-identical duplicate groups were collapsed.

## `bqec-20261009T094853Z-1-001.zip` (Overleaf/Drive export of LaTeX projects)

| archive member | content | now at |
|---|---|---|
| `bqec_spider.zip` | REVTeX draft of the spider-code manuscript + two figures + REVTeX sample files | `literature/notes/spider_code_manuscript_2024/` (samples dropped) |
| `rotation-symmetric_bosonic_code.zip` | 440-line study note on RSB codes (Feb 2024) with 5 reference images | `literature/notes/rotation_symmetric_bosonic_codes/`; images also in `images/` |
| `best_state_for_CAN.zip` | note "Ranking states for phase measurement" + 3 figures | `literature/notes/best_state_for_phase_measurement/`; figures in `figures/legacy/` |
| `crot_gate.zip` | note "How to make a good CROT gate" + 2 figures | `literature/notes/crot_gate_design/`; figures in `figures/legacy/` |
| `phase_measurement.zip` | note "How do you measure the phase of an EM field?" (Wiseman–Killip, beam splitters, IQ mixers) | `misc/notes/phase_measurement/` |
| `CAN.zip` | short reading note on Martin et al. 2020 | `misc/notes/canonical_phase_measurement/` |
| `frame.zip` | note on lab vs drive frame for a driven qubit + 2 plots | `misc/notes/frame/` |
| `damn_oscillator.zip` | homework-style Heisenberg–Langevin expansion | `misc/notes/damn_oscillator/` |
| `figures/` (8 png) | spider-state results whose generating code was lost | `figures/legacy/spider_*_legacy.png` |

## `spider.zip`

Identical to `bqec_spider.zip` plus `apssamp.pdf` (kept as `spider_code_manuscript_2024.pdf`),
`apssamp.bbl`, `apssamp.synctex.gz`, `apssampNotes.bib` (dropped).

## `bQEC-main.zip` (GitHub export of the QuTiP code, `bqec@huikhoon`)

| archive member | content | now at |
|---|---|---|
| `README.md` | one joke line | replaced |
| `utility.py` | curve_fit wrapper | `misc/legacy_notebooks/utility.py` |
| `reproduce.ipynb`, `convergence.ipynb`, `fit_Kerr.ipynb`, `transmon_approx.ipynb`, `phase_meas.ipynb` | see `misc/README.md` | `misc/legacy_notebooks/` |
| `cross_Kerr/*.ipynb` (7), `self_Kerr/*.ipynb` (2), `gallery/*.ipynb` (2) | see `misc/README.md` | `misc/legacy_notebooks/<subdir>/` |
| `cross_Kerr/sab/*.npy` (3), `plot_transparent.png` | four-mode chain Kerr data and its plot | `data/legacy/`, `figures/legacy/four_mode_kerr_chain_legacy.png` |
| `self_Kerr/*.pkl` (3) | perturbative/fitted self-Kerr data | `data/legacy/` |
| `gallery/divison1-4.png` | phase-uncertainty comparison plots | `figures/legacy/epu_comparison_order*_legacy.png` |
| `gallery/crossings/*.png` (26) | level-crossing diagnostics | `misc/legacy_figures/level_crossings/` |
| `gallery/wigner/*.png` (5) | Wigner sketches | `misc/legacy_figures/wigner_gallery/` |
| `__pycache__/`, `.DS_Store`, `.ipynb_checkpoints/` | build artefacts | dropped (`constant.py` itself was not in the archive) |

## What was missing and had to be rewritten

The notebooks that computed the spider state (score operator, eigenvector, coefficients, Wigner,
phase uncertainty, mean photon number vs truncation) were not in any archive — only their eight
output figures were. `spiderqec/precession.py`, `spiderqec/rotation_codes.py`, `spiderqec/knill_laflamme.py`
and notebooks 01–03 were therefore written from the definitions in the manuscript draft and in
Zaw et al. (2022); they reproduce the legacy figures (same coefficients, same P_3^∞ = 0.7094, same
√n growth of n̄) and correct one of them (the sign of the embedded phase uncertainty).
