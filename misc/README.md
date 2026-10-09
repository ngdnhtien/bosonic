# misc/

Material from the original project that does not sit on the story line
(spider state → rotation code → Knill–Laflamme verdict → CROT gate). Kept lightly cleaned:
duplicates removed, hard-coded local paths (`os.chdir('/Users/...')`) commented out,
and the embedded inline images (60 MB of SVG) stripped from the notebook outputs — the printed
numbers are kept, and the figures that mattered are in `figures/legacy/` and `misc/legacy_figures/`.
Nothing else touched. Nothing here is needed to run the notebooks in `notebooks/`.

## legacy_notebooks/ — the original QuTiP notebooks (`bQEC-main`)

| notebook | what it does | relation to the new code |
|---|---|---|
| `reproduce.ipynb` | Reproduces Zhang et al. 2022: RWA cavity–transmon Hamiltonian, level-crossing checks vs truncation, dressed-state labelling, polynomial fit of the dressed ladder, sweep of K, β, σ vs drive. | Rewritten as `spiderqec.transmon.self_kerr_vs_drive` (Notebook 04). |
| `phase_meas.ipynb` | Mean modular phase and embedded phase uncertainty for cat (1–4 legs), binomial, Pegg–Barnett, lazy and squeezed-cat primitives at equal n̄; produced `divison1-4.png`. | Rewritten in `spiderqec.rotation_codes` (Notebook 02). |
| `convergence.ipynb` | Time-dependent (lab-frame, cos φ transmon, Blackman-shaped drive) simulation of a dressed coherent state and Nelder–Mead fit of an effective Kerr Hamiltonian to the final state. | Superseded by the dressed-ladder method; kept for reference. |
| `fit_Kerr.ipynb` | Minimal test of fitting a Kerr Hamiltonian to a simulated evolution. | — |
| `transmon_approx.ipynb` | Compares the full cos φ transmon with its quartic (K = 2) and sextic (K = 3) Duffing truncations by fidelity of 10 µs evolutions (hours-long runs). | Background for the Duffing approximation used in `spiderqec.transmon`. |
| `cross_Kerr/main.ipynb` | Method 1 (eigen-energy fit) and Method 2 (emulated Wigner-tomography experiment) for the self-Kerr. | Method 1 → `spiderqec.transmon`. |
| `cross_Kerr/standalone_Kerr.ipynb` | Clean two-mode version of the eigen-energy fit; K_A = −2.63 kHz at zero drive; truncation study up to (100, 70). | Reproduced in Notebook 04 (−2.62 kHz). |
| `cross_Kerr/two_mode.ipynb` | Three-mode A–c–B lab-frame simulation (dims 20×20×10; the 10 µs run took ~3 days) and an unfinished fit of (δ_a, δ_b, K_a, K_b, χ_ab). | Replaced by `spiderqec.transmon.kerr_matrix` (seconds). |
| `cross_Kerr/trivial.ipynb`, `cross_Kerr/time_dependent_fit.ipynb` | One-mode Kerr-evolution fits (sanity checks of the fitting procedure). | — |
| `cross_Kerr/time_dependentH.ipynb` | Time-step convergence of `mesolve` with a time-dependent drive. | — |
| `cross_Kerr/sab_plot.ipynb` | Plots the four-mode chain results from `data/legacy/SAB_*.npy`. | Re-plotted in Notebook 04. |
| `self_Kerr/cat.ipynb`, `self_Kerr/coherent.ipynb` | Perturbative formula for K_B (from Zhang et al.) vs numerical time-domain fits, for a cat / coherent initial state; produced `data/legacy/*.pkl`. | — |
| `gallery/playground.ipynb`, `gallery/wigner.ipynb` | Wigner-function sketches of coherent, squeezed, cat (N = 2, 3) and kitten states. | — |
| `utility.py` | `fit_function` (scipy curve_fit wrapper) and `are_elements_unique`, imported by the notebooks above. | — |

The legacy notebooks also imported a `constant.py` that was not in the archive (only its `.pyc`); the
constants they need are all defined inline, so nothing is lost.

## legacy_figures/

* `level_crossings/crossings_{dimCavity}_{dimTransmon}.png` — eigenvalue-vs-coupling plots used to choose
  Hilbert-space truncations in `reproduce.ipynb`.
* `wigner_gallery/` — the sketches from `gallery/wigner.ipynb`.

## notes/

| folder | content |
|---|---|
| `frame/` | "Frame in quantum simulation of a driven oscillator": lab frame vs drive frame for a driven qubit, with two QuTiP plots. Tutorial-level. |
| `damn_oscillator/` | Homework-style expansion of a Heisenberg–Langevin equation around one period. Unrelated. |
| `phase_measurement/` | "How do you measure the phase of an electromagnetic field?": POVMs for canonical/heterodyne phase measurement (Wiseman & Killip 1998) and how beam splitters / IQ mixers implement dyne measurements. Background for the ancilla readout, not used in the story. |
| `canonical_phase_measurement/` | Short reading note on Martin et al., Nat. Phys. 2020 (canonical phase measurement with feedback). |

## What was dropped

REVTeX/AIP/AAPM sample files (`aipsamp.tex`, `aapmsamp.tex`, `sorsamp.tex`, `fig_1.eps`, `vid_1*.eps`,
`*Notes.bib`), Overleaf's `frog.jpg`, `.DS_Store`, `__pycache__`, `.ipynb_checkpoints`, and the
`apssamp.synctex.gz` build artefact. Three identical copies of the bosonic-code note images were reduced to one.
