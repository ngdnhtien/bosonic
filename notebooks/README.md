# notebooks/

Four notebooks, to be read in order; together they are the argument of `report/main.pdf`.
They import the package from the parent folder (`sys.path.insert(0, "..")`), write their
plots to `figures/<nn>_<topic>/` and small tables to `data/`.

| # | notebook | runtime | core result |
|---|---|---|---|
| 01 | `01_precession_protocol_and_spider_state.ipynb` | ~3 min | P_3^∞ = 0.7094 > 2/3; spider state has support on \|3m⟩, n̄ ∝ √(truncation) |
| 02 | `02_spider_as_rotation_code.ipynb` | ~2 min | codewords, n̄_1/n̄_0 = 3, fidelity witness F ≥ 1 − δ/ε, phase-uncertainty comparison |
| 03 | `03_knill_laflamme_test.ipynb` | ~10 s | orthogonality exact, non-deformation fails by O(n̄); repair kills the maximal violation |
| 04 | `04_crot_gate_via_cross_kerr.ipynb` | ~2 min | K_A = −2.62 kHz undriven, Kerr-free point, χ_AB enhanced ×4 (untargeted) and ×50–80 at the Qin et al. 2026 resonance, CROT needs \|K\| < 0.7 % \|χ\| |

Re-run everything from the command line with

```
cd notebooks
for nb in 0*.ipynb; do jupyter nbconvert --to notebook --execute --inplace "$nb"; done
```
