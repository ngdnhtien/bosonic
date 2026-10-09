# Bosonic quantum error correction: rotation-symmetric ``Spider'' state and logical gate between two bosonic qubits

This repo summarizes my work during 1 year spent at the Centre for Quantum Technologies. Results were obtained during 2023-2024. With frontier LLMs, those unorganized results are now organized, summarized and written up into a nice little report (see `main.pdf`). Of course there will be some "alien maths" along the way, but essentially the core physics was verified by a human (myself). Below this line, everything is LLM-generated.

## The story in five lines

1. **Bosonic QEC** stores a qubit in one oscillator; the *rotation-symmetric* codes (cat, binomial, …) are built from a primitive state with support on the Fock states |mK⟩.
2. **Zaw et al., PRA 106, 032222 (2022)** found a state that maximally violates a classical bound for a precessing oscillator. It is K-fold rotation symmetric by construction: the **spider state**.
3. We turn it into a rotation code. It can be *certified* from sign measurements alone (no other code can), but …
4. … it **fails the Knill–Laflamme conditions**: its two codewords have photon numbers in a fixed ratio of 3, no reweighting fixes that, and the exact state has infinite energy. A witness, not a code.
5. So the bottleneck is the gates every rotation code needs. We engineer the **CROT gate** from the drive-tuned **cross-Kerr** between two superconducting cavities sharing a transmon, find that the self-Kerr produced by the same drive is the real enemy, and then reproduce the published solution — the single-tone drive placed next to a two-cavity resonance of **Qin et al., arXiv:2609.07076 (2026)** — with the same code.

The full argument, written as a PRA-style article, is `report/main.pdf`.

## Reading order

| step | notebook | what you will see |
|---|---|---|
| 1 | `notebooks/01_precession_protocol_and_spider_state.ipynb` | the protocol, the classical bound, the spider state, its coefficients, Wigner function and divergent energy |
| 2 | `notebooks/02_spider_as_rotation_code.ipynb` | codewords, photon numbers, the certification witness, comparison with cat/binomial primitives |
| 3 | `notebooks/03_knill_laflamme_test.ipynb` | the KL test, the verdict, and why "renormalising" cannot save it |
| 4 | `notebooks/04_crot_gate_via_cross_kerr.ipynb` | drive-induced Kerr, cross-Kerr between two cavities, CROT fidelity, the recovered four-mode study, and the targeted resonance of Qin et al. 2026 |

Each notebook runs in a few minutes on a laptop and calls the Python package `spiderqec/`:

| module | contents |
|---|---|
| `spiderqec/precession.py` | score operator Q_K, classical bound, spider state, convergence |
| `spiderqec/rotation_codes.py` | codewords from a primitive, cat/binomial/Pegg–Barnett primitives, photon number, phase uncertainty, Wigner |
| `spiderqec/knill_laflamme.py` | KL tensor and violation measures |
| `spiderqec/transmon.py` | cavity–transmon Hamiltonians, dressed-state labelling, self-/cross-Kerr, CROT fidelity |

## Folders

```
report/       PRA-style LaTeX article (main.tex, abstract.tex, references.bib, figures/, main.pdf)
notebooks/    the four story notebooks
spiderqec/    the Python package
figures/      every plot the notebooks make, plus figures/legacy/ recovered from the old project
images/       reference circuit diagrams (Grimsmo et al. 2020, CC-BY)
literature/   reading list, merged bibliography, and the original study notes (LaTeX)
data/         small result tables and the recovered four-mode Kerr data
misc/         old notebooks and notes that are not on the story line
```
Every folder has its own `README.md` labelling each file.

## Running

```
pip install -r requirements.txt
jupyter lab notebooks/
```
The report compiles with `cd report && pdflatex main && bibtex main && pdflatex main && pdflatex main`
(REVTeX 4.2 for the journal layout; it falls back to the article class if REVTeX is absent).

## License

MIT — see `LICENSE`.
