# literature/

`references.bib` is the single, de-duplicated bibliography of the project (it is what
`report/main.tex` uses). The papers themselves are not redistributed; DOIs are in the bib file.
`notes/` holds the study notes written during the original project, kept as compilable LaTeX.

## Reading list, in story order

**Why bosonic QEC.** Chuang–Leung–Yamamoto 1997 (first bosonic codes); Gottesman–Kitaev–Preskill 2001
(translation-symmetric codes); Joshi–Noh–Gao 2021 (review of bosonic qubits in circuit QED);
Ofek et al. 2016 and Sivak et al. 2023 (break-even with cat and GKP codes).

**Rotation-symmetric codes.** Grimsmo–Combes–Baragiola 2020 is the framework: codewords from a primitive,
Fock-grid structure, number/rotation distances, the hybrid Steane–Knill gadget built from CROT and a phase
measurement, and the embedded phase uncertainty. Michael et al. 2016 (binomial codes), Mirrahimi et al. 2014
(cat codes), Albert et al. 2018 (performance comparison under loss).

**The spider state.** Zaw–Aw–Lasmar–Scarani 2022 generalise Tsirelson 2006: sign measurements at K equally
spaced times of a uniformly precessing variable give a score bounded by (1 + 1/K)/2 classically; the
quantum optimum (0.709 for K = 3) is reached by a state with support on |mK⟩ — the spider state.

**Error correction criteria.** Knill–Laflamme 1997; Nielsen & Chuang for the textbook statement.

**Phase measurement (for the ancilla readout).** Wiseman–Killip 1998 (adaptive single-shot phase
measurement, POVMs for canonical vs heterodyne); Martin et al. 2020 (canonical phase measurement with
quantum feedback in circuit QED).

**Hardware for the CROT gate.** Koch et al. 2007 (transmon); Blais et al. 2021 (circuit-QED review);
Zhang et al. 2022 — the key hardware paper: an off-resonant drive on the transmon tunes the self- and
cross-Kerr of the cavities, through zero.

**The published outcome of this project.** Qin, My, Copetudo, Kasper, Gao & Ng, *Single-tone
drive-enhanced CROT gate for bosonic quantum error correction*, arXiv:2609.07076 (2026). A single-tone
drive placed next to the two-cavity conversion resonance ω_a = ω_b + ε_10 enhances the cross-Kerr by
more than 60× while the self-Kerrs stay small; a driven side transmon cancels the residual self-Kerr of
cavity a; the Grimsmo gadget is then simulated end-to-end (binomial N=2 data, cat M=1 ancilla) with
CROT fidelities 0.96–0.99 and a ~19 µs error-correction cycle that beats the bare cavity at high loss.
Notebook 04 (Step 6) reproduces the resonance with the dressed-ladder code.

## notes/ — study notes from the original project (author's own)

| folder | what it is | status |
|---|---|---|
| `spider_code_manuscript_2024/` | REVTeX draft "Spider code for single-mode bosonic quantum error detection" (T. D. Nguyen), with its compiled PDF. The seed of `report/main.tex`. | draft; bibliography unfinished |
| `rotation_symmetric_bosonic_codes/` | Study notes on Grimsmo et al.: codewords, Fock/phase-space pictures, Kronecker comb, logical gates, teleportation gadget, modular phase measurement (with reference figures in `images/`). | complete notes, Feb 2024 |
| `best_state_for_phase_measurement/` | "Ranking states for phase measurement": why the CROT gate and phase measurement are needed; mean modular phase and embedded phase uncertainty for cat, binomial, Pegg–Barnett and lazy primitives. | notes with open questions |
| `crot_gate_design/` | "How to make a good CROT gate": derivation of the driven cavity–transmon Hamiltonian, dressed-state labelling, Kerr fitting, convergence checks, and the strategy for arbitrary CROT angles. | notes, contains leftover template text |

Each folder has `main.tex` and `refs.bib`; compile with `pdflatex main && bibtex main && pdflatex main && pdflatex main`
(the manuscript needs REVTeX 4.2).
