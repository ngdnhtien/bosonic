# figures/

Every plot in this folder is produced by one of the notebooks in `notebooks/`
(PNG for viewing, PDF for the report), except the `legacy/` sub-folder, which
holds the figures recovered from the original project files whose generating
code was lost.

## 01_precession/ — Notebook 01
| file | what it shows |
|---|---|
| `score_operator_structure_K3` | Matrix elements of the score operator Q_3 in the Fock basis; only m ≡ n (mod 3) survive (Kronecker comb). |
| `spider_coefficients_K3` | Fock-grid coefficients s_m of the spider state \|S_3⟩; sign pattern (−1)^⌈m/2⌉. |
| `spider_convergence_K3` | Left: score P_3 vs Fock-grid truncation with power-law fit to P_3^∞ = 0.7094. Right: mean photon number ∝ √n (divergent). |
| `spider_coefficient_tail_K3` | Log-log tail of \|s_m\|: power-law decay. |
| `spider_wigner_K3_K7` | Wigner functions of \|S_3⟩ and \|S_7⟩ (the "spider legs"). |

## 02_rotation_code/ — Notebook 02
| file | what it shows |
|---|---|
| `spider_codewords_coefficients_K3` | \|+_L⟩ = \|S_3⟩ and the codewords \|0_L⟩ (even grid), \|1_L⟩ (odd grid). |
| `spider_codewords_wigner_K3` | Wigner functions of \|0_L⟩, \|1_L⟩, \|+_L⟩, \|−_L⟩. |
| `spider_fidelity_witness_K3` | Fidelity with \|S_3⟩ guaranteed by a measured precession score (spectral-gap witness). |
| `embedded_phase_uncertainty_vs_nbar` | Embedded phase uncertainty Δ_K vs mean photon number for cat, binomial, Pegg–Barnett, lazy and spider primitives, K = 1…4. |

## 03_knill_laflamme/ — Notebook 03
| file | what it shows |
|---|---|
| `spider_kl_deformation_vs_truncation` | Photon-number mismatch between \|0_L⟩ and \|1_L⟩ (the KL deformation for one loss) vs truncation; the relative mismatch locks at 0.5. |
| `spider_kl_repair_tradeoff` | Trying to repair the code by reweighting \|1_L⟩: KL deformation for one and two losses (left) vs precession score and fidelity with the spider (right). |

## 04_crot_gate/ — Notebook 04
| file | what it shows |
|---|---|
| `self_kerr_vs_drive_one_cavity` | Drive-induced self-Kerr K, β, σ of one cavity vs drive amplitude (format of Zhang et al. 2022), including the Kerr-free point. |
| `cross_kerr_vs_drive_two_cavities` | Self-Kerr K_A, K_B and cross-Kerr χ_AB of two cavities sharing a driven transmon. |
| `legacy_four_mode_kerr_chain` | Re-plot of the legacy data for the chains A–q_c–B and A–q_c–B–q_b (`data/legacy/`). |
| `crot_fidelity_and_time_vs_drive` | Average CROT(2π/9) fidelity with/without self-Kerr, and gate time, vs drive. |
| `crot_fidelity_vs_residual_self_kerr` | How small the residual self-Kerr must be for a 99 % CROT. |
| `qin2026_targeted_resonance_cross_kerr` | Reproduction of the targeted resonance of Qin et al. (arXiv:2609.07076): K_A, K_B, χ_AB vs drive detuning at \|Ω_d/δ_d\| = 1, the hybridised window, the on-off ratio and \|χ_AB/K_A\|. |

## legacy/ — recovered from the original project (not regenerated)
| file | origin / content |
|---|---|
| `spider_coefficients_K3_legacy.png` | Original `ss_coeffs.png`: coefficients of \|S_3⟩, \|0_L⟩, \|1_L⟩ at truncation 5000. Reproduced by Notebook 02. |
| `spider_wigner_K3_legacy.png`, `spider_wigner_K7_legacy.png` | Original Wigner functions of the spider state. Reproduced by Notebook 01. |
| `spider_score_convergence_K3_legacy.png` | Original score-vs-truncation plot (P_3^∞ = 0.7094). Reproduced by Notebook 01. |
| `spider_mean_photon_vs_truncation_K3/K7_legacy.png` | Original mean-photon plots (∝ √n). Reproduced by Notebook 01. |
| `spider_epu_vs_truncation_K3/K7_legacy.png` | Original embedded-phase-uncertainty plots. **Caution:** they report a negative Δ_K, which is unphysical (a normalisation slip); Notebook 02 gives the corrected values. |
| `epu_comparison_order1…4_legacy.png` | Original `divison1-4.png` from `phase_meas.ipynb`: Δ_K vs n̄ for cat/binomial/Pegg–Barnett/lazy (and squeezed cat for K = 1). Reproduced by Notebook 02. |
| `epu_comparison_cat_lazy_legacy.png`, `epu_coherent_state_legacy.png` | Earlier versions of the same comparison, from the note `best_state_for_phase_measurement`. |
| `four_mode_kerr_chain_legacy.png` | Original `plot_transparent.png` from `sab_plot.ipynb`: K_A, K_B, χ_AB for the A–q_c–B–q_b chain. Re-plotted in Notebook 04. |
| `ode_truncation_error_legacy.png` | Truncation-error study of the QuTiP time-dependent solver (from the CROT design note). |
| `dressed_cat_evolution_wigner_legacy.png` | Wigner snapshots of a dressed cat state evolving under the full cavity–transmon Hamiltonian (from the CROT design note). |
