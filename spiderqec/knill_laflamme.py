"""
Knill-Laflamme (KL) test for single-mode bosonic codes.

A code with codewords {|0_L>, |1_L>} corrects the error set {E_a} iff

    <i_L| E_a^dag E_b |j_L> = C_ab  delta_ij        for all a, b, i, j,

i.e. (i) errors never map one codeword onto the other ("orthogonality",
off-diagonal elements vanish) and (ii) both codewords see the same error
statistics ("non-deformation", the diagonal elements are equal).

For rotation codes and loss errors a^l with l < K, condition (i) is automatic
from the Fock-grid spacing.  Condition (ii) for l = 1 reads
<0_L| n |0_L> = <1_L| n |1_L>: both codewords must have the *same mean photon
number*.  The SPIDER code fails exactly there.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import diags


# ---------------------------------------------------------------------------
# Error operators (dense, Fock basis, dimension N)
# ---------------------------------------------------------------------------
def annihilation(N: int) -> np.ndarray:
    return np.diag(np.sqrt(np.arange(1, N)), k=1)


def number(N: int) -> np.ndarray:
    return np.diag(np.arange(N).astype(float))


def loss_errors(N: int, max_order: int) -> dict[str, np.ndarray]:
    """{I, a, a^2, ..., a^max_order}."""
    a = annihilation(N)
    errs = {"I": np.eye(N)}
    cur = np.eye(N)
    for l in range(1, max_order + 1):
        cur = cur @ a
        errs[f"a^{l}" if l > 1 else "a"] = cur
    return errs


def dephasing_errors(N: int, max_order: int = 1) -> dict[str, np.ndarray]:
    """{n, n^2, ...}  (small-rotation errors exp(i theta n) expanded in theta)."""
    n = number(N)
    errs = {}
    cur = np.eye(N)
    for l in range(1, max_order + 1):
        cur = cur @ n
        errs[f"n^{l}" if l > 1 else "n"] = cur
    return errs


# ---------------------------------------------------------------------------
# KL matrix and violation measures
# ---------------------------------------------------------------------------
def kl_tensor(codewords: list[np.ndarray], errors: dict[str, np.ndarray]):
    """M[a, b, i, j] = <i| E_a^dag E_b |j>.  Returns (names, M)."""
    names = list(errors)
    E = [errors[k] for k in names]
    W = np.array(codewords)  # shape (2, N)
    M = np.zeros((len(E), len(E), len(W), len(W)), dtype=complex)
    for a, Ea in enumerate(E):
        EaW = (Ea @ W.T).T  # E_a |i>
        for b, Eb in enumerate(E):
            EbW = (Eb @ W.T).T
            M[a, b] = EaW.conj() @ EbW.T
    return names, M


def kl_violation(codewords, errors):
    """Summarise how far the code is from satisfying KL for the given error set.

    Returns a dict with
      'orthogonality' : max |<0| E_a^dag E_b |1>|                (should be 0)
      'deformation'   : max |<0|E_a^dag E_b|0> - <1|E_a^dag E_b|1>|  (should be 0)
      'table'         : per-pair (a, b) numbers, for printing.
    """
    names, M = kl_tensor(codewords, errors)
    orth = np.abs(M[:, :, 0, 1])
    deform = np.abs(M[:, :, 0, 0] - M[:, :, 1, 1])
    table = []
    for a, na in enumerate(names):
        for b, nb in enumerate(names):
            if b < a:
                continue
            table.append(
                dict(
                    pair=f"{na}†{nb}",
                    c00=float(np.real(M[a, b, 0, 0])),
                    c11=float(np.real(M[a, b, 1, 1])),
                    off=float(orth[a, b]),
                    deform=float(deform[a, b]),
                )
            )
    return dict(orthogonality=float(orth.max()), deformation=float(deform.max()), table=table, names=names, M=M)


def mean_photon_mismatch(zero: np.ndarray, one: np.ndarray) -> float:
    """<1|n|1> - <0|n|0>: the leading KL deformation for single photon loss."""
    n = np.arange(len(zero))
    return float(np.sum(n * np.abs(one) ** 2) - np.sum(n * np.abs(zero) ** 2))


# ---------------------------------------------------------------------------
# "Can we repair the spider code by reweighting?"  (trade-off helper)
# ---------------------------------------------------------------------------
def reweighted_plus(zero: np.ndarray, one: np.ndarray, lam: float) -> np.ndarray:
    """Return |+'> = (|0_L> + lam |1_L>)/norm.

    Changing the relative weight of the two codewords does NOT change the
    codewords themselves, hence cannot fix a KL deformation; it only changes
    which logical state is the precession-optimal one.  Used in the notebooks to
    make that point quantitatively.
    """
    v = zero + lam * one
    return v / np.linalg.norm(v)


def equalise_mean_photon(zero: np.ndarray, one: np.ndarray, power: float):
    """Deform |1_L> -> N (n+1)^(-power/2) |1_L> (and renormalise) so that its mean
    photon number is pulled towards that of |0_L>.  Returns the deformed |1_L'>.

    A naive 'renormalisation' of the Fock-grid coefficients.  The notebooks show
    that whatever power equalises the photon numbers destroys the maximal
    precession-score property of |+_L>, i.e. the SPIDER state and a KL-correct
    code cannot coexist.
    """
    n = np.arange(len(one))
    v = one * (n + 1.0) ** (-power / 2.0)
    return v / np.linalg.norm(v)
