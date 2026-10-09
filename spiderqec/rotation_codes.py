"""
Rotation-symmetric bosonic (RSB) codes.

An order-K rotation code (Grimsmo, Combes & Baragiola, PRX 10, 011058 (2020))
is a two-dimensional subspace on which the discrete rotation

    Z_L = exp(i pi n / K)

acts as the logical Z.  Consequently

    |0_L> has support on |2mK>,     |1_L> has support on |(2m+1)K>,
    |+_L>, |-_L> have support on every |mK>.

Any normalised state with support on the full K-grid, such as the SPIDER state
|S_K>, can therefore be *declared* to be |+_L>; its even-m and odd-m parts then
define |0_L> and |1_L>.  This module builds codewords from such a primitive,
and from the classical primitives (cat, binomial, Pegg-Barnett, "lazy") used for
comparison, and computes the quantities used to compare codes:

* mean photon number of each codeword and of the code,
* mean modular phase  <e^{i K theta}>  and the embedded phase uncertainty
  Delta_K = 1/|<e^{iK theta}>|^2 - 1   (Grimsmo et al., Eq. (37)-(38)),
* Wigner functions (via QuTiP).
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import gammaln

try:  # QuTiP is only needed for Wigner functions
    import qutip as qt
except ImportError:  # pragma: no cover
    qt = None


# ---------------------------------------------------------------------------
# Projectors and codewords
# ---------------------------------------------------------------------------
def grid_projector(N: int, K: int, ell: int) -> np.ndarray:
    """Diagonal projector Pi^ell_{2K} onto Fock states |2mK + ell>, as an N-vector mask.

    By the Kronecker comb, Pi^ell_{2K} = (1/2K) sum_{j=0}^{2K-1} (e^{-i pi ell/K} Z_L)^j.
    """
    n = np.arange(N)
    return ((n - ell) % (2 * K) == 0).astype(float)


def codewords_from_plus(psi_plus: np.ndarray, K: int) -> tuple[np.ndarray, np.ndarray]:
    """Split a K-grid state |+_L> = sum_m s_m |mK> into normalised |0_L>, |1_L>.

    |0_L> is proportional to the even-m part, |1_L> to the odd-m part.
    Raises if psi_plus has weight outside the K-grid.
    """
    N = len(psi_plus)
    n = np.arange(N)
    off_grid = np.sum(np.abs(psi_plus[n % K != 0]) ** 2)
    if off_grid > 1e-10:
        raise ValueError(f"state has weight {off_grid:.2e} outside the K-grid")
    zero = psi_plus * grid_projector(N, K, 0)
    one = psi_plus * grid_projector(N, K, K)
    return zero / np.linalg.norm(zero), one / np.linalg.norm(one)


def dual_codewords(zero: np.ndarray, one: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    plus = (zero + one) / math.sqrt(2)
    minus = (zero - one) / math.sqrt(2)
    return plus, minus


def code_projector(zero: np.ndarray, one: np.ndarray) -> np.ndarray:
    return np.outer(zero, zero.conj()) + np.outer(one, one.conj())


# ---------------------------------------------------------------------------
# Classical primitives, all returned as normalised |+_K> states of length N
# ---------------------------------------------------------------------------
def coherent(alpha: complex, N: int) -> np.ndarray:
    n = np.arange(N)
    log_mag = -0.5 * abs(alpha) ** 2 + n * np.log(abs(alpha) + 1e-300) - 0.5 * gammaln(n + 1)
    phase = np.exp(1j * n * np.angle(alpha))
    psi = np.exp(log_mag) * phase
    return psi / np.linalg.norm(psi)


def cat_plus(K: int, alpha: float, N: int) -> np.ndarray:
    """|+_K^cat> ~ sum_{j=0}^{K-1} |alpha e^{2 pi i j/K}>  (support on |mK>)."""
    psi = coherent(alpha, N)
    n = np.arange(N)
    psi = psi * (n % K == 0)  # projection onto the K-grid == sum over rotated copies
    return psi / np.linalg.norm(psi)


def binomial_plus(K: int, L: int, N: int) -> np.ndarray:
    """|+_{K,bin}> = sum_{k=0}^{L} sqrt( C(L,k)/2^L ) |kK>   (Michael et al. 2016)."""
    psi = np.zeros(N)
    for k in range(L + 1):
        if k * K < N:
            psi[k * K] = math.sqrt(math.comb(L, k) / 2**L)
    return psi / np.linalg.norm(psi)


def pegg_barnett_plus(K: int, s: int, N: int) -> np.ndarray:
    """Truncated phase state |phi=0, s> projected on the K-grid."""
    psi = np.zeros(N)
    for n in range(0, min(s, N), K):
        psi[n] = 1.0
    return psi / np.linalg.norm(psi)


def lazy_plus(K: int, b: float, N: int) -> np.ndarray:
    """Geometric ("lazy") state sum_m b^{-m} |mK>, |b| > 1."""
    if abs(b) <= 1:
        raise ValueError("|b| must be > 1")
    psi = np.zeros(N)
    m = np.arange(0, N, K)
    psi[m] = b ** (-np.arange(len(m)))
    return psi / np.linalg.norm(psi)


# ---------------------------------------------------------------------------
# Figures of merit
# ---------------------------------------------------------------------------
def mean_photon(psi: np.ndarray) -> float:
    return float(np.real(np.sum(np.arange(len(psi)) * np.abs(psi) ** 2)))


def photon_variance(psi: np.ndarray) -> float:
    n = np.arange(len(psi))
    p = np.abs(psi) ** 2
    return float(np.sum(n**2 * p) - np.sum(n * p) ** 2)


def code_mean_photon(zero: np.ndarray, one: np.ndarray) -> float:
    """n_code = tr(P n)/2 = (n_0 + n_1)/2."""
    return 0.5 * (mean_photon(zero) + mean_photon(one))


def mean_modular_phase(psi_plus: np.ndarray, K: int) -> float:
    """<e^{i K theta}> = sum_m |f_{mK} f_{(m+1)K}|  for a state on the K-grid.

    For a state written as sum_m f_{mK} |mK> this is the modulus of the
    expectation of the "number-translation" operator sum_n |n><n+K|, which is
    how the mean modular phase is defined in Grimsmo et al.  It equals 1 only
    for the (unnormalisable) phase state.  (The 1/2 that appears in some notes
    belongs to a different normalisation of the codewords and must not be used
    with a normalised |+_L>.)
    """
    f = psi_plus[::K]
    return float(np.sum(np.abs(f[:-1] * f[1:])))


def embedded_phase_uncertainty(psi_plus: np.ndarray, K: int) -> float:
    """Delta_K(theta) = 1/|<e^{iK theta}>|^2 - 1  (0 for a perfect phase state)."""
    m = mean_modular_phase(psi_plus, K)
    return 1.0 / m**2 - 1.0


# ---------------------------------------------------------------------------
# Wigner function
# ---------------------------------------------------------------------------
def wigner(psi: np.ndarray, xvec=None, truncate: int | None = None):
    """Wigner function W(x, p) on a grid, via QuTiP.  Returns (xvec, W)."""
    if qt is None:
        raise ImportError("QuTiP is required for Wigner functions")
    if xvec is None:
        xvec = np.linspace(-6, 6, 241)
    v = np.asarray(psi, dtype=complex)
    if truncate is not None:
        v = v[:truncate]
    state = qt.Qobj(v.reshape(-1, 1))
    W = qt.wigner(state, xvec, xvec)
    return xvec, W


def qobj(psi: np.ndarray):
    """Fock-basis vector -> QuTiP ket."""
    if qt is None:
        raise ImportError("QuTiP is required")
    return qt.Qobj(np.asarray(psi, dtype=complex).reshape(-1, 1))
