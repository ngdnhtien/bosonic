"""
The precession protocol and the SPIDER state.

Physics summary
---------------
A harmonic oscillator with H = (p^2 + x^2)/2 (units hbar = omega = 1) makes x(t)
precess uniformly with period T = 2*pi.  In the protocol of Tsirelson (2006) and
Zaw, Aw, Lasmar & Scarani, PRA 106, 032222 (2022) one measures the *sign* of x at
one of K equally spaced times t_k = k T / K (k = 0..K-1), chosen at random, and
estimates

    P_K = (1/K) * sum_k  Pr[ x(t_k) > 0 ].

Any classical oscillator obeys  P_K <= P_K^c = (1 + 1/K)/2  for odd K
(and P_K = 1/2 exactly for even K).  Quantum mechanically

    P_K = <psi| Q_K |psi>,   Q_K = (1/K) sum_k U_k^dag pos(x) U_k,

with pos(x) the projector on x > 0 and U_k = exp(-i n t_k).  The largest
eigenvalue of Q_K exceeds the classical bound (0.709 > 2/3 for K = 3).  Its
eigenvector is the SPIDER state |S_K>.

Because U_k multiplies Fock state |n> by exp(-2 pi i n k/K), the Kronecker comb
(1/K) sum_k exp(2 pi i (m-n) k/K) = delta_{m = n mod K} shows that Q_K is simply
pos(x) with every matrix element between Fock states of *different residue mod K*
erased.  Q_K is therefore block diagonal in the residue classes n mod K, and the
SPIDER state lives in the residue-0 block: it has support only on |mK>.  That is
exactly the statement that |S_K> is K-fold rotation symmetric -- the link to
rotation-symmetric bosonic codes.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.linalg import eigh
from scipy.special import gammaln


# ---------------------------------------------------------------------------
# Harmonic-oscillator wavefunctions at the origin
# ---------------------------------------------------------------------------
def _psi0(n: int) -> float:
    """psi_n(0) for the dimensionless oscillator (x in units sqrt(hbar/m omega)).

    psi_{2j}(0) = pi^{-1/4} (-1)^j sqrt((2j)!) / (2^j j!),  psi_{odd}(0) = 0.
    Computed in log space so that n ~ 1e5 is still fine.
    """
    if n % 2:
        return 0.0
    j = n // 2
    log_abs = -0.25 * math.log(math.pi) + 0.5 * gammaln(n + 1) - j * math.log(2.0) - gammaln(j + 1)
    return (-1) ** j * math.exp(log_abs)


def _dpsi0(n: int) -> float:
    """psi_n'(0).  Using d/dx = (a - a^dag)/sqrt(2):
    psi_n' = ( sqrt(n) psi_{n-1} - sqrt(n+1) psi_{n+1} ) / sqrt(2).
    Non-zero only for odd n."""
    if n % 2 == 0:
        return 0.0
    left = math.sqrt(n) * _psi0(n - 1)
    right = math.sqrt(n + 1) * _psi0(n + 1)
    return (left - right) / math.sqrt(2.0)


def pos_x_block(rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    """<m| pos(x) |n> for m in `rows`, n in `cols` (exact closed form).

    Obtained by integrating the Schrodinger equation from 0 to infinity:

        <m| pos(x) |n> = 1/2                                                  (m = n)
                       = [ psi_n(0) psi_m'(0) - psi_m(0) psi_n'(0) ] / (2 (m - n))   (m != n)

    which vanishes unless m - n is odd.  Only the requested block is built, so
    a residue class of a 10^4-dimensional truncation costs little memory.
    """
    rows = np.asarray(rows)
    cols = np.asarray(cols)
    psi_r = np.array([_psi0(int(n)) for n in rows])
    dpsi_r = np.array([_dpsi0(int(n)) for n in rows])
    psi_c = np.array([_psi0(int(n)) for n in cols])
    dpsi_c = np.array([_dpsi0(int(n)) for n in cols])
    m = rows[:, None]
    n = cols[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        P = (psi_c[None, :] * dpsi_r[:, None] - psi_r[:, None] * dpsi_c[None, :]) / (2.0 * (m - n))
    P[m == n] = 0.5
    return P


def pos_x_matrix(N: int) -> np.ndarray:
    """Full matrix of the projector pos(x) = Theta(x) on the first N Fock states."""
    idx = np.arange(N)
    return pos_x_block(idx, idx)


# ---------------------------------------------------------------------------
# Classical bound and quantum score operator
# ---------------------------------------------------------------------------
def classical_bound(K: int) -> float:
    """Largest score a classical uniformly precessing system can reach."""
    return 0.5 * (1.0 + 1.0 / K) if K % 2 else 0.5


def score_operator(K: int, N: int) -> np.ndarray:
    """Q_K = (1/K) sum_k U_k^dag pos(x) U_k in the Fock basis, truncated at N.

    Thanks to the Kronecker comb this is pos(x) masked to m = n (mod K)."""
    P = pos_x_matrix(N)
    idx = np.arange(N)
    mask = (idx[:, None] - idx[None, :]) % K == 0
    return P * mask


def residue_block(K: int, N: int, residue: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """The block of Q_K acting on Fock states n = residue (mod K).

    Returns (block_matrix, fock_indices)."""
    idx = np.arange(residue, N, K)
    return pos_x_block(idx, idx), idx


@dataclass
class SpiderState:
    """Result of maximising the precession score.

    Attributes
    ----------
    K         : order of the protocol / rotation symmetry
    N         : Fock truncation used
    score     : P_K^infty (largest eigenvalue of Q_K in the chosen block)
    gap       : distance to the second largest eigenvalue of the same block
    psi       : coefficients on |0>, |1>, ..., |N-1>   (real, normalised)
    residue   : Fock residue class mod K in which the state lives
    """

    K: int
    N: int
    score: float
    gap: float
    psi: np.ndarray
    residue: int = 0

    @property
    def grid_coefficients(self) -> np.ndarray:
        """s_m = <mK + residue | S_K>  (the Fock-grid coefficients)."""
        return self.psi[self.residue :: self.K]

    @property
    def mean_photon(self) -> float:
        return float(np.sum(np.arange(self.N) * self.psi**2))


def spider_state(K: int, N: int, residue: int = 0, sign_convention: bool = True) -> SpiderState:
    """Largest-eigenvalue eigenvector of Q_K restricted to Fock states n = residue (mod K).

    For residue = 0 this is the SPIDER state |S_K> of Zaw et al.  The global sign
    is fixed so that the coefficient of |residue> is positive, which reproduces
    the (-1)^{ceil(m/2)} pattern of the coefficients reported in the paper.
    """
    block, idx = residue_block(K, N, residue)
    vals, vecs = eigh(block)
    top = vecs[:, -1]
    if sign_convention and top[0] < 0:
        top = -top
    psi = np.zeros(N)
    psi[idx] = top
    return SpiderState(K=K, N=N, score=float(vals[-1]), gap=float(vals[-1] - vals[-2]), psi=psi, residue=residue)


def best_residue(K: int, N: int) -> int:
    """Which residue class mod K hosts the largest eigenvalue of Q_K (sanity check)."""
    scores = [eigh(residue_block(K, N, r)[0], eigvals_only=True)[-1] for r in range(K)]
    return int(np.argmax(scores))


def score_of_state(psi: np.ndarray, K: int) -> float:
    """<psi| Q_K |psi> for an arbitrary (complex) Fock-basis vector."""
    N = len(psi)
    Q = score_operator(K, N)
    return float(np.real(np.vdot(psi, Q @ psi)))


def sign_probabilities(psi: np.ndarray, K: int) -> np.ndarray:
    """Pr[x(t_k) > 0] for each of the K measurement times separately."""
    N = len(psi)
    P = pos_x_matrix(N)
    n = np.arange(N)
    out = []
    for k in range(K):
        phase = np.exp(-2j * np.pi * n * k / K)
        chi = phase * psi
        out.append(float(np.real(np.vdot(chi, P @ chi))))
    return np.array(out)


def convergence_scan(K: int, N_list, residue: int = 0) -> list[SpiderState]:
    """Spider state for a list of truncations (to study convergence)."""
    return [spider_state(K, int(N), residue) for N in N_list]


def fit_power_law_limit(N_list, scores):
    """Fit  score(N) = P_inf - c N^(-p)  and return (P_inf, c, p)."""
    from scipy.optimize import curve_fit

    N_arr = np.asarray(N_list, float)
    s_arr = np.asarray(scores, float)

    def f(N, Pinf, c, p):
        return Pinf - c * N ** (-p)

    p0 = [s_arr[-1] + 1e-3, 0.1, 0.5]
    popt, _ = curve_fit(f, N_arr, s_arr, p0=p0, maxfev=20000)
    return popt
