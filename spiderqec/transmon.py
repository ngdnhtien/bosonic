"""
Cavity-transmon models, drive-induced Kerr nonlinearities and the CROT gate.

Conventions
-----------
* All frequencies are *ordinary* frequencies in GHz (not angular); a Kerr
  coefficient K in GHz is quoted as K*1e3 MHz or K*1e6 kHz.
* Transmon (Koch et al. 2007) in the Duffing approximation,
      H_c = omega_c c^dag c - (alpha/12) (c + c^dag)^4,
  with omega_c = omega_10 + alpha, E_C = alpha, E_J = omega_c^2/(8 alpha).
* A classical drive  -2 e n V_d cos(omega_d t)  on the transmon charge.  In the
  frame rotating at omega_d and after the rotating-wave approximation (RWA) the
  one-cavity Hamiltonian is (Zhang et al., PRA 105, 022423 (2022))

      H_RWA = -delta_da a^dag a - delta_dc c^dag c - (alpha/2) c^dag c (c^dag c + 1)
              + Omega_d c^dag + Omega_d^* c + g_a a c^dag + g_a^* a^dag c,

  delta_dx = omega_d - omega_x.  The dimensionless knobs used throughout the
  legacy notebooks are

      r1 = delta_a/alpha  (cavity-transmon detuning),
      r2 = g_a/delta_a     (dispersive parameter, kept << 1),
      r3 = delta_d/alpha   (drive-transmon detuning),
      r4 = Omega_d/delta_d (drive amplitude);   plots use r4^2 = |Omega_d/delta_d|^2.

* The cavity inherits a nonlinearity from the transmon.  With the transmon in
  its (dressed) ground state the dressed cavity energies are expanded as

      E(N) = E_0 + delta N + (K/2) N^2 + (beta/6) N^3 + (sigma/24) N^4 + ...

  K is the self-Kerr.  For two cavities A, B sharing a transmon,

      E(N_A, N_B) = ... + chi_AB N_A N_B + ...

  and chi_AB is the cross-Kerr, which generates  CROT(phi) = exp(i phi n_A n_B)
  after a time t = phi / (2 pi chi_AB).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import qutip as qt
from scipy.linalg import eigh


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
@dataclass
class TransmonParams:
    omega10: float = 4.936  # GHz, transmon 0-1 frequency
    alpha: float = 0.168  # GHz, anharmonicity (= E_C)

    @property
    def omega_c(self):
        return self.omega10 + self.alpha

    @property
    def EC(self):
        return self.alpha

    @property
    def EJ(self):
        return self.omega_c**2 / (8 * self.alpha)


@dataclass
class OneCavityParams:
    """Dimensionless knobs of the one-cavity-one-transmon model."""

    r1: float = 9.64  # delta_a / alpha
    r2: float = 0.064  # g_a / delta_a
    r3: float = 0.08  # delta_d / alpha
    r4: float = 0.0  # Omega_d / delta_d
    tp: TransmonParams = None

    def __post_init__(self):
        if self.tp is None:
            self.tp = TransmonParams()

    @property
    def delta_a(self):
        return self.r1 * self.tp.alpha

    @property
    def g_a(self):
        return self.r2 * self.delta_a

    @property
    def delta_d(self):
        return self.r3 * self.tp.alpha

    @property
    def Omega_d(self):
        return self.r4 * self.delta_d

    @property
    def omega_a(self):
        return self.tp.omega10 + self.delta_a

    @property
    def omega_d(self):
        return self.tp.omega10 + self.delta_d


# ---------------------------------------------------------------------------
# Operators and Hamiltonians
# ---------------------------------------------------------------------------
def ops_two_mode(dim_a: int, dim_c: int):
    a = qt.tensor(qt.destroy(dim_a), qt.qeye(dim_c))
    c = qt.tensor(qt.qeye(dim_a), qt.destroy(dim_c))
    return a, c


def h_rwa_one_cavity(dim_a: int, dim_c: int, p: OneCavityParams, coupling: bool = True, drive: bool = True):
    """Drive-frame RWA Hamiltonian of one cavity + driven transmon (GHz)."""
    a, c = ops_two_mode(dim_a, dim_c)
    alpha = p.tp.alpha
    delta_da = p.omega_d - p.omega_a
    delta_dc = p.omega_d - p.tp.omega_c
    nc = c.dag() * c
    H = -delta_da * a.dag() * a - delta_dc * nc - (alpha / 2) * nc * (nc + 1)
    if drive:
        H += p.Omega_d * (c.dag() + c)
    if coupling:
        g = p.g_a
        H += g * (a * c.dag() + a.dag() * c)
    return H


def h_lab_one_cavity(dim_a: int, dim_c: int, p: OneCavityParams, duffing: bool = True):
    """Lab-frame static Hamiltonian (no drive) with either the Duffing (quartic)
    transmon or the full cos(phi) transmon.  Used for cross-checks."""
    a, c = ops_two_mode(dim_a, dim_c)
    tp = p.tp
    if duffing:
        H_anc = tp.omega_c * c.dag() * c - (tp.alpha / 12) * (c + c.dag()) ** 4
    else:
        eta = (8 * tp.EC / tp.EJ) ** 0.25
        n_op = -1j * (c - c.dag()) / (np.sqrt(2) * eta)
        phi_op = eta * (c + c.dag()) / np.sqrt(2)
        H_anc = 4 * tp.EC * n_op**2 - tp.EJ * phi_op.cosm()
    H_int = 1j * p.g_a * (c - c.dag()) * (a + a.dag())
    return p.omega_a * a.dag() * a + H_anc + H_int


# ---------------------------------------------------------------------------
# Dressed-state labelling
# ---------------------------------------------------------------------------
def _eig(H):
    """Dense eigendecomposition; returns (energies, column eigenvectors)."""
    M = H.full() if isinstance(H, qt.Qobj) else np.asarray(H)
    return eigh(M)


def label_by_overlap(ref_vectors: np.ndarray, ref_labels: list, vectors: np.ndarray):
    """Label each column of `vectors` by the label of the reference column with
    largest overlap.  Returns (labels, max_overlap_probabilities)."""
    ov = np.abs(ref_vectors.conj().T @ vectors) ** 2  # (n_ref, n_vec)
    best = np.argmax(ov, axis=0)
    return [ref_labels[b] for b in best], ov[best, np.arange(vectors.shape[1])]


def product_labels(dims):
    """Labels of the product Fock basis in QuTiP's tensor ordering."""
    return [tuple(idx) for idx in np.ndindex(*dims)]


def dressed_ladder(H_uncoupled, H_coupled, dims, ancilla_axis: int = -1, warn_threshold: float = 0.5):
    """Energies of the coupled Hamiltonian labelled by bare occupation numbers.

    1. Diagonalise H_uncoupled (product form) and label its eigenstates with the
       bare Fock labels (max overlap).  This handles a driven transmon, whose
       bare eigenstates are not Fock states.
    2. Diagonalise H_coupled and label its eigenstates by max overlap with the
       labelled uncoupled eigenstates (adiabatic continuation).

    Returns a dict  label -> (energy, state_vector, overlap)  for the states in
    which the ancilla is in its ground state.
    """
    E0, V0 = _eig(H_uncoupled)
    labels0, _ = label_by_overlap(np.eye(V0.shape[0]), product_labels(dims), V0)
    E1, V1 = _eig(H_coupled)
    labels1, ov1 = label_by_overlap(V0, labels0, V1)
    out = {}
    for lab, e, v, o in zip(labels1, E1, V1.T, ov1):
        if lab[ancilla_axis] != 0:
            continue
        if lab in out and out[lab][2] >= o:
            continue  # keep the best-overlap representative
        out[lab] = (float(e), v, float(o))
    return out


# ---------------------------------------------------------------------------
# Kerr extraction
# ---------------------------------------------------------------------------
def fit_nonlinearities(N: np.ndarray, E: np.ndarray, order: int = 4) -> dict:
    """Fit E(N) - E(0) = delta N + (K/2) N^2 + (beta/6) N^3 + (sigma/24) N^4.

    Returns the physical coefficients {delta, K, beta, sigma} (GHz)."""
    N = np.asarray(N, float)
    E = np.asarray(E, float) - E[0]
    A = np.vstack([N**k for k in range(1, order + 1)]).T
    c, *_ = np.linalg.lstsq(A, E, rcond=None)
    names = ["delta", "K", "beta", "sigma", "c5", "c6"]
    facts = [1, 2, 6, 24, 120, 720]
    return {names[k]: float(c[k] * facts[k]) for k in range(order)}


def self_kerr(dim_a: int, dim_c: int, p: OneCavityParams, n_fit: int | None = None, order: int = 4) -> dict:
    """Self-Kerr (and higher nonlinearities) of the cavity for parameters p."""
    H0 = h_rwa_one_cavity(dim_a, dim_c, p, coupling=False, drive=True)
    H1 = h_rwa_one_cavity(dim_a, dim_c, p, coupling=True, drive=True)
    ladder = dressed_ladder(H0, H1, (dim_a, dim_c))
    Ns = sorted(lab[0] for lab in ladder)
    if n_fit is None:
        n_fit = max(4, int(0.6 * dim_a))  # stay away from the truncation edge
    Ns = [n for n in Ns if n < n_fit]
    E = np.array([ladder[(n, 0)][0] for n in Ns])
    res = fit_nonlinearities(np.array(Ns), E, order=order)
    res["min_overlap"] = float(min(ladder[(n, 0)][2] for n in Ns))
    return res


def self_kerr_vs_drive(dim_a, dim_c, r4sq_list, base: OneCavityParams | None = None, **kw):
    """Sweep |Omega_d/delta_d|^2 and return arrays of K, beta, sigma (GHz)."""
    base = base or OneCavityParams()
    rows = []
    for r4sq in r4sq_list:
        p = OneCavityParams(r1=base.r1, r2=base.r2, r3=base.r3, r4=np.sqrt(r4sq), tp=base.tp)
        rows.append(self_kerr(dim_a, dim_c, p, **kw))
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


# ---------------------------------------------------------------------------
# Two cavities sharing one transmon: cross-Kerr
# ---------------------------------------------------------------------------
@dataclass
class TwoCavityParams:
    r1a: float = 9.64  # delta_a/alpha
    r2a: float = 0.064  # g_a/delta_a
    r1b: float = -9.0  # delta_b/alpha
    r2b: float = 0.05  # g_b/delta_b
    r3: float = 0.26  # delta_d/alpha
    r4: float = 0.0  # Omega_d/delta_d
    tp: TransmonParams = None

    def __post_init__(self):
        if self.tp is None:
            self.tp = TransmonParams()

    @property
    def delta_a(self):
        return self.r1a * self.tp.alpha

    @property
    def delta_b(self):
        return self.r1b * self.tp.alpha

    @property
    def g_a(self):
        return self.r2a * self.delta_a

    @property
    def g_b(self):
        return self.r2b * self.delta_b

    @property
    def delta_d(self):
        return self.r3 * self.tp.alpha

    @property
    def Omega_d(self):
        return self.r4 * self.delta_d

    @property
    def omega_a(self):
        return self.tp.omega10 + self.delta_a

    @property
    def omega_b(self):
        return self.tp.omega10 + self.delta_b

    @property
    def omega_d(self):
        return self.tp.omega10 + self.delta_d


def ops_three_mode(dim_a, dim_b, dim_c):
    a = qt.tensor(qt.destroy(dim_a), qt.qeye(dim_b), qt.qeye(dim_c))
    b = qt.tensor(qt.qeye(dim_a), qt.destroy(dim_b), qt.qeye(dim_c))
    c = qt.tensor(qt.qeye(dim_a), qt.qeye(dim_b), qt.destroy(dim_c))
    return a, b, c


def h_rwa_two_cavity(dim_a, dim_b, dim_c, p: TwoCavityParams, coupling=True, drive=True):
    """Drive-frame RWA Hamiltonian of cavities A, B both coupled to one driven transmon."""
    a, b, c = ops_three_mode(dim_a, dim_b, dim_c)
    alpha = p.tp.alpha
    nc = c.dag() * c
    H = (
        -(p.omega_d - p.omega_a) * a.dag() * a
        - (p.omega_d - p.omega_b) * b.dag() * b
        - (p.omega_d - p.tp.omega_c) * nc
        - (alpha / 2) * nc * (nc + 1)
    )
    if drive:
        H += p.Omega_d * (c.dag() + c)
    if coupling:
        H += p.g_a * (a * c.dag() + a.dag() * c) + p.g_b * (b * c.dag() + b.dag() * c)
    return H


def kerr_matrix(dim_a, dim_b, dim_c, p: TwoCavityParams, nmax: int = 2) -> dict:
    """Self- and cross-Kerr from finite differences of the dressed ladder E(N_A, N_B).

        K_A    = E(2,0) - 2E(1,0) + E(0,0)
        K_B    = E(0,2) - 2E(0,1) + E(0,0)
        chi_AB = E(1,1) - E(1,0) - E(0,1) + E(0,0)
    """
    H0 = h_rwa_two_cavity(dim_a, dim_b, dim_c, p, coupling=False)
    H1 = h_rwa_two_cavity(dim_a, dim_b, dim_c, p, coupling=True)
    lad = dressed_ladder(H0, H1, (dim_a, dim_b, dim_c))
    E = lambda i, j: lad[(i, j, 0)][0]
    out = dict(
        K_A=E(2, 0) - 2 * E(1, 0) + E(0, 0),
        K_B=E(0, 2) - 2 * E(0, 1) + E(0, 0),
        chi_AB=E(1, 1) - E(1, 0) - E(0, 1) + E(0, 0),
        min_overlap=min(lad[(i, j, 0)][2] for i in range(nmax + 1) for j in range(nmax + 1)),
    )
    return out


def kerr_matrix_vs_drive(dim_a, dim_b, dim_c, r4sq_list, base: TwoCavityParams | None = None):
    base = base or TwoCavityParams()
    rows = []
    for r4sq in r4sq_list:
        p = TwoCavityParams(base.r1a, base.r2a, base.r1b, base.r2b, base.r3, np.sqrt(r4sq), base.tp)
        rows.append(kerr_matrix(dim_a, dim_b, dim_c, p))
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


# ---------------------------------------------------------------------------
# CROT gate
# ---------------------------------------------------------------------------
def crot(phi: float, dim_a: int, dim_b: int) -> np.ndarray:
    """Ideal CROT(phi) = exp(i phi n_A n_B) as a dense matrix on dim_a x dim_b."""
    na = np.arange(dim_a)[:, None]
    nb = np.arange(dim_b)[None, :]
    return np.diag(np.exp(1j * phi * (na * nb)).ravel())


def kerr_unitary(t: float, chi: float, K_A: float, K_B: float, dim_a: int, dim_b: int) -> np.ndarray:
    """U = exp(-2 pi i t [chi n_A n_B + (K_A/2) n_A^2 + (K_B/2) n_B^2])  (t in ns, rates in GHz).

    Linear terms (delta n) are dropped: they are absorbed into the rotating frame
    of each cavity (equivalently undone by a Z_L rotation)."""
    na = np.arange(dim_a)[:, None]
    nb = np.arange(dim_b)[None, :]
    E = chi * na * nb + 0.5 * K_A * na**2 + 0.5 * K_B * nb**2
    return np.diag(np.exp(-2j * np.pi * t * E).ravel())


def crot_time(phi: float, chi: float) -> float:
    """Time (ns) for the cross-Kerr chi (GHz) to accumulate CROT(phi): phi = -2 pi chi t."""
    return -phi / (2 * np.pi * chi)


def average_gate_fidelity(U_ideal: np.ndarray, U_actual: np.ndarray, basis: np.ndarray) -> float:
    """Average fidelity of U_actual against U_ideal on the subspace spanned by the
    orthonormal columns of `basis` (leakage counts as error):

        F = ( |Tr M|^2 + Tr(M^dag M) ) / ( d (d+1) ),   M = B^dag U_ideal^dag U_actual B.
    """
    M = basis.conj().T @ U_ideal.conj().T @ U_actual @ basis
    d = basis.shape[1]
    return float((abs(np.trace(M)) ** 2 + np.real(np.trace(M.conj().T @ M))) / (d * (d + 1)))


def logical_two_qubit_basis(zero_a, one_a, zero_b, one_b) -> np.ndarray:
    """Columns |00>,|01>,|10>,|11> of two bosonic codes (Kronecker order A x B)."""
    cols = [np.kron(x, y) for x in (zero_a, one_a) for y in (zero_b, one_b)]
    return np.array(cols).T
