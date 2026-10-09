"""
spiderqec
=========

Reusable code behind the *spider code* repository.

Modules
-------
precession      Tsirelson / Zaw et al. precession protocol, the score operator
                Q_K and the SPIDER state (its top eigenvector).
rotation_codes  Rotation-symmetric bosonic codes: codewords from a primitive,
                Fock-grid structure, mean photon number, modular phase
                uncertainty, Wigner functions.
knill_laflamme  Knill-Laflamme test of a code against a set of bosonic errors.
transmon        Cavity-transmon Hamiltonians, dressed-state labelling, fitting
                of drive-induced self-/cross-Kerr, and the CROT gate.
plotting        Small helpers for consistent figures.
"""

from . import precession, rotation_codes, knill_laflamme, transmon, plotting  # noqa: F401

__all__ = ["precession", "rotation_codes", "knill_laflamme", "transmon", "plotting"]
__version__ = "0.1.0"
