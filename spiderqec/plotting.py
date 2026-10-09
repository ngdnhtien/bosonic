"""Small plotting helpers so that all figures in the repository look alike."""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

FIG_DIR = Path(os.environ.get("SPIDERQEC_FIG_DIR", Path(__file__).resolve().parents[1] / "figures"))

# A quiet, colour-blind-friendly cycle (Okabe-Ito).
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]


def use_style():
    mpl.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": 300,
            "font.size": 11,
            "axes.titlesize": 11,
            "axes.labelsize": 11,
            "legend.fontsize": 9,
            "legend.frameon": False,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.prop_cycle": mpl.cycler(color=COLORS),
            "lines.linewidth": 1.6,
        }
    )


def savefig(fig, name: str, subdir: str | None = None, also_report: bool = True):
    """Save `fig` as PNG (and PDF) under figures/[subdir]/name and, optionally, a
    copy of the PDF under report/figures for the LaTeX report."""
    out = FIG_DIR / subdir if subdir else FIG_DIR
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.png", bbox_inches="tight")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    if also_report:
        rep = FIG_DIR.parent / "report" / "figures"
        rep.mkdir(parents=True, exist_ok=True)
        fig.savefig(rep / f"{name}.pdf", bbox_inches="tight")
    return out / f"{name}.png"


def plot_wigner(ax, xvec, W, title: str | None = None, vmax: float | None = None, contour: bool = False):
    """Symmetric red/blue Wigner plot (optionally with a zero contour)."""
    vmax = vmax or float(np.max(np.abs(W)))
    im = ax.pcolormesh(xvec, xvec, W, cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto", rasterized=True)
    if contour:
        ax.contour(xvec, xvec, W, levels=[0.0], colors="k", linewidths=0.4)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$p$")
    if title:
        ax.set_title(title)
    return im


def stem_coefficients(ax, coeffs, K: int, label: str | None = None, color=None, nmax: int | None = None):
    """Stem plot of Fock-grid coefficients s_m against n = mK."""
    n = np.arange(len(coeffs)) * K
    if nmax is not None:
        keep = n <= nmax
        n, coeffs = n[keep], np.asarray(coeffs)[keep]
    ax.axhline(0, color="0.6", lw=0.6)
    ax.stem(n, coeffs, linefmt=(color or COLORS[0]), markerfmt="o", basefmt=" ", label=label)
    ax.set_xlabel(r"Fock number $n$")
