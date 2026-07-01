#!/usr/bin/env python3
"""Phi-alpha doublet potential-shape figure from cf-calc potential grids."""
from __future__ import annotations

import csv
import ctypes
import os
import subprocess
import sys
from pathlib import Path

for site in [
    Path.home() / ".local/lib/python3.8/site-packages",
    Path("/cvmfs/star.sdcc.bnl.gov/star-spack/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/py-numpy-1.22.4-3dkbckoarnjk7exoy6cnapgrgaqhjxbv/lib/python3.8/site-packages"),
    Path("/cvmfs/star.sdcc.bnl.gov/star-spack/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/py-pyparsing-3.0.6-slw2fhcqm64fxvaiwcje4riz45fxzstn/lib/python3.8/site-packages"),
]:
    if site.exists():
        sys.path.insert(0, str(site))

for lib in [
    Path("/cvmfs/singularity.opensciencegrid.org/.images/98/227808bc057d4dabb58d49a8da9a456d6083b0060c0ee6d56433f82dd0297a/cvmfs2/gm2.opensciencegrid.org/prod/external/gcc/v6_4_0/Linux64bit+2.6-2.12/lib64/libgfortran.so.3"),
    Path("/usr/lib64/libgfortran.so.3"),
]:
    if lib.exists():
        try:
            ctypes.CDLL(str(lib), mode=ctypes.RTLD_GLOBAL)
        except OSError:
            pass
        break

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
BUILD = Path(os.environ.get("CF_FITTER_BUILD", str(ROOT / "build")))
CF_CALC = BUILD / "cf-calc"

BARE_CSV = ROOT / "figures" / "phi_alpha_doublet_bare_potential_shapes.csv"
FOLDED_CSV = ROOT / "figures" / "phi_alpha_doublet_folded_potential_shapes.csv"
SUMMARY_CSV = ROOT / "figures" / "phi_alpha_doublet_potential_shapes_summary.csv"
OUT_PDF = ROOT / "figures" / "fig_phi_alpha_doublet_potential_shapes.pdf"
OUT_PNG = ROOT / "figures" / "fig_phi_alpha_doublet_potential_shapes.png"

COLOR_Q = "#333333"
COLOR_QD = "#D55E00"
COLOR_BAND = "#D55E00"
COLORS_BD = {
    "b040": "#E69F00",
    "b055": "#D55E00",
    "b070": "#CC6677",
}


def read_csv(path: Path) -> dict[str, list[float]]:
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    out: dict[str, list[float]] = {h: [] for h in rows[0].keys()}
    for row in rows:
        for key, value in row.items():
            out[key].append(float(value))
    return out


def ensure_potential_csvs() -> None:
    if BARE_CSV.exists() and FOLDED_CSV.exists() and SUMMARY_CSV.exists():
        return
    if not CF_CALC.exists():
        raise SystemExit(f"Missing {CF_CALC}; build cf-calc or generate CSVs first.")
    subprocess.run(
        [str(CF_CALC), "phi-alpha-potentials"],
        cwd=ROOT,
        check=True,
    )


def band_limits(*series: list[float]) -> tuple[list[float], list[float]]:
    lo, hi = [], []
    for vals in zip(*series):
        lo.append(min(vals))
        hi.append(max(vals))
    return lo, hi


def add_caption(fig, text: str) -> None:
    fig.text(0.5, 0.02, text, ha="center", va="bottom", fontsize=8.5, wrap=True)


def make_figure(pdf=None) -> None:
    ensure_potential_csvs()
    bare = read_csv(BARE_CSV)
    folded = read_csv(FOLDED_CSV)

    fig, (ax_bare, ax_fold) = plt.subplots(1, 2, figsize=(12.6, 4.8))

    r = bare["r_fm"]
    ax_bare.plot(r, bare["V_bare_q_MeV"], color=COLOR_Q, lw=2.2, label=r"$V_q(r)$ quartet only")
    ax_bare.plot(
        r,
        bare["V_bare_qd_b040_MeV"],
        color=COLORS_BD["b040"],
        lw=1.6,
        ls="--",
        label=r"$V_c(r)$, $b_d$=0.40 fm",
    )
    ax_bare.plot(
        r,
        bare["V_bare_qd_b055_MeV"],
        color=COLORS_BD["b055"],
        lw=2.0,
        label=r"$V_c(r)$, $b_d$=0.55 fm",
    )
    ax_bare.plot(
        r,
        bare["V_bare_qd_b070_MeV"],
        color=COLORS_BD["b070"],
        lw=1.6,
        ls="--",
        label=r"$V_c(r)$, $b_d$=0.70 fm",
    )
    lo, hi = band_limits(
        bare["V_bare_qd_b040_MeV"],
        bare["V_bare_qd_b055_MeV"],
        bare["V_bare_qd_b070_MeV"],
    )
    ax_bare.fill_between(r, lo, hi, color=COLOR_BAND, alpha=0.15, lw=0)
    ax_bare.set_xlim(0, 3.0)
    ax_bare.set_xlabel(r"$r$ [fm]")
    ax_bare.set_ylabel(r"$V(r)$ [MeV]")
    ax_bare.set_title(r"Bare $\phi$-$N$ central potential")
    ax_bare.grid(alpha=0.18, lw=0.5)
    ax_bare.legend(fontsize=7.0, loc="lower right", frameon=False)

    R = folded["R_fm"]
    ax_fold.plot(R, folded["V_fold_q_MeV"], color=COLOR_Q, lw=2.2, label=r"$V_{\phi\alpha}^q(R)$ quartet / fold")
    ax_fold.plot(
        R,
        folded["V_fold_qd_b040_MeV"],
        color=COLORS_BD["b040"],
        lw=1.6,
        ls="--",
        label=r"$b_d$=0.40 fm",
    )
    ax_fold.plot(
        R,
        folded["V_fold_qd_b055_MeV"],
        color=COLORS_BD["b055"],
        lw=2.0,
        label=r"$b_d$=0.55 fm",
    )
    ax_fold.plot(
        R,
        folded["V_fold_qd_b070_MeV"],
        color=COLORS_BD["b070"],
        lw=1.6,
        ls="--",
        label=r"$b_d$=0.70 fm",
    )
    flo, fhi = band_limits(
        folded["V_fold_qd_b040_MeV"],
        folded["V_fold_qd_b055_MeV"],
        folded["V_fold_qd_b070_MeV"],
    )
    ax_fold.fill_between(R, flo, fhi, color=COLOR_BAND, alpha=0.22, lw=0, label=r"$b_d$ band")
    ax_fold.set_xlim(0, 6.0)
    ax_fold.set_xlabel(r"$R$ [fm]")
    ax_fold.set_ylabel(r"$V_{\phi\alpha}(R)$ [MeV]")
    ax_fold.set_title(r"Folded $\phi\alpha$ potential (HAL fit-B quartet)")
    ax_fold.grid(alpha=0.18, lw=0.5)
    ax_fold.legend(fontsize=7.0, loc="lower right", frameon=False)

    fig.suptitle(
        r"$\phi\alpha$ doublet-shape family: $V_c=(2/3)V_q+(1/3)V_d(b_d)$, $f_0\approx-1.54$ fm",
        y=0.98,
    )
    add_caption(
        fig,
        r"Central q+d curves use the existing Gaussian doublet family ($b_d=0.40,0.55,0.70$ fm); "
        r"not a unique Chizzali doublet reproduction. Folding uses the same $\alpha$ density as "
        r"q_fold / qd_fold scenarios (analytic Gaussian fold).",
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    OUT_PDF.parent.mkdir(exist_ok=True)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=180)
    if pdf is not None:
        pdf.savefig(fig)
    plt.close(fig)


def main() -> None:
    make_figure()
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")
    print(f"wrote {OUT_PNG.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
