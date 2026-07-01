#!/usr/bin/env python3
"""Standalone phi-alpha Chizzali TPE quartet comparison figure."""
from __future__ import annotations

import csv
import ctypes
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

CSV = ROOT / "validation" / "phi_alpha" / "generated_q_chizzali_tpe.csv"
OUT_PDF = ROOT / "figures" / "fig_phi_alpha_chizzali_tpe_quartet_comparison.pdf"
OUT_PNG = ROOT / "figures" / "fig_phi_alpha_chizzali_tpe_quartet_comparison.png"


def read_csv(path: Path) -> dict[str, list[float]]:
    with path.open(newline="") as f:
        lines = [line for line in f if line.strip() and not line.lstrip().startswith("#")]
    rows = list(csv.DictReader(lines))
    out: dict[str, list[float]] = {h: [] for h in rows[0].keys()} if rows else {}
    for row in rows:
        for key, value in row.items():
            out[key].append(float(value))
    return out


def add_caption(fig, text: str):
    fig.text(0.5, 0.02, text, ha="center", va="bottom", fontsize=8.5, wrap=True)


def make_figure(pdf=None):
    if not CSV.exists():
        raise SystemExit(
            f"Missing {CSV.relative_to(ROOT)}. "
            "Run scripts/regenerate_phialpha_chizzali_tpe.sh first."
        )
    data = read_csv(CSV)
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.6), sharey=True)
    for ax, R in zip(axes, [1.0, 3.0, 5.0]):
        rs = f"R{R:.2f}"
        ax.plot(
            data["k_MeV"],
            data[f"KP_q_fold_{rs}"],
            color="#0072B2",
            lw=2.2,
            label="HAL fit-B q/fold KP",
        )
        ax.plot(
            data["k_MeV"],
            data[f"KP_q_chizzali_tpe_fold_{rs}"],
            color="#D55E00",
            lw=2.2,
            label="Chizzali t/a=12 TPE q/fold KP",
        )
        ax.plot(
            data["k_MeV"],
            data[f"KP_q_chizzali_tpe_nofold_{rs}"],
            color="#D55E00",
            lw=1.2,
            ls=":",
            alpha=0.75,
            label="Chizzali TPE q/no-fold KP",
        )
        ax.plot(
            data["k_MeV"],
            data[f"LL_q_chizzali_tpe_fold_{rs}"],
            color="#D55E00",
            lw=1.1,
            ls="--",
            alpha=0.8,
            label="Chizzali TPE q/fold LL",
        )
        ax.axhline(1.0, color="0.55", lw=0.8)
        ax.set_xlim(0, 200)
        ax.set_xlabel(r"$k^*$ [MeV/$c$]")
        ax.set_title(rf"$R={R:g}$ fm")
        ax.grid(alpha=0.18, lw=0.5)
    axes[0].set_ylabel(r"$C(k^*)$")
    axes[0].set_ylim(0.15, 1.15)
    axes[0].legend(fontsize=7.0, loc="lower right", frameon=False)
    fig.suptitle(
        r"$\phi\alpha$: HAL fit-B vs Chizzali $t/a{=}12$ TPE quartet folding",
        y=0.98,
    )
    add_caption(
        fig,
        "Solid: KP. Orange dashed: LL for Chizzali TPE folded cell. "
        "Numeric folding of the bare Chizzali TPE quartet; no doublet.",
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    OUT_PDF.parent.mkdir(exist_ok=True)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=180)
    if pdf is not None:
        pdf.savefig(fig)
    plt.close(fig)


def main():
    make_figure()
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")
    print(f"wrote {OUT_PNG.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
