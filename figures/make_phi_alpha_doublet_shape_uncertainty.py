#!/usr/bin/env python3
"""Phi-alpha doublet-shape uncertainty band figure and CSV tables."""
from __future__ import annotations

import csv
import ctypes
import math
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
PHIALPHA = ROOT / "cf_phialpha_folded.csv"
OUT_PDF = ROOT / "figures" / "fig_phi_alpha_doublet_shape_uncertainty.pdf"
OUT_PNG = ROOT / "figures" / "fig_phi_alpha_doublet_shape_uncertainty.png"
OUT_CURVES = ROOT / "figures" / "phi_alpha_doublet_shape_uncertainty.csv"
OUT_SUMMARY = ROOT / "figures" / "phi_alpha_doublet_shape_uncertainty_summary.csv"

RADII_FM = (1.0, 3.0, 5.0)
SUMMARY_K = (2.0, 10.0, 50.0, 100.0)
BAND_SCENARIOS = ("qd_b040_fold", "qd_fold", "qd_b070_fold")

COLOR_Q = "#333333"
COLOR_QD = "#D55E00"
COLOR_BAND = "#D55E00"


def read_csv(path: Path) -> dict[str, list[float]]:
    with path.open(newline="") as f:
        lines = [
            line
            for line in f
            if line.strip() and not line.lstrip().startswith("#")
        ]
    rows = list(csv.DictReader(lines))
    out: dict[str, list[float]] = {h: [] for h in rows[0].keys()} if rows else {}
    for row in rows:
        for key, value in row.items():
            if value is None or str(value).strip() == "":
                out[key].append(float("nan"))
            else:
                out[key].append(float(value))
    return out


def interp_linear(xs: list[float], ys: list[float], x: float) -> float:
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]:
            t = (x - xs[i]) / (xs[i + 1] - xs[i])
            return ys[i] + t * (ys[i + 1] - ys[i])
    return ys[-1]


def kp_col(scenario: str, R: float) -> str:
    return f"KP_{scenario}_R{R:.2f}"


def ll_col(scenario: str, R: float) -> str:
    return f"LL_{scenario}_R{R:.2f}"


def band_limits_at_k(
    data: dict[str, list[float]], scenarios: tuple[str, ...], R: float, k: float
) -> tuple[float, float]:
    vals = [
        interp_linear(data["k_MeV"], data[kp_col(s, R)], k) for s in scenarios
    ]
    return min(vals), max(vals)


def curve_row(
    data: dict[str, list[float]], R: float, k: float
) -> dict[str, float]:
    c_q = interp_linear(data["k_MeV"], data[kp_col("q_fold", R)], k)
    c_b040 = interp_linear(data["k_MeV"], data[kp_col("qd_b040_fold", R)], k)
    c_qd = interp_linear(data["k_MeV"], data[kp_col("qd_fold", R)], k)
    c_b070 = interp_linear(data["k_MeV"], data[kp_col("qd_b070_fold", R)], k)
    band_min = min(c_b040, c_qd, c_b070)
    band_max = max(c_b040, c_qd, c_b070)
    return {
        "kstar_MeV_c": k,
        "R_fm": R,
        "C_KP_q_fold": c_q,
        "C_KP_qd_b040_fold": c_b040,
        "C_KP_qd_fold": c_qd,
        "C_KP_qd_b070_fold": c_b070,
        "C_KP_qd_band_min": band_min,
        "C_KP_qd_band_max": band_max,
        "C_LL_qd_fold": interp_linear(data["k_MeV"], data[ll_col("qd_fold", R)], k),
        "C_KP_qd_b040_nofold": interp_linear(
            data["k_MeV"], data[kp_col("qd_b040_nofold", R)], k
        ),
        "C_KP_qd_nofold": interp_linear(
            data["k_MeV"], data[kp_col("qd_nofold", R)], k
        ),
        "C_KP_qd_b070_nofold": interp_linear(
            data["k_MeV"], data[kp_col("qd_b070_nofold", R)], k
        ),
    }


def relative_band_width_percent(c_qd: float, band_min: float, band_max: float) -> float:
    denom = abs(c_qd - 1.0)
    if denom < 1e-6:
        return float("nan")
    return 100.0 * (band_max - band_min) / denom


def write_curve_csv(data: dict[str, list[float]], path: Path) -> None:
    fieldnames = [
        "kstar_MeV_c",
        "R_fm",
        "C_KP_q_fold",
        "C_KP_qd_b040_fold",
        "C_KP_qd_fold",
        "C_KP_qd_b070_fold",
        "C_KP_qd_band_min",
        "C_KP_qd_band_max",
        "C_LL_qd_fold",
        "C_KP_qd_b040_nofold",
        "C_KP_qd_nofold",
        "C_KP_qd_b070_nofold",
    ]
    rows: list[dict[str, float]] = []
    for R in RADII_FM:
        for k in data["k_MeV"]:
            rows.append(curve_row(data, R, k))
    path.parent.mkdir(exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: f"{row[key]:.8g}" for key in fieldnames})


def write_summary_csv(data: dict[str, list[float]], path: Path) -> list[dict[str, float]]:
    fieldnames = [
        "R_fm",
        "kstar_MeV_c",
        "C_q_fold",
        "C_qd_b040_fold",
        "C_qd_fold",
        "C_qd_b070_fold",
        "C_qd_band_min",
        "C_qd_band_max",
        "Delta_C_qd_band",
        "Delta_C_qd_minus_q",
        "relative_band_width_percent",
    ]
    rows: list[dict[str, float]] = []
    for R in RADII_FM:
        for k in SUMMARY_K:
            row = curve_row(data, R, k)
            delta_band = row["C_KP_qd_band_max"] - row["C_KP_qd_band_min"]
            rows.append(
                {
                    "R_fm": R,
                    "kstar_MeV_c": k,
                    "C_q_fold": row["C_KP_q_fold"],
                    "C_qd_b040_fold": row["C_KP_qd_b040_fold"],
                    "C_qd_fold": row["C_KP_qd_fold"],
                    "C_qd_b070_fold": row["C_KP_qd_b070_fold"],
                    "C_qd_band_min": row["C_KP_qd_band_min"],
                    "C_qd_band_max": row["C_KP_qd_band_max"],
                    "Delta_C_qd_band": delta_band,
                    "Delta_C_qd_minus_q": row["C_KP_qd_fold"] - row["C_KP_q_fold"],
                    "relative_band_width_percent": relative_band_width_percent(
                        row["C_KP_qd_fold"],
                        row["C_KP_qd_band_min"],
                        row["C_KP_qd_band_max"],
                    ),
                }
            )
    path.parent.mkdir(exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            out = {}
            for key in fieldnames:
                val = row[key]
                out[key] = "" if math.isnan(val) else f"{val:.8g}"
            writer.writerow(out)
    return rows


def add_caption(fig, text: str) -> None:
    fig.text(0.5, 0.02, text, ha="center", va="bottom", fontsize=8.5, wrap=True)


def make_figure(data: dict[str, list[float]], pdf=None) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.6), sharey=True)
    for ax, R in zip(axes, RADII_FM):
        k = data["k_MeV"]
        band_cols = [kp_col(s, R) for s in BAND_SCENARIOS]
        lo, hi = [], []
        for vals in zip(*(data[c] for c in band_cols)):
            lo.append(min(vals))
            hi.append(max(vals))
        ax.fill_between(k, lo, hi, color=COLOR_BAND, alpha=0.22, lw=0, label=r"$b_d$=0.40–0.70 fm band")
        ax.plot(k, data[kp_col("q_fold", R)], color=COLOR_Q, lw=2.2, label="KP quartet only / fold")
        ax.plot(k, data[kp_col("qd_fold", R)], color=COLOR_QD, lw=2.4, label=r"KP q+d / fold ($b_d$=0.55 fm)")
        ax.plot(
            k,
            data[ll_col("qd_fold", R)],
            color=COLOR_QD,
            lw=1.0,
            ls="--",
            alpha=0.45,
            label="LL q+d / fold (diagnostic)",
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
        r"$\phi\alpha$: Gaussian doublet-shape uncertainty ($f_0 \approx -1.54$ fm)",
        y=0.98,
    )
    add_caption(
        fig,
        "Band spans KP folded q+d curves for $b_d=0.40$, $0.55$, $0.70$ fm (V0 tuned to "
        "the same scattering length). Not a unique Chizzali doublet reproduction; "
        "$d_0=0.39$ fm is not forced. Quartet baseline: HAL fit-B / fold.",
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    OUT_PDF.parent.mkdir(exist_ok=True)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=180)
    if pdf is not None:
        pdf.savefig(fig)
    plt.close(fig)


def main() -> list[dict[str, float]]:
    if not PHIALPHA.exists():
        raise SystemExit(f"Missing {PHIALPHA.relative_to(ROOT)}")
    data = read_csv(PHIALPHA)
    write_curve_csv(data, OUT_CURVES)
    summary_rows = write_summary_csv(data, OUT_SUMMARY)
    make_figure(data)
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")
    print(f"wrote {OUT_PNG.relative_to(ROOT)}")
    print(f"wrote {OUT_CURVES.relative_to(ROOT)}")
    print(f"wrote {OUT_SUMMARY.relative_to(ROOT)}")
    return summary_rows


if __name__ == "__main__":
    main()
