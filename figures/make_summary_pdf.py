#!/usr/bin/env python3
"""Build the cf-fitter multi-page summary PDF from in-repository CSV outputs."""
from __future__ import annotations

import csv
import ctypes
import math
import sys
from pathlib import Path

# STAR/CVMFS environments can put an old Python-2 pyparsing ahead of the
# Python-3 stack, and NumPy may need libgfortran.so.3. Keep these fallbacks
# local to plotting so the PDF can be regenerated with one script invocation.
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
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "figures"))
OUT = ROOT / "figures" / "cf_fitter_summary_revised.pdf"
PHIP_REF = ROOT / "validation" / "phi_p"
PHIP_MODELS = PHIP_REF / "generated_phi_proton_spinavg_R1p08.csv"
PHIP_GAUSS = PHIP_REF / "generated_ALICEgauss_KP.csv"
PHIP_HAL_KP = PHIP_REF / "generated_HALquartetFitB_KP_Rband.csv"
PHIP_CHIZZALI_KP = PHIP_REF / "generated_HALquartetChizzaliTa12TPE_KP_Rband.csv"
PHIP_CHIZZALI_RSCAN = PHIP_REF / "generated_HALquartetChizzaliTa12TPE_KP_Rscan.csv"
CHIZZALI_REF = ROOT / "validation" / "chizzali_reference"
CHIZZALI_C3HALF = CHIZZALI_REF / "chizzali_fig2_c3half_screenshot_digitized_v0.csv"
PHIALPHA = ROOT / "cf_phialpha_folded.csv"
PHIALPHA_SUMMARY = ROOT / "figures" / "phi_alpha_potential_summary.csv"

COLORS = {
    "alice": "#000000",
    "hal": "#0072B2",
    "bound": "#D55E00",
    "gauss": "#CC79A7",
    "q": "#0072B2",
    "qd": "#D55E00",
    "fold": "#009E73",
}


def read_csv(path: Path) -> dict[str, list[float]]:
    with path.open(newline="") as f:
        lines = [line for line in f if line.strip() and not line.lstrip().startswith("#") and not line.lstrip().startswith('"#')]
    rows = list(csv.DictReader(lines))
    out: dict[str, list[float]] = {h: [] for h in rows[0].keys()} if rows else {}
    for row in rows:
        for key, value in row.items():
            if value is None or value.strip() == "":
                out[key].append(float("nan"))
            else:
                out[key].append(float(value))
    return out


def digitized_points(path: Path) -> dict[str, list[float]] | None:
    data = read_csv(path)
    if not data or "k_MeV" not in data or not data["k_MeV"]:
        return None
    if all(math.isnan(v) for v in data.get("C", [])):
        return None
    return data


def nearest_value(data: dict[str, list[float]], col: str, k: float) -> float:
    idx = min(range(len(data["k_MeV"])), key=lambda i: abs(data["k_MeV"][i] - k))
    return data[col][idx]


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


def max_delta(model: dict[str, list[float]], model_col: str, ref: dict[str, list[float]], ref_col: str, kmax: float) -> float:
    worst = 0.0
    for k, rval in zip(ref["k_MeV"], ref[ref_col]):
        if k <= kmax:
            worst = max(worst, abs(nearest_value(model, model_col, k) - rval))
    return worst


def add_caption(fig, text: str):
    fig.text(0.5, 0.025, text, ha="center", va="bottom", fontsize=9, wrap=True)


def dawson(x: float) -> float:
    nmax = 6
    h = 0.4
    a1, a2, a3 = 2.0 / 3.0, 0.4, 2.0 / 7.0
    coeff = [0.0] + [math.exp(-((2.0 * i - 1.0) * h) ** 2) for i in range(1, nmax + 1)]
    ax = abs(x)
    if ax < 0.2:
        x2 = x * x
        return x * (1.0 - a1 * x2 * (1.0 - a2 * x2 * (1.0 - a3 * x2)))
    n0 = 2 * int(0.5 * ax / h + 0.5)
    xp = ax - n0 * h
    e1 = math.exp(2.0 * xp * h)
    e2 = e1 * e1
    d1 = n0 + 1.0
    d2 = d1 - 2.0
    total = 0.0
    for i in range(1, nmax + 1):
        total += coeff[i] * (e1 / d1 + 1.0 / (d2 * e1))
        d1 += 2.0
        d2 -= 2.0
        e1 *= e2
    ans = (1.0 / math.sqrt(math.pi)) * math.exp(-xp * xp) * total
    return ans if x >= 0.0 else -ans


def F1(z: float) -> float:
    if z < 1e-8:
        return 1.0
    return dawson(z) / z


def F2(z: float) -> float:
    if z < 1e-8:
        return z
    return (1.0 - math.exp(-z * z)) / z


def ll_c(k_mev: float, R: float, f0: complex, d0: float) -> float:
    hbarc = 197.3269804
    k = max(k_mev / hbarc, 1e-6)
    amp = 1.0 / (1.0 / f0 + 0.5 * d0 * k * k - 1j * k)
    z = 2.0 * k * R
    t1 = 0.5 * abs(amp) ** 2 / (R * R) * (1.0 - d0 / (2.0 * math.sqrt(math.pi) * R))
    t2 = 2.0 * amp.real / (math.sqrt(math.pi) * R) * F1(z)
    t3 = -amp.imag / R * F2(z)
    return 1.0 + t1 + t2 + t3


def build_phi_p_debug_table() -> dict[str, list[float]]:
    ref = read_csv(PHIP_REF / "phi_p_ALICE_LL_reference.csv")
    repro = read_csv(PHIP_REF / "generated_ALICEfig_LL.csv")
    kgrid = [float(k) for k in range(1, 301)]
    out = {
        "kstar_MeV_c": kgrid,
        "C_alice_LL_R108": [],
        "C_alice_LL_R103": [],
        "C_alice_LL_R113": [],
        "C_cf_fitter_repro": [],
        "C_doublet_chizzali_literal_ERE": [],
        "C_doublet_chizzali_signflip_test": [],
        "C_quartet_HAL": [],
        "C_spinavg_chizzali_literal_ERE_HAL": [],
        "C_spinavg_doublet_free_HAL": [],
    }
    for k in kgrid:
        out["C_alice_LL_R108"].append(nearest_value(ref, "C_LL_R1.08", k))
        out["C_alice_LL_R103"].append(nearest_value(ref, "C_LL_R1.03", k))
        out["C_alice_LL_R113"].append(nearest_value(ref, "C_LL_R1.13", k))
        out["C_cf_fitter_repro"].append(nearest_value(repro, "LL_ALICEfig_R1.08", k))
        c_d = ll_c(k, 1.08, complex(-1.54, 0.0), 0.39)
        c_d_flip = ll_c(k, 1.08, complex(+1.54, 0.0), 0.39)
        c_q = ll_c(k, 1.08, complex(+1.43, 0.0), 2.36)
        out["C_doublet_chizzali_literal_ERE"].append(c_d)
        out["C_doublet_chizzali_signflip_test"].append(c_d_flip)
        out["C_quartet_HAL"].append(c_q)
        out["C_spinavg_chizzali_literal_ERE_HAL"].append((1.0 / 3.0) * c_d + (2.0 / 3.0) * c_q)
        out["C_spinavg_doublet_free_HAL"].append((1.0 / 3.0) * 1.0 + (2.0 / 3.0) * c_q)
    write_debug_csv(out)
    return out


def write_debug_csv(data: dict[str, list[float]]):
    path = ROOT / "figures" / "phi_p_spin_average_debug.csv"
    keys = list(data.keys())
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for row in zip(*(data[k] for k in keys)):
            w.writerow([f"{v:.8g}" for v in row])


def band_limits(data: dict[str, list[float]], cols: list[str]) -> tuple[list[float], list[float]]:
    lo, hi = [], []
    for vals in zip(*(data[c] for c in cols)):
        lo.append(min(vals))
        hi.append(max(vals))
    return lo, hi


def figure_phi_p_baseline(pdf: PdfPages | None = None):
    ref = read_csv(PHIP_REF / "phi_p_ALICE_LL_reference.csv")
    reproduced = read_csv(PHIP_REF / "generated_ALICEfig_LL.csv")
    dmax = max_delta(reproduced, "LL_ALICEfig_R1.08", ref, "C_LL_R1.08", 200.0)

    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    ax.fill_between(ref["k_MeV"], ref["C_LL_R1.03"], ref["C_LL_R1.13"], color="0.82", label="ALICE LL source band R=1.03-1.13 fm")
    ax.plot(ref["k_MeV"], ref["C_LL_R1.08"], color=COLORS["alice"], lw=2.3, label="ALICE LL fit, R=1.08 fm")
    ax.plot(reproduced["k_MeV"], reproduced["LL_ALICEfig_R1.08"], color="#D55E00", ls="--", lw=2.0, label="cf-fitter ALICE reproduction")
    ax.annotate(r"ALICE: $f_0=0.85+0.16i$ fm, $d_0=7.85$ fm", xy=(0.04, 0.92), xycoords="axes fraction", fontsize=10)
    ax.annotate(rf"max $|\Delta C|$ = {dmax:.2g} for $k^*<200$ MeV/$c$", xy=(0.04, 0.86), xycoords="axes fraction", fontsize=9)
    ax.axhline(1.0, color="0.55", lw=0.8)
    ax.set_xlim(0, 300)
    ax.set_ylim(0.75, 1.75)
    ax.set_xlabel(r"$k^*$ [MeV/$c$]")
    ax.set_ylabel(r"$C(k^*)$")
    ax.set_title(r"Figure 1. p-$\phi$ ALICE LL baseline validation at $R=1.08$ fm")
    ax.legend(fontsize=8.5, loc="upper right", frameon=False)
    add_caption(fig, "This page validates the empirical ALICE LL baseline only. No microscopic spin-channel model curves are shown here.")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(ROOT / "figures" / "fig01_alice_baseline_validation.pdf")
    fig.savefig(ROOT / "figures" / "fig01_alice_baseline_validation.png", dpi=180)
    if pdf is not None:
        pdf.savefig(fig)
    plt.close(fig)


def figure_phi_p_diagnostic(pdf: PdfPages | None = None):
    data = build_phi_p_debug_table()
    k = data["kstar_MeV_c"]
    ref = read_csv(PHIP_REF / "phi_p_ALICE_LL_reference.csv")
    dmax = max(abs(a - b) for kk, a, b in zip(data["kstar_MeV_c"], data["C_cf_fitter_repro"], data["C_alice_LL_R108"]) if kk <= 200.0)

    fig, ax = plt.subplots(figsize=(8.8, 5.5))
    ax.fill_between(ref["k_MeV"], ref["C_LL_R1.03"], ref["C_LL_R1.13"], color="0.88", label="ALICE LL source band")
    ax.plot(ref["k_MeV"], ref["C_LL_R1.08"], color="black", lw=2.2, label="ALICE empirical LL fit")
    ax.plot(k, data["C_doublet_chizzali_literal_ERE"], color="#009E73", lw=1.7, label=r"doublet LL/ERE using Chizzali $f_0,d_0$")
    ax.plot(k, data["C_quartet_HAL"], color=COLORS["hal"], lw=1.7, label=r"quartet LL/ERE using HAL $f_0,d_0$")
    ax.plot(k, data["C_spinavg_doublet_free_HAL"], color="#56B4E9", lw=1.6, ls=":", label=r"spin avg LL/ERE: free doublet + HAL quartet")
    ax.plot(k, data["C_spinavg_chizzali_literal_ERE_HAL"], color=COLORS["bound"], lw=2.2, label=r"spin avg LL/ERE: Chizzali doublet params + HAL quartet")
    ax.axhline(1.0, color="0.55", lw=0.8)
    ax.set_xlim(0, 300)
    ax.set_ylim(0.0, 3.05)
    ax.set_xlabel(r"$k^*$ [MeV/$c$]")
    ax.set_ylabel(r"$C(k^*)$")
    ax.set_title(r"Figure 2. p-$\phi$ LL/ERE spin-channel diagnostic at $R=1.08$ fm")
    ax.legend(fontsize=7.7, loc="upper right", frameon=False)
    ibin = min(range(len(k)), key=lambda i: abs(k[i] - 2.0))
    text = (rf"$k^*={k[ibin]:.0f}$ MeV/$c$: ALICE={data['C_alice_LL_R108'][ibin]:.3f}, "
            rf"doublet={data['C_doublet_chizzali_literal_ERE'][ibin]:.3f}, quartet={data['C_quartet_HAL'][ibin]:.3f}, "
            rf"spin avg={data['C_spinavg_chizzali_literal_ERE_HAL'][ibin]:.3f}")
    ax.annotate(text, xy=(0.56, 0.63), xycoords="axes fraction", fontsize=8.2)
    add_caption(fig, f"LL/ERE diagnostic only: Chizzali et al. use finite-range CATS potentials, not this direct scattering-parameter LL/ERE construction. cf-fitter ALICE reproduction remains max|Delta C|={dmax:.2g}.")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(ROOT / "figures" / "fig02_phi_p_spin_average_diagnostic.pdf")
    fig.savefig(ROOT / "figures" / "fig02_phi_p_spin_average_diagnostic.png", dpi=180)
    if pdf is not None:
        pdf.savefig(fig)
    plt.close(fig)


def write_table_csv(path: Path, data: dict[str, list[float]]):
    keys = list(data.keys())
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for row in zip(*(data[k] for k in keys)):
            w.writerow([f"{v:.8g}" for v in row])


def write_phi_p_kp_debug_csv(data: dict[str, list[float]]):
    write_table_csv(ROOT / "figures" / "phi_p_HAL_quartet_KP_debug.csv", data)
    write_table_csv(ROOT / "figures" / "phi_p_HAL_quartet_chizzali_C3half_comparison.csv", data)


def quartet_from_weighted_kp(c_weighted: float) -> float:
    return 1.0 + (c_weighted - 1.0) / (2.0 / 3.0)


PHI_P_SOURCE_R_SCAN_GRID = [0.90, 1.03, 1.08, 1.13, 1.20, 1.40, 1.60, 1.80, 2.00]
PHI_P_SOURCE_R_PLOT = [1.08, 1.40, 1.60, 1.80, 2.00]
PHI_P_SOURCE_R_NOMINAL = (1.03, 1.13)

PHI_P_KSTAR_SCALE_ALPHA_REQUIRED = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
PHI_P_KSTAR_SCALE_ALPHA_EXTRA = [1.1, 1.2, 1.3, 1.4, 1.6]
PHI_P_KSTAR_SCALE_ALPHA_ALL = sorted(set(PHI_P_KSTAR_SCALE_ALPHA_REQUIRED + PHI_P_KSTAR_SCALE_ALPHA_EXTRA))
PHI_P_KSTAR_SCALE_R_PRIMARY = 1.08
PHI_P_KSTAR_SCALE_R_BEST = 1.60
PHI_P_KSTAR_SCALE_R_GRID = [1.08, 1.40, 1.60, 1.80, 2.00]
PHI_P_KSTAR_SCALE_ALPHA_GRID = [0.75, 1.0, 1.25, 1.5, 2.0]


def source_moments_fm(R: float) -> tuple[float, float]:
    return 4.0 * R / math.sqrt(math.pi), math.sqrt(6.0) * R


def kp_col_chizzali_tpe(R: float) -> str:
    return f"KP_HALquartetChizzaliTa12TPE_R{R:.2f}"


def debug_col_chizzali_tpe_quartet(R: float) -> str:
    return "C_ChizzaliTa12TPE_KP_R" + f"{R:.2f}".replace(".", "p")


def quartet_kp_series(kp_data: dict[str, list[float]], R: float) -> list[float]:
    col = kp_col_chizzali_tpe(R)
    return [quartet_from_weighted_kp(v) for v in kp_data[col]]


def quartet_kp_at_k(kp_data: dict[str, list[float]], R: float, k: float) -> float:
    ks = kp_data["k_MeV"]
    ys = quartet_kp_series(kp_data, R)
    return interp_linear(ks, ys, k)


def chi2_like_to_chizzali_center(
    kp_data: dict[str, list[float]],
    chiz: dict[str, list[float]],
    R: float,
    kmax: float = 120.0,
    alpha: float = 1.0,
) -> float:
    diffs = []
    for k, c_ref in zip(chiz["kstar_MeV_c"], chiz["C3half_center"]):
        if k <= 0.0 or k > kmax:
            continue
        c_kp = quartet_kp_at_k(kp_data, R, alpha * k)
        diffs.append((c_kp - c_ref) ** 2)
    if not diffs:
        return float("nan")
    return sum(diffs) / len(diffs)


def alpha_col_suffix(alpha: float) -> str:
    s = f"{alpha:.2f}".rstrip("0").rstrip(".")
    return s.replace(".", "p")


def scaled_quartet_kp_series(kp_data: dict[str, list[float]], R: float, alpha: float) -> list[float]:
    ks = kp_data["k_MeV"]
    return [quartet_kp_at_k(kp_data, R, alpha * k) for k in ks]


def build_phi_p_kstar_scale_diagnostic_tables() -> tuple[dict[str, list[float]], list[dict[str, float]]]:
    kp = read_csv(PHIP_CHIZZALI_RSCAN)
    chiz = read_csv(CHIZZALI_C3HALF)
    kgrid = kp["k_MeV"]
    out: dict[str, list[float]] = {
        "kstar_MeV_c": kgrid,
        "C_chizzali_C3half_center_v0": [],
        "C_chizzali_C3half_lower_v0": [],
        "C_chizzali_C3half_upper_v0": [],
    }
    for alpha in PHI_P_KSTAR_SCALE_ALPHA_ALL:
        out[f"C_R1p08_alpha{alpha_col_suffix(alpha)}"] = scaled_quartet_kp_series(
            kp, PHI_P_KSTAR_SCALE_R_PRIMARY, alpha)
    out["C_R1p60_alpha1p0"] = scaled_quartet_kp_series(kp, PHI_P_KSTAR_SCALE_R_BEST, 1.0)

    for k in kgrid:
        out["C_chizzali_C3half_center_v0"].append(
            interp_linear(chiz["kstar_MeV_c"], chiz["C3half_center"], k))
        out["C_chizzali_C3half_lower_v0"].append(
            interp_linear(chiz["kstar_MeV_c"], chiz["C3half_lower"], k))
        out["C_chizzali_C3half_upper_v0"].append(
            interp_linear(chiz["kstar_MeV_c"], chiz["C3half_upper"], k))

    summary_rows: list[dict[str, float]] = []
    seen_ra: set[tuple[float, float]] = set()

    def append_summary_row(R: float, alpha: float) -> None:
        key = (R, alpha)
        if key in seen_ra:
            return
        seen_ra.add(key)
        summary_rows.append({
            "R_fm": R,
            "alpha": alpha,
            "chi2_like_to_chizzali_C3half_v0_0ltkle120": chi2_like_to_chizzali_center(
                kp, chiz, R, alpha=alpha),
            "C_scaled_at_k2": quartet_kp_at_k(kp, R, alpha * 2.0),
            "C_scaled_at_k10": quartet_kp_at_k(kp, R, alpha * 10.0),
            "C_scaled_at_k50": quartet_kp_at_k(kp, R, alpha * 50.0),
        })

    for alpha in PHI_P_KSTAR_SCALE_ALPHA_ALL:
        append_summary_row(PHI_P_KSTAR_SCALE_R_PRIMARY, alpha)
    append_summary_row(PHI_P_KSTAR_SCALE_R_BEST, 1.0)
    for R in PHI_P_KSTAR_SCALE_R_GRID:
        for alpha in PHI_P_KSTAR_SCALE_ALPHA_GRID:
            append_summary_row(R, alpha)

    write_table_csv(ROOT / "figures" / "phi_p_kstar_scale_diagnostic.csv", out)
    with (ROOT / "figures" / "phi_p_kstar_scale_diagnostic_summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        w.writeheader()
        for row in summary_rows:
            w.writerow({k: f"{v:.8g}" for k, v in row.items()})
    return out, summary_rows


def build_phi_p_source_R_scan_tables() -> tuple[dict[str, list[float]], list[dict[str, float]]]:
    kp = read_csv(PHIP_CHIZZALI_RSCAN)
    chiz = read_csv(CHIZZALI_C3HALF)
    kgrid = kp["k_MeV"]
    out: dict[str, list[float]] = {
        "kstar_MeV_c": kgrid,
        "C_chizzali_C3half_digitized_center_v0": [],
        "C_chizzali_C3half_digitized_lower_v0": [],
        "C_chizzali_C3half_digitized_upper_v0": [],
    }
    for R in PHI_P_SOURCE_R_SCAN_GRID:
        out[debug_col_chizzali_tpe_quartet(R)] = quartet_kp_series(kp, R)
    for k in kgrid:
        out["C_chizzali_C3half_digitized_center_v0"].append(
            interp_linear(chiz["kstar_MeV_c"], chiz["C3half_center"], k))
        out["C_chizzali_C3half_digitized_lower_v0"].append(
            interp_linear(chiz["kstar_MeV_c"], chiz["C3half_lower"], k))
        out["C_chizzali_C3half_digitized_upper_v0"].append(
            interp_linear(chiz["kstar_MeV_c"], chiz["C3half_upper"], k))

    summary_rows: list[dict[str, float]] = []
    for R in PHI_P_SOURCE_R_SCAN_GRID:
        mean_r, rms_r = source_moments_fm(R)
        summary_rows.append({
            "R_fm": R,
            "source_mean_r_fm": mean_r,
            "source_rms_r_fm": rms_r,
            "C_k2": quartet_kp_at_k(kp, R, 2.0),
            "C_k10": quartet_kp_at_k(kp, R, 10.0),
            "C_k50": quartet_kp_at_k(kp, R, 50.0),
            "chi2_like_to_chizzali_C3half_v0_0to120": chi2_like_to_chizzali_center(kp, chiz, R),
        })

    write_table_csv(ROOT / "figures" / "phi_p_source_R_scan_debug.csv", out)
    with (ROOT / "figures" / "phi_p_source_R_scan_summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        w.writeheader()
        for row in summary_rows:
            w.writerow({k: f"{v:.8g}" for k, v in row.items()})
    return out, summary_rows


def build_phi_p_chizzali_ta12_tpe_debug_table() -> dict[str, list[float]]:
    ref = read_csv(PHIP_REF / "phi_p_ALICE_LL_reference.csv")
    kp_fitb = read_csv(PHIP_HAL_KP)
    kp_chiz = read_csv(PHIP_CHIZZALI_KP)
    chiz = read_csv(CHIZZALI_C3HALF)
    kgrid = kp_chiz["k_MeV"]
    out = {
        "kstar_MeV_c": kgrid,
        "C_alice_empirical_LL_R108": [],
        "C_quartet_HAL_fitB_t14_LL_ERE": [],
        "C_quartet_HAL_fitB_t14_KP_R108": [],
        "C_quartet_HAL_ChizzaliTa12TPE_KP_R103": [],
        "C_quartet_HAL_ChizzaliTa12TPE_KP_R108": [],
        "C_quartet_HAL_ChizzaliTa12TPE_KP_R113": [],
        "C_chizzali_C3half_digitized_lower_v0": [],
        "C_chizzali_C3half_digitized_center_v0": [],
        "C_chizzali_C3half_digitized_upper_v0": [],
        "C_spinavg_free_doublet_plus_fitB_t14_KP": [],
        "C_spinavg_free_doublet_plus_ChizzaliTa12TPE_KP": [],
    }
    for k in kgrid:
        c_ll = ll_c(k, 1.08, complex(+1.43, 0.0), 2.36)
        c_fitb_kp = quartet_from_weighted_kp(nearest_value(kp_fitb, "KP_HALquartetFitB_R1.08", k))
        c_chiz_103 = quartet_from_weighted_kp(nearest_value(kp_chiz, "KP_HALquartetChizzaliTa12TPE_R1.03", k))
        c_chiz_108 = quartet_from_weighted_kp(nearest_value(kp_chiz, "KP_HALquartetChizzaliTa12TPE_R1.08", k))
        c_chiz_113 = quartet_from_weighted_kp(nearest_value(kp_chiz, "KP_HALquartetChizzaliTa12TPE_R1.13", k))
        out["C_alice_empirical_LL_R108"].append(nearest_value(ref, "C_LL_R1.08", k))
        out["C_quartet_HAL_fitB_t14_LL_ERE"].append(c_ll)
        out["C_quartet_HAL_fitB_t14_KP_R108"].append(c_fitb_kp)
        out["C_quartet_HAL_ChizzaliTa12TPE_KP_R103"].append(c_chiz_103)
        out["C_quartet_HAL_ChizzaliTa12TPE_KP_R108"].append(c_chiz_108)
        out["C_quartet_HAL_ChizzaliTa12TPE_KP_R113"].append(c_chiz_113)
        out["C_chizzali_C3half_digitized_lower_v0"].append(interp_linear(chiz["kstar_MeV_c"], chiz["C3half_lower"], k))
        out["C_chizzali_C3half_digitized_center_v0"].append(interp_linear(chiz["kstar_MeV_c"], chiz["C3half_center"], k))
        out["C_chizzali_C3half_digitized_upper_v0"].append(interp_linear(chiz["kstar_MeV_c"], chiz["C3half_upper"], k))
        out["C_spinavg_free_doublet_plus_fitB_t14_KP"].append(nearest_value(kp_fitb, "KP_HALquartetFitB_R1.08", k))
        out["C_spinavg_free_doublet_plus_ChizzaliTa12TPE_KP"].append(nearest_value(kp_chiz, "KP_HALquartetChizzaliTa12TPE_R1.08", k))
    write_table_csv(ROOT / "figures" / "phi_p_chizzali_ta12_TPE_quartet_debug.csv", out)
    return out


def build_phi_p_kp_debug_table() -> dict[str, list[float]]:
    ref = read_csv(PHIP_REF / "phi_p_ALICE_LL_reference.csv")
    kp = read_csv(PHIP_HAL_KP)
    chiz = read_csv(CHIZZALI_C3HALF)
    kgrid = kp["k_MeV"]
    out = {
        "kstar_MeV_c": kgrid,
        "C_alice_empirical_LL_R108": [],
        "C_quartet_HAL_LL_ERE": [],
        "C_quartet_HAL_KP_Numerov_R103": [],
        "C_quartet_HAL_KP_Numerov_R108": [],
        "C_quartet_HAL_KP_Numerov_R113": [],
        "C_spinavg_free_doublet_plus_HAL_quartet_LL_ERE": [],
        "C_spinavg_free_doublet_plus_HAL_quartet_KP_R108": [],
        "C_chizzali_C3half_lower_digitized_v0": [],
        "C_chizzali_C3half_center_digitized_v0": [],
        "C_chizzali_C3half_upper_digitized_v0": [],
        "Delta_quartet_HAL_LL_ERE_minus_chizzali_center": [],
        "Delta_quartet_HAL_KP_Numerov_minus_chizzali_center": [],
    }
    for k in kgrid:
        c_alice = nearest_value(ref, "C_LL_R1.08", k)
        c_ll = ll_c(k, 1.08, complex(+1.43, 0.0), 2.36)
        c_kp_spinavg_103 = nearest_value(kp, "KP_HALquartetFitB_R1.03", k)
        c_kp_spinavg_108 = nearest_value(kp, "KP_HALquartetFitB_R1.08", k)
        c_kp_spinavg_113 = nearest_value(kp, "KP_HALquartetFitB_R1.13", k)
        c_kp_103 = quartet_from_weighted_kp(c_kp_spinavg_103)
        c_kp_108 = quartet_from_weighted_kp(c_kp_spinavg_108)
        c_kp_113 = quartet_from_weighted_kp(c_kp_spinavg_113)
        out["C_alice_empirical_LL_R108"].append(c_alice)
        out["C_quartet_HAL_LL_ERE"].append(c_ll)
        out["C_quartet_HAL_KP_Numerov_R103"].append(c_kp_103)
        out["C_quartet_HAL_KP_Numerov_R108"].append(c_kp_108)
        out["C_quartet_HAL_KP_Numerov_R113"].append(c_kp_113)
        c_chiz_lo = interp_linear(chiz["kstar_MeV_c"], chiz["C3half_lower"], k)
        c_chiz = interp_linear(chiz["kstar_MeV_c"], chiz["C3half_center"], k)
        c_chiz_hi = interp_linear(chiz["kstar_MeV_c"], chiz["C3half_upper"], k)
        out["C_spinavg_free_doublet_plus_HAL_quartet_LL_ERE"].append((1.0 / 3.0) + (2.0 / 3.0) * c_ll)
        out["C_spinavg_free_doublet_plus_HAL_quartet_KP_R108"].append(c_kp_spinavg_108)
        out["C_chizzali_C3half_lower_digitized_v0"].append(c_chiz_lo)
        out["C_chizzali_C3half_center_digitized_v0"].append(c_chiz)
        out["C_chizzali_C3half_upper_digitized_v0"].append(c_chiz_hi)
        out["Delta_quartet_HAL_LL_ERE_minus_chizzali_center"].append(c_ll - c_chiz)
        out["Delta_quartet_HAL_KP_Numerov_minus_chizzali_center"].append(c_kp_108 - c_chiz)
    write_phi_p_kp_debug_csv(out)
    return out


def figure_phi_p_kp_diagnostic(pdf: PdfPages | None = None):
    data = build_phi_p_kp_debug_table()
    chiz_data = build_phi_p_chizzali_ta12_tpe_debug_table()
    chiz = read_csv(CHIZZALI_C3HALF)
    k = data["kstar_MeV_c"]
    fig, ax = plt.subplots(figsize=(8.8, 5.5))
    kp_lo, kp_hi = band_limits(data, ["C_quartet_HAL_KP_Numerov_R103", "C_quartet_HAL_KP_Numerov_R113"])
    ax.fill_between(chiz["kstar_MeV_c"], chiz["C3half_lower"], chiz["C3half_upper"], color="0.55", alpha=0.22, label=r"Chizzali Fig.2 $C_{3/2}$, screenshot digitized v0")
    ax.plot(chiz["kstar_MeV_c"], chiz["C3half_center"], color="0.25", lw=1.5, marker="o", ms=3.0, label=r"Chizzali $C_{3/2}$ center, qualitative")
    ax.fill_between(k, kp_lo, kp_hi, color="#E69F00", alpha=0.18, label=r"HAL fit-B t/a=14 KP source band $R=1.03$-$1.13$ fm")
    ax.plot(k, data["C_alice_empirical_LL_R108"], color="black", lw=2.2, label=r"ALICE empirical LL fit, $R=1.08$ fm")
    ax.plot(k, data["C_quartet_HAL_LL_ERE"], color=COLORS["hal"], lw=1.9, ls="--", label=r"HAL quartet t/a=14 fit-B LL/ERE")
    ax.plot(k, data["C_quartet_HAL_KP_Numerov_R108"], color="#D55E00", lw=2.2, label=r"HAL quartet t/a=14 fit-B KP/Numerov")
    ax.plot(k, chiz_data["C_quartet_HAL_ChizzaliTa12TPE_KP_R108"], color="#0072B2", lw=2.2, label=r"HAL quartet Chizzali t/a=12 TPE KP/Numerov")
    ax.plot(k, chiz_data["C_spinavg_free_doublet_plus_ChizzaliTa12TPE_KP"], color="#56B4E9", lw=1.4, ls=":", label=r"spin avg diagnostic: free doublet + Chizzali t/a=12 quartet KP")
    ax.plot(k, data["C_spinavg_free_doublet_plus_HAL_quartet_KP_R108"], color="#CC79A7", lw=1.8, ls="-.", label=r"spin avg diagnostic: free doublet + fit-B quartet KP")
    ax.axhline(1.0, color="0.55", lw=0.8)
    ax.set_xlim(0, 300)
    ax.set_ylim(0.75, 3.15)
    ax.set_xlabel(r"$k^*$ [MeV/$c$]")
    ax.set_ylabel(r"$C(k^*)$")
    ax.set_title(r"Figure 3. p-$\phi$ HAL quartet comparison with Chizzali $C_{3/2}$ reference")
    ax.legend(fontsize=6.2, loc="upper right", frameon=False)
    ibin = min(range(len(k)), key=lambda i: abs(k[i] - 2.0))
    text = (rf"$k^*={k[ibin]:.0f}$ MeV/$c$: fit-B KP={data['C_quartet_HAL_KP_Numerov_R108'][ibin]:.3f}, "
            rf"Chizzali TPE KP={chiz_data['C_quartet_HAL_ChizzaliTa12TPE_KP_R108'][ibin]:.3f}, "
            rf"Chizzali ref center={chiz_data['C_chizzali_C3half_digitized_center_v0'][ibin]:.3f}")
    ax.annotate(text, xy=(0.50, 0.55), xycoords="axes fraction", fontsize=8.0)
    ax.annotate(r"Chizzali reference is screenshot-digitized qualitative v0.", xy=(0.50, 0.49), xycoords="axes fraction", fontsize=8.0)
    ax.annotate(r"Spin averages are diagnostics only; Chizzali full doublet+quartet CATS model is not implemented.", xy=(0.50, 0.43), xycoords="axes fraction", fontsize=8.0)
    add_caption(fig, "Chizzali C3/2 reference is a rough screenshot digitization for qualitative comparison only. Chizzali t/a=12 TPE uses arXiv:2212.12690 quartet parameters (beta=1, gamma=0); fit-B uses t/a=14 three-Gaussian HAL input. Neither includes Chizzali fitted doublet KP/CATS.")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(ROOT / "figures" / "fig03_phi_p_HAL_quartet_KP_diagnostic.pdf")
    fig.savefig(ROOT / "figures" / "fig03_phi_p_HAL_quartet_KP_diagnostic.png", dpi=180)
    if pdf is not None:
        pdf.savefig(fig)
    plt.close(fig)


def figure_phi_p_source_R_scan(pdf: PdfPages | None = None):
    data, summary = build_phi_p_source_R_scan_tables()
    chiz = read_csv(CHIZZALI_C3HALF)
    k = data["kstar_MeV_c"]
    best = min(summary, key=lambda row: row["chi2_like_to_chizzali_C3half_v0_0to120"])
    best_R = best["R_fm"]

    fig, ax = plt.subplots(figsize=(8.8, 5.5))
    ax.fill_between(
        chiz["kstar_MeV_c"], chiz["C3half_lower"], chiz["C3half_upper"],
        color="0.55", alpha=0.22, label=r"Chizzali Fig.2 $C_{3/2}$, screenshot digitized v0",
    )
    nominal_lo = data[debug_col_chizzali_tpe_quartet(PHI_P_SOURCE_R_NOMINAL[0])]
    nominal_hi = data[debug_col_chizzali_tpe_quartet(PHI_P_SOURCE_R_NOMINAL[1])]
    ax.fill_between(k, nominal_lo, nominal_hi, color="#0072B2", alpha=0.10,
                    label=r"Chizzali t/a=12 TPE KP nominal band, $R=1.03$-$1.13$ fm")

    cmap = plt.cm.viridis
    n_plot = len(PHI_P_SOURCE_R_PLOT)
    for i, R in enumerate(PHI_P_SOURCE_R_PLOT):
        color = cmap(0.15 + 0.75 * i / max(n_plot - 1, 1))
        col = debug_col_chizzali_tpe_quartet(R)
        ax.plot(k, data[col], color=color, lw=2.0 if R == 1.08 else 1.8,
                label=rf"Chizzali t/a=12 TPE KP, $R={R:.2f}$ fm")

    ax.axhline(1.0, color="0.55", lw=0.8)
    ax.set_xlim(0, 300)
    ax.set_ylim(0.75, 3.15)
    ax.set_xlabel(r"$k^*$ [MeV/$c$]")
    ax.set_ylabel(r"$C(k^*)$")
    ax.set_title(r"Figure 4. p-$\phi$ source-radius scan for HAL quartet KP/Numerov")
    ax.legend(fontsize=6.4, loc="upper right", frameon=False)
    ibin = min(range(len(k)), key=lambda i: abs(k[i] - 2.0))
    text = (rf"$k^*={k[ibin]:.0f}$ MeV/$c$: $R=1.08$ KP={data[debug_col_chizzali_tpe_quartet(1.08)][ibin]:.3f}, "
            rf"Chizzali ref center={data['C_chizzali_C3half_digitized_center_v0'][ibin]:.3f}, "
            rf"qualitative best $R$={best_R:.2f} fm")
    ax.annotate(text, xy=(0.48, 0.58), xycoords="axes fraction", fontsize=8.0)
    ax.annotate(
        r"Source scan uses $S(r)=\exp[-r^2/(4R^2)]/(4\pi R^2)^{3/2}$. "
        r"Chizzali reference is screenshot-digitized qualitative v0.",
        xy=(0.48, 0.52), xycoords="axes fraction", fontsize=8.0,
    )
    ax.annotate(
        r"Qualitative effective-$R$ matching excludes $k^*=0$ because the KP grid starts at $k^*=1$ MeV/$c$.",
        xy=(0.48, 0.46), xycoords="axes fraction", fontsize=8.0,
    )
    add_caption(
        fig,
        "Pure quartet KP/Numerov curves from HALquartetChizzaliTa12TPE. "
        "Effective-R metric is a qualitative chi2-like mean over Chizzali digitized center points with "
        "0 < k* <= 120 MeV/c; it is not a quantitative Chizzali reproduction.",
    )
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(ROOT / "figures" / "fig04_phi_p_source_R_scan.pdf")
    fig.savefig(ROOT / "figures" / "fig04_phi_p_source_R_scan.png", dpi=180)
    if pdf is not None:
        pdf.savefig(fig)
    plt.close(fig)
    return summary


def figure_phi_p_kstar_scale_diagnostic(pdf: PdfPages | None = None):
    data, summary = build_phi_p_kstar_scale_diagnostic_tables()
    chiz = read_csv(CHIZZALI_C3HALF)
    k = data["kstar_MeV_c"]
    r108_rows = [row for row in summary if row["R_fm"] == PHI_P_KSTAR_SCALE_R_PRIMARY]
    best_alpha = min(r108_rows, key=lambda row: row["chi2_like_to_chizzali_C3half_v0_0ltkle120"])

    fig, ax = plt.subplots(figsize=(8.8, 5.5))
    ax.fill_between(
        chiz["kstar_MeV_c"], chiz["C3half_lower"], chiz["C3half_upper"],
        color="0.55", alpha=0.22, label=r"Chizzali Fig.2 $C_{3/2}$, screenshot digitized v0",
    )
    ax.plot(
        chiz["kstar_MeV_c"], chiz["C3half_center"],
        color="0.25", lw=1.5, marker="o", ms=3.0,
    )

    plot_alphas = [1.0, 1.25, 1.5, 2.0]
    colors = ["#0072B2", "#009E73", "#D55E00", "#CC79A7"]
    for alpha, color in zip(plot_alphas, colors):
        col = f"C_R1p08_alpha{alpha_col_suffix(alpha)}"
        lw = 2.2 if alpha == 1.0 else 1.8
        ax.plot(k, data[col], color=color, lw=lw, label=rf"$R=1.08$, $\alpha={alpha:g}$")

    ax.plot(
        k, data["C_R1p60_alpha1p0"], color="#56B4E9", lw=1.8, ls="--",
        label=r"$R=1.60$, $\alpha=1.0$ source-scan best grid",
    )

    ax.axhline(1.0, color="0.55", lw=0.8)
    ax.set_xlim(0, 300)
    ax.set_ylim(0.75, 3.15)
    ax.set_xlabel(r"$k^*$ [MeV/$c$]")
    ax.set_ylabel(r"$C(k^*)$")
    ax.set_title(r"Figure 5. p-$\phi$ $k^*$ scale diagnostic for Chizzali $C_{3/2}$ reference")
    ax.legend(fontsize=6.4, loc="upper right", frameon=False)
    ibin = min(range(len(k)), key=lambda i: abs(k[i] - 2.0))
    text = (rf"$k^*={k[ibin]:.0f}$ MeV/$c$: $R=1.08$, $\alpha=1$ scaled={data['C_R1p08_alpha1'][ibin]:.3f}, "
            rf"Chizzali ref center={data['C_chizzali_C3half_center_v0'][ibin]:.3f}, "
            rf"best $\alpha$ at $R=1.08$={best_alpha['alpha']:g}")
    ax.annotate(text, xy=(0.48, 0.58), xycoords="axes fraction", fontsize=8.0)
    ax.annotate(
        r"Diagnostic only: $C_{\mathrm{scaled}}(k^*)=C_{\mathrm{KP}}(\alpha k^*)$. "
        r"$\alpha=2$ tests a possible $k^*$ vs $q=2k^*$ convention mismatch.",
        xy=(0.48, 0.52), xycoords="axes fraction", fontsize=8.0,
    )
    ax.annotate(
        r"Chizzali reference is screenshot-digitized qualitative v0. No physical correction is applied.",
        xy=(0.48, 0.46), xycoords="axes fraction", fontsize=8.0,
    )
    add_caption(
        fig,
        "Pure quartet KP/Numerov from HALquartetChizzaliTa12TPE. "
        "Horizontal-axis rescaling is a diagnostic interpolation only; "
        "chi2-like metric is mean over Chizzali native points with 0 < k* <= 120 MeV/c.",
    )
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(ROOT / "figures" / "fig05_phi_p_kstar_scale_diagnostic.pdf")
    fig.savefig(ROOT / "figures" / "fig05_phi_p_kstar_scale_diagnostic.png", dpi=180)
    if pdf is not None:
        pdf.savefig(fig)
    plt.close(fig)
    return summary


def shade_qd_band(ax, data, base: str, R: float, color: str, alpha: float = 0.16):
    rs = f"R{R:.2f}"
    lo, hi = band_limits(data, [f"KP_qd_b040_{base}_{rs}", f"KP_qd_{base}_{rs}", f"KP_qd_b070_{base}_{rs}"])
    ax.fill_between(data["k_MeV"], lo, hi, color=color, alpha=alpha, lw=0)


def figure_phi_alpha_matrix(pdf: PdfPages):
    data = read_csv(PHIALPHA)
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.6), sharey=True)
    for ax, R in zip(axes, [1.0, 3.0, 5.0]):
        rs = f"R{R:.2f}"
        shade_qd_band(ax, data, "nofold", R, COLORS["qd"], 0.13)
        shade_qd_band(ax, data, "fold", R, COLORS["qd"], 0.20)
        ax.plot(data["k_MeV"], data[f"KP_q_nofold_{rs}"], color=COLORS["q"], lw=1.3, alpha=0.72, label="KP quartet only, no fold")
        ax.plot(data["k_MeV"], data[f"KP_q_fold_{rs}"], color=COLORS["q"], lw=2.3, label="KP quartet only, fold")
        ax.plot(data["k_MeV"], data[f"KP_qd_nofold_{rs}"], color=COLORS["qd"], lw=1.3, alpha=0.80, label="KP q+d, no fold")
        ax.plot(data["k_MeV"], data[f"KP_qd_fold_{rs}"], color=COLORS["qd"], lw=2.3, label="KP q+d, fold")
        for name, label, color in [("q_fold", "LL quartet only, fold", COLORS["q"]), ("qd_fold", "LL q+d, fold", COLORS["qd"])]:
            ax.plot(data["k_MeV"], data[f"LL_{name}_{rs}"], color=color, lw=1.4, ls="--", alpha=0.95, label=label)
        ax.axhline(1.0, color="0.55", lw=0.8)
        ax.set_xlim(0, 200)
        ax.set_xlabel(r"$k^*$ [MeV/$c$]")
        ax.set_title(rf"$R={R:g}$ fm")
        ax.grid(alpha=0.18, lw=0.5)
    axes[0].set_ylabel(r"$C(k^*)$")
    axes[0].set_ylim(0.15, 1.15)
    axes[0].legend(fontsize=7.0, loc="lower right", frameon=False)
    fig.suptitle(r"Figure 6. $\phi\alpha$ folded-potential scenario matrix", y=0.98)
    add_caption(fig, "Solid lines are KP central curves; orange shading is the doublet range band from b_d=0.4-0.7 fm. Dashed lines are LL for folded cells. The doublet volume integral is now constrained to the phi-N range, removing the previous near-threshold pathology.")
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    pdf.savefig(fig)
    plt.close(fig)

def figure_etminan(pdf: PdfPages):
    data = read_csv(PHIALPHA)
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.5), sharey=True)
    note = (
        "Quartet-only/fold is the Etminan-like configuration in this framework; "
        "orange band scans the doublet range."
    )

    for ax, R in zip(axes, [1.0, 3.0, 5.0]):
        rs = f"R{R:.2f}"
        shade_qd_band(ax, data, "fold", R, COLORS["qd"], 0.22)
        ax.plot(data["k_MeV"], data[f"KP_q_fold_{rs}"], color=COLORS["q"], lw=2.2, label="Etminan configuration: q/fold KP")
        ax.plot(data["k_MeV"], data[f"KP_qd_fold_{rs}"], color=COLORS["qd"], lw=2.2, label="This work: q+d/fold KP")
        ax.plot(data["k_MeV"], data[f"LL_q_fold_{rs}"], color=COLORS["q"], lw=1.3, ls="--", label="LL q/fold")
        ax.plot(data["k_MeV"], data[f"LL_qd_fold_{rs}"], color=COLORS["qd"], lw=1.3, ls="--", label="LL q+d/fold")
        ax.axhline(1.0, color="0.55", lw=0.8)
        ax.set_xlim(0, 200)
        ax.set_ylim(0.15, 1.15)
        ax.set_xlabel(r"$k^*$ [MeV/$c$]")
        ax.set_title(rf"$R={R:g}$ fm")
        ax.grid(alpha=0.18, lw=0.5)
    axes[0].set_ylabel(r"$C(k^*)$")
    axes[0].legend(fontsize=7.3, frameon=False, loc="lower right")
    fig.suptitle(r"Figure 7. $\phi\alpha$: Etminan configuration vs doublet increment", y=0.98)
    add_caption(fig, note)
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    pdf.savefig(fig)
    plt.close(fig)


def figure_table(pdf: PdfPages):
    if not PHIALPHA_SUMMARY.exists():
        raise SystemExit(
            "Missing phi-alpha potential summary CSV: "
            f"{PHIALPHA_SUMMARY.relative_to(ROOT)}\n"
            "Generate with: build/cf-calc phi-alpha-summary "
            "or scripts/regenerate_phialpha_pipeline.sh"
        )

    with PHIALPHA_SUMMARY.open(newline="") as f:
        summary_rows = list(csv.DictReader(f))
    primary_order = ["q_nofold", "q_fold", "qd_nofold", "qd_fold"]
    by_name = {row["scenario"]: row for row in summary_rows}
    rows = []
    for name in primary_order:
        if name not in by_name:
            raise SystemExit(f"phi_alpha_potential_summary.csv missing scenario row: {name}")
        row = by_name[name]
        label = (row.get("table_label") or "").strip() or name
        be_raw = (row.get("BE_MeV") or "").strip()
        rows.append([
            label,
            f"{float(row['f0_fm']):+.3f}",
            f"{float(row['d0_fm']):.3f}",
            f"{float(be_raw):.2f}" if be_raw else "",
            f"{float(row['V_at_0_MeV']):.1f}",
            f"{float(row['volume_integral_MeV_fm3']):.1f}",
        ])

    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.axis("off")
    ax.set_title(r"Figure 8. $\phi\alpha$ scenario summary", fontsize=14, pad=16)
    table = ax.table(
        cellText=rows,
        colLabels=["cell", r"$f_0$ [fm]", r"$d_0$ [fm]", "BE [MeV]", r"$V(0)$ [MeV]", r"$\int V d^3r$ [MeV fm$^3$]"],
        cellLoc="center",
        colLoc="center",
        loc="center",
        colWidths=[0.30, 0.11, 0.11, 0.11, 0.14, 0.20],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.45)
    for (r, _c), cell in table.get_celld().items():
        if r == 0:
            cell.set_text_props(weight="bold")
            cell.set_facecolor("#EAEAEA")
    note = (
        "Inputs: HAL QCD quartet fit-B; central doublet is Gaussian b=0.55 fm, V0=571.5 MeV "
        "(f0=-1.539 fm, d0=0.582 fm, ERE BE~30.1 MeV, integral=-529.5 MeV fm^3). "
        "Doublet band scans b=0.4-0.7 fm retuned to f0=-1.54; d0 remains above Chizzali's 0.39 fm because a single Gaussian cannot reach that short effective range. "
        "He-4 single-Gaussian matter density uses rms=1.56 fm, A=4."
    )
    add_caption(fig, note)
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))
    pdf.savefig(fig)
    plt.close(fig)


def figure_phi_alpha_chizzali_tpe(pdf: PdfPages):
    chizzali_csv = ROOT / "validation" / "phi_alpha" / "generated_q_chizzali_tpe.csv"
    if not chizzali_csv.exists():
        print("skip Chizzali TPE phi-alpha page: missing validation/phi_alpha/generated_q_chizzali_tpe.csv")
        return
    sys.path.insert(0, str(ROOT / "figures"))
    from make_phi_alpha_chizzali_tpe_figure import make_figure as chizzali_tpe_figure

    chizzali_tpe_figure(pdf)


def figure_phi_alpha_doublet_shape_uncertainty(pdf: PdfPages):
    if not PHIALPHA.exists():
        print("skip doublet-shape uncertainty page: missing cf_phialpha_folded.csv")
        return
    sys.path.insert(0, str(ROOT / "figures"))
    from make_phi_alpha_doublet_shape_uncertainty import make_figure as doublet_uncertainty_figure

    doublet_uncertainty_figure(read_csv(PHIALPHA), pdf)


def figure_phi_alpha_doublet_potential_shapes(pdf: PdfPages):
    sys.path.insert(0, str(ROOT / "figures"))
    from make_phi_alpha_doublet_potential_shapes import make_figure as potential_shapes_figure

    potential_shapes_figure(pdf)


def main():
    required = [
        PHIP_REF / "phi_p_ALICE_LL_reference.csv",
        PHIP_REF / "generated_ALICEfig_LL.csv",
        PHIP_HAL_KP,
        PHIP_CHIZZALI_KP,
        PHIP_CHIZZALI_RSCAN,
        CHIZZALI_C3HALF,
        PHIALPHA,
        PHIALPHA_SUMMARY,
    ]
    missing = [str(p.relative_to(ROOT)) for p in required if not p.exists()]
    if missing:
        raise SystemExit("Missing required input CSVs:\n" + "\n".join(missing))
    OUT.parent.mkdir(exist_ok=True)
    with PdfPages(OUT) as pdf:
        figure_phi_p_baseline(pdf)
        figure_phi_p_diagnostic(pdf)
        figure_phi_p_kp_diagnostic(pdf)
        figure_phi_p_source_R_scan(pdf)
        figure_phi_p_kstar_scale_diagnostic(pdf)
        figure_phi_alpha_matrix(pdf)
        figure_etminan(pdf)
        figure_table(pdf)
        figure_phi_alpha_chizzali_tpe(pdf)
        figure_phi_alpha_doublet_shape_uncertainty(pdf)
        figure_phi_alpha_doublet_potential_shapes(pdf)
    legacy = ROOT / "figures" / "cf_fitter_summary.pdf"
    legacy.write_bytes(OUT.read_bytes())
    print(f"wrote {OUT.relative_to(ROOT)}")
    print(f"wrote {legacy.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
