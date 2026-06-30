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
OUT = ROOT / "figures" / "cf_fitter_summary_revised.pdf"
PHIP_REF = ROOT / "validation" / "phi_p"
PHIP_MODELS = PHIP_REF / "generated_phi_proton_spinavg_R1p08.csv"
PHIP_GAUSS = PHIP_REF / "generated_ALICEgauss_KP.csv"
PHIALPHA = ROOT / "cf_phialpha_folded.csv"
ETMINAN = ROOT / "validation" / "etminan_phialpha_reference.csv"

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
    ax.plot(ref["k_MeV"], ref["C_LL_R1.08"], color="black", lw=2.2, label="ALICE empirical LL")
    ax.plot(k, data["C_doublet_chizzali_literal_ERE"], color="#009E73", lw=1.7, label=r"doublet literal Chizzali ERE: $f_0=-1.54$, $d_0=0.39$")
    ax.plot(k, data["C_quartet_HAL"], color=COLORS["hal"], lw=1.7, label=r"quartet HAL ERE: $f_0=+1.43$, $d_0=2.36$")
    ax.plot(k, data["C_spinavg_doublet_free_HAL"], color="#56B4E9", lw=1.6, ls=":", label=r"spin avg: free doublet + HAL quartet")
    ax.plot(k, data["C_spinavg_chizzali_literal_ERE_HAL"], color=COLORS["bound"], lw=2.2, label=r"spin avg: Chizzali doublet + HAL quartet")
    ax.axhline(1.0, color="0.55", lw=0.8)
    ax.set_xlim(0, 300)
    ax.set_ylim(0.45, 2.75)
    ax.set_xlabel(r"$k^*$ [MeV/$c$]")
    ax.set_ylabel(r"$C(k^*)$")
    ax.set_title(r"Figure 2. p-$\phi$ spin-channel diagnostic at $R=1.08$ fm")
    ax.legend(fontsize=7.7, loc="upper right", frameon=False)
    first = 0
    text = (rf"first bin $k^*={k[first]:.0f}$ MeV/$c$: ALICE={data['C_alice_LL_R108'][first]:.3f}, "
            rf"doublet={data['C_doublet_chizzali_literal_ERE'][first]:.3f}, quartet={data['C_quartet_HAL'][first]:.3f}, "
            rf"spin avg={data['C_spinavg_chizzali_literal_ERE_HAL'][first]:.3f}")
    ax.annotate(text, xy=(0.04, 0.06), xycoords="axes fraction", fontsize=8.2)
    add_caption(fig, f"The p-phi doublet uses literal Chizzali ERE parameters, not the phi-alpha folding Gaussian. cf-fitter ALICE reproduction remains max|Delta C|={dmax:.2g}; the microscopic spin average still overshoots the empirical ALICE LL curve.")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(ROOT / "figures" / "fig02_phi_p_spin_average_diagnostic.pdf")
    fig.savefig(ROOT / "figures" / "fig02_phi_p_spin_average_diagnostic.png", dpi=180)
    if pdf is not None:
        pdf.savefig(fig)
    plt.close(fig)


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
    fig.suptitle(r"Figure 3. $\phi\alpha$ folded-potential scenario matrix", y=0.98)
    add_caption(fig, "Solid lines are KP central curves; orange shading is the doublet range band from b_d=0.4-0.7 fm. Dashed lines are LL for folded cells. The doublet volume integral is now constrained to the phi-N range, removing the previous near-threshold pathology.")
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    pdf.savefig(fig)
    plt.close(fig)

def figure_etminan(pdf: PdfPages):
    data = read_csv(PHIALPHA)
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.5), sharey=True)
    note = "Literal digitization of Etminan's published curve: pending. Quartet-only/fold is the Etminan configuration regenerated in this framework; orange band scans the doublet range."

    et_data = None
    et_col = None
    if ETMINAN.exists():
        et_data = read_csv(ETMINAN)
        et_col = next((c for c in et_data if c != "k_MeV"), None)
        note = "Etminan reference CSV was found and plotted where a matching/common curve column is available."

    for ax, R in zip(axes, [1.0, 3.0, 5.0]):
        rs = f"R{R:.2f}"
        if et_data and et_col:
            ax.plot(et_data["k_MeV"], et_data[et_col], color="0.25", lw=2.0, label="Etminan published reference")
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
    fig.suptitle(r"Figure 4. $\phi\alpha$: Etminan configuration vs doublet increment", y=0.98)
    add_caption(fig, note)
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    pdf.savefig(fig)
    plt.close(fig)

def folded_terms(phi_terms, rms=1.56, A=4.0, fold=True):
    br = rms * math.sqrt(2.0 / 3.0)
    rho0 = A / (math.pi ** 1.5 * br ** 3)
    if not fold:
        G = sum(a * math.pi ** 1.5 * b ** 3 for a, b in phi_terms)
        return [(G * rho0, br)]
    terms = []
    for a, b in phi_terms:
        B2 = b * b + br * br
        amp = a * rho0 * (math.pi * b * b * br * br / B2) ** 1.5
        terms.append((amp, math.sqrt(B2)))
    return terms


def volume_integral(terms):
    return sum(a * math.pi ** 1.5 * b ** 3 for a, b in terms)


def figure_table(pdf: PdfPages):
    q = [(-371.0, 0.15), (-50.0, 0.66), (-31.0, 1.09)]
    d = [(-571.5, 0.55)]
    qd = [(2.0 / 3.0 * a, b) for a, b in q] + [(1.0 / 3.0 * a, b) for a, b in d]
    rows = []
    params = {
        "q_nofold": (-1.916, 0.905, "17.34"),
        "q_fold": (-3.029, 1.380, "6.29"),
        "qd_nofold": (-1.531, 0.669, "22.59"),
        "qd_fold": (-2.111, 0.979, "13.53"),
    }
    labels = {
        "q_nofold": "quartet only / no fold",
        "q_fold": "quartet only / fold",
        "qd_nofold": "quartet+doublet / no fold",
        "qd_fold": "quartet+doublet / fold",
    }
    for key, phi_terms, fold in [
        ("q_nofold", q, False),
        ("q_fold", q, True),
        ("qd_nofold", qd, False),
        ("qd_fold", qd, True),
    ]:
        terms = folded_terms(phi_terms, fold=fold)
        f0, d0, be = params[key]
        rows.append([labels[key], f"{f0:+.3f}", f"{d0:.3f}", be, f"{sum(a for a, _ in terms):.1f}", f"{volume_integral(terms):.1f}"])

    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.axis("off")
    ax.set_title(r"Figure 5. $\phi\alpha$ scenario summary", fontsize=14, pad=16)
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


def main():
    required = [PHIP_REF / "phi_p_ALICE_LL_reference.csv", PHIP_REF / "generated_ALICEfig_LL.csv", PHIALPHA]
    missing = [str(p.relative_to(ROOT)) for p in required if not p.exists()]
    if missing:
        raise SystemExit("Missing required input CSVs:\n" + "\n".join(missing))
    OUT.parent.mkdir(exist_ok=True)
    with PdfPages(OUT) as pdf:
        figure_phi_p_baseline(pdf)
        figure_phi_p_diagnostic(pdf)
        figure_phi_alpha_matrix(pdf)
        figure_etminan(pdf)
        figure_table(pdf)
    legacy = ROOT / "figures" / "cf_fitter_summary.pdf"
    legacy.write_bytes(OUT.read_bytes())
    print(f"wrote {OUT.relative_to(ROOT)}")
    print(f"wrote {legacy.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
