#!/usr/bin/env python3
"""Build phi-alpha ERE / potential summary for doublet-shape scenarios.

Uses cf-calc phi-alpha-summary (folded HAL/Chizzali-motivated scenarios) and
cf-calc make-cf (toy Gaussian Vcrit scan reference). Does not retune potentials.
"""
from __future__ import annotations

import csv
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "build" / "cf-calc"
SUMMARY_IN = ROOT / "figures" / "phi_alpha_potential_summary.csv"
OUT_CSV = ROOT / "figures" / "phi_alpha_current_doublet_ere_summary.csv"
OUT_MD = ROOT / "PHI_ALPHA_CURRENT_DOUBLET_ERE_SUMMARY.md"

PRIMARY_FOLDED = [
    "q_fold",
    "qd_b040_fold",
    "qd_fold",
    "qd_b070_fold",
]
DIAGNOSTIC_NOFOLD = [
    "q_nofold",
    "qd_b040_nofold",
    "qd_nofold",
    "qd_b070_nofold",
]
TOY_SCENARIOS = ["weak", "strong", "bound"]

# Documented toy Gaussian depths (b=2 fm, Vcrit scan); matches docs/writeup.md.
TOY_V0_MEV = {"weak": 6.5, "strong": 15.0, "bound": 21.9}
TOY_B_FM = 2.0
MU_REDUCED_MEV = 800.5

FIELDNAMES = [
    "scenario",
    "scenario_group",
    "fold",
    "doublet_b_fm",
    "V0_MeV",
    "f0_fa_fm",
    "d0_fa_fm",
    "BE_MeV",
    "V_at_0_MeV",
    "volume_integral_MeV_fm3",
    "mu_reduced_MeV",
    "notes",
]


def run_cf_calc(*args: str) -> str:
    cmd = [str(BUILD), *args]
    env = dict(**{k: v for k, v in __import__("os").environ.items() if k != "LD_LIBRARY_PATH"})
    proc = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True, check=True)
    return proc.stdout + proc.stderr


def refresh_potential_summary() -> None:
    if not BUILD.is_file():
        raise SystemExit(f"missing {BUILD}; build cf-calc first")
    run_cf_calc(
        "phi-alpha-summary",
        "--output",
        str(SUMMARY_IN),
    )


def read_potential_summary() -> dict[str, dict[str, str]]:
    if not SUMMARY_IN.is_file():
        refresh_potential_summary()
    with SUMMARY_IN.open(newline="") as f:
        rows = list(csv.DictReader(f))
    return {row["scenario"]: row for row in rows}


def toy_gaussian_volume_integral(v0_mev: float, b_fm: float) -> float:
    import math

    return -v0_mev * (math.pi ** 1.5) * (b_fm ** 3)


def parse_toy_make_cf() -> dict[str, dict[str, float | str]]:
    text = run_cf_calc(
        "make-cf",
        "--channel",
        "phi_alpha",
        "--scenarios",
        ",".join(TOY_SCENARIOS),
        "--radii",
        "3.0",
    )
    pat = re.compile(
        r"scenario\s+(\w+)\s+mu=[\d.]+\s+MeV\s+\[\w+\]\s+"
        r"f0=([+-][\d.]+)[+-][\d.]+i\s+d0=([\d.]+)"
        r"(?:\s+BE~([\d.]+)\s+MeV)?"
    )
    out: dict[str, dict[str, float | str]] = {}
    for m in pat.finditer(text):
        name = m.group(1)
        f0 = float(m.group(2))
        d0 = float(m.group(3))
        be = float(m.group(4)) if m.group(4) else ""
        v0 = TOY_V0_MEV[name]
        out[name] = {
            "f0_fa_fm": f0,
            "d0_fa_fm": d0,
            "BE_MeV": be,
            "V_at_0_MeV": -v0,
            "volume_integral_MeV_fm3": toy_gaussian_volume_integral(v0, TOY_B_FM),
            "V0_MeV": v0,
        }
    if set(out) != set(TOY_SCENARIOS):
        raise RuntimeError(f"failed to parse toy scenarios from make-cf output:\n{text}")
    return out


def row_from_summary(
    scenario: str,
    group: str,
    src: dict[str, str],
    notes: str = "",
) -> dict[str, str | float]:
    be = src.get("BE_MeV", "").strip()
    return {
        "scenario": scenario,
        "scenario_group": group,
        "fold": src["fold"],
        "doublet_b_fm": src.get("doublet_b_fm", ""),
        "V0_MeV": src.get("V0_MeV", ""),
        "f0_fa_fm": float(src["f0_fm"]),
        "d0_fa_fm": float(src["d0_fm"]),
        "BE_MeV": float(be) if be else "",
        "V_at_0_MeV": float(src["V_at_0_MeV"]),
        "volume_integral_MeV_fm3": float(src["volume_integral_MeV_fm3"]),
        "mu_reduced_MeV": MU_REDUCED_MEV,
        "notes": notes or src.get("notes", ""),
    }


def row_from_toy(name: str, data: dict[str, float | str]) -> dict[str, str | float]:
    labels = {
        "weak": "toy Gaussian scan: 0.40 Vcrit, unbound",
        "strong": "toy Gaussian scan: 0.92 Vcrit, near-unitary unbound",
        "bound": "toy Gaussian scan: 1.35 Vcrit, shallow bound",
    }
    return {
        "scenario": name,
        "scenario_group": "toy_gaussian_scan",
        "fold": "n/a",
        "doublet_b_fm": "",
        "V0_MeV": data["V0_MeV"],
        "f0_fa_fm": data["f0_fa_fm"],
        "d0_fa_fm": data["d0_fa_fm"],
        "BE_MeV": data["BE_MeV"],
        "V_at_0_MeV": data["V_at_0_MeV"],
        "volume_integral_MeV_fm3": round(float(data["volume_integral_MeV_fm3"]), 1),
        "mu_reduced_MeV": MU_REDUCED_MEV,
        "notes": labels[name],
    }


def write_csv(rows: list[dict[str, str | float]]) -> None:
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def fmt(v: object, nd: int = 3) -> str:
    if v == "" or v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.{nd}f}"
    return str(v)


def write_markdown(rows: list[dict[str, str | float]]) -> None:
    folded = [r for r in rows if r["scenario_group"] == "current_folded"]
    nofold = [r for r in rows if r["scenario_group"] == "current_nofold"]
    toy = [r for r in rows if r["scenario_group"] == "toy_gaussian_scan"]

    def table(rs: list[dict[str, str | float]]) -> list[str]:
        lines = [
            "| scenario | f₀ [fm] | d₀ [fm] | BE [MeV] | V(0) [MeV] | ∫V d³r [MeV fm³] |",
            "|----------|---------|---------|----------|------------|-------------------|",
        ]
        for r in rs:
            lines.append(
                f"| `{r['scenario']}` | {fmt(r['f0_fa_fm'])} | {fmt(r['d0_fa_fm'])} "
                f"| {fmt(r['BE_MeV'], 2)} | {fmt(r['V_at_0_MeV'], 1)} "
                f"| {fmt(r['volume_integral_MeV_fm3'], 1)} |"
            )
        return lines

    f0_fold = [float(r["f0_fa_fm"]) for r in folded]
    be_fold = [float(r["BE_MeV"]) for r in folded if r["BE_MeV"] != ""]
    toy_by = {r["scenario"]: r for r in toy}

    text = f"""# Phi-alpha current doublet ERE summary

Date: 2026-07-01

## Purpose

Low-energy and potential-shape summary for the **current HAL/Chizzali-motivated φα folded scenarios** used in the doublet-shape uncertainty study. Compared against the older **toy single-Gaussian φα strength scan** (`weak`, `strong`, `bound`).

Femtoscopy convention: **f₀ > 0** = attractive and **unbound**; **f₀ < 0** = **bound** state exists. This is **not** a unique Chizzali doublet reproduction — the folded q+d family propagates Gaussian doublet-shape uncertainty at fixed φ-N scattering length.

Reduced mass μ = **{MU_REDUCED_MEV} MeV** (same as φα KP/Numerov).

## Current folded scenarios (primary)

{chr(10).join(table(folded))}

Doublet depths **V₀** and ranges **b_d** come from `figures/phi_alpha_potential_summary.csv` (via `cf-calc phi-alpha-summary`).

## Diagnostic no-fold scenarios

{chr(10).join(table(nofold))}

These use the same φ-N central potentials without α folding; **f₀** and **BE** are not directly comparable to the folded φα correlation curves.

## Toy Gaussian φα scan (reference)

| scenario | V₀ [MeV] | b [fm] | f₀ [fm] | d₀ [fm] | BE [MeV] |
|----------|----------|--------|---------|---------|----------|
| `weak` | {fmt(toy_by['weak']['V0_MeV'], 1)} | {TOY_B_FM} | {fmt(toy_by['weak']['f0_fa_fm'])} | {fmt(toy_by['weak']['d0_fa_fm'])} | — |
| `strong` (near-unitary) | {fmt(toy_by['strong']['V0_MeV'], 1)} | {TOY_B_FM} | {fmt(toy_by['strong']['f0_fa_fm'])} | {fmt(toy_by['strong']['d0_fa_fm'])} | — |
| `bound` | {fmt(toy_by['bound']['V0_MeV'], 1)} | {TOY_B_FM} | {fmt(toy_by['bound']['f0_fa_fm'])} | {fmt(toy_by['bound']['d0_fa_fm'])} | {fmt(toy_by['bound']['BE_MeV'], 2)} |

Values from `cf-calc make-cf --channel phi_alpha --scenarios weak,strong,bound` (Numerov on the toy Gaussian potential).

## Comparison: pole region vs current folded q+d

**Main question:** Are the current folded φα q+d scenarios close to the toy near-unitary / shallow-bound pole region?

**Answer: No — they are far from the pole region.**

1. **Sign / binding sector.** All four folded scenarios have **f₀ < 0** with **BE ≈ {min(be_fold):.2f}–{max(be_fold):.2f} MeV**. The toy near-unitary case (`strong`) has **f₀ = +{toy_by['strong']['f0_fa_fm']} fm** (still **unbound**). The current folded family sits firmly on the **bound** side, not straddling the threshold.

2. **Distance from unitarity.** The toy scan brackets the pole with **f₀ = +{toy_by['strong']['f0_fa_fm']} fm** (large positive, unbound) and **f₀ = {toy_by['bound']['f0_fa_fm']} fm** (shallow bound, BE ≈ {fmt(toy_by['bound']['BE_MeV'], 2)} MeV). Current folded **f₀** spans **{min(f0_fold):.3f} to {max(f0_fold):.3f} fm** — much smaller |f₀| and **much deeper binding** than the toy shallow-bound point.

3. **Effective range.** Folded **d₀ ≈ 0.73–1.38 fm** (short, folded φ-N structure). Toy scan uses **d₀ ≈ 2.4–5.5 fm** on a single **b = 2 fm** Gaussian — a different potential class altogether.

4. **Interpretation for C(k*).** The toy scan explored correlation enhancement/depletion near **unitarity** (large |f₀|, threshold BE). The current folded HAL/Chizzali-motivated scenarios are **moderately deep bound states** with **modest |f₀| ~ 2–3 fm**. The doublet-shape band in C(k*) is therefore a **propagated input-shape uncertainty within a bound-state regime**, not sensitivity near the φα scattering pole.

## Files and commands

```bash
build/cf-calc phi-alpha-summary
python3 figures/make_phi_alpha_current_doublet_ere_summary.py
```

Outputs:

- `{OUT_CSV.relative_to(ROOT)}` — machine-readable table
- `{OUT_MD.name}` — this note

Source of truth for folded potentials: `cf-calc phi-alpha-summary` → `figures/phi_alpha_potential_summary.csv` (same potentials as `cf_phialpha_folded.csv` KP curves).
"""
    OUT_MD.write_text(text)


def main() -> int:
    refresh_potential_summary()
    summary = read_potential_summary()
    rows: list[dict[str, str | float]] = []

    for name in PRIMARY_FOLDED:
        rows.append(row_from_summary(name, "current_folded", summary[name]))
    for name in DIAGNOSTIC_NOFOLD:
        rows.append(
            row_from_summary(
                name,
                "current_nofold",
                summary[name],
                notes="bare phi-N central; no alpha fold",
            )
        )
    toy = parse_toy_make_cf()
    for name in TOY_SCENARIOS:
        rows.append(row_from_toy(name, toy[name]))

    write_csv(rows)
    write_markdown(rows)
    print(f"wrote {OUT_CSV}")
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
