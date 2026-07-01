#!/usr/bin/env python3
"""Axis-calibrated digitization of Etminan Fig. 4 HAL QCD rms=1.56 fm KP curves.

Source: arXiv:2410.22756, Fig. 4 panels (a,c,e) left column (HAL QCD).
Curve: green dashed line, rms=1.56 fm, Koonin-Pratt (KP) formula.
"""
from __future__ import annotations

import csv
import statistics
from collections import defaultdict
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
FIG4_HIRES = Path(__file__).resolve().parent / "fig4_hires.png"
OUT_CSV = ROOT / "validation" / "etminan_phialpha_reference.csv"

GREEN = (143, 237, 143)


def _near_green(r: int, g: int, b: int) -> bool:
    return (
        abs(r - GREEN[0]) <= 10
        and abs(g - GREEN[1]) <= 15
        and abs(b - GREEN[2]) <= 10
    )


def _extract_curve(
    px,
    *,
    xl: int,
    xr: int,
    yt: int,
    yb: int,
    c_min: float,
    c_max: float,
    k_max: float = 100.0,
) -> dict[int, float]:
    bins: dict[int, list[int]] = defaultdict(list)
    for y in range(yt, yb + 1):
        for x in range(xl, xr + 1):
            if _near_green(*px[x, y]):
                k = round((x - xl) / (xr - xl) * k_max)
                bins[k].append(y)
    out: dict[int, float] = {}
    for k, ys in sorted(bins.items()):
        y = statistics.median(ys)
        out[k] = c_max - (y - yt) / (yb - yt) * (c_max - c_min)
    return out


def _interp(k: float, curve: dict[int, float], *, k_hi_cap: int | None, plateau: float) -> float:
    if not curve:
        return plateau
    ks = sorted(curve)
    if k_hi_cap is not None and k >= k_hi_cap:
        return plateau
    if k <= ks[0]:
        return curve[ks[0]]
    if k >= ks[-1]:
        return curve[ks[-1]]
    lo = max(kk for kk in ks if kk <= k)
    hi = min(kk for kk in ks if kk >= k)
    if lo == hi:
        return curve[lo]
    t = (k - lo) / (hi - lo)
    return curve[lo] * (1.0 - t) + curve[hi] * t


def main() -> None:
    if not FIG4_HIRES.exists():
        raise SystemExit(f"Missing {FIG4_HIRES}; run PDF extract first.")

    img = Image.open(FIG4_HIRES).convert("RGB")
    px = img.load()
    xl, xr = 460, 1119

    curves = {
        "C_R1.00": _extract_curve(
            px, xl=xl, xr=xr, yt=370, yb=871, c_min=0.0, c_max=1.5
        ),
        "C_R3.00": _extract_curve(
            px, xl=xl, xr=xr, yt=1150, yb=1547, c_min=0.2, c_max=1.1
        ),
        "C_R5.00": _extract_curve(
            px, xl=xl, xr=xr, yt=1880, yb=2270, c_min=0.5, c_max=1.1
        ),
    }

    k_grid = (
        [0.0, 1.0, 2.0, 3.0, 5.0]
        + [float(k) for k in range(10, 101, 5)]
        + [110.0, 120.0, 140.0, 160.0, 180.0, 200.0]
    )

    rows: list[dict[str, float]] = []
    for k in k_grid:
        row: dict[str, float] = {"k_MeV": k}
        row["C_R1.00"] = min(_interp(k, curves["C_R1.00"], k_hi_cap=None, plateau=1.0), 1.15)
        row["C_R3.00"] = min(_interp(k, curves["C_R3.00"], k_hi_cap=48, plateau=1.0), 1.0)
        row["C_R5.00"] = min(_interp(k, curves["C_R5.00"], k_hi_cap=20, plateau=1.0), 1.0)
        rows.append(row)

    header = [
        "# Etminan phi-alpha published curve digitization",
        "# Source: F. Etminan, arXiv:2410.22756 (2024), Fig. 4 panels (a,c,e) HAL QCD column",
        "# Curve: HAL QCD single-folding, alpha rms=1.56 fm, KP formula (green dashed line)",
        "# Method: axis-calibrated color-specific extraction from fig4_hires.png (Tier A/B hybrid)",
        "# Method_in_paper: KP (Koonin-Pratt); compared to LL markers in same figure",
        "# R_source_fm: 1, 3, 5 (Gaussian S(r) per Eq. 5 in paper)",
        "# Caveat: qualitative overlay only; HAL Filikhin t/a=14 WS params, CATS/KP in paper",
        "#         vs HAL fit-B three-Gaussian + analytic fold in cf-fitter q_fold",
        "k_MeV,C_R1.00,C_R3.00,C_R5.00",
    ]

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="") as f:
        for line in header:
            f.write(line + "\n")
        writer = csv.DictWriter(
            f, fieldnames=["k_MeV", "C_R1.00", "C_R3.00", "C_R5.00"]
        )
        for row in rows:
            writer.writerow(
                {
                    "k_MeV": f"{row['k_MeV']:.1f}",
                    "C_R1.00": f"{row['C_R1.00']:.4f}",
                    "C_R3.00": f"{row['C_R3.00']:.4f}",
                    "C_R5.00": f"{row['C_R5.00']:.4f}",
                }
            )

    print(f"Wrote {OUT_CSV} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
