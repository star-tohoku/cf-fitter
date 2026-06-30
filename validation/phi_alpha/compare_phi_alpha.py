#!/usr/bin/env python3
"""Compare regenerated phi-alpha curves against checked-in reference points."""
import argparse
import csv
from pathlib import Path
from typing import Dict, List


def read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as f:
        lines = [
            line for line in f
            if line.strip() and not line.lstrip().startswith("#")
        ]
    return list(csv.DictReader(lines))


def nearest(rows: List[Dict[str, str]], k_mev: float) -> Dict[str, str]:
    return min(rows, key=lambda row: abs(float(row["k_MeV"]) - k_mev))


def main() -> int:
    root = Path(__file__).resolve().parent
    repo = root.parents[1]

    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--folded-csv",
        default=str(repo / "cf_phialpha_folded.csv"),
        help="regenerated folded CSV from cf-calc",
    )
    ap.add_argument(
        "--reference",
        default=str(root / "reference_points.csv"),
        help="checked-in reference points",
    )
    ap.add_argument("--kp-tol", type=float, default=5e-3)
    ap.add_argument("--ll-tol", type=float, default=5e-3)
    args = ap.parse_args()

    folded_path = Path(args.folded_csv)
    ref_path = Path(args.reference)
    if not folded_path.exists():
        print(f"error: missing folded CSV: {folded_path}")
        return 1
    if not ref_path.exists():
        print(f"error: missing reference CSV: {ref_path}")
        return 1

    got_rows = read_rows(folded_path)
    ref_rows = read_rows(ref_path)

    total_fail = 0
    for ref in ref_rows:
        scenario = ref["scenario"]
        R = float(ref["R_fm"])
        rs = f"R{R:.2f}"
        kp_col = f"KP_{scenario}_{rs}"
        ll_col = f"LL_{scenario}_{rs}"
        if kp_col not in got_rows[0]:
            print(f"error: folded CSV missing column {kp_col}")
            return 1

        k = float(ref["k_MeV"])
        got = nearest(got_rows, k)
        kp_delta = abs(float(got[kp_col]) - float(ref["C_KP"]))
        ll_delta = abs(float(got[ll_col]) - float(ref["C_LL"]))
        print(
            f"{scenario} R={R:g} fm k*={k:g} MeV/c: "
            f"|dKP|={kp_delta:.6g} (tol {args.kp_tol:g}), "
            f"|dLL|={ll_delta:.6g} (tol {args.ll_tol:g})"
        )
        if kp_delta > args.kp_tol or ll_delta > args.ll_tol:
            total_fail += 1

    if total_fail:
        print(f"compare_phi_alpha: {total_fail} point(s) outside tolerance")
        return 1
    print("compare_phi_alpha: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
