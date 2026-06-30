#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path


def read_csv(path):
    with open(path, newline='') as f:
        lines = [line for line in f if not line.lstrip().startswith('#') and not line.lstrip().startswith('\"#')]
    return list(csv.DictReader(lines))


def nearest(rows, k):
    return min(rows, key=lambda r: abs(float(r['k_MeV']) - k))


def max_delta(got_path, got_col, ref_path, ref_col, kmax):
    got = read_csv(got_path)
    ref = read_csv(ref_path)
    worst = (0.0, 0.0, 0.0, 0.0)
    for rr in ref:
        k = float(rr['k_MeV'])
        if k > kmax:
            continue
        gr = nearest(got, k)
        delta = abs(float(gr[got_col]) - float(rr[ref_col]))
        if delta > worst[0]:
            worst = (delta, k, float(gr[got_col]), float(rr[ref_col]))
    return worst


def fig2_status(path):
    rows = read_csv(path)
    filled = []
    for r in rows:
        values = [v.strip() for v in r.values() if v is not None]
        if any(values):
            filled.append(r)
    if not filled:
        return 'Fig.2 digitized-data comparison: pending (template is empty)'
    return f'Fig.2 digitized-data comparison: {len(filled)} points present; chi2 harness not yet wired to model columns'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ll-csv')
    ap.add_argument('--kp-csv')
    ap.add_argument('--root', default=str(Path(__file__).resolve().parent))
    ap.add_argument('--kmax', type=float, default=200.0)
    args = ap.parse_args()
    root = Path(args.root)
    if args.ll_csv:
        d, k, got, ref = max_delta(args.ll_csv, 'LL_ALICEfig_R1.08', root / 'phi_p_ALICE_LL_reference.csv', 'C_LL_R1.08', args.kmax)
        print(f'LL ALICE parameter reproduction: max|Delta C|={d:.6g} at k={k:.2f} MeV/c (got={got:.6f}, ref={ref:.6f})')
    if args.kp_csv:
        d, k, got, ref = max_delta(args.kp_csv, 'KP_ALICEgauss_R1.08', root / 'phi_p_ALICE_gauss_KP_reference.csv', 'C_KP_R1.08', args.kmax)
        print(f'KP ALICE Gaussian reproduction: max|Delta C|={d:.6g} at k={k:.2f} MeV/c (got={got:.6f}, ref={ref:.6f})')
    print(fig2_status(root / 'phi_p_ALICE_Fig2_data_template.csv'))


if __name__ == '__main__':
    main()
