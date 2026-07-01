#!/usr/bin/env bash
# Regenerate Chizzali TPE phi-alpha quartet comparison artifacts.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="${CF_FITTER_BUILD:-$ROOT/build}"
CF_CALC="$BUILD/cf-calc"

if [[ ! -x "$CF_CALC" ]]; then
  echo "error: cf-calc not found at $CF_CALC" >&2
  exit 1
fi

CURVES_CMD=(
  "$CF_CALC" make-cf
  --channel phi_alpha
  --scenarios q_fold,q_chizzali_tpe_nofold,q_chizzali_tpe_fold
  --radii 1,3,5
  --output "$ROOT/validation/phi_alpha/generated_q_chizzali_tpe.csv"
)
printf 'running: %s\n' "${CURVES_CMD[*]}"
cd "$ROOT"
"${CURVES_CMD[@]}"

SUMMARY_CMD=(
  "$CF_CALC" phi-alpha-summary
  --scenarios q_chizzali_tpe_nofold,q_chizzali_tpe_fold
  --output "$ROOT/figures/phi_alpha_chizzali_tpe_potential_summary.csv"
)
printf 'running: %s\n' "${SUMMARY_CMD[*]}"
"${SUMMARY_CMD[@]}"

FIG_CMD=(python3 "$ROOT/figures/make_phi_alpha_chizzali_tpe_figure.py")
printf 'running: %s\n' "${FIG_CMD[*]}"
"${FIG_CMD[@]}"

echo "Chizzali TPE phi-alpha artifacts regenerated"
