#!/usr/bin/env bash
# Regenerate phi-alpha pipeline artifacts (folded CSV + potential summary CSV).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="${CF_FITTER_BUILD:-$ROOT/build}"
CF_CALC="$BUILD/cf-calc"

if [[ ! -x "$CF_CALC" ]]; then
  echo "error: cf-calc not found at $CF_CALC" >&2
  exit 1
fi

"$ROOT/scripts/regenerate_phialpha_folded.sh"

SUMMARY_CMD=(
  "$CF_CALC" phi-alpha-summary
  --output "$ROOT/figures/phi_alpha_potential_summary.csv"
)
printf 'running: %s\n' "${SUMMARY_CMD[*]}"
cd "$ROOT"
"${SUMMARY_CMD[@]}"

echo "phi-alpha pipeline artifacts regenerated"
