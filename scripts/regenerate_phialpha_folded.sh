#!/usr/bin/env bash
# Regenerate cf_phialpha_folded.csv from the eight primary folded phi-alpha scenarios.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="${CF_FITTER_BUILD:-$ROOT/build}"
CF_CALC="$BUILD/cf-calc"
OUT="$ROOT/cf_phialpha_folded.csv"

SCENARIOS="q_nofold,q_fold,qd_nofold,qd_fold,qd_b040_nofold,qd_b040_fold,qd_b070_nofold,qd_b070_fold"

if [[ ! -x "$CF_CALC" ]]; then
  echo "error: cf-calc not found at $CF_CALC" >&2
  echo "Build first: mkdir -p build && cd build && cmake .. && cmake --build ." >&2
  exit 1
fi

CMD=(
  "$CF_CALC" make-cf
  --channel phi_alpha
  --scenarios "$SCENARIOS"
  --radii 1,3,5
  --output "$OUT"
)

printf 'running: %s\n' "${CMD[*]}"
cd "$ROOT"
"${CMD[@]}"

echo "wrote $OUT"
