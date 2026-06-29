#!/bin/bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
BUILD_DIR="${ROOT_DIR}/build"
export LD_LIBRARY_PATH="${BUILD_DIR}:${LD_LIBRARY_PATH:-}"
"${BUILD_DIR}/cf-fit" --demo --mode both --channel phi_proton --scenario HAL "$@"
