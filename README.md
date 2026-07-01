# cf-fitter

Femtoscopic correlation function \(C(k^*)\) framework: **Lednický–Lyuboshitz (LL)** and **Koonin–Pratt (KP)** models with YAML-driven channels/scenarios, calc tools, and ROOT fitting.

**Location:** `/star/u/oura/gpfs/others/cf-fitter`

## Quick start

```bash
cd /star/u/oura/gpfs/others/cf-fitter
mkdir -p build && cd build
cmake ..
cmake --build .
./cf-calc validate
./cf-calc make-cf --channel phi_proton
./cf-calc make-cf --channel phi_alpha
```

`cf-calc validate` runs the square-well Numerov benchmark and the LL/KP consistency benchmark.

## Fit (ROOT)

```bash
cd build
export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH
./cf-fit --demo --mode both --channel phi_proton --scenario HAL
./cf-fit --demo --mode ll    # LL only
./cf-fit --demo --mode kp    # KP only
./cf-fit --input data.root --hist hCF --mode both --channel phi_proton --scenario HAL
```

Or: `scripts/run_fit_demo.sh`

Legacy ROOT macro `macro/fitCf.C` shells out to `cf-fit` (set `CF_FITTER_BUILD=build`).

## Conventions

- Femto ERE: \(k^*\cot\delta = 1/f_0 + d_0 k^{*2}/2\), with **\(f_0 = -a_0\)** (nuclear scattering length).
- Units: fm, MeV/\(c\) on histogram axis (internally converted with \(\hbar c = 197.327\) MeV·fm).
- Fit function: \(C_{\rm fit} = N(1+b_1 k^*)[1+\lambda(C_{\rm model}-1)]\).

## Layout

| Path | Role |
|------|------|
| `include/femto/FemtoModels.hpp` | Core LL/KP physics (ROOT-free) |
| `config/channels/` | Pair channel definitions (φ–p, φ–α, …) |
| `config/scenarios/` | Interaction presets |
| `cf-calc` | `validate` (square-well + LL/KP), `make-cf` |
| `libCfFitter.so` | ROOT fit API (`FitSession`) |
| `cf-fit` | CLI fitter (`--mode ll|kp|both`) |
| `macro/fitCf.C` | Thin wrapper calling `cf-fit` |

## Adding a channel

1. Add `config/channels/<pair>.yaml` (masses, spin channels, weights).
2. Add `config/scenarios/<pair>_<name>.yaml` (potentials or ERE per spin).
3. Run `cf-calc make-cf --channel <pair>`.

See [docs/ROADMAP.md](docs/ROADMAP.md) for folding (Phase 5) and coupling (Phase 6).

## Documentation

- [docs/writeup.md](docs/writeup.md) — theory and validation (from original prototype)
- [docs/ROADMAP.md](docs/ROADMAP.md) — future extensions
