# cf-fitter roadmap

## Implemented (Phase 1–4)

- YAML channel/scenario configuration (φ–p, φ–α)
- `cf-calc validate` — square-well solver benchmark + LL/KP cross-check
- `cf-calc make-cf` — CSV generation
- `libCfFitter.so` + `fitCf.C` — fit modes `ll`, `kp`, `both`
- Extension hooks: `PotentialBuilder::folded` (stub), `CouplingScheme` enum

## Phase 5 — Potential folding

\[
V_{\phi\alpha}(r) = \int d^3r'\, V_{\phi N}(|\mathbf{r}-\mathbf{r}'|)\,\rho_\alpha(\mathbf{r}')
\]

- Implement `PotentialBuilder::folded` in `src/FoldingPotential.cpp`
- Add `config/densities/` (e.g. α ground-state density)
- Scenario YAML: `potential.type: folded`
- Validation against thin-shell or known limits

## Phase 6 — Channel coupling

| Stage | Description |
|-------|-------------|
| **6a** | Effective single channel + complex \(f_0\) (partially via `coupling: effective`) |
| **6b** | `CoupledRadialSolver` — 2×2 real potential matrix, shared \(\mu\) |
| **6c** | Open channels / complex potentials / R-matrix |

Current default: `coupling: independent` (incoherent sum over spin channels).

## Other future work

- Coulomb for charged pairs
- Momentum-resolution smearing in fit (`smearModel` utility)
- star-analyzer integration via git submodule + `TH1D` fit macro
- Multi-channel LL fit with all \(f_{0,i}, d_{0,i}\) free (φ–p Set B full fit)
