
---

## cf-fitter package (YAML configuration)

The standalone package at `/star/u/oura/gpfs/others/cf-fitter` extends this prototype with:

- `config/channels/*.yaml` — pair metadata and spin-channel weights
- `config/scenarios/*.yaml` — Gaussian potentials, ERE presets, or `ll_only` effective channels
- `cf-calc validate` / `cf-calc make-cf` — batch theory curves
- `cf-fit --mode ll|kp|both` — ROOT fits with optional overlay

See the package `README.md` and `docs/ROADMAP.md` for build instructions and future folding/coupling plans.
