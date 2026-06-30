# phi-p literature comparison package (ALICE PRL 127, 172301 / arXiv:2105.05578)

Reference data for validating the cf-fitter phi-p correlation function against the
published ALICE measurement. Hand this folder to the remote agent together with the
planning-document prompt.

## Files
| file | role |
|------|------|
| `phi_p_ALICE_params.yaml` | All extracted ALICE numbers, conventions, provenance, and two drop-in scenarios. |
| `phi_p_ALICE_LL_reference.csv` | **Comparison (2)** target: ALICE's headline Lednicky-Lyuboshits curve, computed from the published f0, d0 at the ALICE source size. This *is* the Fig. 2 red band. |
| `phi_p_ALICE_gauss_KP_reference.csv` | Independent cross-check: Koonin-Pratt curve from ALICE's fitted Gaussian potential (exercises the Numerov path, not just LL). |
| `phi_p_ALICE_Fig2_data_template.csv` | **Comparison (1)** target: empty template for the measured Fig. 2 data points (must be digitized; not on HEPData). |

## How the two comparisons work

**(2) Reproduce-from-parameters (primary, exact).**
Feed cf-fitter the published ALICE inputs and check the output lands on the reference curve:
- f0 = 0.85 + 0.16i fm, d0 = 7.85 fm, Gaussian source r_eff = 1.08 fm, lambda = 1 (genuine CF).
- The reference columns (`C_LL_R1.08`, plus 1.03/1.13 for the +-0.05 fm source band, and an
  `Imf0_0` elastic variant) were generated with the *same* LL formula as
  `include/femto/FemtoModels.hpp` and verified against `cf-calc` to ~1e-7. So the agent's
  `cf-calc make-cf` (corrected ALICEfig scenario, R = 1.08) must reproduce `C_LL_R1.08`
  to numerical precision. Success criterion: max |Delta C| < 1e-3 over k* in [0, 200] MeV/c.

**(1) Compare-to-data (digitized).**
Overlay the agent's computed curve on the measured Fig. 2 points. Because ALICE did not
deposit this CF on HEPData, the points must be digitized from Fig. 2 (template provided).
Report a chi2/point of the model vs the digitized points within k* < 200 MeV/c.

## Physics caveats (get these right or the overlay is meaningless)
- **k\*** convention matches the framework (k* = half the relative momentum, pair rest frame).
  No factor-2 / q = 2k* rescaling.
- **Source size is fixed to ALICE's value R = 1.08 fm.** Do not use the repo's default R = 1.2/3.0.
- **lambda = 1** for the genuine CF (Fig. 2 already has the lambda/feed-down/purity divided out).
- The repo's current `phi_proton_ALICEeff.yaml` sets **d0 = 0.0, which is wrong**; ALICE's
  d0 = 7.85 fm. The corrected scenario is in `phi_p_ALICE_params.yaml` (ALICEfig).
- **Large d0 breaks the LL approximation at small R.** With d0 = 7.85-15 fm and R ~ 1.08 fm,
  the LL term [1 - d0/(2 sqrt(pi) R)] can go negative and LL diverges from the exact KP
  result. The KP (Numerov) curve is the trustworthy one in this regime; ALICE itself notes
  the d0/r0 asymptotic correction matters here. Expect LL(ERE) and KP(potential) to agree
  only roughly at low k*; this is physics, not a bug.
- **Im(f0) is consistent with zero.** Keep the 0.16i central value but treat the elastic
  (`Imf0_0`) variant as an equally valid reference band.
