# Etminan φα reference data

Digitized published φα correlation functions for PDF overlay against cf-fitter `q_fold`.

## Canonical file

```text
validation/etminan_phialpha_reference.csv
```

Plot scripts load this via `figures/etminan_overlay.py` (Figure 7 and Chizzali TPE last page).

## Source

| Field | Value |
|-------|-------|
| Citation | F. Etminan, *Exploring the φ-α interaction via femtoscopic study*, arXiv:2410.22756 (2024) |
| Figure | Fig. 4, left column (HAL QCD): panels (a) R=1 fm, (c) R=3 fm, (e) R=5 fm |
| Curve | Green dashed line — HAL QCD single-folding, **α rms = 1.56 fm**, **KP formula** |
| Source sizes | Gaussian S(r) with R = 1, 3, 5 fm (Eq. 5 in paper) |
| k range | 0–100 MeV/c digitized; 110–200 MeV/c set to plateau where curve saturates |
| Screenshot | `etminan_source_figure.png` (cropped HAL column from `fig4_hires.png`) |
| PDF | `etminan_arxiv_2410.22756.pdf` |

## Digitization

- **Tier:** A/B hybrid — axis-calibrated, color-specific extraction of RGB (143,237,143) from 300 DPI PDF render.
- **Script:** `digitize_etminan_fig4.py` (reproducible regeneration).
- **Not publication-grade:** dashed-line sparsity and y-axis saturation above C≈1.0 add ~0.01–0.03 uncertainty.

## Convention caveats (cf-fitter `q_fold` vs Etminan)

Etminan's HAL curve uses Wood–Saxon parameters from Filikhin et al. (HAL QCD t/a=14, quartet channel) folded with a Gaussian α density (rms = 1.56 fm), computed with **KP via CATS** in Fig. 3 and **KP lines** in Fig. 4.

In cf-fitter, the Etminan-like scenario is **`q_fold`**: HAL **fit-B three-Gaussian** quartet potential + **analytic** Gaussian fold (A=4, rms=1.56 fm), KP solver in this codebase.

Similar physics setup; differences in potential shape, folding implementation, and mass inputs imply **qualitative** overlay only.

## Sanity check vs `KP_q_fold_R3.00` (`cf_phialpha_folded.csv`)

| k* [MeV/c] | Etminan C | q_fold KP | ΔC |
|------------|-----------|-----------|-----|
| ≈2 | 0.358 | 0.306 | +0.052 |
| ≈10 | 0.441 | 0.337 | +0.104 |
| ≈50 | 1.000 | 0.778 | +0.222 |

Etminan's R=3 curve rises faster toward unity; low-k agreement is reasonable, mid-k diverges as the published curve approaches C→1 earlier.

## Regenerate

```bash
python3 validation/etminan_phialpha/digitize_etminan_fig4.py
python3 figures/make_summary_pdf.py
python3 figures/make_phi_alpha_chizzali_tpe_figure.py
```
