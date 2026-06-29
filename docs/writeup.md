# Femtoscopic Correlation-Function Fitting: Lednický–Lyuboshitz and Koonin–Pratt Implementations

*Implementation note and validation summary*

> **Abstract.** This note documents a self-contained C++ framework for computing and fitting two-particle femtoscopic correlation functions $C(k^*)$ using two complementary formalisms: the analytic Lednický–Lyuboshitz (LL) model and the numerical Koonin–Pratt (KP) source-integral with wave functions from a Numerov solver. The two methods are cross-checked against each other and against an analytic square-well benchmark. The framework is then applied to generate $\phi p$ and $\phi\alpha$ correlation functions for several interaction scenarios. The purpose of this document is to state precisely *what* is computed, *how* each ingredient is defined, and *how reliable* the results are, so that the methodology can be defended to others.

---

## 1. What is computed

The observable is the femtoscopic correlation function

$$
C(k^*) = \int d^3r \; S_{12}(r)\,\bigl|\Psi^{(-)}_{k^*}(\mathbf{r})\bigr|^2 ,
$$

where $k^*$ is half the pair relative momentum in the pair rest frame, $S_{12}(r)$ the pair emission source, and $\Psi^{(-)}$ the relative wave function carrying the final-state interaction. Both formalisms compute this same quantity under the $s$-wave approximation; they differ only in how $|\Psi|^2$ is obtained.

**Sign / unit conventions.** All low-energy parameters use the *femtoscopy* effective-range expansion (ERE)

$$
k^*\cot\delta_0 = \frac{1}{f_0} + \frac{1}{2}d_0\,k^{*2} ,
\tag{ERE}
$$

which is related to the nuclear convention $k^*\cot\delta_0 = -1/a_0 + \tfrac{1}{2}r_e\,k^{*2}$ by $f_0=-a_0$, $d_0=r_e$. Lengths in fm, momenta in fm⁻¹ (converted from MeV/$c$ with $\hbar c = 197.327$ MeV·fm), potentials in MeV.

### 1.1 Definitions of the low-energy parameters

**$s$-wave phase shift $\delta_0$.** For a finite-range potential, the $\ell=0$ reduced radial wave function $u(r)=r\,\psi(r)$ obeys $u'' + [k^{*2} - 2\mu V/(\hbar c)^2]\,u = 0$ and, beyond the range of $V$, becomes a free wave shifted in phase:

$$
u(r) \;\xrightarrow{\;r\to\infty\;}\; \text{const}\times \sin\!\bigl(k^* r + \delta_0(k^*)\bigr).
$$

$\delta_0(k^*)$ is the *$s$-wave phase shift*: the shift of the asymptotic wave relative to the non-interacting case ($\delta_0=0$ when $V=0$). It is the single real number that fully encodes elastic $s$-wave scattering at momentum $k^*$. The $s$-wave scattering amplitude is $f(k^*) = (k^*\cot\delta_0 - ik^*)^{-1} = e^{i\delta_0}\sin\delta_0/k^*$, and in the code $\delta_0$ is what the Numerov matching returns.

**Effective-range expansion.** At low energy $k^*\cot\delta_0$ is analytic in $k^{*2}$ and is expanded as Eq. (ERE). The first two coefficients are the scattering length and the effective range; higher terms are neglected here.

**Scattering length $a_0$ (nuclear convention).** The leading, zero-energy coefficient,

$$
-\frac{1}{a_0} = \lim_{k^*\to0} k^*\cot\delta_0
\quad\Longleftrightarrow\quad
a_0 = -\lim_{k^*\to0}\frac{\tan\delta_0}{k^*}.
$$

Geometrically, the zero-energy wave function extrapolates to a straight line $u_0(r)\propto(r-a_0)$ outside the range, so $a_0$ is the intercept of that line on the $r$-axis. The threshold cross section is $\sigma_0 = 4\pi a_0^2$. Sign content in this convention:

- $a_0>0$ small — repulsion (hard-sphere limit $a_0\to R$);
- $a_0<0$ — attraction with **no** bound state;
- $a_0\to+\infty$ (large positive) — attraction with a shallow bound state.

**Effective range $r_e$.** The next coefficient,

$$
r_e = 2\int_0^\infty \bigl[\,v_0^2(r) - u_0^2(r)\,\bigr]\,dr ,
$$

where $u_0$ is the true zero-energy solution and $v_0(r)=1-r/a_0$ its asymptotic straight-line form, both normalised to $u_0(0)=v_0(0)=1$. It sets the energy dependence of $k^*\cot\delta_0$ and is of the order of the spatial range of the interaction (hence "effective range").

**Femtoscopy parameters $f_0,d_0$.** The Lednický–Lyuboshitz formula is written with the opposite overall sign in the leading term, Eq. (ERE), i.e. the scattering amplitude reads

$$
f(k^*) = \left(\frac{1}{f_0} + \frac{1}{2}d_0\,k^{*2} - ik^*\right)^{-1}.
$$

Thus $f_0$ is the *scattering length in the femtoscopy convention* ($f_0=-a_0$) and $d_0$ the effective range ($d_0=r_e$). The sign content is therefore reversed relative to $a_0$:

- $f_0>0$ — attraction without a bound state;
- $f_0<0$ small — repulsion;
- $f_0\to-\infty$ (large negative) — a bound state appears.

This is the source of the most common cross-comparison error and is the reason the two scenarios in the $\phi\alpha$ study with $f_0=+24.5$ fm (near-unitary, unbound) and $f_0=-8.53$ fm (bound) sit on opposite sides of the threshold. In the code $f_0$ is complex, so $\mathrm{Im}\,f_0>0$ additionally encodes flux loss to inelastic channels (absorption).

---

## 2. Implementation

The code is organised in three layers, all in a single header (`FemtoModels.hpp`, ROOT-independent) plus driver programs.

### 2.1 Layer 1 — physics models

**(a) Lednický–Lyuboshitz (analytic).** For each spin channel with scattering amplitude $f(k^*) = (1/f_0 + \tfrac{1}{2}d_0 k^{*2} - ik^*)^{-1}$ and a Gaussian source of radius $R$,

$$
C(k^*) = 1 + \alpha_{\rm qs}\,e^{-4k^{*2}R^2}
+ \sum_{\rm ch} w_{\rm ch}\,g_{\rm sym}
\left[\frac{|f|^2}{2R^2}\!\left(1-\frac{d_0}{2\sqrt\pi R}\right)
+ \frac{2\,\mathrm{Re}f}{\sqrt\pi R}F_1(2k^*R)
- \frac{\mathrm{Im}f}{R}F_2(2k^*R)\right],
\tag{LL}
$$

with $F_1(z)=\int_0^1 e^{z^2(x^2-1)}dx = D(z)/z$ (Dawson function), $F_2(z)=(1-e^{-z^2})/z$. Implementation details:

- $f_0$ is stored as `std::complex`, so open/absorptive channels (e.g. an effective $K^-p$ or $\phi p$ amplitude) are handled directly via $\mathrm{Im}\,f_0$.
- Identical pairs set $g_{\rm sym}=2$ (symmetrised $s$-wave) and a quantum-statistics term $\alpha_{\rm qs}$ ($+1$ for identical spin-$0$ bosons, $-\tfrac{1}{2}$ for an identical spin-$\tfrac{1}{2}$ pair, $0$ for non-identical pairs).
- The Dawson function is implemented from scratch (Rybicki algorithm) since it is not in the standard library.

The known limitation of (LL) is the factor $(1 - d_0/2\sqrt\pi R)$, which becomes unphysical for $R\lesssim d_0/2\sqrt\pi$; this is exactly where KP is required.

**(b) Koonin–Pratt (numerical).**

$$
C(k^*) = 1 + \alpha_{\rm qs}\,e^{-4k^{*2}R^2}
+ \sum_{\rm ch} w_{\rm ch}\,g_{\rm sym}
\int_0^\infty 4\pi r^2\,S(r;R)\,
\bigl[\,|\psi_0(k^*,r)|^2 - j_0(k^* r)^2\,\bigr]\,dr .
\tag{KP}
$$

The relative wave function is obtained by solving the radial Schrödinger equation for a local potential $V(r)$,

$$
u''(r) = \Bigl[\tfrac{2\mu}{(\hbar c)^2}V(r) - k^{*2}\Bigr]u(r),\qquad u(0)=0,
$$

with the **Numerov** method. The regular solution is started from $u(h)=h$, $u(2h)=r(1+w r^2/6)$, propagated outward, and matched to $A\sin k^* r + B\cos k^* r$ at two points beyond the potential range to extract the phase shift and normalisation ($\psi_0\to\sin(k^* r+\delta_0)/k^* r$ asymptotically).

- **Performance:** $|\psi_0|^2$ does not depend on $R$, so the quantity $|\psi_0|^2 - j_0^2$ is computed once per $k^*$ and cached on the $r$-grid; only the Simpson integral over $S(r;R)$ is repeated during a fit. This is the same design principle as the CATS package.
- The source $S(r;R)$ is a normalised Gaussian by default but is a `std::function`, so a resonance-halo or transport-model source can be substituted without touching the integrator.
- `scatteringParams()` reads $(f_0,d_0)$ back out of the solved phase shifts via Eq. (ERE) at two low momenta, providing the bridge that lets the **same** interaction be fed to both LL and KP.

### 2.2 Layer 2 — experimental effects

The fit function applied to data is

$$
C_{\rm fit}(k^*) = N\,(1 + b_1 k^*)\,\bigl[\,1 + \lambda\,(C_{\rm model}(k^*) - 1)\,\bigr],
$$

with normalisation $N$, linear non-femtoscopic baseline $b_1$, and $\lambda$ the pair purity × primary fraction. A momentum-resolution smearing utility (folding a model curve with an MC response matrix $M(k^*_{\rm true}, k^*_{\rm rec})$) is provided.

### 2.3 Layer 3 — fitting

ROOT `TF1` wrappers (`fitCorrelation.C`) expose two modes:

- **LL mode**: $(R, f_0, d_0, \lambda, N, b_1)$ free — model-independent extraction of scattering parameters from data.
- **KP mode**: interaction *fixed* by a chosen potential, $(R, \lambda, N, b_1)$ free — for testing a specific theory (lattice / chiral / phenomenological) against data.

### 2.4 How the potential is specified

The potential is a sum of Gaussians,

$$
V(r) = \sum_i V_i\,\exp\!\bigl(-r^2/b_i^2\bigr),
$$

where $b_i$ is the $1/e$ width (not the standard deviation; $b=\sqrt2\,\sigma$). In the present study a single Gaussian is used. The depth $V_0$ is **not** taken from first principles: a bisection routine (`tuneV0`) adjusts $V_0$ so that the solved $f_0$ matches a target value, with $b$ chosen to place $d_0$ in the desired range. This makes the potential a phenomenological device reproducing a chosen $(f_0,d_0)$, which is sufficient because $C(k^*)$ is controlled by $(f_0, d_0, R, \text{source})$. It is **not** a microscopic $\phi$–nucleus potential (see Section 5).

---

## 3. Validation

### 3.1 Solver vs. analytic square well

For an attractive square well ($V_0=-40$ MeV, range $1$ fm, $\mu=488.6$ MeV, the $\phi p$ reduced mass) the numerical phase shift matches the closed-form result $\delta_0 = -k^* b + \arctan[(k^*/K)\tan Kb]$:

| $k^*$ [fm⁻¹] | $\delta_{\rm num}$ | $\delta_{\rm analytic}$ | rel. diff |
|:---:|:---:|:---:|:---:|
| 0.02 | $+0.01117$ | $+0.01121$ | $0.4\%$ |
| 0.05 | $+0.02789$ | $+0.02789$ | $0.0\%$ |
| 0.10 | $+0.05553$ | $+0.05485$ | $1.2\%$ |

The extracted scattering length, $f_0=0.5586$ fm, agrees with the analytic value $-a_0=0.5610$ fm to $0.5\%$. (At $k^*=0.30$ fm⁻¹ the two-point ERE extraction starts to deviate, as expected, since the fit assumes Eq. (ERE) holds; the low-momentum region used for $(f_0,d_0)$ is unaffected.)

**Bug found and fixed during validation.** The initial Numerov start used a crude second point that admitted a small admixture of the irregular ($\propto 1/r$) solution, distorting the normalisation. Enforcing the regular boundary condition $u\propto r$ near the origin removed it; the square-well benchmark above is the post-fix result.

### 3.2 KP vs. LL consistency

The cross-check uses the simplest fully-specified system — a pair of **identical spin-0 bosons** — so that spin weights and channel structure are unambiguous (single channel, weight 1; symmetrised $s$-wave, $g_{\rm sym}=2$; statistics term $+e^{-4k^{*2}R^2}$). Feeding the **same** Gaussian potential ($f_0=0.696$ fm, $d_0=2.892$ fm extracted) to both formalisms and comparing $C(k^*)$:

| $R$ [fm] | 1.2 | 2.0 | 3.0 | 5.0 |
|:---|:---:|:---:|:---:|:---:|
| max $\|C_{\rm KP}-C_{\rm LL}\|/C_{\rm LL}$ | $0.57\%$ | $0.15\%$ | $0.05\%$ | $0.01\%$ |

The two independent methods agree to well below the percent level for $R\gtrsim2$ fm, and to $<1\%$ even at $R=1.2$ fm. The residual difference is the combined effect of the ERE truncation (LL) and the finite numerical grid (KP), and it grows as $R$ shrinks — consistent with the expected breakdown of the LL analytic form at small source size. (Note that the statistics term is added identically in both formalisms, so it cancels in the difference; this comparison therefore tests the interaction part, with the $g_{\rm sym}=2$ symmetrisation.) This agreement, together with the analytic benchmark, is the basis for trusting both implementations.

---

## 4. Application: $\phi p$ and $\phi\alpha$

Both pairs are charge-neutral, so no Coulomb interaction is needed and the present $s$-wave code applies directly. Reduced masses: $\mu_{\phi p}=488.6$ MeV, $\mu_{\phi\alpha}=800.5$ MeV.

**$\phi p$** ($^2S_{1/2}$ weight $1/3$, $^4S_{3/2}$ weight $2/3$):

- **Set A** (HAL QCD-like): attraction in $^4S_{3/2}$ only, tuned to $f_0=1.43$ fm, $d_0=2.35$ fm.
- **Set B**: Set A plus a bound $^2S_{1/2}$ ($f_0=-4.24$ fm, $B\approx3$ MeV from the ERE pole), producing a characteristic dip below unity near $k^*\sim130$–$200$ MeV/$c$.
- **Set C**: effective single channel with complex $f_0=0.85+0.16i$ fm (absorption; LL only).

**$\phi\alpha$** (single channel): depths set as fractions of the critical depth $V_{\rm crit}\approx16.3$ MeV (the depth at which the first bound state appears):

| scenario | $V_0$ [MeV] | $f_0$ [fm] | $d_0$ [fm] | $B$ (ERE) |
|:---|:---:|:---:|:---:|:---:|
| weak ($0.40\,V_{\rm crit}$) | $-6.5$ | $+1.54$ | $5.52$ | — |
| near-unitary ($0.92\,V_{\rm crit}$) | $-15.0$ | $+24.5$ | $3.05$ | — |
| bound ($1.35\,V_{\rm crit}$) | $-21.9$ | $-8.53$ | $2.38$ | $\sim0.5$ MeV |

![phi-p correlation functions](cf_phip.png)

*Figure 1. $\phi p$ correlation functions (solid: LL, dashed: KP) for the three interaction scenarios at $R=1.2$ and $3.0$ fm. The bound-$^2S_{1/2}$ scenario (red) shows the dip below unity. LL and KP overlap for these parameters.*

![phi-alpha correlation functions](cf_phialpha.png)

*Figure 2. $\phi\alpha$ correlation functions (log scale; solid: LL, dashed: KP) at $R=1.2, 2.5, 5.0$ fm. The bound scenario (red) inverts from enhancement to depletion as $R$ grows — the classic signature used to diagnose a bound state. At $R=1.2$ fm the bound case shows a visible $\sim50\%$ LL/KP discrepancy, demonstrating the LL breakdown when $R\lesssim d_0$ and $|f_0|>R$; the two agree for $R\gtrsim2.5$ fm.*

---

## 5. Reliability and limitations

**What is solid.**

- The $s$-wave solver is verified against an analytic square well ($\lesssim1\%$ in $\delta_0$ and $f_0$ at low $k^*$).
- The two independent formalisms (analytic LL, numerical KP) agree at the sub-percent level for $R\gtrsim2$ fm, providing strong internal consistency.
- The qualitative physics is correctly reproduced: enhancement for attraction, $C(0)$ growing toward unitarity, and the enhancement→depletion inversion with $R$ for a bound state.

**What is approximate or assumed.**

- **$s$-wave only**; no higher partial waves or $d$-wave mixing.
- **Coulomb** is not implemented. This is fine for the neutral $\phi p$ and $\phi\alpha$ pairs studied here, but would be required for charged pairs ($K^-p$, etc.).
- **KP uses real potentials only**; absorption is currently captured only through a complex $f_0$ in LL. A coupled-channel or complex-potential extension is needed for genuinely inelastic systems.
- **The potential is phenomenological**, tuned to reproduce a target $(f_0,d_0)$. In particular the $\phi\alpha$ potential is *not* folded from a $\phi N$ interaction over the ⁴He density; the "two-body to many-body" logic is therefore not yet built in. A folding option $V_{\phi\alpha}(r)=\int d^3r'\,V_{\phi N}(|\mathbf{r}-\mathbf{r}'|)\rho_\alpha(r')$ is the natural next step.
- **Literature input values** (HAL QCD $\phi N$ parameters, the ALICE-like effective $f_0$) are approximate and *must be checked against the primary sources* before any quantitative comparison.
- In the demonstration fits $\lambda$ is left free; a real analysis fixes it from purity/feed-down MC and adds residual-correlation transformations.

**One-line summary for a questioner.** *Two independent calculations of the same observable (an analytic FSI formula and a direct numerical solution of the Schrödinger equation) are cross-validated against each other and against an exactly solvable test case; they agree where they should and diverge only in the regime where the analytic approximation is known to fail. The interaction inputs are phenomenological parametrisations of chosen scattering parameters, not microscopic potentials.*
