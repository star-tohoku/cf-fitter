#ifndef FEMTO_MODELS_HPP
#define FEMTO_MODELS_HPP
// ============================================================================
// FemtoModels.hpp
//   Femtoscopic correlation function models:
//     (1) Lednicky-Lyuboshitz (LL) analytical model
//     (2) Koonin-Pratt (KP) formula with numerical wave functions (Numerov)
//   Header-only, ROOT-independent. Units: r [fm], k [fm^-1], V,mu [MeV].
//   Convention: f(k) = ( 1/f0 + d0 k^2 / 2 - i k )^{-1}   [femtoscopy sign!]
//               (nuclear convention: k cot(delta) = -1/a0 + re k^2/2 -> f0 = -a0)
// ============================================================================
#include <cmath>
#include <complex>
#include <vector>
#include <map>
#include <functional>
#include <memory>
#include <stdexcept>

namespace femto {

constexpr double HBARC = 197.3269804;   // MeV fm
constexpr double SQRTPI = 1.7724538509055160;

// ---------------------------------------------------------------------------
// Dawson function D(z) = exp(-z^2) * int_0^z exp(t^2) dt   (Rybicki algorithm)
// ---------------------------------------------------------------------------
inline double dawson(double x) {
    const int    NMAX = 6;
    const double H = 0.4, A1 = 2.0/3.0, A2 = 0.4, A3 = 2.0/7.0;
    static double c[NMAX + 1];
    static bool init = false;
    if (!init) {
        for (int i = 1; i <= NMAX; ++i) {
            double t = (2.0*i - 1.0)*H;
            c[i] = std::exp(-t*t);
        }
        init = true;
    }
    double ax = std::fabs(x);
    if (ax < 0.2) {                       // small-x Maclaurin series
        double x2 = x*x;
        return x*(1.0 - A1*x2*(1.0 - A2*x2*(1.0 - A3*x2)));
    }
    int    n0 = 2*int(0.5*ax/H + 0.5);
    double xp = ax - n0*H;
    double e1 = std::exp(2.0*xp*H), e2 = e1*e1;
    double d1 = n0 + 1.0, d2 = d1 - 2.0, sum = 0.0;
    for (int i = 1; i <= NMAX; ++i, d1 += 2.0, d2 -= 2.0, e1 *= e2)
        sum += c[i]*(e1/d1 + 1.0/(d2*e1));
    double ans = (1.0/SQRTPI)*std::exp(-xp*xp)*sum;
    return (x >= 0.0) ? ans : -ans;
}

// LL auxiliary functions
inline double F1(double z) {              // = D(z)/z
    if (z < 1e-8) return 1.0;
    return dawson(z)/z;
}
inline double F2(double z) {              // = (1 - exp(-z^2))/z
    if (z < 1e-8) return z;
    return (1.0 - std::exp(-z*z))/z;
}

// ---------------------------------------------------------------------------
// Effective-range scattering amplitude (f0 may be complex: open channels)
// ---------------------------------------------------------------------------
inline std::complex<double> f_ERE(double k, std::complex<double> f0, double d0) {
    std::complex<double> inv = 1.0/f0 + 0.5*d0*k*k - std::complex<double>(0.0, k);
    return 1.0/inv;
}

// ---------------------------------------------------------------------------
// Spin channel: statistical weight + ERE parameters
// ---------------------------------------------------------------------------
struct SpinChannel {
    double weight;                 // rho_S (e.g. singlet 1/4, triplet 3/4)
    std::complex<double> f0;       // [fm]  femtoscopy convention
    double d0;                     // [fm]
    bool   interacting = true;     // e.g. LL-LL triplet: Pauli-blocked s-wave
};

// ===========================================================================
// (1) Lednicky-Lyuboshitz model
// ===========================================================================
//   C(k) = 1 + qs_amplitude * exp(-4 k^2 R^2)
//        + sum_ch  w * g_sym * [ |f|^2/(2R^2) (1 - d0/(2 sqrt(pi) R))
//                               + 2 Re f /(sqrt(pi) R) F1(2kR)
//                               - Im f / R           F2(2kR) ]
//   g_sym = 2 for identical pairs (symmetrized s-wave), 1 otherwise.
//   qs_amplitude: 0 (non-identical), -1/2 (identical spin-1/2 pair),
//                 +1 (identical spin-0 bosons), ...
// ===========================================================================
class LLModel {
public:
    std::vector<SpinChannel> channels;
    bool   identical    = false;
    double qsAmplitude  = 0.0;

    double C(double k, double R) const {
        if (k < 1e-6) k = 1e-6;
        double c = 1.0 + qsAmplitude*std::exp(-4.0*k*k*R*R);
        const double gsym = identical ? 2.0 : 1.0;
        const double z = 2.0*k*R;
        for (const auto& ch : channels) {
            if (!ch.interacting) continue;
            std::complex<double> f = f_ERE(k, ch.f0, ch.d0);
            double t1 = 0.5*std::norm(f)/(R*R)*(1.0 - ch.d0/(2.0*SQRTPI*R));
            double t2 = 2.0*f.real()/(SQRTPI*R)*F1(z);
            double t3 = -f.imag()/R*F2(z);
            c += ch.weight*gsym*(t1 + t2 + t3);
        }
        return c;
    }
};

// Pre-baked constructors -----------------------------------------------------
// phi-p: non-identical pair (phi spin-1, p spin-1/2) -> total spin 1/2 or 3/2.
//   doublet 2S_{1/2} weight 1/3, quartet 4S_{3/2} weight 2/3.
inline LLModel makeLL_phip(std::complex<double> f0_doublet, double d0_doublet,
                           std::complex<double> f0_quartet, double d0_quartet) {
    LLModel m;                              // non-identical, two spin channels
    m.identical   = false;
    m.qsAmplitude = 0.0;
    m.channels = { {1.0/3.0, f0_doublet, d0_doublet, true},
                   {2.0/3.0, f0_quartet, d0_quartet, true} };
    return m;
}
// phi-alpha: alpha is spin-0, phi spin-1 -> single spin-1 channel.
inline LLModel makeLL_phialpha(std::complex<double> f0, double d0) {
    LLModel m;                              // non-identical, single channel
    m.identical   = false;
    m.qsAmplitude = 0.0;
    m.channels    = { {1.0, f0, d0, true} };
    return m;
}
// Generic single channel (e.g. an effective amplitude; complex f0 = absorption).
inline LLModel makeLL_singleChannel(std::complex<double> f0, double d0) {
    LLModel m;
    m.channels = { {1.0, f0, d0, true} };
    return m;
}
// Generic identical spin-0 boson pair (used by the LL<->KP consistency check).
inline LLModel makeLL_identicalBosons(std::complex<double> f0, double d0) {
    LLModel m;
    m.identical   = true;                   // symmetrised s-wave (gsym = 2)
    m.qsAmplitude = +1.0;                   // spin-0 boson statistics
    m.channels    = { {1.0, f0, d0, true} };
    return m;
}

// ===========================================================================
// (2) Koonin-Pratt model
// ===========================================================================
// Radial Schroedinger solver (s-wave, real local potential), Numerov method.
//   u'' = [ 2 mu V(r)/(hbarc)^2 - k^2 ] u ,  u(0)=0
//   Matching to A sin(kr) + B cos(kr) at two points beyond the range.
// ---------------------------------------------------------------------------
class RadialSolverS {
public:
    RadialSolverS(std::function<double(double)> V_MeV,
                  double mu_MeV,
                  double rMax = 60.0, double h = 0.01, double rMatch = 12.0)
        : V_(std::move(V_MeV)), mu_(mu_MeV),
          rMax_(rMax), h_(h), rMatch_(rMatch) {}

    // |psi_0(k,r)|^2 sampled on the internal r-grid; also returns delta(k).
    // Normalization: psi_0 -> sin(kr+delta)/(kr) asymptotically.
    struct Wave {
        std::vector<double> psi0sq;   // |psi_0|^2 on grid r_i = (i+1)*h
        double delta;                 // phase shift [rad]
    };

    int    nGrid() const { return int(rMax_/h_); }
    double rAt(int i) const { return (i + 1)*h_; }
    double hStep() const { return h_; }

    Wave solve(double k) const {
        const int N = nGrid();
        std::vector<double> u(N);
        const double fac = 2.0*mu_/(HBARC*HBARC);
        auto w = [&](double r){ return fac*V_(r) - k*k; };

        // Numerov: u'' = w u ;  near origin u(r) ~ r + O(r^3) for regular V
        u[0] = rAt(0);                               // u(h)  = h
        u[1] = rAt(1)*(1.0 + rAt(1)*rAt(1)*w(rAt(1))/6.0);  // u(2h) ~ r(1 + w r^2/6)
        const double h2_12 = h_*h_/12.0;
        for (int i = 1; i < N - 1; ++i) {
            double wm = w(rAt(i-1)), w0 = w(rAt(i)), wp = w(rAt(i+1));
            u[i+1] = ( 2.0*u[i]*(1.0 + 5.0*h2_12*w0)
                     - u[i-1]*(1.0 - h2_12*wm) ) / (1.0 - h2_12*wp);
        }

        // Match at r1, r2 (> potential range), quarter wavelength apart
        int i1 = int(rMatch_/h_) - 1;
        int i2 = i1 + std::max(1, int((M_PI/(2.0*k))/h_));
        if (i2 >= N) i2 = N - 1;
        double r1 = rAt(i1), r2 = rAt(i2);
        double det = std::sin(k*(r1 - r2));
        if (std::fabs(det) < 1e-12) { i2 = (i1 + N - 1)/2; r2 = rAt(i2);
                                      det = std::sin(k*(r1 - r2)); }
        double A = ( u[i1]*std::cos(k*r2) - u[i2]*std::cos(k*r1) )/det;
        double B = ( u[i2]*std::sin(k*r1) - u[i1]*std::sin(k*r2) )/det;
        double C0 = std::hypot(A, B);
        double delta = std::atan2(B, A);

        Wave out;
        out.delta = delta;
        out.psi0sq.resize(N);
        for (int i = 0; i < N; ++i) {
            double x = u[i]/(C0*k*rAt(i));
            out.psi0sq[i] = x*x;
        }
        return out;
    }

    // Low-energy parameters from k cot(delta) = 1/f0 + d0 k^2/2 (femto conv.)
    void scatteringParams(double& f0, double& d0,
                          double k1 = 0.02, double k2 = 0.08) const {
        double y1 = k1/std::tan(solve(k1).delta);
        double y2 = k2/std::tan(solve(k2).delta);
        d0 = 2.0*(y2 - y1)/(k2*k2 - k1*k1);
        f0 = 1.0/(y1 - 0.5*d0*k1*k1);
    }

private:
    std::function<double(double)> V_;
    double mu_, rMax_, h_, rMatch_;
};

// ---------------------------------------------------------------------------
// KP correlation function with Gaussian (or user) source.
//   C(k;R) = 1 + qsAmp * exp(-4k^2R^2)
//          + sum_ch w * g_sym * Int 4 pi r^2 S(r) [ |psi_0|^2 - j0(kr)^2 ] dr
// Wave functions cached per k (independent of R)  ->  fast R-fits.
// ---------------------------------------------------------------------------
class KPModel {
public:
    struct Channel {
        double weight;
        std::shared_ptr<RadialSolverS> solver;   // nullptr = non-interacting
    };

    std::vector<Channel> channels;
    bool   identical   = false;
    double qsAmplitude = 0.0;

    // Gaussian source by default; replaceable (resonance halo, EPOS, ...)
    // S(r) must be normalized: Int 4 pi r^2 S(r) dr = 1, R passed separately.
    std::function<double(double r, double R)> source =
        [](double r, double R) {
            double n = std::pow(4.0*M_PI*R*R, -1.5);
            return n*std::exp(-r*r/(4.0*R*R));
        };

    double C(double k, double R) const {
        if (k < 1e-4) k = 1e-4;
        double c = 1.0 + qsAmplitude*std::exp(-4.0*k*k*R*R);
        const double gsym = identical ? 2.0 : 1.0;
        for (std::size_t ic = 0; ic < channels.size(); ++ic) {
            const auto& ch = channels[ic];
            if (!ch.solver) continue;
            const auto& dpsi = cached(ic, k);          // |psi0|^2 - j0^2 on grid
            const double h = ch.solver->hStep();
            const int    N = (int)dpsi.size();
            // Simpson integration of 4 pi r^2 S(r;R) * dpsi(r)
            double sum = 0.0;
            for (int i = 0; i < N; ++i) {
                double r = ch.solver->rAt(i);
                double f = 4.0*M_PI*r*r*source(r, R)*dpsi[i];
                double wgt = (i == 0 || i == N-1) ? 1.0 : (i % 2 ? 4.0 : 2.0);
                sum += wgt*f;
            }
            c += ch.weight*gsym*(h/3.0)*sum;
        }
        return c;
    }

    void clearCache() const { cache_.clear(); }

private:
    // cache key: (channel, k rounded to 1e-6 fm^-1)
    mutable std::map<std::pair<std::size_t,long long>,
                     std::vector<double>> cache_;

    const std::vector<double>& cached(std::size_t ic, double k) const {
        long long kk = (long long)std::llround(k*1e6);
        auto key = std::make_pair(ic, kk);
        auto it = cache_.find(key);
        if (it != cache_.end()) return it->second;

        const auto& sol = *channels[ic].solver;
        auto wave = sol.solve(k);
        std::vector<double> dpsi(wave.psi0sq.size());
        for (std::size_t i = 0; i < dpsi.size(); ++i) {
            double kr = k*sol.rAt(i);
            double j0 = (kr < 1e-8) ? 1.0 : std::sin(kr)/kr;
            dpsi[i] = wave.psi0sq[i] - j0*j0;
        }
        return cache_.emplace(key, std::move(dpsi)).first->second;
    }
};

// Common potentials ----------------------------------------------------------
inline std::function<double(double)>
gaussianPotential(std::vector<std::pair<double,double>> terms /*(V_i MeV, b_i fm)*/) {
    return [terms](double r) {
        double v = 0.0;
        for (const auto& t : terms) v += t.first*std::exp(-r*r/(t.second*t.second));
        return v;
    };
}

} // namespace femto
#endif
