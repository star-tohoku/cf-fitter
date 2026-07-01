#include "femto/PotentialBuilder.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace femto {
namespace {

double densityWidth(const GaussianDensitySpec& density) {
    if (density.b_fm > 0.0) return density.b_fm;
    if (density.rms_fm <= 0.0) throw std::runtime_error("Gaussian density requires rms_fm or b_fm");
    return density.rms_fm * std::sqrt(2.0 / 3.0);
}

double densityNorm(double A, double b) {
    return A / (std::pow(M_PI, 1.5) * b * b * b);
}

double rhoAt(double s, double rho0, double b_rho) {
    return rho0 * std::exp(-s * s / (b_rho * b_rho));
}

double simpsonRadial(const std::function<double(double)>& integrand,
                     double rMax, double dr) {
    if (rMax <= 0.0 || dr <= 0.0) return 0.0;
    const int n = std::max(2, int(std::ceil(rMax / dr)));
    const double h = rMax / n;
    double sum = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double r = i * h;
        const double f = integrand(r);
        const double w = (i == 0 || i == n) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
        sum += w * f;
    }
    return sum * h / 3.0;
}

double angularMuIntegral(double R, double s, const std::function<double(double)>& VN, int nMu) {
    if (nMu < 4) nMu = 4;
    const double dmu = 2.0 / nMu;
    double sum = 0.0;
    for (int i = 0; i <= nMu; ++i) {
        const double mu = -1.0 + i * dmu;
        const double rad2 = R * R + s * s - 2.0 * R * s * mu;
        const double rad = rad2 > 0.0 ? std::sqrt(rad2) : 0.0;
        const double w = (i == 0 || i == nMu) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
        sum += w * VN(rad);
    }
    return sum * dmu / 3.0;
}

double foldAtR(double R, double rho0, double b_rho,
               const std::function<double(double)>& VN,
               const NumericFoldingOptions& opts) {
    if (R <= 0.0) {
        return 4.0 * M_PI * simpsonRadial(
            [&](double s) { return s * s * rhoAt(s, rho0, b_rho) * VN(s); },
            opts.sMax_fm, opts.ds_fm);
    }

    double sum = 0.0;
    const int nS = std::max(2, int(std::ceil(opts.sMax_fm / opts.ds_fm)));
    const double ds = opts.sMax_fm / nS;
    for (int is = 0; is <= nS; ++is) {
        const double s = is * ds;
        const double angular = angularMuIntegral(R, s, VN, opts.nMu);
        const double f = s * s * rhoAt(s, rho0, b_rho) * angular;
        const double w = (is == 0 || is == nS) ? 1.0 : ((is % 2) ? 4.0 : 2.0);
        sum += w * f;
    }
    return 2.0 * M_PI * sum * ds / 3.0;
}

std::function<double(double)> makeTabulatedPotential(std::vector<double> rGrid,
                                                     std::vector<double> vGrid) {
    return [rGrid = std::move(rGrid), vGrid = std::move(vGrid)](double R) -> double {
        if (R <= 0.0) return vGrid.empty() ? 0.0 : vGrid.front();
        if (R >= rGrid.back()) return 0.0;
        const auto it = std::lower_bound(rGrid.begin(), rGrid.end(), R);
        const std::size_t j = std::size_t(it - rGrid.begin());
        if (j == 0) return vGrid.front();
        const double r0 = rGrid[j - 1];
        const double r1 = rGrid[j];
        const double v0 = vGrid[j - 1];
        const double v1 = vGrid[j];
        const double t = (R - r0) / (r1 - r0);
        return v0 + t * (v1 - v0);
    };
}

} // namespace

double PotentialBuilder::phiNVolumeIntegral(const std::function<double(double)>& V2body,
                                          const NumericFoldingOptions& opts) {
    return 4.0 * M_PI * simpsonRadial(
        [&](double r) { return r * r * V2body(r); },
        opts.volumeRMax_fm, opts.volumeDr_fm);
}

double PotentialBuilder::potentialVolumeIntegral(const std::function<double(double)>& V,
                                                 double rMax_fm, double dr_fm) {
    return 4.0 * M_PI * simpsonRadial(
        [&](double r) { return r * r * V(r); },
        rMax_fm, dr_fm);
}

std::function<double(double)> PotentialBuilder::foldedPhiN(
    const std::function<double(double)>& V2body,
    const GaussianDensitySpec& density,
    bool fold,
    const NumericFoldingOptions& opts) {
    const double b_rho = densityWidth(density);
    const double rho0 = densityNorm(density.A, b_rho);

    if (!fold) {
        const double G = phiNVolumeIntegral(V2body, opts);
        return [G, rho0, b_rho](double R) {
            return G * rhoAt(R, rho0, b_rho);
        };
    }

    std::vector<double> rGrid;
    std::vector<double> vGrid;
    const int nR = std::max(2, int(std::ceil(opts.rGridMax_fm / opts.dr_fm)));
    rGrid.reserve(nR + 1);
    vGrid.reserve(nR + 1);
    for (int iR = 0; iR <= nR; ++iR) {
        const double R = iR * opts.rGridMax_fm / nR;
        rGrid.push_back(R);
        vGrid.push_back(foldAtR(R, rho0, b_rho, V2body, opts));
    }
    return makeTabulatedPotential(std::move(rGrid), std::move(vGrid));
}

} // namespace femto
