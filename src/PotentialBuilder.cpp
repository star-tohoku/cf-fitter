#include "femto/PotentialBuilder.hpp"
#include "femto/FemtoModels.hpp"
#include <cmath>
#include <stdexcept>

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

} // namespace

std::function<double(double)> PotentialBuilder::gaussian(const GaussianSpec& spec) {
    std::vector<std::pair<double, double>> terms;
    for (const auto& t : spec.terms)
        terms.emplace_back(t.V_MeV, t.b_fm);
    return gaussianPotential(terms);
}

GaussianSpec PotentialBuilder::foldedGaussianSpec(const FoldedPotentialSpec& spec) {
    if (spec.phiNCentral.terms.empty())
        throw std::runtime_error("folded potential requires phiN central Gaussian terms");

    const double brho = densityWidth(spec.density);
    const double rho0 = densityNorm(spec.density.A, brho);

    GaussianSpec out;
    if (!spec.fold) {
        double Gc = 0.0;
        for (const auto& t : spec.phiNCentral.terms) {
            if (t.b_fm <= 0.0) throw std::runtime_error("Gaussian width must be positive");
            Gc += t.V_MeV * std::pow(M_PI, 1.5) * t.b_fm * t.b_fm * t.b_fm;
        }
        out.terms.push_back({Gc * rho0, brho});
        return out;
    }

    for (const auto& t : spec.phiNCentral.terms) {
        if (t.b_fm <= 0.0) throw std::runtime_error("Gaussian width must be positive");
        const double b2 = t.b_fm * t.b_fm;
        const double br2 = brho * brho;
        const double B2 = b2 + br2;
        const double B = std::sqrt(B2);
        const double amp = t.V_MeV * rho0 * std::pow(M_PI * b2 * br2 / B2, 1.5);
        out.terms.push_back({amp, B});
    }
    return out;
}

std::function<double(double)> PotentialBuilder::foldedGaussian(const FoldedPotentialSpec& spec) {
    return gaussian(foldedGaussianSpec(spec));
}

double PotentialBuilder::gaussianSpecV0(const GaussianSpec& spec) {
    double v = 0.0;
    for (const auto& t : spec.terms) v += t.V_MeV;
    return v;
}

double PotentialBuilder::gaussianSpecVolumeIntegral(const GaussianSpec& spec) {
    double integral = 0.0;
    for (const auto& t : spec.terms)
        integral += t.V_MeV * std::pow(M_PI, 1.5) * std::pow(t.b_fm, 3);
    return integral;
}

std::function<double(double)> PotentialBuilder::folded(
    const std::function<double(double)>& V2body,
    const std::function<double(double)>& density,
    const FoldingOptions&) {
    (void)density;
    return V2body;
}

} // namespace femto
