#include "femto/PotentialBuilder.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace {

constexpr double kTol = 1e-6;

bool near(double a, double b, double tol = kTol) {
    return std::fabs(a - b) <= tol;
}

double densityWidth(double rms_fm) {
    return rms_fm * std::sqrt(2.0 / 3.0);
}

double densityNorm(double A, double b_rho) {
    return A / (std::pow(M_PI, 1.5) * b_rho * b_rho * b_rho);
}

int testSingleGaussianFolded() {
    using namespace femto;

    const double a = -100.0;
    const double b = 0.50;
    const double A = 4.0;
    const double rms = 2.0;
    const double b_rho = densityWidth(rms);
    const double rho0 = densityNorm(A, b_rho);
    const double B2 = b * b + b_rho * b_rho;
    const double B = std::sqrt(B2);
    const double expectedAmp = a * rho0 * std::pow(M_PI * b * b * b_rho * b_rho / B2, 1.5);

    FoldedPotentialSpec spec;
    spec.fold = true;
    spec.density.A = A;
    spec.density.rms_fm = rms;
    spec.phiNCentral.terms.push_back({a, b});

    const GaussianSpec folded = PotentialBuilder::foldedGaussianSpec(spec);
    if (folded.terms.size() != 1) {
        std::printf("FAIL single folded: expected 1 term, got %zu\n", folded.terms.size());
        return 1;
    }
    if (!near(folded.terms[0].V_MeV, expectedAmp) || !near(folded.terms[0].b_fm, B)) {
        std::printf("FAIL single folded: amp=%.8f (exp %.8f), width=%.8f (exp %.8f)\n",
                    folded.terms[0].V_MeV, expectedAmp, folded.terms[0].b_fm, B);
        return 1;
    }

    auto V = PotentialBuilder::foldedGaussian(spec);
    const double r = 0.75;
    const double expectedV = expectedAmp * std::exp(-r * r / B2);
    if (!near(V(r), expectedV, 1e-4)) {
        std::printf("FAIL single folded V(%.2f): got %.8f expected %.8f\n", r, V(r), expectedV);
        return 1;
    }
    return 0;
}

int testSingleGaussianNoFold() {
    using namespace femto;

    const double a = -80.0;
    const double b = 0.40;
    const double A = 4.0;
    const double rms = 1.56;
    const double b_rho = densityWidth(rms);
    const double rho0 = densityNorm(A, b_rho);
    const double G = a * std::pow(M_PI, 1.5) * b * b * b;
    const double expectedAmp = G * rho0;

    FoldedPotentialSpec spec;
    spec.fold = false;
    spec.density.A = A;
    spec.density.rms_fm = rms;
    spec.phiNCentral.terms.push_back({a, b});

    const GaussianSpec folded = PotentialBuilder::foldedGaussianSpec(spec);
    if (folded.terms.size() != 1) {
        std::printf("FAIL no-fold: expected 1 term, got %zu\n", folded.terms.size());
        return 1;
    }
    if (!near(folded.terms[0].V_MeV, expectedAmp) || !near(folded.terms[0].b_fm, b_rho)) {
        std::printf("FAIL no-fold: amp=%.8f (exp %.8f), width=%.8f (exp %.8f)\n",
                    folded.terms[0].V_MeV, expectedAmp, folded.terms[0].b_fm, b_rho);
        return 1;
    }

    const double r = 0.5;
    auto V = PotentialBuilder::foldedGaussian(spec);
    const double expectedV = expectedAmp * std::exp(-r * r / (b_rho * b_rho));
    if (!near(V(r), expectedV, 1e-4)) {
        std::printf("FAIL no-fold V(%.2f): got %.8f expected %.8f\n", r, V(r), expectedV);
        return 1;
    }
    return 0;
}

int testHalQuartetFolded() {
    using namespace femto;

    FoldedPotentialSpec spec;
    spec.fold = true;
    spec.density.A = 4.0;
    spec.density.rms_fm = 1.56;
    spec.phiNCentral.terms = {
        {-371.0, 0.15},
        { -50.0, 0.66},
        { -31.0, 1.09},
    };

    const GaussianSpec folded = PotentialBuilder::foldedGaussianSpec(spec);
    if (folded.terms.size() != 3) {
        std::printf("FAIL HAL fold: expected 3 terms, got %zu\n", folded.terms.size());
        return 1;
    }

    const double b_rho = densityWidth(1.56);
    const double rho0 = densityNorm(4.0, b_rho);
    for (std::size_t i = 0; i < spec.phiNCentral.terms.size(); ++i) {
        const auto& t = spec.phiNCentral.terms[i];
        const double B2 = t.b_fm * t.b_fm + b_rho * b_rho;
        const double B = std::sqrt(B2);
        const double expectedAmp =
            t.V_MeV * rho0 * std::pow(M_PI * t.b_fm * t.b_fm * b_rho * b_rho / B2, 1.5);
        if (!near(folded.terms[i].V_MeV, expectedAmp, 1e-4) ||
            !near(folded.terms[i].b_fm, B, 1e-6)) {
            std::printf("FAIL HAL fold term %zu\n", i);
            return 1;
        }
    }

    auto V = PotentialBuilder::foldedGaussian(spec);
    if (!near(V(0.0), PotentialBuilder::gaussianSpecV0(folded), 1e-6)) {
        std::printf("FAIL HAL fold V(0)\n");
        return 1;
    }
    return 0;
}

int testQdCentralFolded() {
    using namespace femto;

    FoldedPotentialSpec spec;
    spec.fold = true;
    spec.density.A = 4.0;
    spec.density.rms_fm = 1.56;
    const double V0 = 571.5;
    spec.phiNCentral.terms = {
        {-371.0 * 2.0 / 3.0, 0.15},
        { -50.0 * 2.0 / 3.0, 0.66},
        { -31.0 * 2.0 / 3.0, 1.09},
        {-V0 / 3.0, 0.55},
    };

    const GaussianSpec folded = PotentialBuilder::foldedGaussianSpec(spec);
    if (folded.terms.size() != 4) {
        std::printf("FAIL q+d fold: expected 4 terms, got %zu\n", folded.terms.size());
        return 1;
    }

    auto V = PotentialBuilder::foldedGaussian(spec);
    const double vol = PotentialBuilder::gaussianSpecVolumeIntegral(folded);
    if (!(vol < 0.0) || std::fabs(vol) > 5000.0) {
        std::printf("FAIL q+d fold volume integral: %.3f\n", vol);
        return 1;
    }
    if (V(0.0) >= 0.0) {
        std::printf("FAIL q+d fold V(0) should be negative, got %.3f\n", V(0.0));
        return 1;
    }
    return 0;
}

} // namespace

int main() {
    int failures = 0;
    failures += testSingleGaussianFolded();
    failures += testSingleGaussianNoFold();
    failures += testHalQuartetFolded();
    failures += testQdCentralFolded();

    if (failures) {
        std::printf("test_folded_gaussian: %d failures\n", failures);
        return EXIT_FAILURE;
    }
    std::printf("test_folded_gaussian: OK\n");
    return EXIT_SUCCESS;
}
