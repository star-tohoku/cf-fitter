#include "femto/FemtoModels.hpp"
#include "femto/PotentialBuilder.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>

namespace {

bool relFail(double got, double ref, double relTol = 1e-3) {
    if (std::fabs(ref) < 1e-9) return std::fabs(got - ref) > 1e-6;
    return std::fabs((got - ref) / ref) > relTol;
}

bool mixedFail(double got, double ref, double relTol = 5e-3, double absTol = 0.05) {
    if (std::fabs(ref) < absTol) return std::fabs(got - ref) > 6e-4;
    return relFail(got, ref, relTol);
}

femto::NumericFoldingOptions testFoldingOpts() {
    femto::NumericFoldingOptions opts;
    opts.dr_fm = 0.02;
    opts.ds_fm = 0.02;
    opts.nMu = 96;
    opts.sMax_fm = 20.0;
    opts.volumeRMax_fm = 20.0;
    opts.volumeDr_fm = 0.02;
    return opts;
}

int testGaussianFolded() {
    using namespace femto;

    const double a = -120.0;
    const double b = 0.45;
    GaussianDensitySpec density;
    density.A = 4.0;
    density.rms_fm = 1.56;

    FoldedPotentialSpec analyticSpec;
    analyticSpec.fold = true;
    analyticSpec.density = density;
    analyticSpec.phiNCentral.terms.push_back({a, b});

    const auto opts = testFoldingOpts();
    auto VN = gaussianPotential({{a, b}});
    auto Vnum = PotentialBuilder::foldedPhiN(VN, density, true, opts);
    auto Vanalytic = PotentialBuilder::foldedGaussian(analyticSpec);

    int failures = 0;
    for (double R : {0.0, 0.5, 1.0, 2.0, 4.0}) {
        const double vn = Vnum(R);
        const double va = Vanalytic(R);
        if (!mixedFail(vn, va)) continue;
        std::printf("FAIL Gaussian fold R=%.1f: numeric=%.6f analytic=%.6f\n", R, vn, va);
        ++failures;
    }
    return failures;
}

int testGaussianNoFold() {
    using namespace femto;

    const double a = -90.0;
    const double b = 0.55;
    GaussianDensitySpec density;
    density.A = 4.0;
    density.rms_fm = 1.56;

    FoldedPotentialSpec analyticSpec;
    analyticSpec.fold = false;
    analyticSpec.density = density;
    analyticSpec.phiNCentral.terms.push_back({a, b});

    const auto opts = testFoldingOpts();
    auto VN = gaussianPotential({{a, b}});
    const double Gnum = PotentialBuilder::phiNVolumeIntegral(VN, opts);
    const double Ganalytic = a * std::pow(M_PI, 1.5) * b * b * b;
    if (relFail(Gnum, Ganalytic)) {
        std::printf("FAIL no-fold G: numeric=%.6f analytic=%.6f\n", Gnum, Ganalytic);
        return 1;
    }

    auto Vnum = PotentialBuilder::foldedPhiN(VN, density, false, opts);
    auto Vanalytic = PotentialBuilder::foldedGaussian(analyticSpec);
    int failures = 0;
    for (double R : {0.0, 0.8, 1.5, 3.0}) {
        const double vn = Vnum(R);
        const double va = Vanalytic(R);
        if (relFail(vn, va)) {
            std::printf("FAIL Gaussian no-fold R=%.1f: numeric=%.6f analytic=%.6f\n", R, vn, va);
            ++failures;
        }
    }
    return failures;
}

int testHalFitBFolded() {
    using namespace femto;

    GaussianDensitySpec density;
    density.A = 4.0;
    density.rms_fm = 1.56;

    FoldedPotentialSpec fs;
    fs.fold = true;
    fs.density = density;
    fs.phiNCentral.terms = {
        {-371.0, 0.15},
        {-50.0, 0.66},
        {-31.0, 1.09},
    };

    const auto opts = testFoldingOpts();
    auto VN = gaussianPotential({{-371.0, 0.15}, {-50.0, 0.66}, {-31.0, 1.09}});
    auto Vnum = PotentialBuilder::foldedPhiN(VN, density, true, opts);
    auto Vanalytic = PotentialBuilder::foldedGaussian(fs);

    int failures = 0;
    if (relFail(Vnum(0.0), Vanalytic(0.0), 5e-3)) {
        std::printf("FAIL HAL fold V(0): numeric=%.4f analytic=%.4f\n", Vnum(0.0), Vanalytic(0.0));
        ++failures;
    }
    const double volNum = PotentialBuilder::potentialVolumeIntegral(Vnum, 20.0, 0.02);
    const double volAnalytic = PotentialBuilder::gaussianSpecVolumeIntegral(
        PotentialBuilder::foldedGaussianSpec(fs));
    if (relFail(volNum, volAnalytic, 2e-2)) {
        std::printf("FAIL HAL fold volume: numeric=%.4f analytic=%.4f\n", volNum, volAnalytic);
        ++failures;
    }
    return failures;
}

int testHalQFoldCorrelations() {
    using namespace femto;

    const double mPhi = 1019.461;
    const double mAlpha = 3727.379;
    const double mu = mPhi * mAlpha / (mPhi + mAlpha);

    GaussianDensitySpec density;
    density.A = 4.0;
    density.rms_fm = 1.56;
    auto VN = gaussianPotential({{-371.0, 0.15}, {-50.0, 0.66}, {-31.0, 1.09}});

    const auto opts = testFoldingOpts();
    FoldedPotentialSpec fs;
    fs.fold = true;
    fs.density = density;
    fs.phiNCentral.terms = {
        {-371.0, 0.15},
        {-50.0, 0.66},
        {-31.0, 1.09},
    };

    auto solverNum = std::make_shared<RadialSolverS>(
        PotentialBuilder::foldedPhiN(VN, density, true, opts), mu);
    auto solverAnalytic = std::make_shared<RadialSolverS>(
        PotentialBuilder::foldedGaussian(fs), mu);

    double f0n = 0, d0n = 0, f0a = 0, d0a = 0;
    solverNum->scatteringParams(f0n, d0n);
    solverAnalytic->scatteringParams(f0a, d0a);

    LLModel llNum;
    llNum.channels = {{1.0, {f0n, 0.0}, d0n, true}};
    LLModel llAnalytic;
    llAnalytic.channels = {{1.0, {f0a, 0.0}, d0a, true}};

    KPModel kpNum;
    kpNum.channels = {{1.0, solverNum}};
    KPModel kpAnalytic;
    kpAnalytic.channels = {{1.0, solverAnalytic}};

    int failures = 0;
    const double R = 3.0;
    for (double kMeV : {2.0, 10.0, 50.0}) {
        const double k = kMeV / HBARC;
        const double cNum = kpNum.C(k, R);
        const double cAnalytic = kpAnalytic.C(k, R);
        if (relFail(cNum, cAnalytic, 2e-2)) {
            std::printf("FAIL q_fold KP k*=%.0f: numeric=%.6f analytic=%.6f\n",
                        kMeV, cNum, cAnalytic);
            ++failures;
        }
        const double lNum = llNum.C(k, R);
        const double lAnalytic = llAnalytic.C(k, R);
        if (relFail(lNum, lAnalytic, 2e-2)) {
            std::printf("FAIL q_fold LL k*=%.0f: numeric=%.6f analytic=%.6f\n",
                        kMeV, lNum, lAnalytic);
            ++failures;
        }
    }
    return failures;
}

int testChizzaliTpeVolumeStability() {
    using namespace femto;

    auto V = halQuartetChizzaliTa12TPEPotential();
    int failures = 0;
    double prev = 0.0;
    bool havePrev = false;
    for (double rMax : {8.0, 10.0, 12.0, 16.0}) {
        NumericFoldingOptions opts;
        opts.volumeRMax_fm = rMax;
        opts.volumeDr_fm = 0.02;
        const double G = PotentialBuilder::phiNVolumeIntegral(V, opts);
        if (havePrev && relFail(G, prev, 5e-3)) {
            std::printf("FAIL TPE G stability rMax=%.0f: G=%.4f prev=%.4f\n", rMax, G, prev);
            ++failures;
        }
        prev = G;
        havePrev = true;
    }
    return failures;
}

} // namespace

int main() {
    int failures = 0;
    failures += testGaussianFolded();
    failures += testGaussianNoFold();
    failures += testHalFitBFolded();
    failures += testHalQFoldCorrelations();
    failures += testChizzaliTpeVolumeStability();

    if (failures) {
        std::printf("test_numeric_folding: %d failures\n", failures);
        return EXIT_FAILURE;
    }
    std::printf("test_numeric_folding: OK\n");
    return EXIT_SUCCESS;
}
