#include "calc_common.hpp"
#include "square_well_validation.hpp"
#include "femto/ChannelRegistry.hpp"
#include "femto/FemtoModels.hpp"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>

std::string getConfigRoot(const std::vector<std::string>& args) {
    for (std::size_t i = 0; i + 1 < args.size(); ++i)
        if (args[i] == "--config" || args[i] == "-c") return args[i + 1];
    return femto::defaultConfigRoot();
}

int runValidate(const std::vector<std::string>&) {
    using namespace femto;

    int fail = 0;

    {
        const auto sq = square_well_validation::run();
        std::printf("Square-well solver benchmark: V0 = -40 MeV, b = 1 fm, mu = 488.6 MeV\n");
        std::printf("  k*[fm^-1]   delta_num   delta_analytic   rel.diff[%%]\n");
        for (const auto& p : sq.phasePoints) {
            std::printf("      %4.2f     %+9.5f      %+9.5f      %6.2f%s\n",
                        p.k, p.deltaNumeric, p.deltaAnalytic, p.relDiffPercent,
                        p.strict ? "" : "  diagnostic");
        }
        std::printf("  f0_numeric = %.4f fm, f0_analytic = %.4f fm, diff = %.2f%%\n",
                    sq.f0Numeric, sq.f0Analytic, sq.f0RelDiffPercent);
        std::printf("  d0_numeric = %.4f fm\n", sq.d0Numeric);
        if (!sq.passed) ++fail;
    }

    std::printf("\nLL/KP consistency benchmark\n");

    const double mPhi = 1019.461, mP = 938.272;
    const double mu = mPhi * mP / (mPhi + mP);

    auto V = gaussianPotential({{-40.0, 1.0}});
    auto solver = std::make_shared<RadialSolverS>(V, mu);

    double f0 = 0, d0 = 0;
    solver->scatteringParams(f0, d0);
    std::printf("Test system: identical spin-0 bosons, mu = %.2f MeV\n", mu);
    std::printf("Extracted ERE params (femto convention): f0 = %.4f fm, d0 = %.4f fm\n",
                f0, d0);

    LLModel ll = makeLL_identicalBosons({f0, 0.0}, d0);

    KPModel kp;
    kp.identical = true;
    kp.qsAmplitude = +1.0;
    kp.channels = {{1.0, solver}};

    for (double R : {1.2, 2.0, 3.0, 5.0}) {
        std::printf("\n  R = %.1f fm\n  k*[MeV/c]   C_LL      C_KP      diff[%%]\n", R);
        for (double kMeV : {10., 20., 40., 60., 80., 120., 160., 200.}) {
            double k = kMeV / HBARC;
            double cll = ll.C(k, R);
            double ckp = kp.C(k, R);
            double diff = 100.0 * (ckp - cll) / cll;
            std::printf("  %8.0f   %.5f   %.5f   %+6.2f\n", kMeV, cll, ckp, diff);
            if (R >= 2.0 && std::fabs(diff) > 2.0) ++fail;
        }
    }
    return fail > 0 ? 1 : 0;
}
