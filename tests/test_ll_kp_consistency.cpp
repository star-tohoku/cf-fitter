#include "femto/FemtoModels.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>

int main() {
    using namespace femto;

    const double mPhi = 1019.461, mP = 938.272;
    const double mu = mPhi * mP / (mPhi + mP);

    auto V = gaussianPotential({{-40.0, 1.0}});
    auto solver = std::make_shared<RadialSolverS>(V, mu);

    double f0 = 0, d0 = 0;
    solver->scatteringParams(f0, d0);

    LLModel ll = makeLL_identicalBosons({f0, 0.0}, d0);
    KPModel kp;
    kp.identical = true;
    kp.qsAmplitude = +1.0;
    kp.channels = {{1.0, solver}};

    int failures = 0;
    for (double R : {2.0, 3.0, 5.0}) {
        for (double kMeV : {10., 40., 80., 120., 200.}) {
            double k = kMeV / HBARC;
            double cll = ll.C(k, R);
            double ckp = kp.C(k, R);
            double rel = std::fabs(100.0 * (ckp - cll) / cll);
            if (rel > 1.0) {
                std::printf("FAIL R=%.1f k=%.0f MeV diff=%.2f%%\n", R, kMeV, rel);
                ++failures;
            }
        }
    }

    if (failures) {
        std::printf("test_ll_kp_consistency: %d failures\n", failures);
        return EXIT_FAILURE;
    }
    std::printf("test_ll_kp_consistency: OK\n");
    return EXIT_SUCCESS;
}
