#include "femto/FemtoModels.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>

int main() {
    using namespace femto;

    struct Ref { double r; double v; };
    const Ref refs[] = {
        {0.10, -347.118}, {0.20, -141.886}, {0.30, -82.123},
        {0.50, -51.481},  {0.80, -30.573},  {1.00, -18.410},
        {1.20, -10.455},  {1.50,  -4.514},  {2.00,  -1.265},
    };

    auto V = halQuartetChizzaliTa12TPEPotential();
    int failures = 0;
    constexpr double tol = 0.02;  // MeV

    if (std::fabs(V(0.0) - (-537.0)) > 1e-6) {
        std::printf("FAIL r=0: V=%.6f (expected -537 from Gaussians, TPE=0)\n", V(0.0));
        ++failures;
    }

    for (const auto& ref : refs) {
        double got = V(ref.r);
        double diff = std::fabs(got - ref.v);
        if (diff > tol) {
            std::printf("FAIL r=%.2f fm: V=%.6f ref=%.3f diff=%.4f MeV\n",
                        ref.r, got, ref.v, diff);
            ++failures;
        }
    }

    if (failures) {
        std::printf("test_chizzali_ta12_tpe_potential: %d failures\n", failures);
        return EXIT_FAILURE;
    }
    std::printf("test_chizzali_ta12_tpe_potential: OK\n");
    return EXIT_SUCCESS;
}
