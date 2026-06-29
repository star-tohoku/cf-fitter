#include "square_well_validation.hpp"
#include <cstdio>
#include <cstdlib>

int main() {
    const auto result = square_well_validation::run();
    int failures = 0;

    for (const auto& point : result.phasePoints) {
        if (!point.strict) continue;
        if (point.relDiffPercent > point.tolerancePercent) {
            std::printf("FAIL k=%.2f fm^-1 delta_num=%+.5f delta_analytic=%+.5f "
                        "diff=%.2f%% tolerance=%.2f%%\n",
                        point.k, point.deltaNumeric, point.deltaAnalytic,
                        point.relDiffPercent, point.tolerancePercent);
            ++failures;
        }
    }

    if (result.f0RelDiffPercent > 1.0) {
        std::printf("FAIL f0_numeric=%.4f fm f0_analytic=%.4f fm diff=%.2f%%\n",
                    result.f0Numeric, result.f0Analytic, result.f0RelDiffPercent);
        ++failures;
    }

    if (failures) {
        std::printf("test_square_well_solver: %d failures\n", failures);
        return EXIT_FAILURE;
    }

    std::printf("test_square_well_solver: OK "
                "(f0_numeric=%.4f fm, f0_analytic=%.4f fm, diff=%.2f%%)\n",
                result.f0Numeric, result.f0Analytic, result.f0RelDiffPercent);
    return EXIT_SUCCESS;
}
