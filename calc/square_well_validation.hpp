#pragma once

#include "femto/FemtoModels.hpp"
#include <cmath>
#include <vector>

namespace square_well_validation {

struct PhasePoint {
    double k;
    double deltaNumeric;
    double deltaAnalytic;
    double relDiffPercent;
    double tolerancePercent;
    bool strict;
};

struct Result {
    std::vector<PhasePoint> phasePoints;
    double f0Numeric;
    double f0Analytic;
    double d0Numeric;
    double f0RelDiffPercent;
    bool passed;
};

inline double relDiffPercent(double numeric, double expected) {
    return 100.0 * std::fabs(numeric - expected) / std::fabs(expected);
}

inline Result run() {
    using namespace femto;

    constexpr double mu = 488.6;
    constexpr double v0Abs = 40.0;
    constexpr double radius = 1.0;

    auto Vsq = [](double r) { return r < radius ? -v0Abs : 0.0; };
    RadialSolverS solver(Vsq, mu, 60.0, 0.002, 12.0);

    const double K = std::sqrt(2.0 * mu * v0Abs) / HBARC;

    Result out;
    out.passed = true;
    for (double k : {0.02, 0.05, 0.10, 0.30}) {
        const double deltaAnalytic =
            -k * radius + std::atan((k / K) * std::tan(K * radius));
        const double deltaNumeric = solver.solve(k).delta;
        const bool strict = k < 0.30;
        const double tolerance = (k <= 0.05) ? 1.0 : (strict ? 2.0 : 0.0);
        const double rel = relDiffPercent(deltaNumeric, deltaAnalytic);
        out.phasePoints.push_back({k, deltaNumeric, deltaAnalytic, rel, tolerance, strict});
        if (strict && rel > tolerance) out.passed = false;
    }

    solver.scatteringParams(out.f0Numeric, out.d0Numeric);
    const double a0Analytic = radius - std::tan(K * radius) / K;
    out.f0Analytic = -a0Analytic;
    out.f0RelDiffPercent = relDiffPercent(out.f0Numeric, out.f0Analytic);
    if (out.f0RelDiffPercent > 1.0) out.passed = false;

    return out;
}

} // namespace square_well_validation
