#include "femto/PotentialUtils.hpp"
#include "femto/FemtoModels.hpp"

namespace femto {

double tuneV0(double targetF0, double b_fm, double mu_MeV, double& f0_out, double& d0_out) {
    auto f0of = [&](double V0) {
        auto s = RadialSolverS(gaussianPotential({{-V0, b_fm}}), mu_MeV);
        double f0 = 0, d0 = 0;
        s.scatteringParams(f0, d0);
        return f0;
    };
    double Vlo = 0.1, Vhi = 1.0;
    while (f0of(Vhi) > 0 && f0of(Vhi) < targetF0) {
        Vlo = Vhi;
        Vhi *= 1.5;
    }
    for (int it = 0; it < 60; ++it) {
        double Vm = 0.5 * (Vlo + Vhi);
        double f = f0of(Vm);
        if (f > 0 && f < targetF0) Vlo = Vm;
        else Vhi = Vm;
    }
    double V0 = 0.5 * (Vlo + Vhi);
    auto s = RadialSolverS(gaussianPotential({{-V0, b_fm}}), mu_MeV);
    s.scatteringParams(f0_out, d0_out);
    return V0;
}

double bindingEnergyERE(double f0, double d0, double mu_MeV) {
    double disc = 1.0 + 2.0 * d0 / f0;
    if (disc < 0) return -1;
    double kap = (1.0 - std::sqrt(disc)) / d0;
    if (kap <= 0) kap = (1.0 + std::sqrt(disc)) / d0;
    return kap * kap * HBARC * HBARC / (2.0 * mu_MeV);
}

} // namespace femto
