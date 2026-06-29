#ifndef FEMTO_POTENTIAL_UTILS_HPP
#define FEMTO_POTENTIAL_UTILS_HPP

namespace femto {

// Tune attractive Gaussian depth V0 (>0) so extracted f0 matches targetF0 (>0).
double tuneV0(double targetF0, double b_fm, double mu_MeV,
              double& f0_out, double& d0_out);

double bindingEnergyERE(double f0, double d0, double mu_MeV);

} // namespace femto

#endif
