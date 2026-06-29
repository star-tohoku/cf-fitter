#include "femto/PotentialBuilder.hpp"
#include "femto/FemtoModels.hpp"
#include <stdexcept>

namespace femto {

std::function<double(double)> PotentialBuilder::gaussian(const GaussianSpec& spec) {
    std::vector<std::pair<double, double>> terms;
    for (const auto& t : spec.terms)
        terms.emplace_back(t.V_MeV, t.b_fm);
    return gaussianPotential(terms);
}

std::function<double(double)> PotentialBuilder::folded(
    const std::function<double(double)>& V2body,
    const std::function<double(double)>& density,
    const FoldingOptions&) {
    (void)density;
    // Phase 5 stub: pass through 2-body potential until folding is implemented.
    return V2body;
}

} // namespace femto
