#ifndef FEMTO_POTENTIAL_BUILDER_HPP
#define FEMTO_POTENTIAL_BUILDER_HPP

#include <functional>
#include <vector>

namespace femto {

struct GaussianTerm {
    double V_MeV = 0.0;
    double b_fm = 1.0;
};

struct GaussianSpec {
    std::vector<GaussianTerm> terms;
};

struct GaussianDensitySpec {
    double A = 4.0;
    double rms_fm = 1.56;
    double b_fm = 0.0;
};

struct FoldedPotentialSpec {
    GaussianSpec phiNCentral;
    GaussianDensitySpec density;
    bool fold = true;
};

struct FoldingOptions {
    double rMax_fm = 20.0;
    double dr_fm = 0.02;
};

class PotentialBuilder {
public:
    static std::function<double(double)> gaussian(const GaussianSpec& spec);
    static GaussianSpec foldedGaussianSpec(const FoldedPotentialSpec& spec);
    static std::function<double(double)> foldedGaussian(const FoldedPotentialSpec& spec);
    // Legacy numeric-folding hook; prefer foldedGaussian for Gaussian-sum inputs.
    static std::function<double(double)> folded(
        const std::function<double(double)>& V2body,
        const std::function<double(double)>& density,
        const FoldingOptions& opts = {});
};

} // namespace femto

#endif
