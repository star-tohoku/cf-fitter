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

struct NumericFoldingOptions {
    double rGridMax_fm = 20.0;
    double dr_fm = 0.02;
    double sMax_fm = 16.0;
    double ds_fm = 0.02;
    int nMu = 80;
    double volumeRMax_fm = 20.0;
    double volumeDr_fm = 0.02;
};

class PotentialBuilder {
public:
    static std::function<double(double)> gaussian(const GaussianSpec& spec);
    static GaussianSpec foldedGaussianSpec(const FoldedPotentialSpec& spec);
    static std::function<double(double)> foldedGaussian(const FoldedPotentialSpec& spec);
    static double gaussianSpecV0(const GaussianSpec& spec);
    static double gaussianSpecVolumeIntegral(const GaussianSpec& spec);
    // Legacy hook kept for API compatibility.
    static std::function<double(double)> folded(
        const std::function<double(double)>& V2body,
        const std::function<double(double)>& density,
        const FoldingOptions& opts = {});
    static std::function<double(double)> foldedPhiN(
        const std::function<double(double)>& V2body,
        const GaussianDensitySpec& density,
        bool fold,
        const NumericFoldingOptions& opts = {});
    static double phiNVolumeIntegral(
        const std::function<double(double)>& V2body,
        const NumericFoldingOptions& opts = {});
    static double potentialVolumeIntegral(
        const std::function<double(double)>& V,
        double rMax_fm, double dr_fm);
};

} // namespace femto

#endif
