#ifndef FEMTO_CF_FITTER_ROOT_HPP
#define FEMTO_CF_FITTER_ROOT_HPP

#include "femto/FemtoModels.hpp"
#include "femto/Scenario.hpp"
#include <memory>
#include <string>
#include <vector>

class TF1;
class TH1D;

namespace femto {

enum class FitMode { LL, KP, Both };

struct FitOptions {
    FitMode mode = FitMode::Both;
    std::string channel = "phi_proton";
    std::string scenario = "HAL";
    std::string configRoot;
    double kFitMax_MeV = 250.0;
    bool fixLambda = false;
    double lambdaFixed = 0.6;
    std::vector<std::string> activeSpinChannels;
};

struct FitResultSummary {
    double R = 0, R_err = 0;
    double f0_re = 0, f0_re_err = 0;
    double f0_im = 0, f0_im_err = 0;
    double d0 = 0, d0_err = 0;
    double lambda = 0, lambda_err = 0;
    double N = 0, b1 = 0;
    double chi2_ndf_ll = -1;
    double chi2_ndf_kp = -1;
    bool hasLL = false;
    bool hasKP = false;
};

struct FitResult {
    TF1* fLL = nullptr;
    TF1* fKP = nullptr;
    FitResultSummary summary;
};

class FitSession {
public:
    explicit FitSession(FitOptions opts);

    FitResult fit(TH1D* hCF) const;
    TH1D* makeDemoHistogram() const;

private:
    FitOptions opts_;
};

std::string fitModeToString(FitMode m);
FitMode parseFitMode(const std::string& s);

} // namespace femto

#endif
