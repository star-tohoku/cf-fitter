#include "femto/CfFitterRoot.hpp"
#include "femto/ChannelRegistry.hpp"
#include "femto/ModelFactory.hpp"
#include "femto/Scenario.hpp"
#include "TF1.h"
#include "TH1D.h"
#include "TRandom3.h"
#include <cmath>
#include <stdexcept>

namespace femto {
namespace {

struct LLFitFunc {
    bool identical = false;
    double qsAmp = 0.0;
    double weight = 1.0;
    double operator()(double* x, double* p) const {
        const double kMeV = x[0];
        const double k = kMeV / HBARC;
        const double R = p[0], lam = p[4], N = p[5], b1 = p[6];
        LLModel m;
        m.identical = identical;
        m.qsAmplitude = qsAmp;
        m.channels = {{weight, {p[1], p[2]}, p[3], true}};
        const double c = m.C(k, R);
        return N * (1.0 + b1 * kMeV) * (1.0 + lam * (c - 1.0));
    }
};

struct KPFitFunc {
    std::shared_ptr<KPModel> model;
    double operator()(double* x, double* p) const {
        const double kMeV = x[0];
        const double k = kMeV / HBARC;
        const double R = p[0], lam = p[1], N = p[2], b1 = p[3];
        const double c = model->C(k, R);
        return N * (1.0 + b1 * kMeV) * (1.0 + lam * (c - 1.0));
    }
};

double spinWeight(const ChannelSpec& ch, const std::vector<std::string>& active,
                  std::string& chosenSpin) {
    if (!active.empty()) {
        for (const auto& name : active) {
            const SpinChannelSpec* sc = findSpinChannel(ch, name);
            if (sc) {
                chosenSpin = name;
                return sc->weight;
            }
        }
        throw std::runtime_error("active spin channel not found");
    }
    if (ch.name == "phi_proton") {
        chosenSpin = "quartet_4S32";
        return 2.0 / 3.0;
    }
    if (!ch.spinChannels.empty()) {
        chosenSpin = ch.spinChannels[0].name;
        return ch.spinChannels[0].weight;
    }
    return 1.0;
}

} // namespace

std::string fitModeToString(FitMode m) {
    switch (m) {
    case FitMode::LL: return "ll";
    case FitMode::KP: return "kp";
    case FitMode::Both: return "both";
    }
    return "both";
}

FitMode parseFitMode(const std::string& s) {
    if (s == "ll" || s == "LL") return FitMode::LL;
    if (s == "kp" || s == "KP") return FitMode::KP;
    return FitMode::Both;
}

FitSession::FitSession(FitOptions opts) : opts_(std::move(opts)) {
    if (opts_.configRoot.empty()) opts_.configRoot = defaultConfigRoot();
}

TH1D* FitSession::makeDemoHistogram() const {
    ChannelRegistry registry(opts_.configRoot);
    const ChannelSpec& ch = registry.get(opts_.channel);
    Scenario sc = loadScenarioByName(opts_.channel, opts_.scenario, opts_.configRoot);
    BuiltModels models = buildModels(sc, ch);
    if (!models.kp) throw std::runtime_error("demo requires KP scenario");

    const double Rtrue = 2.5;
    const double lamTrue = opts_.fixLambda ? opts_.lambdaFixed : 0.6;

    auto* hCF = new TH1D("hCF", "demo CF;k* [MeV/c];C(k*)", 60, 0, 300);
    TRandom3 rng(42);
    for (int i = 1; i <= hCF->GetNbinsX(); ++i) {
        double kMeV = hCF->GetBinCenter(i);
        double c = 1.0 + lamTrue * (models.kp->C(kMeV / HBARC, Rtrue) - 1.0);
        double err = 0.01 + 0.02 * std::exp(-kMeV / 80.0);
        hCF->SetBinContent(i, c + rng.Gaus(0, err));
        hCF->SetBinError(i, err);
    }
    return hCF;
}

FitResult FitSession::fit(TH1D* hCF) const {
    if (!hCF) throw std::runtime_error("null histogram");

    ChannelRegistry registry(opts_.configRoot);
    const ChannelSpec& ch = registry.get(opts_.channel);
    Scenario sc = loadScenarioByName(opts_.channel, opts_.scenario, opts_.configRoot);
    BuiltModels models = buildModels(sc, ch);

    FitResult result;
    const double kFitMax = opts_.kFitMax_MeV;

    if (opts_.mode == FitMode::LL || opts_.mode == FitMode::Both) {
        std::string spinName;
        const double weight = spinWeight(ch, opts_.activeSpinChannels, spinName);
        LLFitFunc llf{ch.identical, ch.qsAmplitude, weight};
        auto* fLL = new TF1("fLL", llf, 0.0, kFitMax, 7);
        fLL->SetParNames("R", "Re_f0", "Im_f0", "d0", "lambda", "N", "b1");
        fLL->SetParameters(2.0, 0.5, 0.0, 3.0,
                           opts_.fixLambda ? opts_.lambdaFixed : 0.6, 1.0, 0.0);
        fLL->FixParameter(2, 0.0);
        if (opts_.fixLambda) fLL->FixParameter(4, opts_.lambdaFixed);
        fLL->SetParLimits(0, 0.5, 8.0);
        fLL->SetParLimits(4, 0.0, 1.0);
        fLL->SetNpx(600);
        hCF->Fit(fLL, opts_.mode == FitMode::Both ? "RME0" : "RME0");
        result.fLL = fLL;
        result.summary.hasLL = true;
        result.summary.R = fLL->GetParameter(0);
        result.summary.R_err = fLL->GetParError(0);
        result.summary.f0_re = fLL->GetParameter(1);
        result.summary.f0_re_err = fLL->GetParError(1);
        result.summary.f0_im = fLL->GetParameter(2);
        result.summary.f0_im_err = fLL->GetParError(2);
        result.summary.d0 = fLL->GetParameter(3);
        result.summary.d0_err = fLL->GetParError(3);
        result.summary.lambda = fLL->GetParameter(4);
        result.summary.lambda_err = fLL->GetParError(4);
        result.summary.N = fLL->GetParameter(5);
        result.summary.b1 = fLL->GetParameter(6);
        if (fLL->GetNDF() > 0) result.summary.chi2_ndf_ll = fLL->GetChisquare() / fLL->GetNDF();
    }

    if (opts_.mode == FitMode::KP || opts_.mode == FitMode::Both) {
        if (!models.kp) throw std::runtime_error("KP fit requires scenario with potential");
        KPFitFunc kpf{models.kp};
        auto* fKP = new TF1("fKP", kpf, 0.0, kFitMax, 4);
        fKP->SetParNames("R", "lambda", "N", "b1");
        fKP->SetParameters(2.0, opts_.fixLambda ? opts_.lambdaFixed : 0.6, 1.0, 0.0);
        if (opts_.fixLambda) fKP->FixParameter(1, opts_.lambdaFixed);
        fKP->SetParLimits(0, 0.5, 8.0);
        fKP->SetParLimits(1, 0.0, 1.0);
        fKP->SetNpx(600);
        hCF->Fit(fKP, opts_.mode == FitMode::Both ? "RME0+" : "RME0");
        models.kp->clearCache();
        result.fKP = fKP;
        result.summary.hasKP = true;
        if (result.summary.R == 0) {
            result.summary.R = fKP->GetParameter(0);
            result.summary.R_err = fKP->GetParError(0);
        }
        result.summary.lambda = fKP->GetParameter(1);
        result.summary.lambda_err = fKP->GetParError(1);
        result.summary.N = fKP->GetParameter(2);
        result.summary.b1 = fKP->GetParameter(3);
        if (fKP->GetNDF() > 0) result.summary.chi2_ndf_kp = fKP->GetChisquare() / fKP->GetNDF();
    }

    return result;
}

} // namespace femto
