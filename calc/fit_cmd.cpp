#include "femto/CfFitterRoot.hpp"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>

namespace {

void usage() {
    std::printf("cf-fit — femtoscopic CF fitter\n\n");
    std::printf("Usage:\n");
    std::printf("  cf-fit [--demo] [--mode ll|kp|both] [--channel NAME] [--scenario NAME]\n");
    std::printf("         [--config DIR] [--input FILE] [--hist NAME] [--output PDF]\n");
}

} // namespace

int main(int argc, char** argv) {
    femto::FitOptions opts;
    std::string input;
    std::string histName = "hCF";
    std::string output = "cf_fit.pdf";
    bool demo = false;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--demo") demo = true;
        else if (a == "--mode" && i + 1 < argc) opts.mode = femto::parseFitMode(argv[++i]);
        else if (a == "--channel" && i + 1 < argc) opts.channel = argv[++i];
        else if (a == "--scenario" && i + 1 < argc) opts.scenario = argv[++i];
        else if ((a == "--config" || a == "-c") && i + 1 < argc) opts.configRoot = argv[++i];
        else if ((a == "--input" || a == "-i") && i + 1 < argc) input = argv[++i];
        else if (a == "--hist" && i + 1 < argc) histName = argv[++i];
        else if ((a == "--output" || a == "-o") && i + 1 < argc) output = argv[++i];
        else if (a == "--help" || a == "-h") {
            usage();
            return 0;
        } else {
            std::cerr << "unknown argument: " << a << std::endl;
            usage();
            return 1;
        }
    }

    try {
        femto::FitSession session(opts);
        TH1D* hCF = nullptr;
        bool owned = false;
        TFile* f = nullptr;

        if (demo || input.empty()) {
            hCF = session.makeDemoHistogram();
            owned = true;
            std::printf("demo data: channel=%s scenario=%s mode=%s\n", opts.channel.c_str(),
                        opts.scenario.c_str(), femto::fitModeToString(opts.mode).c_str());
        } else {
            f = TFile::Open(input.c_str());
            if (!f || f->IsZombie()) throw std::runtime_error("cannot open " + input);
            hCF = dynamic_cast<TH1D*>(f->Get(histName.c_str()));
            if (!hCF) throw std::runtime_error("histogram not found: " + histName);
        }

        femto::FitResult res = session.fit(hCF);

        if (res.summary.hasLL) {
            std::printf("LL: R=%.3f+-%.3f fm  f0=%.3f+-%.3f fm  d0=%.2f+-%.2f fm  chi2/ndf=%.2f\n",
                        res.summary.R, res.summary.R_err, res.summary.f0_re, res.summary.f0_re_err,
                        res.summary.d0, res.summary.d0_err, res.summary.chi2_ndf_ll);
        }
        if (res.summary.hasKP) {
            std::printf("KP: R=%.3f+-%.3f fm  lambda=%.3f+-%.3f  chi2/ndf=%.2f\n",
                        res.summary.R, res.summary.R_err, res.summary.lambda, res.summary.lambda_err,
                        res.summary.chi2_ndf_kp);
        }

        auto* c1 = new TCanvas("c1", "CF fit", 700, 550);
        hCF->SetMinimum(0.4);
        hCF->SetMarkerStyle(20);
        hCF->Draw("E");
        auto* leg = new TLegend(0.5, 0.18, 0.88, 0.42);
        leg->AddEntry(hCF, "data", "pe");
        if (res.fLL) {
            res.fLL->SetLineColor(kRed);
            res.fLL->Draw("same");
            leg->AddEntry(res.fLL, "LL fit", "l");
        }
        if (res.fKP) {
            res.fKP->SetLineColor(kBlue);
            res.fKP->SetLineStyle(2);
            res.fKP->Draw("same");
            leg->AddEntry(res.fKP, "KP fit", "l");
        }
        leg->Draw();
        c1->SaveAs(output.c_str());
        std::printf("wrote %s\n", output.c_str());

        if (owned) delete hCF;
        if (f) f->Close();
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << std::endl;
        return 1;
    }
    return 0;
}
