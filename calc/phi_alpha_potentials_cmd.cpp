#include "calc_common.hpp"
#include "femto/ChannelRegistry.hpp"
#include "femto/PotentialBuilder.hpp"
#include "femto/Scenario.hpp"
#include <cmath>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

const femto::SpinInteractionSpec* foldedSpin(const femto::Scenario& sc) {
    for (const auto& si : sc.channels) {
        if (si.potentialType == femto::PotentialType::Folded)
            return &si;
    }
    return nullptr;
}

void writeBareCsv(const std::string& path,
                  const std::function<double(double)>& vq,
                  const std::function<double(double)>& vb040,
                  const std::function<double(double)>& vb055,
                  const std::function<double(double)>& vb070) {
    FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + path);
    std::fprintf(f, "r_fm,V_bare_q_MeV,V_bare_qd_b040_MeV,V_bare_qd_b055_MeV,V_bare_qd_b070_MeV\n");
    for (int i = 0; i <= 300; ++i) {
        const double r = 0.01 * static_cast<double>(i);
        std::fprintf(f, "%.2f,%.8g,%.8g,%.8g,%.8g\n",
                     r, vq(r), vb040(r), vb055(r), vb070(r));
    }
    std::fclose(f);
}

void writeFoldedCsv(const std::string& path,
                    const std::function<double(double)>& vq,
                    const std::function<double(double)>& vf040,
                    const std::function<double(double)>& vf055,
                    const std::function<double(double)>& vf070) {
    FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + path);
    std::fprintf(f, "R_fm,V_fold_q_MeV,V_fold_qd_b040_MeV,V_fold_qd_b055_MeV,V_fold_qd_b070_MeV\n");
    for (int i = 0; i <= 300; ++i) {
        const double R = 0.02 * static_cast<double>(i);
        std::fprintf(f, "%.2f,%.8g,%.8g,%.8g,%.8g\n",
                     R, vq(R), vf040(R), vf055(R), vf070(R));
    }
    std::fclose(f);
}

void writeSummaryRow(FILE* f,
                     const std::string& scenario,
                     const std::string& potType,
                     double b_d,
                     const std::function<double(double)>& V,
                     const femto::FoldedPotentialSummary& summary) {
    std::fprintf(f, "%s,%s,", scenario.c_str(), potType.c_str());
    if (b_d >= 0.0)
        std::fprintf(f, "%.2f,", b_d);
    else
        std::fprintf(f, ",");
    std::fprintf(f, "%.4f,%.4f,%.4f,%.4f,%.2f,%.3f,%.3f,",
                 V(0.0), V(0.5), V(1.0), V(2.0),
                 summary.volume_integral_MeV_fm3,
                 summary.f0_fm, summary.d0_fm);
    if (summary.hasBE)
        std::fprintf(f, "%.2f\n", summary.BE_MeV);
    else
        std::fprintf(f, "\n");
}

} // namespace

int runPhiAlphaPotentials(const std::vector<std::string>& args) {
    std::string configRoot = getConfigRoot(args);
    std::string bareOut = "figures/phi_alpha_doublet_bare_potential_shapes.csv";
    std::string foldedOut = "figures/phi_alpha_doublet_folded_potential_shapes.csv";
    std::string summaryOut = "figures/phi_alpha_doublet_potential_shapes_summary.csv";

    for (std::size_t i = 0; i + 1 < args.size(); ++i) {
        if (args[i] == "--bare-output") bareOut = args[i + 1];
        if (args[i] == "--folded-output") foldedOut = args[i + 1];
        if (args[i] == "--summary-output") summaryOut = args[i + 1];
    }

    femto::ChannelRegistry registry(configRoot);
    const femto::ChannelSpec& ch = registry.get("phi_alpha");

    femto::Scenario sc_q = femto::loadScenarioByName("phi_alpha", "q_fold", configRoot);
    femto::Scenario sc_b040 = femto::loadScenarioByName("phi_alpha", "qd_b040_fold", configRoot);
    femto::Scenario sc_b055 = femto::loadScenarioByName("phi_alpha", "qd_fold", configRoot);
    femto::Scenario sc_b070 = femto::loadScenarioByName("phi_alpha", "qd_b070_fold", configRoot);

    const auto* fold_q = foldedSpin(sc_q);
    const auto* fold_b040 = foldedSpin(sc_b040);
    const auto* fold_b055 = foldedSpin(sc_b055);
    const auto* fold_b070 = foldedSpin(sc_b070);
    if (!fold_q || !fold_b040 || !fold_b055 || !fold_b070)
        throw std::runtime_error("missing folded phi-alpha scenario");

    auto vBareQ = femto::barePhiNCentralPotential(fold_q->folded.central);
    auto vBare040 = femto::barePhiNCentralPotential(fold_b040->folded.central);
    auto vBare055 = femto::barePhiNCentralPotential(fold_b055->folded.central);
    auto vBare070 = femto::barePhiNCentralPotential(fold_b070->folded.central);

    auto vFoldQ = femto::phiAlphaFoldedPotential(sc_q);
    auto vFold040 = femto::phiAlphaFoldedPotential(sc_b040);
    auto vFold055 = femto::phiAlphaFoldedPotential(sc_b055);
    auto vFold070 = femto::phiAlphaFoldedPotential(sc_b070);

    writeBareCsv(bareOut, vBareQ, vBare040, vBare055, vBare070);
    writeFoldedCsv(foldedOut, vFoldQ, vFold040, vFold055, vFold070);

    FILE* f = std::fopen(summaryOut.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + summaryOut);
    std::fprintf(f,
                 "scenario,potential_type,doublet_b_fm,V_at_0_MeV,V_at_0p5_MeV,"
                 "V_at_1_MeV,V_at_2_MeV,volume_integral_MeV_fm3,f0_fm,d0_fm,BE_MeV\n");

    auto sumQ = femto::summarizePhiAlphaPotential(sc_q, ch);
    auto sum040 = femto::summarizePhiAlphaPotential(sc_b040, ch);
    auto sum055 = femto::summarizePhiAlphaPotential(sc_b055, ch);
    auto sum070 = femto::summarizePhiAlphaPotential(sc_b070, ch);

    writeSummaryRow(f, "q_fold", "bare_phiN", -1.0, vBareQ, sumQ);
    writeSummaryRow(f, "qd_b040_fold", "bare_phiN", 0.40, vBare040, sum040);
    writeSummaryRow(f, "qd_fold", "bare_phiN", 0.55, vBare055, sum055);
    writeSummaryRow(f, "qd_b070_fold", "bare_phiN", 0.70, vBare070, sum070);
    writeSummaryRow(f, "q_fold", "folded_phiAlpha", -1.0, vFoldQ, sumQ);
    writeSummaryRow(f, "qd_b040_fold", "folded_phiAlpha", 0.40, vFold040, sum040);
    writeSummaryRow(f, "qd_fold", "folded_phiAlpha", 0.55, vFold055, sum055);
    writeSummaryRow(f, "qd_b070_fold", "folded_phiAlpha", 0.70, vFold070, sum070);
    std::fclose(f);

    std::printf("wrote %s\n", bareOut.c_str());
    std::printf("wrote %s\n", foldedOut.c_str());
    std::printf("wrote %s\n", summaryOut.c_str());
    return 0;
}
