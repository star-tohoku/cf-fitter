#include "calc_common.hpp"
#include "femto/ChannelRegistry.hpp"
#include "femto/Scenario.hpp"
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<std::string> splitComma(const std::string& s) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) out.push_back(item);
    }
    return out;
}

const char* tableLabelForScenario(const std::string& name) {
    if (name == "q_nofold") return "quartet only / no fold";
    if (name == "q_fold") return "quartet only / fold";
    if (name == "qd_nofold") return "quartet+doublet / no fold";
    if (name == "q_chizzali_tpe_nofold") return "Chizzali TPE quartet / no fold";
    if (name == "q_chizzali_tpe_fold") return "Chizzali TPE quartet / fold";
    return "";
}

void writeField(FILE* f, double value, bool present, int precision = 3) {
    if (!present) std::fprintf(f, ",");
    else std::fprintf(f, ",%.*f", precision, value);
}

} // namespace

int runPhiAlphaSummary(const std::vector<std::string>& args) {
    std::string configRoot = getConfigRoot(args);
    std::string outFile = "figures/phi_alpha_potential_summary.csv";
    std::string scenarios =
        "q_nofold,q_fold,qd_nofold,qd_fold,"
        "qd_b040_nofold,qd_b040_fold,qd_b070_nofold,qd_b070_fold";

    for (std::size_t i = 0; i + 1 < args.size(); ++i) {
        if (args[i] == "--output" || args[i] == "-o") outFile = args[i + 1];
        if (args[i] == "--scenarios") scenarios = args[i + 1];
    }

    femto::ChannelRegistry registry(configRoot);
    const femto::ChannelSpec& ch = registry.get("phi_alpha");
    const auto scenNames = splitComma(scenarios);
    if (scenNames.empty())
        throw std::runtime_error("no scenarios provided for phi-alpha summary");

    FILE* f = std::fopen(outFile.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + outFile);

    std::fprintf(f,
                 "scenario,fold,doublet_b_fm,V0_MeV,f0_fm,d0_fm,BE_MeV,"
                 "V_at_0_MeV,volume_integral_MeV_fm3,table_label,notes\n");

    for (const auto& sn : scenNames) {
        femto::Scenario sc = femto::loadScenarioByName("phi_alpha", sn, configRoot);
        femto::FoldedPotentialSummary summary = femto::summarizePhiAlphaPotential(sc, ch);

        std::fprintf(f, "%s,%s", sc.name.c_str(), summary.fold ? "true" : "false");
        if (summary.hasDoublet)
            std::fprintf(f, ",%.2f,%.1f", summary.doublet_b_fm, summary.doublet_V0_MeV);
        else
            std::fprintf(f, ",,");

        std::fprintf(f, ",%.3f,%.3f", summary.f0_fm, summary.d0_fm);
        writeField(f, summary.BE_MeV, summary.hasBE, 2);
        std::fprintf(f, ",%.1f,%.1f", summary.V_at_0_MeV, summary.volume_integral_MeV_fm3);

        const char* label = tableLabelForScenario(sc.name);
        std::fprintf(f, ",%s", label);
        if (!summary.notes.empty())
            std::fprintf(f, ",%s", summary.notes.c_str());
        else if (*label == '\0')
            std::fprintf(f, ",doublet range band member");
        else
            std::fprintf(f, ",");
        std::fprintf(f, "\n");
    }

    std::fclose(f);
    std::printf("wrote %s\n", outFile.c_str());
    return 0;
}
