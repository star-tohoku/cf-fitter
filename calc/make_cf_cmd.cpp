#include "calc_common.hpp"
#include "femto/ChannelRegistry.hpp"
#include "femto/ModelFactory.hpp"
#include "femto/PotentialUtils.hpp"
#include "femto/Scenario.hpp"
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Curve {
    std::string name;
    femto::LLModel ll;
    std::shared_ptr<femto::KPModel> kp;
};

std::vector<std::string> splitComma(const std::string& s) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) out.push_back(item);
    }
    return out;
}

void writeCSV(const char* fname, const std::vector<Curve>& curves,
              const std::vector<double>& Rs, double kMax, int nk) {
    FILE* f = std::fopen(fname, "w");
    if (!f) throw std::runtime_error(std::string("cannot write ") + fname);
    std::fprintf(f, "k_MeV");
    for (double R : Rs)
        for (const auto& c : curves) {
            std::fprintf(f, ",LL_%s_R%.1f", c.name.c_str(), R);
            if (c.kp) std::fprintf(f, ",KP_%s_R%.1f", c.name.c_str(), R);
        }
    std::fprintf(f, "\n");
    for (int i = 1; i <= nk; ++i) {
        double kMeV = kMax * i / nk;
        double k = kMeV / femto::HBARC;
        std::fprintf(f, "%.2f", kMeV);
        for (double R : Rs)
            for (const auto& c : curves) {
                std::fprintf(f, ",%.6f", c.ll.C(k, R));
                if (c.kp) std::fprintf(f, ",%.6f", c.kp->C(k, R));
            }
        std::fprintf(f, "\n");
    }
    std::fclose(f);
    std::printf("wrote %s\n", fname);
}

void printScenarioInfo(const femto::Scenario& sc, const femto::ChannelSpec& ch,
                       const femto::BuiltModels& models) {
    const double mu = femto::reducedMass(ch);
    std::printf("  scenario %-14s mu=%.1f MeV", sc.name.c_str(), mu);
    for (const auto& spin : sc.resolve(ch, mu)) {
        if (!spin.interacting) continue;
        std::printf("  [%s] f0=%+.3f%+.3fi d0=%.3f", spin.spin.c_str(), spin.f0_re,
                    spin.f0_im, spin.d0_fm);
        if (spin.f0_re < 0) {
            double be = femto::bindingEnergyERE(spin.f0_re, spin.d0_fm, mu);
            if (be > 0) std::printf(" BE~%.2f MeV", be);
        }
    }
    std::printf("  %s\n", models.llOnly ? "(LL only)" : "");
}

} // namespace

int runMakeCf(const std::vector<std::string>& args) {
    std::string channel = "phi_proton";
    std::string scenarios;
    std::string configRoot = getConfigRoot(args);
    std::string outFile;

    for (std::size_t i = 0; i + 1 < args.size(); ++i) {
        if (args[i] == "--channel") channel = args[i + 1];
        if (args[i] == "--scenarios") scenarios = args[i + 1];
        if (args[i] == "--output" || args[i] == "-o") outFile = args[i + 1];
    }

    femto::ChannelRegistry registry(configRoot);
    const femto::ChannelSpec& ch = registry.get(channel);

    std::vector<std::string> scenNames;
    if (scenarios.empty()) scenNames = femto::listScenariosForChannel(channel, configRoot);
    else scenNames = splitComma(scenarios);

    if (scenNames.empty())
        throw std::runtime_error("no scenarios for channel " + channel);

    std::printf("== %s ==\n", channel.c_str());
    std::vector<Curve> curves;
    for (const auto& sn : scenNames) {
        femto::Scenario sc = femto::loadScenarioByName(channel, sn, configRoot);
        femto::BuiltModels models = femto::buildModels(sc, ch);
        printScenarioInfo(sc, ch, models);
        Curve c;
        c.name = sc.name;
        c.ll = models.ll;
        c.kp = models.kp;
        curves.push_back(std::move(c));
    }

    std::vector<double> Rs;
    double kMax = 300.0;
    int nk = 150;
    if (channel == "phi_alpha") {
        Rs = {1.2, 2.5, 5.0};
        kMax = 200.0;
    } else {
        Rs = {1.2, 3.0};
    }

    if (outFile.empty()) {
        if (channel == "phi_proton") outFile = "cf_phip.csv";
        else if (channel == "phi_alpha") outFile = "cf_phialpha.csv";
        else outFile = "cf_" + channel + ".csv";
    }

    writeCSV(outFile.c_str(), curves, Rs, kMax, nk);
    return 0;
}
