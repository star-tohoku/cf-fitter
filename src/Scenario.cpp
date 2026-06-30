#include "femto/Scenario.hpp"
#include "femto/ChannelRegistry.hpp"
#include "femto/PotentialBuilder.hpp"
#include "femto/PotentialUtils.hpp"
#include "femto/FemtoModels.hpp"
#include "femto/YamlReader.hpp"
#include <dirent.h>
#include <cmath>
#include <stdexcept>

namespace femto {
namespace {

PotentialType parsePotentialType(const std::string& s) {
    if (s == "gaussian") return PotentialType::Gaussian;
    if (s == "ere") return PotentialType::Ere;
    if (s == "folded") return PotentialType::Folded;
    if (s == "hal_quartet_fitB") return PotentialType::HalQuartetFitB;
    if (s == "hal_quartet_chizzali_ta12_tpe") return PotentialType::HalQuartetChizzaliTa12TPE;
    return PotentialType::Gaussian;
}

ScenarioMode parseScenarioMode(const std::string& s) {
    if (s == "ll_only") return ScenarioMode::LlOnly;
    return ScenarioMode::KpAndLl;
}

CouplingScheme parseCoupling(const std::string& s) {
    if (s == "effective") return CouplingScheme::EffectiveSingleChannel;
    if (s == "coupled") return CouplingScheme::CoupledRadial;
    return CouplingScheme::Independent;
}


FoldedScenarioSpec parseFolded(const YamlNode& pot) {
    FoldedScenarioSpec f;
    f.fold = yamlBool(pot, "fold", true);

    auto densityIt = pot.children.find("density");
    if (densityIt != pot.children.end()) {
        f.density.A = yamlDouble(densityIt->second, "A", f.density.A);
        f.density.rms_fm = yamlDouble(densityIt->second, "rms_fm", f.density.rms_fm);
        f.density.b_fm = yamlDouble(densityIt->second, "b_fm", f.density.b_fm);
    }

    auto centralIt = pot.children.find("phiN_central");
    if (centralIt != pot.children.end()) {
        const std::string mode = yamlScalar(centralIt->second, "mode", "quartet_only");
        if (mode == "quartet_plus_doublet" || mode == "spin_average") {
            f.central.quartetWeight = 2.0 / 3.0;
            f.central.doubletWeight = 1.0 / 3.0;
        } else if (mode == "quartet_only") {
            f.central.quartetWeight = 1.0;
            f.central.doubletWeight = 0.0;
        }
        f.central.quartetWeight = yamlDouble(centralIt->second, "quartet_weight", f.central.quartetWeight);
        f.central.doubletWeight = yamlDouble(centralIt->second, "doublet_weight", f.central.doubletWeight);
        f.central.doublet_b_fm = yamlDouble(centralIt->second, "doublet_b_fm", f.central.doublet_b_fm);
        f.central.doublet_target_f0 = yamlDouble(centralIt->second, "doublet_target_f0", f.central.doublet_target_f0);
    }
    return f;
}

GaussianPotentialSpec parseGaussian(const YamlNode& pot) {
    GaussianPotentialSpec g;
    g.V0_MeV = yamlDouble(pot, "V0_MeV", 0.0);
    g.b_fm = yamlDouble(pot, "b_fm", 1.0);
    g.VcritFraction = yamlDouble(pot, "VcritFraction", 0.0);
    auto tuneIt = pot.children.find("tune");
    if (tuneIt != pot.children.end()) {
        g.useTune = true;
        g.targetF0 = yamlDouble(tuneIt->second, "targetF0", 1.0);
        g.b_fm = yamlDouble(tuneIt->second, "b_fm", g.b_fm);
    }
    return g;
}

Scenario parseScenarioNode(const YamlNode& root) {
    Scenario sc;
    sc.channelName = yamlScalar(root, "channel");
    sc.name = yamlScalar(root, "name");
    sc.mode = parseScenarioMode(yamlScalar(root, "mode", "kp_and_ll"));
    sc.coupling = parseCoupling(yamlScalar(root, "coupling", "independent"));
    if (sc.coupling == CouplingScheme::CoupledRadial)
        throw std::runtime_error("coupling scheme not implemented in Phase 1-4: " + sc.name);

    auto chIt = root.children.find("channels");
    if (chIt == root.children.end() || !chIt->second.isSequence)
        throw std::runtime_error("scenario missing channels list: " + sc.name);

    for (const auto& node : chIt->second.sequence) {
        SpinInteractionSpec si;
        si.spin = yamlScalar(node, "spin");
        si.interacting = yamlBool(node, "interacting", true);
        auto potIt = node.children.find("potential");
        if (potIt != node.children.end()) {
            std::string ptype = yamlScalar(potIt->second, "type", "gaussian");
            si.potentialType = parsePotentialType(ptype);
            if (si.potentialType == PotentialType::Gaussian)
                si.gaussian = parseGaussian(potIt->second);
            else if (si.potentialType == PotentialType::Folded)
                si.folded = parseFolded(potIt->second);
            else if (si.potentialType == PotentialType::Ere) {
                si.ere.f0_re = yamlDouble(potIt->second, "f0_re", 0.0);
                si.ere.f0_im = yamlDouble(potIt->second, "f0_im", 0.0);
                si.ere.d0_fm = yamlDouble(potIt->second, "d0_fm", 0.0);
            }
        }
        auto ereIt = node.children.find("ere");
        if (ereIt != node.children.end()) {
            si.potentialType = PotentialType::Ere;
            si.ere.f0_re = yamlDouble(ereIt->second, "f0_re", 0.0);
            si.ere.f0_im = yamlDouble(ereIt->second, "f0_im", 0.0);
            si.ere.d0_fm = yamlDouble(ereIt->second, "d0_fm", 0.0);
        }
        sc.channels.push_back(si);
    }
    return sc;
}


double doubletV0ForB_fm(double b_fm) {
    if (std::fabs(b_fm - 0.40) < 1e-9) return 931.5;
    if (std::fabs(b_fm - 0.55) < 1e-9) return 571.5;
    if (std::fabs(b_fm - 0.70) < 1e-9) return 413.5;
    throw std::runtime_error("custom doublet folded tuning is not implemented; use doublet_b_fm: 0.40, 0.55, or 0.70");
}

GaussianSpec buildPhiNCentral(const PhiNCentralSpec& central) {
    GaussianSpec spec;
    const GaussianTerm quartet[] = {
        {-371.0, 0.15},
        { -50.0, 0.66},
        { -31.0, 1.09},
    };
    for (const auto& t : quartet) {
        if (central.quartetWeight != 0.0)
            spec.terms.push_back({central.quartetWeight * t.V_MeV, t.b_fm});
    }

    if (central.doubletWeight != 0.0) {
        if (std::fabs(central.doublet_target_f0 + 1.54) > 1e-9)
            throw std::runtime_error("custom doublet target is not implemented; use target -1.54");
        const double V0 = doubletV0ForB_fm(central.doublet_b_fm);
        spec.terms.push_back({-central.doubletWeight * V0, central.doublet_b_fm});
    }
    return spec;
}

} // namespace

Scenario loadScenario(const std::string& path) {
    return parseScenarioNode(loadYamlFile(path));
}

Scenario loadScenarioByName(const std::string& channel, const std::string& scenarioName,
                            const std::string& configRoot) {
    const std::string root = configRoot.empty() ? defaultConfigRoot() : configRoot;
    std::string path = root + "/scenarios/" + channel + "_" + scenarioName + ".yaml";
    Scenario sc = loadScenario(path);
    if (sc.channelName != channel)
        throw std::runtime_error("scenario channel mismatch: " + path);
    return sc;
}

std::vector<std::string> listScenariosForChannel(const std::string& channel,
                                                 const std::string& configRoot) {
    const std::string root = configRoot.empty() ? defaultConfigRoot() : configRoot;
    std::string dir = root + "/scenarios";
    std::string prefix = channel + "_";
    std::vector<std::string> out;
    DIR* d = opendir(dir.c_str());
    if (!d) return out;
    struct dirent* ent;
    while ((ent = readdir(d)) != nullptr) {
        std::string name = ent->d_name;
        if (name.size() <= prefix.size() + 5) continue;
        if (name.compare(0, prefix.size(), prefix) != 0) continue;
        if (name.substr(name.size() - 5) != ".yaml") continue;
        out.push_back(name.substr(prefix.size(), name.size() - prefix.size() - 5));
    }
    closedir(d);
    return out;
}

std::vector<Scenario::ResolvedSpin> Scenario::resolve(const ChannelSpec& channel,
                                                    double mu_MeV) const {
    std::vector<ResolvedSpin> out;
    double Vcrit = 0;
    bool haveVcrit = false;

  for (const auto& si : channels) {
        const SpinChannelSpec* sc = findSpinChannel(channel, si.spin);
        if (!sc) throw std::runtime_error("unknown spin channel: " + si.spin);
        ResolvedSpin rs;
        rs.spin = si.spin;
        rs.weight = sc->weight;
        rs.interacting = si.interacting && sc->defaultInteracting;

        if (si.potentialType == PotentialType::Ere) {
            rs.f0_re = si.ere.f0_re;
            rs.f0_im = si.ere.f0_im;
            rs.d0_fm = si.ere.d0_fm;
        } else if (si.potentialType == PotentialType::Gaussian) {
            double V0 = si.gaussian.V0_MeV;
            double b = si.gaussian.b_fm;
            double f0 = 0, d0 = 0;
            if (si.gaussian.useTune) {
                V0 = tuneV0(si.gaussian.targetF0, b, mu_MeV, f0, d0);
                rs.f0_re = f0;
                rs.d0_fm = d0;
            } else if (si.gaussian.VcritFraction > 0) {
                if (!haveVcrit) {
                    double f0c = 0, d0c = 0;
                    Vcrit = tuneV0(500.0, b, mu_MeV, f0c, d0c);
                    haveVcrit = true;
                }
                V0 = si.gaussian.VcritFraction * Vcrit;
            }
            rs.solver = std::make_shared<RadialSolverS>(
                gaussianPotential({{-V0, b}}), mu_MeV);
            if (!si.gaussian.useTune)
                rs.solver->scatteringParams(rs.f0_re, rs.d0_fm);
        } else if (si.potentialType == PotentialType::Folded) {
            FoldedPotentialSpec fs;
            fs.phiNCentral = buildPhiNCentral(si.folded.central);
            fs.density = si.folded.density;
            fs.fold = si.folded.fold;
            rs.solver = std::make_shared<RadialSolverS>(
                PotentialBuilder::foldedGaussian(fs), mu_MeV);
            rs.solver->scatteringParams(rs.f0_re, rs.d0_fm);
        } else if (si.potentialType == PotentialType::HalQuartetFitB) {
            rs.solver = std::make_shared<RadialSolverS>(
                gaussianPotential({{-371.0, 0.15}, {-50.0, 0.66}, {-31.0, 1.09}}),
                mu_MeV);
            rs.solver->scatteringParams(rs.f0_re, rs.d0_fm);
        } else if (si.potentialType == PotentialType::HalQuartetChizzaliTa12TPE) {
            rs.solver = std::make_shared<RadialSolverS>(
                halQuartetChizzaliTa12TPEPotential(), mu_MeV);
            rs.solver->scatteringParams(rs.f0_re, rs.d0_fm);
        }
        out.push_back(std::move(rs));
    }
    return out;
}

FoldedPotentialSummary summarizeFoldedScenario(const Scenario& scenario,
                                             const ChannelSpec& channel) {
    const SpinInteractionSpec* foldedSpin = nullptr;
    for (const auto& si : scenario.channels) {
        if (si.potentialType == PotentialType::Folded) {
            foldedSpin = &si;
            break;
        }
    }
    if (!foldedSpin)
        throw std::runtime_error("scenario has no folded potential: " + scenario.name);

    FoldedPotentialSummary out;
    out.fold = foldedSpin->folded.fold;
    if (foldedSpin->folded.central.doubletWeight != 0.0) {
        out.hasDoublet = true;
        out.doublet_b_fm = foldedSpin->folded.central.doublet_b_fm;
        out.doublet_V0_MeV = doubletV0ForB_fm(out.doublet_b_fm);
    }

    FoldedPotentialSpec fs;
    fs.phiNCentral = buildPhiNCentral(foldedSpin->folded.central);
    fs.density = foldedSpin->folded.density;
    fs.fold = foldedSpin->folded.fold;
    const GaussianSpec folded = PotentialBuilder::foldedGaussianSpec(fs);
    out.V_at_0_MeV = PotentialBuilder::gaussianSpecV0(folded);
    out.volume_integral_MeV_fm3 = PotentialBuilder::gaussianSpecVolumeIntegral(folded);

    const double mu = reducedMass(channel);
    auto resolved = scenario.resolve(channel, mu);
    for (const auto& rs : resolved) {
        if (!rs.interacting) continue;
        out.f0_fm = rs.f0_re;
        out.d0_fm = rs.d0_fm;
        if (rs.f0_re < 0.0) {
            out.BE_MeV = bindingEnergyERE(rs.f0_re, rs.d0_fm, mu);
            out.hasBE = out.BE_MeV > 0.0;
        }
        break;
    }
    return out;
}

} // namespace femto
