#ifndef FEMTO_SCENARIO_HPP
#define FEMTO_SCENARIO_HPP

#include "femto/ChannelSpec.hpp"
#include "femto/PotentialBuilder.hpp"
#include <complex>
#include <memory>
#include <string>
#include <vector>

namespace femto {

class RadialSolverS;

enum class ScenarioMode { KpAndLl, LlOnly };

enum class CouplingScheme {
    Independent,
    EffectiveSingleChannel,
    CoupledRadial
};

enum class PotentialType { Gaussian, Ere, Folded, HalQuartetFitB, HalQuartetChizzaliTa12TPE };

struct GaussianPotentialSpec {
    double V0_MeV = 0.0;
    double b_fm = 1.0;
    bool useTune = false;
    double targetF0 = 0.0;
    double VcritFraction = 0.0;
};

struct EreSpec {
    double f0_re = 0.0;
    double f0_im = 0.0;
    double d0_fm = 0.0;
};

struct PhiNCentralSpec {
    double quartetWeight = 1.0;
    double doubletWeight = 0.0;
    double doublet_b_fm = 0.55;
    double doublet_target_f0 = -1.54;
};

struct FoldedScenarioSpec {
    PhiNCentralSpec central;
    GaussianDensitySpec density;
    bool fold = true;
};

struct SpinInteractionSpec {
    std::string spin;
    bool interacting = true;
    PotentialType potentialType = PotentialType::Gaussian;
    GaussianPotentialSpec gaussian;
    EreSpec ere;
    FoldedScenarioSpec folded;
};

struct Scenario {
    std::string channelName;
    std::string name;
    ScenarioMode mode = ScenarioMode::KpAndLl;
    CouplingScheme coupling = CouplingScheme::Independent;
    std::vector<SpinInteractionSpec> channels;

    struct ResolvedSpin {
        std::string spin;
        double weight = 0.0;
        bool interacting = true;
        double f0_re = 0.0;
        double f0_im = 0.0;
        double d0_fm = 0.0;
        std::shared_ptr<RadialSolverS> solver;
    };

    std::vector<ResolvedSpin> resolve(const ChannelSpec& channel, double mu_MeV) const;
};

Scenario loadScenario(const std::string& path);
Scenario loadScenarioByName(const std::string& channel, const std::string& scenarioName,
                            const std::string& configRoot = "");

std::vector<std::string> listScenariosForChannel(const std::string& channel,
                                                 const std::string& configRoot = "");

struct FoldedPotentialSummary {
    bool fold = false;
    bool hasDoublet = false;
    double doublet_b_fm = 0.0;
    double doublet_V0_MeV = 0.0;
    double f0_fm = 0.0;
    double d0_fm = 0.0;
    double BE_MeV = 0.0;
    bool hasBE = false;
    double V_at_0_MeV = 0.0;
    double volume_integral_MeV_fm3 = 0.0;
    std::string notes;
};

FoldedPotentialSummary summarizeFoldedScenario(const Scenario& scenario,
                                               const ChannelSpec& channel);

} // namespace femto

#endif
