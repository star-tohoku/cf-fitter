#include "femto/ModelFactory.hpp"
#include <stdexcept>

namespace femto {

BuiltModels buildModels(const Scenario& scenario, const ChannelSpec& channel) {
    BuiltModels out;
    const double mu = reducedMass(channel);
    auto resolved = scenario.resolve(channel, mu);

    out.ll.identical = channel.identical;
    out.ll.qsAmplitude = channel.qsAmplitude;
    out.llOnly = (scenario.mode == ScenarioMode::LlOnly);

    if (scenario.coupling == CouplingScheme::EffectiveSingleChannel) {
        double wsum = 0;
        for (const auto& rs : resolved) {
            if (!rs.interacting) continue;
            wsum += rs.weight;
        }
        if (wsum <= 0) throw std::runtime_error("no interacting channel in scenario");
        double f0_re = 0, f0_im = 0, d0 = 0;
        for (const auto& rs : resolved) {
            if (!rs.interacting) continue;
            double w = rs.weight / wsum;
            f0_re += w * rs.f0_re;
            f0_im += w * rs.f0_im;
            d0 += w * rs.d0_fm;
        }
        out.ll.channels.push_back({1.0, {f0_re, f0_im}, d0, true});
        if (!out.llOnly) {
            out.kp = std::make_shared<KPModel>();
            out.kp->identical = channel.identical;
            out.kp->qsAmplitude = channel.qsAmplitude;
            for (const auto& rs : resolved) {
                if (!rs.interacting || !rs.solver) continue;
                out.kp->channels.push_back({rs.weight, rs.solver});
            }
        }
        return out;
    }

    for (const auto& rs : resolved) {
        out.ll.channels.push_back(
            {rs.weight, {rs.f0_re, rs.f0_im}, rs.d0_fm, rs.interacting});
    }

    if (!out.llOnly) {
        out.kp = std::make_shared<KPModel>();
        out.kp->identical = channel.identical;
        out.kp->qsAmplitude = channel.qsAmplitude;
        for (const auto& rs : resolved) {
            if (!rs.interacting || !rs.solver) continue;
            out.kp->channels.push_back({rs.weight, rs.solver});
        }
        if (out.kp->channels.empty()) {
            out.kp.reset();
            out.llOnly = true;
        }
    }
    return out;
}

} // namespace femto
