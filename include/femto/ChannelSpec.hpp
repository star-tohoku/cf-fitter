#ifndef FEMTO_CHANNEL_SPEC_HPP
#define FEMTO_CHANNEL_SPEC_HPP

#include <string>
#include <vector>

namespace femto {

struct SpinChannelSpec {
    std::string name;
    double weight = 1.0;
    bool defaultInteracting = true;
};

struct ChannelSpec {
    std::string name;
    double massA_MeV = 0.0;
    double massB_MeV = 0.0;
    bool identical = false;
    double qsAmplitude = 0.0;
    std::vector<SpinChannelSpec> spinChannels;
};

double reducedMass(const ChannelSpec& ch);

const SpinChannelSpec* findSpinChannel(const ChannelSpec& ch, const std::string& name);

} // namespace femto

#endif
