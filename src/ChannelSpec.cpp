#include "femto/ChannelSpec.hpp"
#include <stdexcept>

namespace femto {

double reducedMass(const ChannelSpec& ch) {
    return ch.massA_MeV * ch.massB_MeV / (ch.massA_MeV + ch.massB_MeV);
}

const SpinChannelSpec* findSpinChannel(const ChannelSpec& ch, const std::string& name) {
    for (const auto& sc : ch.spinChannels)
        if (sc.name == name) return &sc;
    return nullptr;
}

} // namespace femto
