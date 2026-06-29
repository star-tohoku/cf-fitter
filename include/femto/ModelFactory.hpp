#ifndef FEMTO_MODEL_FACTORY_HPP
#define FEMTO_MODEL_FACTORY_HPP

#include "femto/FemtoModels.hpp"
#include "femto/Scenario.hpp"
#include "femto/ChannelSpec.hpp"
#include <memory>

namespace femto {

enum class ModelKind { LL, KP };

struct BuiltModels {
    LLModel ll;
    std::shared_ptr<KPModel> kp;
    bool llOnly = false;
};

BuiltModels buildModels(const Scenario& scenario, const ChannelSpec& channel);

} // namespace femto

#endif
