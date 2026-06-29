#ifndef FEMTO_CHANNEL_REGISTRY_HPP
#define FEMTO_CHANNEL_REGISTRY_HPP

#include "femto/ChannelSpec.hpp"
#include <string>
#include <vector>

namespace femto {

class ChannelRegistry {
public:
    explicit ChannelRegistry(const std::string& configRoot = "");

    const ChannelSpec& get(const std::string& name) const;
    std::vector<std::string> listChannels() const;

private:
    std::vector<ChannelSpec> channels_;
    void loadFromDirectory(const std::string& dir);
};

std::string defaultConfigRoot();

} // namespace femto

#endif
