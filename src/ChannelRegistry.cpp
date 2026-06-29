#include "femto/ChannelRegistry.hpp"
#include "femto/YamlReader.hpp"
#include <cstdlib>
#include <dirent.h>
#include <stdexcept>
#include <sys/stat.h>

namespace femto {

std::string defaultConfigRoot() {
#ifdef CF_FITTER_CONFIG_DIR
    return CF_FITTER_CONFIG_DIR;
#else
    return "config";
#endif
}

namespace {

ChannelSpec parseChannelFile(const std::string& path) {
    YamlNode root = loadYamlFile(path);
    ChannelSpec ch;
    ch.name = yamlScalar(root, "name");
    ch.massA_MeV = yamlDouble(root, "massA_MeV");
    ch.massB_MeV = yamlDouble(root, "massB_MeV");
    ch.identical = yamlBool(root, "identical", false);
    ch.qsAmplitude = yamlDouble(root, "qsAmplitude", 0.0);

    auto it = root.children.find("spinChannels");
    if (it != root.children.end() && it->second.isSequence) {
        for (const auto& node : it->second.sequence) {
            SpinChannelSpec sc;
            sc.name = yamlScalar(node, "name");
            sc.weight = yamlDouble(node, "weight", 1.0);
            sc.defaultInteracting = yamlBool(node, "defaultInteracting", true);
            ch.spinChannels.push_back(sc);
        }
    }
    if (ch.name.empty()) throw std::runtime_error("channel yaml missing name: " + path);
    return ch;
}

} // namespace

ChannelRegistry::ChannelRegistry(const std::string& configRoot) {
    const std::string root = configRoot.empty() ? defaultConfigRoot() : configRoot;
    loadFromDirectory(root + "/channels");
}

void ChannelRegistry::loadFromDirectory(const std::string& dir) {
    DIR* d = opendir(dir.c_str());
    if (!d) throw std::runtime_error("cannot open channel config dir: " + dir);
    struct dirent* ent;
    while ((ent = readdir(d)) != nullptr) {
        std::string name = ent->d_name;
        if (name.size() < 5 || name.substr(name.size() - 5) != ".yaml") continue;
        channels_.push_back(parseChannelFile(dir + "/" + name));
    }
    closedir(d);
    if (channels_.empty()) throw std::runtime_error("no channel yaml in " + dir);
}

const ChannelSpec& ChannelRegistry::get(const std::string& name) const {
    for (const auto& ch : channels_)
        if (ch.name == name) return ch;
    throw std::runtime_error("unknown channel: " + name);
}

std::vector<std::string> ChannelRegistry::listChannels() const {
    std::vector<std::string> out;
    for (const auto& ch : channels_) out.push_back(ch.name);
    return out;
}

} // namespace femto
