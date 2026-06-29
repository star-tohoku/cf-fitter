#ifndef FEMTO_YAML_READER_HPP
#define FEMTO_YAML_READER_HPP

#include <map>
#include <string>
#include <vector>

namespace femto {

struct YamlNode {
    std::string scalar;
    std::map<std::string, YamlNode> children;
    std::vector<YamlNode> sequence;
    bool isSequence = false;
};

YamlNode loadYamlFile(const std::string& path);
std::string yamlScalar(const YamlNode& n, const std::string& key, const std::string& def = "");
double yamlDouble(const YamlNode& n, const std::string& key, double def = 0.0);
bool yamlBool(const YamlNode& n, const std::string& key, bool def = false);

} // namespace femto

#endif
