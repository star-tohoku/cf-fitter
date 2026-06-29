#include "femto/YamlReader.hpp"
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace femto {
namespace {

std::string trim(const std::string& s) {
    std::size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    std::size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

std::string stripComment(const std::string& line) {
    bool inQuote = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        if (line[i] == '"') inQuote = !inQuote;
        if (!inQuote && line[i] == '#') return line.substr(0, i);
    }
    return line;
}

int countIndent(const std::string& line) {
    int n = 0;
    for (char c : line) {
        if (c == ' ') ++n;
        else if (c == '\t') n += 2;
        else break;
    }
    return n;
}

bool parseScalarPair(const std::string& line, std::string& key, std::string& val) {
    auto pos = line.find(':');
    if (pos == std::string::npos) return false;
    key = trim(line.substr(0, pos));
    val = trim(line.substr(pos + 1));
    if (!val.empty() && val.front() == '"' && val.back() == '"')
        val = val.substr(1, val.size() - 2);
    return !key.empty();
}

void parseBlock(const std::vector<std::string>& lines, std::size_t& i, int baseIndent,
                YamlNode& out) {
    while (i < lines.size()) {
        std::string raw = stripComment(lines[i]);
        if (trim(raw).empty()) { ++i; continue; }
        int ind = countIndent(raw);
        if (ind < baseIndent) break;
        std::string content = trim(raw.substr(ind));
        if (content.empty()) { ++i; continue; }

        if (content[0] == '-') {
            YamlNode item;
            std::string rest = trim(content.substr(1));
            if (!rest.empty()) {
                std::string k, v;
                if (parseScalarPair(rest, k, v)) item.children[k].scalar = v;
                else item.scalar = rest;
            }
            ++i;
            if (i < lines.size() && countIndent(lines[i]) > ind) {
                parseBlock(lines, i, ind + 2, item);
            }
            out.isSequence = true;
            out.sequence.push_back(std::move(item));
            continue;
        }

        std::string key, val;
        if (!parseScalarPair(content, key, val)) { ++i; continue; }
        YamlNode& child = out.children[key];
        if (!val.empty()) {
            child.scalar = val;
            ++i;
            continue;
        }
        ++i;
        if (i < lines.size() && countIndent(lines[i]) > ind)
            parseBlock(lines, i, ind + 2, child);
    }
}

} // namespace

YamlNode loadYamlFile(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open yaml: " + path);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) lines.push_back(line);
    YamlNode root;
    std::size_t i = 0;
    parseBlock(lines, i, 0, root);
    return root;
}

std::string yamlScalar(const YamlNode& n, const std::string& key, const std::string& def) {
    auto it = n.children.find(key);
    if (it == n.children.end()) return def;
    return it->second.scalar;
}

double yamlDouble(const YamlNode& n, const std::string& key, double def) {
    auto s = yamlScalar(n, key, "");
    if (s.empty()) return def;
    return std::stod(s);
}

bool yamlBool(const YamlNode& n, const std::string& key, bool def) {
    auto s = yamlScalar(n, key, "");
    if (s.empty()) return def;
    if (s == "true" || s == "yes" || s == "1") return true;
    if (s == "false" || s == "no" || s == "0") return false;
    return def;
}

} // namespace femto
