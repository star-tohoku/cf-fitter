#include "calc_common.hpp"
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

namespace {

void usage() {
    std::printf("cf-calc — femtoscopic correlation function calculator\n\n");
    std::printf("Usage:\n");
    std::printf("  cf-calc validate [--config DIR]\n");
    std::printf("  cf-calc make-cf --channel NAME [--scenarios A,B] [--output FILE] [--config DIR]\n");
    std::printf("\nExamples:\n");
    std::printf("  cf-calc validate\n");
    std::printf("  cf-calc make-cf --channel phi_proton --scenarios HAL,HALplusBound,ALICEeff\n");
    std::printf("  cf-calc make-cf --channel phi_alpha --scenarios weak,strong,bound\n");
}

} // namespace

int main(int argc, char** argv) {
    std::vector<std::string> args;
    for (int i = 1; i < argc; ++i) args.emplace_back(argv[i]);
    if (args.empty()) {
        usage();
        return 1;
    }

    const std::string cmd = args[0];
    args.erase(args.begin());

    try {
        if (cmd == "validate") return runValidate(args);
        if (cmd == "make-cf") return runMakeCf(args);
        if (cmd == "help" || cmd == "-h" || cmd == "--help") {
            usage();
            return 0;
        }
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << std::endl;
        return 1;
    }

    std::cerr << "unknown command: " << cmd << std::endl;
    usage();
    return 1;
}
