#ifndef MCLUSTER_CONFIG_H
#define MCLUSTER_CONFIG_H

#include <string>
#include <vector>

// Configuration for McLuster IC generator
// Parsed from [mcluster] section in TOML config
struct MclusterConfig {
    // Required: user must specify N or M (not both, M takes precedence)
    int N = 0;              // Number of stars (0 = use M instead)
    double M = 0.0;         // Total mass in Msun (0 = use N instead)

    // Optional with McLuster defaults
    int P = 0;              // Density profile: -1=none, 0=Plummer, 1=King, 2=Subr, 3=EFF
    double R = 0.8;         // Half-mass radius [pc]
    int f = 1;              // IMF: 0=single mass, 1=Kroupa, 2=user-defined
    double Z = 0.02;        // Metallicity [0.0001-0.03, solar=0.02]
    double b = 0.0;         // Binary fraction [0.0-1.0]
    double e = 0.0;         // Stellar evolution epoch [Myr]

    // ABYSS-specific
    bool generate_only = false;  // If true, exit after IC generation

    // Track if [mcluster] section was present in config
    bool has_mcluster_section = false;
};

// Valid parameter names for unknown key detection
const std::vector<std::string> MCLUSTER_VALID_PARAMS = {
    "N", "M", "P", "R", "f", "Z", "b", "e", "generate_only"
};

#endif // MCLUSTER_CONFIG_H
