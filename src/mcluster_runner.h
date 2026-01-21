#ifndef MCLUSTER_RUNNER_H
#define MCLUSTER_RUNNER_H

#include <string>
#include <vector>
#include "mcluster_config.h"

// Default McLuster binary path (symlink created by build system)
const std::string MCLUSTER_BINARY = "./src/mcluster";

// Base filename for McLuster output (McLuster adds .txt extension)
const std::string MCLUSTER_OUTPUT_BASE = "mcluster_ic";

// Result from subprocess execution
struct RunResult {
    int exit_code;
    std::string stdout_content;
    std::string stderr_content;
    bool success;
};

// Build McLuster command-line arguments from configuration
// Converts MclusterConfig struct to vector of CLI arguments
// M takes precedence over N when both are specified
std::vector<std::string> buildMclusterArgs(const MclusterConfig& config);

// Execute McLuster subprocess with output capture
// Uses fork/exec pattern for proper exit code handling and stderr capture
// Returns RunResult with exit_code, captured stdout/stderr, and success flag
RunResult runMclusterSubprocess(const std::string& binary_path,
                                 const std::vector<std::string>& args);

// Validate McLuster output file
// Checks: file exists, non-empty, has expected header, contains data lines
// If expected_count > 0, verifies line count matches
// If expected_count = 0, just checks file exists and has data
// Returns true on success, prints diagnostic messages
bool validateMclusterOutput(const std::string& output_file, int expected_count);

// Transform McLuster output to ABYSS format
// McLuster format: mass x y z vx vy vz [extras...]
// ABYSS format: x y z vx vy vz mass
// Skips header/comment lines, writes with scientific notation (8 decimal places)
// Returns true on success
bool transformMclusterOutput(const std::string& mcluster_file,
                              const std::string& abyss_file);

#endif // MCLUSTER_RUNNER_H
