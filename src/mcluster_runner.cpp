#include "mcluster_runner.h"

#include <unistd.h>      // fork, exec, pipe, dup2, read, close
#include <sys/wait.h>    // waitpid, WIFEXITED, WEXITSTATUS
#include <sys/stat.h>    // stat
#include <fstream>       // ifstream, ofstream
#include <sstream>       // ostringstream
#include <iostream>      // cout, cerr
#include <iomanip>       // setprecision, scientific
#include <cstring>       // strerror

std::vector<std::string> buildMclusterArgs(const MclusterConfig& config) {
    std::vector<std::string> args;

    // M or N (M takes precedence per Phase 26 decision)
    if (config.M > 0.0) {
        args.push_back("-M");
        std::ostringstream oss;
        oss << config.M;
        args.push_back(oss.str());
    } else {
        args.push_back("-N");
        std::ostringstream oss;
        oss << config.N;
        args.push_back(oss.str());
    }

    // Density profile
    {
        args.push_back("-P");
        std::ostringstream oss;
        oss << config.P;
        args.push_back(oss.str());
    }

    // Half-mass radius
    {
        args.push_back("-R");
        std::ostringstream oss;
        oss << config.R;
        args.push_back(oss.str());
    }

    // IMF
    {
        args.push_back("-f");
        std::ostringstream oss;
        oss << config.f;
        args.push_back(oss.str());
    }

    // Metallicity
    {
        args.push_back("-Z");
        std::ostringstream oss;
        oss << config.Z;
        args.push_back(oss.str());
    }

    // Binary fraction
    {
        args.push_back("-b");
        std::ostringstream oss;
        oss << config.b;
        args.push_back(oss.str());
    }

    // Epoch
    {
        args.push_back("-e");
        std::ostringstream oss;
        oss << config.e;
        args.push_back(oss.str());
    }

    // Output format: ASCII table with astrophysical units
    args.push_back("-C");
    args.push_back("3");   // Table of stars format

    args.push_back("-u");
    args.push_back("1");   // Astrophysical units (Msun, pc, km/s)

    // Output filename (base name, McLuster adds .txt)
    args.push_back("-o");
    args.push_back(MCLUSTER_OUTPUT_BASE);

    return args;
}

RunResult runMclusterSubprocess(const std::string& binary_path,
                                 const std::vector<std::string>& args) {
    RunResult result;
    result.exit_code = -1;
    result.success = false;

    int stdout_pipe[2];
    int stderr_pipe[2];

    // Create pipes for output capture
    if (pipe(stdout_pipe) < 0) {
        result.stderr_content = "Failed to create stdout pipe: ";
        result.stderr_content += strerror(errno);
        return result;
    }

    if (pipe(stderr_pipe) < 0) {
        close(stdout_pipe[0]);
        close(stdout_pipe[1]);
        result.stderr_content = "Failed to create stderr pipe: ";
        result.stderr_content += strerror(errno);
        return result;
    }

    pid_t pid = fork();
    if (pid < 0) {
        // Fork failed
        close(stdout_pipe[0]);
        close(stdout_pipe[1]);
        close(stderr_pipe[0]);
        close(stderr_pipe[1]);
        result.stderr_content = "Fork failed: ";
        result.stderr_content += strerror(errno);
        return result;
    }

    if (pid == 0) {
        // Child process
        // Close read ends of pipes
        close(stdout_pipe[0]);
        close(stderr_pipe[0]);

        // Redirect stdout and stderr to pipes
        dup2(stdout_pipe[1], STDOUT_FILENO);
        dup2(stderr_pipe[1], STDERR_FILENO);

        // Close write ends (now duplicated)
        close(stdout_pipe[1]);
        close(stderr_pipe[1]);

        // Build argv array for execv
        std::vector<char*> argv;
        argv.push_back(const_cast<char*>(binary_path.c_str()));
        for (size_t i = 0; i < args.size(); ++i) {
            argv.push_back(const_cast<char*>(args[i].c_str()));
        }
        argv.push_back(nullptr);

        // Execute McLuster
        execv(binary_path.c_str(), argv.data());

        // If execv returns, it failed
        _exit(127);
    }

    // Parent process
    // Close write ends of pipes
    close(stdout_pipe[1]);
    close(stderr_pipe[1]);

    // Read from pipes
    char buffer[4096];
    ssize_t bytes_read;

    // Read stdout
    while ((bytes_read = read(stdout_pipe[0], buffer, sizeof(buffer) - 1)) > 0) {
        buffer[bytes_read] = '\0';
        result.stdout_content += buffer;
    }
    close(stdout_pipe[0]);

    // Read stderr
    while ((bytes_read = read(stderr_pipe[0], buffer, sizeof(buffer) - 1)) > 0) {
        buffer[bytes_read] = '\0';
        result.stderr_content += buffer;
    }
    close(stderr_pipe[0]);

    // Wait for child process to complete
    int status;
    waitpid(pid, &status, 0);

    if (WIFEXITED(status)) {
        result.exit_code = WEXITSTATUS(status);
        result.success = (result.exit_code == 0);
    } else {
        result.exit_code = -1;
        result.success = false;
        if (result.stderr_content.empty()) {
            result.stderr_content = "Child process did not exit normally";
        }
    }

    return result;
}

bool validateMclusterOutput(const std::string& output_file, int expected_count) {
    // Check file exists using stat
    struct stat st;
    if (stat(output_file.c_str(), &st) != 0) {
        std::cerr << "McLuster output file not found: " << output_file << std::endl;
        return false;
    }

    // Check file size > 0
    if (st.st_size == 0) {
        std::cerr << "McLuster output file is empty: " << output_file << std::endl;
        return false;
    }

    // Open and check content
    std::ifstream in(output_file.c_str());
    if (!in) {
        std::cerr << "Failed to open McLuster output file: " << output_file << std::endl;
        return false;
    }

    // Check first line contains expected header
    std::string first_line;
    if (!std::getline(in, first_line)) {
        std::cerr << "Failed to read first line from McLuster output" << std::endl;
        return false;
    }

    if (first_line.find("#Mass") == std::string::npos) {
        std::cerr << "McLuster output missing expected header (expected #Mass...)" << std::endl;
        return false;
    }

    // Count data lines (non-empty, non-comment)
    int line_count = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] != '#') {
            ++line_count;
        }
    }

    if (line_count == 0) {
        std::cerr << "McLuster output has no data lines" << std::endl;
        return false;
    }

    // If expected_count > 0, verify line count matches
    if (expected_count > 0 && line_count != expected_count) {
        std::cerr << "McLuster output has " << line_count << " stars, expected "
                  << expected_count << std::endl;
        return false;
    }

    std::cout << "McLuster generated " << line_count << " stars" << std::endl;
    return true;
}

bool transformMclusterOutput(const std::string& mcluster_file,
                              const std::string& abyss_file) {
    std::ifstream in(mcluster_file.c_str());
    if (!in) {
        std::cerr << "Failed to open McLuster output file: " << mcluster_file << std::endl;
        return false;
    }

    std::ofstream out(abyss_file.c_str());
    if (!out) {
        std::cerr << "Failed to create ABYSS IC file: " << abyss_file << std::endl;
        return false;
    }

    std::string line;
    int lines_written = 0;

    while (std::getline(in, line)) {
        // Skip empty lines and header/comment lines (starting with #)
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Parse McLuster format: mass x y z vx vy vz [extras...]
        std::istringstream iss(line);
        double mass, x, y, z, vx, vy, vz;

        if (!(iss >> mass >> x >> y >> z >> vx >> vy >> vz)) {
            std::cerr << "Warning: Failed to parse line: " << line << std::endl;
            continue;  // Skip malformed lines
        }

        // Write in ABYSS format: x y z vx vy vz mass
        out << std::scientific << std::setprecision(8)
            << x << " " << y << " " << z << " "
            << vx << " " << vy << " " << vz << " "
            << mass << "\n";

        ++lines_written;
    }

    in.close();
    out.close();

    if (lines_written == 0) {
        std::cerr << "No data written to ABYSS IC file" << std::endl;
        return false;
    }

    if (!out.good() && !out.eof()) {
        std::cerr << "Error writing ABYSS IC file" << std::endl;
        return false;
    }

    std::cout << "Transformed " << lines_written << " particles to ABYSS format: "
              << abyss_file << std::endl;
    return true;
}
