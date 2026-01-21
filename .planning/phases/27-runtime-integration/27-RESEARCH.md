# Phase 27: Runtime Integration - Research

**Researched:** 2026-01-20
**Domain:** Subprocess management, McLuster orchestration, POSIX process control
**Confidence:** HIGH

## Summary

This phase integrates McLuster subprocess execution into ABYSS's main.cpp, enabling automatic initial condition (IC) generation from TOML configuration. The core challenge is straightforward: detect `[mcluster]` section presence, spawn McLuster binary with correct arguments, capture output, validate results, and either exit (generate_only mode) or proceed to simulation.

The existing ABYSS codebase provides all necessary infrastructure: mcluster_config struct with `has_mcluster_section` flag from Phase 26, working McLuster binary at `src/mcluster` symlink, and clear data flow in main.cpp between config parsing and data reading. The subprocess spawning requires standard POSIX functions (fork/exec or popen) which are already available through `<unistd.h>` included in multiple source files.

**Primary recommendation:** Use `fork()/exec()` with pipe-based output capture for McLuster subprocess. This approach provides proper exit code handling, stderr capture for error diagnostics, and seamless integration with MPI (must spawn McLuster only on ROOT process before MPI data distribution).

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `<unistd.h>` | POSIX | fork(), exec(), pipe(), dup2() | Already included in ABYSS, universal availability |
| `<sys/wait.h>` | POSIX | waitpid(), WEXITSTATUS macros | Process completion and exit code |
| `<sys/stat.h>` | POSIX | stat(), file existence checks | Already used in read_write.cpp |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `<fstream>` | Standard | File I/O for output validation | Already used throughout codebase |
| `<sstream>` | Standard | Command-line argument building | Already used in read_parameter_file.cpp |
| `<cstdlib>` | Standard | exit(), EXIT_FAILURE | Error handling |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| fork/exec | popen() | Simpler but no stderr capture, no direct exit code access |
| fork/exec | posix_spawn() | More portable to some systems but more complex API |
| fork/exec | system() | Simpler but shell injection risk, no output capture |

**Installation:**
```bash
# No additional installation needed - all POSIX headers already available
```

## Architecture Patterns

### Recommended Project Structure
```
src/
├── main.cpp                 # Add McLuster orchestration before readData()
├── mcluster_runner.cpp      # NEW: McLuster subprocess logic
├── mcluster_runner.h        # NEW: Function declarations
├── mcluster_config.h        # EXISTS: Config struct (from Phase 26)
├── read_parameter_file.cpp  # EXISTS: Config parsing (from Phase 26)
└── read_write.cpp           # Existing IC file reading (unchanged)
```

### Pattern 1: Fork-Exec with Pipe Output Capture
**What:** Spawn subprocess, capture stdout/stderr, wait for completion
**When to use:** Need both output capture and proper exit code handling
**Example:**
```cpp
// Source: Standard POSIX pattern
#include <unistd.h>
#include <sys/wait.h>
#include <cstdio>
#include <string>
#include <vector>

struct RunResult {
    int exit_code;
    std::string stdout_content;
    std::string stderr_content;
    bool success;
};

RunResult runMclusterSubprocess(const std::string& binary_path,
                                 const std::vector<std::string>& args) {
    RunResult result;
    int stdout_pipe[2], stderr_pipe[2];

    if (pipe(stdout_pipe) < 0 || pipe(stderr_pipe) < 0) {
        result.success = false;
        result.stderr_content = "Failed to create pipes";
        return result;
    }

    pid_t pid = fork();
    if (pid < 0) {
        result.success = false;
        result.stderr_content = "Fork failed";
        return result;
    }

    if (pid == 0) {
        // Child process
        close(stdout_pipe[0]);  // Close read end
        close(stderr_pipe[0]);
        dup2(stdout_pipe[1], STDOUT_FILENO);
        dup2(stderr_pipe[1], STDERR_FILENO);
        close(stdout_pipe[1]);
        close(stderr_pipe[1]);

        // Build argv
        std::vector<char*> argv;
        argv.push_back(const_cast<char*>(binary_path.c_str()));
        for (const auto& arg : args) {
            argv.push_back(const_cast<char*>(arg.c_str()));
        }
        argv.push_back(nullptr);

        execv(binary_path.c_str(), argv.data());
        _exit(127);  // execv failed
    }

    // Parent process
    close(stdout_pipe[1]);  // Close write ends
    close(stderr_pipe[1]);

    // Read output
    char buffer[4096];
    ssize_t bytes;
    while ((bytes = read(stdout_pipe[0], buffer, sizeof(buffer)-1)) > 0) {
        buffer[bytes] = '\0';
        result.stdout_content += buffer;
    }
    while ((bytes = read(stderr_pipe[0], buffer, sizeof(buffer)-1)) > 0) {
        buffer[bytes] = '\0';
        result.stderr_content += buffer;
    }
    close(stdout_pipe[0]);
    close(stderr_pipe[0]);

    int status;
    waitpid(pid, &status, 0);
    result.exit_code = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
    result.success = (result.exit_code == 0);

    return result;
}
```

### Pattern 2: MPI-Aware Subprocess Execution
**What:** Execute subprocess only on ROOT rank, broadcast results
**When to use:** MPI programs needing subprocess output available to all ranks
**Example:**
```cpp
// Source: ABYSS pattern from main.cpp line 62
void runMclusterIfConfigured() {
    if (!mcluster_config.has_mcluster_section) {
        return;  // No mcluster section - use existing IC file
    }

    // Only ROOT process runs McLuster
    if (my_rank == ROOT) {
        std::cout << "Generating initial conditions with McLuster..." << std::endl;

        RunResult result = runMclusterSubprocess(
            "./src/mcluster",
            buildMclusterArgs(mcluster_config)
        );

        if (!result.success) {
            std::cerr << "McLuster failed with exit code: " << result.exit_code << std::endl;
            std::cerr << result.stderr_content << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Validate output file exists
        if (!validateMclusterOutput(mcluster_config)) {
            std::cerr << "McLuster output validation failed" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (mcluster_config.generate_only) {
            std::cout << "generate_only=true: Exiting after IC generation" << std::endl;
            MPI_Finalize();
            exit(0);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);  // Wait for ROOT to complete
}
```

### Pattern 3: McLuster Argument Building
**What:** Convert MclusterConfig to command-line arguments
**When to use:** Translating parsed config to subprocess arguments
**Example:**
```cpp
// Source: McLuster README parameters
std::vector<std::string> buildMclusterArgs(const MclusterConfig& config) {
    std::vector<std::string> args;

    // N or M (M takes precedence per CONTEXT.md from Phase 26)
    if (config.M > 0.0) {
        args.push_back("-M");
        args.push_back(std::to_string(config.M));
    } else {
        args.push_back("-N");
        args.push_back(std::to_string(config.N));
    }

    // Density profile
    args.push_back("-P");
    args.push_back(std::to_string(config.P));

    // Half-mass radius
    args.push_back("-R");
    args.push_back(std::to_string(config.R));

    // IMF
    args.push_back("-f");
    args.push_back(std::to_string(config.f));

    // Metallicity
    args.push_back("-Z");
    args.push_back(std::to_string(config.Z));

    // Binary fraction
    args.push_back("-b");
    args.push_back(std::to_string(config.b));

    // Epoch
    args.push_back("-e");
    args.push_back(std::to_string(config.e));

    // Output format: ASCII table with astrophysical units
    args.push_back("-C");
    args.push_back("3");   // Table of stars format

    args.push_back("-u");
    args.push_back("1");   // Astrophysical units (Msun, pc, km/s)

    // Output filename (base name, McLuster adds .txt)
    args.push_back("-o");
    args.push_back("mcluster_ic");

    return args;
}
```

### Anti-Patterns to Avoid
- **System() calls:** Don't use `system()` - vulnerable to shell injection, poor error handling
- **Ignoring exit codes:** Always check subprocess exit code and handle failures
- **Blocking all MPI ranks:** Only ROOT should run subprocess; others barrier-wait
- **Leaving zombie processes:** Always waitpid() after fork to reap child
- **Hardcoded paths:** Use relative path from executable or config-derived paths

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Command execution | Custom shell wrapper | fork/exec | Shell injection, quoting issues |
| Output capture | Temp files | pipe() + read() | Atomicity, cleanup complexity |
| Process completion | sleep/poll loops | waitpid() | Race conditions, resource waste |
| File format conversion | Custom parser | Column reordering on read | Single responsibility |

**Key insight:** McLuster output format differs from ABYSS expected format (column order). Rather than post-processing McLuster output, modify `readData()` to handle McLuster format OR use a simple file transformation.

## McLuster Output Format Analysis

### McLuster Output (with -C 3 -u 1)
```
#Mass_[Msun] x_[pc] y_[pc] z_[pc] vx_[km/s] vy_[km/s] vz_[km/s] ...
0.593310     0.29   -0.41  0.62   0.05      0.11      -0.07    ...
```
Columns: `mass x y z vx vy vz [stellar_evolution_data...]`

### ABYSS Expected Format (readData)
```
x y z vx vy vz mass
```
Columns: `x y z vx vy vz mass` (7 columns, position/velocity first)

### Resolution Strategy
**Option A:** Transform McLuster output file to ABYSS format
- Simple column reordering: `awk '{print $2,$3,$4,$5,$6,$7,$1}' mcluster_ic.txt > ic.dat`
- Or do it in C++ after McLuster completes

**Option B:** Add McLuster format support to readData()
- Detect header line starting with `#Mass`
- Parse with different column order
- More flexible but modifies existing stable code

**Recommendation:** Option A - Transform output file. Keeps readData() unchanged, single transformation point.

## Common Pitfalls

### Pitfall 1: MPI Process Interference
**What goes wrong:** All MPI ranks try to spawn McLuster, causing race conditions and multiple IC files
**Why it happens:** MPI programs run same code on all ranks
**How to avoid:** Guard subprocess call with `if (my_rank == ROOT)` and MPI_Barrier after
**Warning signs:** Multiple mcluster_ic.txt files, corrupted output, MPI hangs

### Pitfall 2: Fork in MPI Program
**What goes wrong:** fork() after MPI_Init can cause undefined behavior on some MPI implementations
**Why it happens:** MPI internal state may not be fork-safe
**How to avoid:** Run McLuster before heavy MPI initialization (after MPI_Init is OK, but before data distribution). Alternatively, use `MPI_Comm_spawn` but that's overkill here.
**Warning signs:** Random hangs, corrupted MPI state, SIGBUS errors

### Pitfall 3: Missing Output Validation
**What goes wrong:** McLuster "succeeds" but produces empty or malformed output
**Why it happens:** Wrong parameters, filesystem issues, partial writes
**How to avoid:** After subprocess completes: check file exists, check file size > 0, check first line matches expected format
**Warning signs:** readData() fails cryptically, simulation starts with 0 particles

### Pitfall 4: Exit Code Masking
**What goes wrong:** McLuster fails but ABYSS doesn't detect it
**Why it happens:** `system()` returns shell exit code, not program exit code; or WEXITSTATUS not used
**How to avoid:** Use `WIFEXITED(status) && WEXITSTATUS(status) == 0` to check success
**Warning signs:** Silent failures, missing IC file not reported

### Pitfall 5: Filename Collision
**What goes wrong:** McLuster overwrites user's existing IC file
**Why it happens:** Output filename same as Filename in config
**How to avoid:** Use dedicated name like `mcluster_ic.txt`, or generate in temp directory
**Warning signs:** User complaints about lost data

### Pitfall 6: generate_only vs Normal Flow
**What goes wrong:** Simulation proceeds even when generate_only=true
**Why it happens:** Early exit path not properly implemented
**How to avoid:** Check `generate_only` flag after McLuster completes, call `MPI_Finalize()` and `exit(0)` on ROOT, have workers detect and exit cleanly
**Warning signs:** Unnecessary simulation time, user confusion

## Code Examples

### Main Integration Point
```cpp
// Source: Proposed modification to main.cpp, after readParameterFile()
// Insert between lines 59 and 61 of current main.cpp

// After config parsing, before readData
readParameterFile();

// Phase 27: McLuster IC generation (if configured)
if (mcluster_config.has_mcluster_section) {
    if (my_rank == ROOT) {
        if (!runMcluster(mcluster_config)) {
            std::cerr << "McLuster IC generation failed" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (mcluster_config.generate_only) {
            std::cout << "IC generation complete (generate_only=true). Exiting." << std::endl;
            // Clean exit for all ranks
        }
    }

    // Broadcast generate_only decision to all ranks
    int should_exit = mcluster_config.generate_only ? 1 : 0;
    MPI_Bcast(&should_exit, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    if (should_exit) {
        // Close MPI windows that were allocated
        MPI_Win_free(&win);
        MPI_Win_free(&win2);
        MPI_Win_free(&win3);
        MPI_Win_free(&win4);
        particle_data.deallocate_shared();
        MPI_Comm_free(&shared_comm);
        MPI_Type_free(&queue_type_mpi);
        MPI_Type_free(&iparticle_type_mpi);
        MPI_Type_free(&jparticle_type_mpi);
        MPI_Finalize();
        return 0;
    }
}

// Write Particles (existing code)
if (my_rank == ROOT && !readData())
    fprintf(stderr, "Read Data Failed!\n");
```

### Output File Transformation
```cpp
// Transform McLuster output to ABYSS format
bool transformMclusterOutput(const std::string& mcluster_file,
                              const std::string& abyss_file) {
    std::ifstream in(mcluster_file);
    std::ofstream out(abyss_file);

    if (!in || !out) return false;

    std::string line;
    while (std::getline(in, line)) {
        // Skip header and comment lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        double mass, x, y, z, vx, vy, vz;

        // McLuster format: mass x y z vx vy vz ...
        if (!(iss >> mass >> x >> y >> z >> vx >> vy >> vz)) {
            continue;  // Skip malformed lines
        }

        // ABYSS format: x y z vx vy vz mass
        out << std::scientific << std::setprecision(8)
            << x << " " << y << " " << z << " "
            << vx << " " << vy << " " << vz << " "
            << mass << "\n";
    }

    return out.good();
}
```

### Validation Function
```cpp
// Validate McLuster completed successfully
bool validateMclusterOutput(const std::string& output_file) {
    struct stat st;
    if (stat(output_file.c_str(), &st) != 0) {
        std::cerr << "McLuster output file not found: " << output_file << std::endl;
        return false;
    }

    if (st.st_size == 0) {
        std::cerr << "McLuster output file is empty: " << output_file << std::endl;
        return false;
    }

    // Check file has expected header
    std::ifstream in(output_file);
    std::string first_line;
    std::getline(in, first_line);

    if (first_line.find("#Mass") == std::string::npos) {
        std::cerr << "McLuster output missing expected header" << std::endl;
        return false;
    }

    // Count data lines
    int line_count = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] != '#') line_count++;
    }

    if (line_count == 0) {
        std::cerr << "McLuster output has no data lines" << std::endl;
        return false;
    }

    std::cout << "McLuster generated " << line_count << " stars" << std::endl;
    return true;
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Manual IC file creation | Automated McLuster integration | Phase 27 | Simplified workflow |
| Shell scripts calling McLuster | Direct subprocess from C++ | Phase 27 | Better error handling |
| External format conversion | Inline transformation | Phase 27 | Atomic operation |

**Deprecated/outdated:**
- McLuster -C 0/1 (Nbody6/4 formats) - ABYSS uses ASCII table format (-C 3)

## Execution Flow Diagram

```
main()
  |
  v
DefaultGlobal()
initializeMPI()
  |
  v
Parser()
readParameterFile()   <-- Phase 26: parses [mcluster] section
  |
  v
[has_mcluster_section?]
  |
  +--NO--> readData() --> simulation continues normally
  |
  +--YES (ROOT only):
        |
        v
      buildMclusterArgs()
        |
        v
      runMclusterSubprocess()  <-- Phase 27
        |
        v
      validateMclusterOutput()
        |
        v
      transformMclusterOutput()
        |
        v
      Update fname to transformed file
        |
        v
      [generate_only?]
         |
         +--YES--> MPI_Bcast(exit) --> clean exit all ranks
         |
         +--NO --> MPI_Barrier()
                     |
                     v
                   readData() --> simulation continues
```

## Open Questions

Things that couldn't be fully resolved:

1. **McLuster binary location**
   - What we know: Symlink at `src/mcluster` -> `../mcluster/mcluster_sse`
   - What's unclear: Should this be configurable? Absolute vs relative path?
   - Recommendation: Use relative path from ABYSS executable (`./src/mcluster`) for now; add config option later if needed.

2. **Temporary file cleanup**
   - What we know: McLuster creates `.txt` and `.info` files
   - What's unclear: Should these be cleaned up automatically? Kept for diagnostics?
   - Recommendation: Keep both files (useful for debugging); add cleanup option later.

3. **Error message detail level**
   - What we know: McLuster prints verbose output to stdout
   - What's unclear: How much to display to user on error?
   - Recommendation: Show last 20 lines of stderr on failure; full output available in captured string.

## Sources

### Primary (HIGH confidence)
- `/gpfs/home/vjl4366/pkg/ABYSS/src/main.cpp` - Current main execution flow
- `/gpfs/home/vjl4366/pkg/ABYSS/src/mcluster_config.h` - Config struct from Phase 26
- `/gpfs/home/vjl4366/pkg/ABYSS/src/read_parameter_file.cpp` - Config parsing from Phase 26
- `/gpfs/home/vjl4366/pkg/ABYSS/src/read_write.cpp` - Existing IC file reading (readData)
- `/gpfs/home/vjl4366/pkg/ABYSS/mcluster/README` - McLuster parameter documentation
- McLuster test runs - Verified output format with -C 3 -u 0/1 options

### Secondary (MEDIUM confidence)
- [POSIX popen() specification](https://pubs.opengroup.org/onlinepubs/9699919799/functions/popen.html) - Standard subprocess API
- [POSIX fork/exec documentation](https://man7.org/linux/man-pages/man3/popen.3.html) - Process creation patterns
- [subprocess.h library](https://github.com/sheredom/subprocess.h) - Reference implementation patterns

### Tertiary (LOW confidence)
- WebSearch results on MPI + fork compatibility - General guidance, implementation-dependent

## Metadata

**Confidence breakdown:**
- Subprocess pattern: HIGH - Standard POSIX, verified on Linux
- McLuster integration: HIGH - Tested McLuster binary, known output format
- MPI compatibility: MEDIUM - fork() after MPI_Init generally safe but implementation-dependent
- Output transformation: HIGH - Simple column reordering, well-understood

**Research date:** 2026-01-20
**Valid until:** Indefinite (stable POSIX APIs, McLuster binary interface unlikely to change)
