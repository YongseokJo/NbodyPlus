# ABYSS Technical Concerns

## Technical Debt

### 1. Incomplete SEVN Binary Integration

**Location**: `src/FewBody/fb_integration.cpp:36-64`

**Issue**: SEVN binary evolution (RLOF handling) is marked as "not tested yet":
```cpp
#ifdef SEVN_BINARY // this code is not tested yet!!!!!
```

**Impact**: Binary stellar evolution may not work correctly when enabled.

### 2. Dual Configuration Formats

**Issue**: Codebase supports both TOML (`config.toml`) and legacy plain-text (`config.txt`) formats.

**Files Affected**:
- `src/read_parameter_file.cpp`
- `tests/test1/config.toml` (new)
- `tests/test_1/config.txt` (legacy)

**Impact**: Maintenance burden; potential for format-specific bugs.

### 3. SDAR Field Name Aliasing

**Location**: `src/FewBody/group.h:13-29`

**Issue**: Uses macro defines/undefs to alias particle field names for SDAR compatibility:
```cpp
#define Mass mass
#define Position position
// ... include SDAR headers ...
#undef Mass
```

**Impact**: Fragile; prone to namespace pollution if headers are reordered.

### 4. Hard-coded GPU Parameters

**Location**: `src/def.h:94-99`

**Issue**: GPU tuning parameters are compile-time constants:
```cpp
#define BATCH_SIZE      64         // Particles per thread
#define GRID_DIM_Y      32         // Blocks
#define NNB_PER_BLOCK   128        // Neighbors per block
```

**Impact**: Requires recompilation to tune for different GPU architectures.

## Known Issues

### 1. TODO/Query Markers in Code

Multiple unresolved items marked in code:

```cpp
// (Query) by EW 2025.1.6
// Eunwoo: this should be fixed later
// commented out by EW 2025.9.9 // not tested yet!!!!!
```

**Locations**:
- `src/FewBody/fb_integration.cpp`
- `src/particle.h`

### 2. README Mentions Ongoing Work

**Location**: `README.md:58-62`

```
I am rewriting the I/O of this repo. I am going to use TOML for input
and hdf5 for output (time series).

refer to test/test1/config and output.

revise the codes to do this.
```

**Impact**: Documentation may be out of sync with current implementation.

## Performance Considerations

### 1. Neighbor List Memory

**Location**: `src/def.h:6-7`

```cpp
#define MAX_NUM_PARTICLE 200000
#define MAX_NUM_NEIGHBOR 1000  // 10000 -> 2000 modified by EW 2025.1.11
```

**Concern**: Pre-allocated neighbor arrays (200K * 1K = 200M ints for neighbors).

### 2. Multi-GPU Load Balancing

**Location**: `src/cuda/cuda_acceleration.cu`

**Concern**: Current implementation splits work by particle count but may not balance GPU load evenly for non-uniform particle distributions.

### 3. Regular Time Finding

**Location**: `src/root_routines.cpp:186-215`

**Issue**: `updateNextRegTime()` iterates all particles each step:
```cpp
for (int i=0; i<=last_particle_index; i++) {
    // find minimum regular time
}
```

**Concern**: O(N) per step; could use heap/priority queue for O(log N).

## Security Considerations

### 1. File Path Handling

**Issue**: Command-line paths (`-c config.txt`) and file I/O don't appear to validate paths.

**Impact**: Low risk for scientific computing context, but could be hardened.

### 2. Output File Overwrites

**Issue**: Output files are created without checking for existing files.

## Code Quality

### 1. Inconsistent Function Naming

Mixed conventions:
- `RootRoutines()` (PascalCase)
- `compute_acceleration_irr()` (snake_case)
- `check_new_group_v4()` (snake_case with version suffix)

### 2. Large Functions

**Example**: `Merge()` in `src/FewBody/fb_integration.cpp` is 300+ lines handling multiple merger types.

**Suggestion**: Could be refactored into smaller, type-specific functions.

### 3. Debug Output Cleanup

**Issue**: Various `fprintf(stdout, ...)` debug statements throughout code.

**Locations**:
- `src/FewBody/fb_integration.cpp` (merger debugging)
- `src/stellar_evolution.cpp`

## Missing Features

### 1. Formal Test Suite

No automated test runner or CI/CD pipeline. Tests are run manually.

### 2. Checkpoint/Restart

**Status**: Configuration option exists but implementation may be incomplete:
```toml
[restart]
Enabled = false
# CheckpointFile = "output/checkpoint_0.h5"
```

### 3. Documentation

- API documentation is minimal
- In-code comments vary in completeness
- SDAR library has Doxygen docs but may be outdated

## Recommendations

1. **Add unit tests** for critical physics routines (energy calculation, force computation)
2. **Consolidate configuration** to TOML-only, deprecate plain-text format
3. **Address TODO/Query markers** systematically
4. **Profile and optimize** regular time finding for large N
5. **Add runtime GPU parameter tuning** for different architectures
6. **Implement proper checkpoint/restart** for long-running simulations
