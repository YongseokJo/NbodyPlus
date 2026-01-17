# ABYSS Code Conventions

## Naming Conventions

### Variables and Functions

| Type | Convention | Examples |
|------|------------|----------|
| Local variables | `snake_case` | `time_step`, `num_neighbors`, `global_time` |
| Global variables | `snake_case` | `my_rank`, `num_processors`, `particles` |
| Functions | `PascalCase` or `snake_case` | `RootRoutines()`, `compute_acceleration_irr()` |
| Constants | `SCREAMING_SNAKE_CASE` | `MAX_NUM_PARTICLE`, `TIME_UNIT` |
| Types | `snake_case_t` or `PascalCase` | `cuda_real_t`, `ull_t`, `Particle` |

### Files

| Type | Convention | Examples |
|------|------------|----------|
| C++ source | `snake_case.cpp` | `root_routines.cpp`, `mpi_routines.cpp` |
| Headers | `snake_case.h` | `global.h`, `particle.h` |
| CUDA | `snake_case.cu` | `cuda_acceleration.cu` |
| Python | `snake_case.py` | `analyze_energy.py` |

## Code Style

### Indentation
- **Tabs** for indentation (visible in source files)
- Spaces for alignment within lines

### Braces
```cpp
// K&R style for control structures
if (condition) {
    // code
} else {
    // code
}

// Functions
void MyFunction() {
    // code
}
```

### Includes Order
1. Conditional feature headers (`#ifdef SEVN` includes)
2. Standard library
3. Project headers
4. MPI/CUDA headers

```cpp
#ifdef SEVN
#include "sevn.h"
#endif

#include <iostream>
#include <vector>
#include "def.h"
#include "particle.h"
#include <mpi.h>
#ifdef CUDA
#include <cuda_runtime.h>
#endif
```

## Preprocessor Patterns

### Feature Toggles
```cpp
#ifdef CUDA
    // CUDA-specific code
#endif

#ifdef SEVN
    // Stellar evolution code
#endif

#ifdef FEWBODY
    // Few-body dynamics code
#endif
```

### Debug Guards
```cpp
#ifdef DEBUG
    fprintf(stderr, "Debug: value = %d\n", value);
#endif
```

## Particle Access Patterns

### Direct Access
```cpp
Particle* ptcl = &particles[index];
ptcl->position[0] = x;
ptcl->velocity[dim] = v;
```

### Dimension Loops
```cpp
for (int dim = 0; dim < DIM; dim++) {
    ptcl->position[dim] = 0.0;
}
```

### Acceleration Array Access
```cpp
// acc_total[dimension][derivative_order]
// Order: 0=acc, 1=jerk, 2=snap, 3=crackle
ptcl->acc_total[dim][0]  // acceleration
ptcl->acc_total[dim][1]  // jerk
```

## Unit System

### Code Units (from `src/def.h`)
```cpp
#define TIME_UNIT     1e10      // in 1e10 yr
#define POSITION_UNIT 4.0       // in 4 pc
#define VELOCITY_UNIT 4e-10     // in 4e-10 pc/yr
#define MASS_UNIT     0.0001424198  // Msun where G = 1
```

### Conversion Pattern
```cpp
// To physical units
double position_pc = position * POSITION_UNIT;
double mass_msun = mass * MASS_UNIT;
double velocity_kms = velocity * velocity_unit / yr * pc / 1e5;
```

## MPI Communication Patterns

### Shared Memory Windows
```cpp
extern MPI_Win win;
extern Particle *particles;  // Shared across all ranks
```

### Root/Worker Dispatch
```cpp
if (my_rank == ROOT) {
    RootRoutines();
} else {
    WorkerRoutines();
}
```

## Error Handling

### Assertions
```cpp
assert(ptcl->is_active);
assert(time > 0);
```

### Runtime Checks
```cpp
if (!readData())
    fprintf(stderr, "Read Data Failed!\n");
```

### CUDA Error Checking
```cpp
cudaDeviceSynchronize();
// Errors typically caught via cudaGetLastError() in debug mode
```

## Output Patterns

### File Handles
```cpp
extern FILE* bin_output_file;
extern FILE* merger_output_file;
extern FILE* worker_output_file;
```

### Logging Style
```cpp
fprintf(merger_output_file, "GW driven merger happens!!! (PID: %d, PID: %d)\n", p1->pid, p2->pid);
fflush(merger_output_file);
```

## SDAR Compatibility

### Macro-Based Field Aliases
```cpp
// In group.h - temporary aliases for SDAR headers
#define Mass mass
#define Position position
#define Velocity velocity
#define PID pid
// ... include SDAR headers ...
#undef Mass
#undef Position
#undef Velocity
#undef PID
```

### SDAR Accessor Methods
```cpp
// Particle provides SDAR-compatible accessors
double* getPos() { return position; }
double* getVel() { return velocity; }
```

## Comments

### Author/Date Attribution
```cpp
// made 2024.08.12 by Eunwoo Chung
// modified by EW 2025.1.11
```

### Section Headers
```cpp
// ============================================================================
// Time variables
// ============================================================================
```

### TODO/Query Markers
```cpp
// (Query) by EW 2025.1.6
// Eunwoo: this should be fixed later
```
