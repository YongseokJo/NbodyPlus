# ABYSS Directory Structure

## Root Layout

```
ABYSS/
├── src/                    # Main source code
├── SDAR/                   # Few-body integration library (submodule/vendored)
├── tools/                  # Python analysis scripts
├── tests/                  # Test cases and initial conditions
├── workflow/               # Build and run automation
├── Makefile                # Root makefile (delegates to src/)
├── README.md               # Project documentation
└── .planning/              # GSD planning documents
```

## Source Code (`src/`)

```
src/
├── main.cpp                 # Entry point, MPI init
├── def.h                    # Global constants, units, GPU parameters
├── global.h                 # Global variable declarations
├── global_state.h           # Shared simulation state struct
├── particle.h               # Particle data structure
├── Makefile                 # Build configuration
│
├── # Core Routines
├── root_routines.cpp        # Root process main loop
├── worker_routines.cpp      # Worker process routines
├── initialization_routines.cpp  # Particle/simulation setup
├── mpi_routines.cpp         # MPI initialization and communication
│
├── # Time Integration
├── regular_routines.cpp     # Regular (far-field) integration
├── irregular_routines.cpp   # Irregular (near-field) integration
├── timestep_routines.cpp    # Timestep calculation
│
├── # Particle Operations
├── Particle/
│   ├── initialize.cpp       # Particle initialization
│   ├── update_particle.cpp  # Position/velocity updates
│   └── compute_acceleration.cpp  # CPU force calculation
│
├── # Few-Body Dynamics
├── FewBody/
│   ├── group.h              # Few-body group structure
│   ├── fb_check.cpp         # Group detection logic
│   ├── fb_initialization.cpp  # Group setup
│   ├── fb_integration.cpp   # SDAR integration, PN corrections, mergers
│   ├── fb_termination.cpp   # Group dissolution
│   ├── group_acceleration.cpp  # Group force calculation
│   ├── ar_interaction.hpp   # SDAR interaction functor
│   └── ar_perturber.hpp     # SDAR perturber interface
│
├── # GPU Code
├── cuda/
│   ├── cuda_defs.h          # GPU type definitions
│   ├── cuda_global.h        # GPU global variables
│   ├── cuda_functions.h     # GPU function declarations
│   ├── cuda_routines.h      # GPU routine declarations
│   ├── cuda_kernels.h       # Kernel declarations
│   ├── cuda_acceleration.cu # Multi-GPU force calculation
│   ├── cuda_kernels.cu      # CUDA kernels
│   ├── cuda_routines.cpp    # GPU initialization/cleanup
│   └── calculate_regular_acceleration.cpp  # GPU dispatch
│
├── # I/O
├── read_write.cpp           # HDF5 output
├── read_parameter_file.cpp  # TOML/config parsing
├── parser.cpp               # Command-line argument parsing
│
├── # Other
├── stellar_evolution.cpp    # SEVN integration (when USE_SEVN=1)
├── queue_scheduler.h        # Task queue for workers
├── queue.h                  # Queue data structure
├── worker.h                 # Worker process state
├── skip_list.h              # Skip list data structure
├── performance.h            # Performance tracking
├── profiler.h               # Profiling macros
└── default_global.cpp       # Default global variable values
```

## SDAR Library (`SDAR/`)

Vendored few-body regularization library:

```
SDAR/
├── src/
│   ├── AR/                  # Algorithmic Regularization
│   │   ├── symplectic_integrator.h  # Time-transformed integrator
│   │   ├── force.h          # Force computation
│   │   ├── information.h    # Binary tree info
│   │   ├── slow_down.h      # Slowdown factor
│   │   └── profile.h        # AR profiling
│   ├── Hermite/             # Hermite integrator components
│   ├── Common/              # Shared utilities
│   │   ├── Float.h          # Floating point types
│   │   ├── binary_tree.h    # Binary tree structure
│   │   ├── matrix.h         # Matrix operations
│   │   └── particle_group.h # Particle group management
│   └── README.md
├── sample/                  # Example implementations
│   ├── AR/                  # AR examples
│   ├── Hermite/             # Hermite examples
│   └── Kepler/              # Kepler solver examples
└── docs/                    # Doxygen documentation
```

## Tools (`tools/`)

Python analysis scripts:

```
tools/
├── analyze_energy.py        # Energy conservation analysis
├── analyze_profiling.py     # Performance profiling analysis
├── read_hdf5_output.py      # HDF5 output reader
└── summarize_run.py         # Run summary generator
```

## Tests (`tests/`)

```
tests/
├── test1/                   # Primary test case
│   ├── config.toml          # TOML configuration
│   └── nbody.dat            # Initial conditions
├── test_1/, test_2/         # Legacy test cases
├── test_10/, test_20/       # Larger test cases
└── ICs/                     # Initial condition files
    ├── c1e4.dat             # 10^4 particle IC
    ├── c1e5.dat             # 10^5 particle IC
    └── core_enzo_*.dat      # Enzo-generated ICs
```

## Workflow (`workflow/`)

Build and run automation:

```
workflow/
├── config.sh                # Configuration (cluster-specific)
├── README.md                # Workflow documentation
├── bin/                     # Automation scripts
│   ├── build.sh             # Compile script
│   ├── run.sh               # Run script
│   └── analyze.sh           # Analysis script
└── runs/                    # Run output directories
    └── run_YYYYMMDD_HHMMSS/
        ├── config.sh        # Run-specific config
        └── work/            # Working directory
```

## Key File Locations

| Purpose | Path |
|---------|------|
| Main entry | `src/main.cpp` |
| Particle definition | `src/particle.h` |
| Constants/units | `src/def.h` |
| Build config | `src/Makefile` |
| Cluster config | `workflow/config.sh` |
| Test configuration | `tests/test1/config.toml` |
| Energy analysis | `tools/analyze_energy.py` |
| GPU kernels | `src/cuda/cuda_kernels.cu` |
| Few-body integration | `src/FewBody/fb_integration.cpp` |
