## Current File Structure

### Main src/ Directory
- `main.cpp` - Main entry point
- `default_global.cpp` - Default global variable definitions
- `global.h` - Global variable declarations
- `global_state.h` - GlobalVariable struct definition
- `def.h` - Constants and type definitions
- `particle.h` - Particle struct
- `parser.cpp` - Command line parser
- `read_parameter_file.cpp` - Configuration file reader
- `read_write.cpp` - I/O routines
- `mpi_routines.cpp` - MPI communication
- `initialization_routines.cpp` - Simulation initialization
- `irregular_routines.cpp` - Irregular force routines
- `regular_routines.cpp` - Regular force routines
- `root_routines.cpp` - Root process routines
- `worker_routines.cpp` - Worker process routines
- `timestep_routines.cpp` - Time step calculation
- `stellar_evolution.cpp` - Stellar evolution (SEVN)
- `queue_scheduler.h` - Queue scheduler class
- `queue.h` - Task queue definitions
- `worker.h` - Worker struct
- `skip_list.h` - Skip list data structure
- `performance.h` - Performance tracking
- `profiler.h` - Profiler
- `toml.hpp` - TOML parser (third-party)

### FewBody/ Directory
- `fb_check.cpp` - Few-body group checking
- `fb_initialization.cpp` - Few-body initialization
- `fb_integration.cpp` - Few-body integration
- `fb_termination.cpp` - Few-body termination
- `group_acceleration.cpp` - Group acceleration calculation
- `group.h` - Group struct
- `ar_interaction.hpp` - SDAR interaction
- `ar_perturber.hpp` - SDAR perturber

### Particle/ Directory
- `initialize.cpp` - Particle initialization
- `update_particle.cpp` - Particle update routines
- `compute_acceleration.cpp` - Acceleration computation

### cuda/ Directory
- `cuda_acceleration.cu` - CUDA acceleration kernels
- `cuda_kernels.cu` - CUDA kernels
- `cuda_kernels.h` - CUDA kernel declarations
- `cuda_defs.h` - CUDA type definitions
- `cuda_functions.h` - CUDA function declarations
- `cuda_global.h` - CUDA global variables
- `cuda_routines.cpp` - CUDA routine wrappers
- `cuda_routines.cu` - CUDA routines
- `cuda_routines.h` - CUDA routine declarations
- `calculate_acceleration_all.cpp` - Full acceleration calculation
- `calculate_regular_acceleration.cpp` - Regular acceleration calculation
