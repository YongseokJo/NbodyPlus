
// ============================================================================
// Array size limits
// ============================================================================
#define MAX_NUM_PARTICLE 200000
#define MAX_NUM_NEIGHBOR 1000  // 10000 -> 2000 modified by EW 2025.1.11
#define MAX_NEIGHBOR_RADIUS 0.2  // pc, can be set in config file later

// ============================================================================
// Particle types (aligned with SEVN by EW 2025.7.5)
// ============================================================================
#define STAR_TYPE_NO_FEEDBACK       0
#define STAR_TYPE_MAIN_SEQUENCE     1
#define STAR_TYPE_TERMINAL_MS       2
#define STAR_TYPE_SHELL_H_BURN      3
#define STAR_TYPE_CORE_HE_BURN      4
#define STAR_TYPE_TERMINAL_HE       5
#define STAR_TYPE_SHELL_HE_BURN     6
#define STAR_TYPE_EMPTY             7  // Not stored in particle_type
#define STAR_TYPE_REMNANT           8  // Not stored in particle_type
#define STAR_TYPE_WD_HE             9
#define STAR_TYPE_WD_CO             10
#define STAR_TYPE_WD_ONE            11
#define STAR_TYPE_NS_ECSN           12
#define STAR_TYPE_NS_CCSN           13
#define STAR_TYPE_BH                14
#define STAR_TYPE_MASSIVE_BH        15

// ============================================================================
// Simulation dimensions and integration order
// ============================================================================
#define DIM 3
#define HERMITE_ORDER 4
#define MIN_LEVEL_BUFFER 30

// ============================================================================
// Custom type definitions
// ============================================================================
#define CUDA_FLOAT
#ifdef CUDA_FLOAT
typedef float cuda_real_t;
#else
typedef double cuda_real_t;
#endif
typedef unsigned long long ull_t;

// ============================================================================
// Physical constants (cgs units)
// ============================================================================
#define PC_TO_CM   3.08567758149137e18
#define YR_TO_SEC  3.1536e7
#define MSUN_TO_G  1.98847e33

// ============================================================================
// Unit conversion factors
// [pc, yr, Msun] = [position_unit, time_unit, mass_unit] * [code_unit]
// ============================================================================
#define TIME_UNIT     1e10      // in 1e10 yr
#define POSITION_UNIT 4.0       // in 4 pc
#define VELOCITY_UNIT 4e-10     // in 4e-10 pc/yr
#define MASS_UNIT     0.0001424198  // Msun where G = 1

// ---------------------------------------------------------------------------
// Backward-compat aliases (legacy lowercase names used across the codebase)
// ---------------------------------------------------------------------------
#ifndef time_unit
#define time_unit TIME_UNIT
#endif
#ifndef position_unit
#define position_unit POSITION_UNIT
#endif
#ifndef velocity_unit
#define velocity_unit VELOCITY_UNIT
#endif
#ifndef mass_unit
#define mass_unit MASS_UNIT
#endif

// Legacy unit scalars used in velocity conversions
#ifndef yr
#define yr YR_TO_SEC
#endif
#ifndef pc
#define pc PC_TO_CM
#endif

// Legacy remnant threshold name
#ifndef REMNANT
#define REMNANT STAR_TYPE_REMNANT
#endif

// ============================================================================
// GPU parameters (tuned for NVIDIA A100)
// ============================================================================
#define NBODY_MAX       100000000  // Maximum particles on GPU
#define BATCH_SIZE      64         // Particles per thread with shared memory
#define GRID_DIM_Y      32         // Blocks calculate NNB/GRID_DIM_Y particles
#define NNB_PER_BLOCK   128        // Neighbors per block

