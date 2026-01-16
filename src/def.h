
#define MaxNumParticle 200000
#define MaxNumNeighbor 1000 // 10000 -> 2000 modified by EW 2025.1.11
#define MaxNeighborRadius 0.2 // Let's set MaxNeighborRadius as 0.2 pc // Newly set by EW 2025.8.18 // This can be set in config file later.

// SDAR // Now configured via TOML config file (RSearch, TSearch global variables) by EW 2025.8.28
// Legacy macros - kept for reference only
// #define RSEARCH 2.5e-4 // 1e-4 // pc (default in config)
// #define TSEARCH 1e-6 // 1e-6 // Myr (default in config)

// Particle Type in line with SEVN by EW 2025.7.5
#define NO_FEEDBACK_STAR            0
#define MAIN_SEQUENCE               1
#define TERMINAL_MAIN_SEQUENCE      2
#define SHELL_H_BURNING             3
#define CORE_HE_BURNING             4
#define TERMINAL_CORE_HE_BURNING    5
#define SHELL_HE_BURNING            6
#define EMPTY                       7 // Useful, but not directly stored in ptcl->ParticleType
#define REMNANT                     8 // Useful, but not directly stored in ptcl->ParticleType
#define WHITE_DWARF_HE              9
#define WHITE_DWARF_CO              10
#define WHITE_DWARF_ONE             11
#define NEUTRON_STAR_ECSN           12
#define NEUTRON_STAR_CCSN           13
#define BLACK_HOLE                  14
#define MASSIVE_BLACK_HOLE          15

#define Dim 3
#define HERMITE_ORDER 4
#define MIN_LEVEL_BUFFER 30

// Custom type definitions
#define CUDA_FLOAT
#ifdef CUDA_FLOAT
typedef float CUDA_REAL;
#else
typedef double CUDA_REAL;
#endif
typedef unsigned long long ULL;

// Physical units in cgs
#define pc 3.08567758149137e18
#define yr 3.1536e7
#define Msun 1.98847e33

// Conversion factors ([pc, yr, Msun] = [position_unit, time_unit, mass_unit] * [code unit])
#define time_unit 1e10 // in 1e10 yr
#define position_unit 4. // in 4 pc
#define velocity_unit 4e-10 // in 4e-10 pc/yr
#define mass_unit 0.0001424198  // Msun in the unit that G = 1.

// GPU related parameters
#define nbodymax 100000000 //100000000 for node14
#define BatchSize 64 // 64. each thread calculates BatchSize particles with a single shared memory
#define GridDimY 32 // 32 original //  each block calcuates NNB/GridDimY particles
#define NNB_per_block 128 //256 original
//#define BatchSize 32 // each thread calculates BatchSize particles
//#define GridDimY 16 // each block calcuates NNB/GridDimY particles
//#define NNB_per_block 128