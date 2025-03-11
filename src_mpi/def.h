#define FAIL    -1
#define SUCCESS  1


#define NumberOfTask 20

#define MaxNumberOfParticle 200000
#define MaxNumberOfCommunication 10000


//#define FixNumNeighbor 100 // 500 -> 100 modified by EW 2025.1.11
#define MaxNumNeighbor 1000 // 10000 -> 2000 modified by EW 2025.1.11
//#define ACRadius 0.013 // 0.11 -> 0.05 modified by EW 2025.1.11


// SDAR
#define RSEARCH 2e-4 // 1e-4 // pc
#define TSEARCH 1e-6 // 1e-6 // Myr

#define NormalStar 1
#define Blackhole 32
#define SingleStar 64


#define MIN_LEVEL_BUFFER 30

#define Dim 3
//#define eta 0.01
#define HERMITE_ORDER 4


typedef unsigned long long ULL;



#define mag(a) (a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
#define mag0(a) (a[0][0]*a[0][0]+a[1][0]*a[1][0]+a[2][0]*a[2][0])
#define dist(a,b) std::sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))



// Physical units in cgs
#define pc 3.08567758149137e18
#define yr 3.1536e7
#define Msun 1.98847e33

// Code unit in pc, yr, Msun
#define time_unit 1e10 // in 1e10 yr
#define position_unit 4. // in 4 pc
#define velocity_unit 4e-10 // in 4e-10 pc/yr
//#define mass_unit 256e-20  // 256e-20 Msun in the unit that G = 1.
#define mass_unit 0.0001424198  // Msun in the unit that G = 1.
// Physical constants
// #define G_cgs 6.67430e-8 // Eunwoo: crash with SEVN
// #define G // pc, yr, Msun // Eunwoo: crash with SEVN

#define RCAST(a)  static_cast<double>(a)
#define ABS(a) static_cast<double>(std::abs(a))
#define MIN(a,b) std::min(RCAST(a),RCAST(b))

#define CUDA_FLOAT
#ifdef CUDA_FLOAT
#define CUDA_REAL float
#else
#define CUDA_REAL double
#endif

#define nbodymax 100000000 //100000000 for node14
#define BatchSize 64 // 64. each thread calculates BatchSize particles with a single shared memory
#define GridDimY 32 // 32 original //  each block calcuates NNB/GridDimY particles
#define NNB_per_block 128 //256 original
#define MultiGPU // for multi-gpu
//#define BatchSize 32 // each thread calculates BatchSize particles
//#define GridDimY 16 // each block calcuates NNB/GridDimY particles
//#define NNB_per_block 128
