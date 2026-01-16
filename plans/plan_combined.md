# ABYSS Combined Optimization Plan (Summary)

## Scope and constraints (high confidence)
- Do not modify any code under SDAR/ or SEVN paths.
- Strict numerical precision required (avoid fast-math).

## Where the plans agree (highest confidence)
1. **CUDA kernel tuning and reduction optimization**
   - Target kernel launch parameters and inner-loop unroll for the main GPU force kernel in [src/cuda/cuda_kernels.cu](src/cuda/cuda_kernels.cu).
   - Parallelize the reduction in the regular-force reduction kernel (currently serial over n) in [src/cuda/cuda_kernels.cu](src/cuda/cuda_kernels.cu).

2. **GPU/CPU overlap and buffer reuse**
   - Reuse staging buffers and avoid per-step vector allocations in [src/cuda/CalculateRegularAcceleration.cpp](src/cuda/CalculateRegularAcceleration.cpp).
   - Increase async overlap (H2D/D2H vs compute) in [src/cuda/cuda_my_acceleration.cu](src/cuda/cuda_my_acceleration.cu).

3. **MPI non-blocking task distribution**
   - Replace blocking sends with non-blocking in task assignment in [src/MPIRoutines.cpp](src/MPIRoutines.cpp).

4. **OpenMP for data-parallel CPU loops**
   - Add OpenMP pragmas to per-particle update loops in [src/IrregularRoutines.cpp](src/IrregularRoutines.cpp) and relevant regular-update loops in [src/RegularRoutines.cpp](src/RegularRoutines.cpp) if thread-safety holds.

5. **Build flags for performance**
   - Use -O3, -DNDEBUG, -flto, and enable OpenMP in release builds in [src/Makefile](src/Makefile).

## Conflicts and risk notes (medium confidence)
- **Fast-math:** One plan recommends CUDA fast math; this conflicts with strict precision requirements. Treat fast-math as out of scope unless explicitly approved.
- **SoA refactor:** One plan proposes an SoA layout for Jparticles. This can be high-impact but high-effort and invasive. It should be considered a phase-3+ refactor after profiling and kernel-tuning wins are confirmed.

## Reliability notes (based on codebase alignment)
- The cited files and functions exist in the repository and appear to match the plans’ targets:
  - GPU kernels: [src/cuda/cuda_kernels.cu](src/cuda/cuda_kernels.cu)
  - GPU staging and overlap: [src/cuda/cuda_my_acceleration.cu](src/cuda/cuda_my_acceleration.cu)
  - GPU data gathering: [src/cuda/CalculateRegularAcceleration.cpp](src/cuda/CalculateRegularAcceleration.cpp)
  - MPI dispatch: [src/MPIRoutines.cpp](src/MPIRoutines.cpp)
  - CPU update loops: [src/IrregularRoutines.cpp](src/IrregularRoutines.cpp), [src/RegularRoutines.cpp](src/RegularRoutines.cpp)
- The recommendations are **optimization proposals** rather than verified measurements. Actual gains depend on hardware, problem size, and current compile/runtime settings.
- Any MPI, OpenMP, or CUDA changes should be validated with existing tests and profiling.

## Consolidated phased plan

### Phase 0: Baseline (high confidence)
- Capture timing baselines with existing profiling or logs.
- Record GPU kernel timing before changes.

### Phase 1: Low-risk performance wins (high confidence)
1. MPI non-blocking task assignment in [src/MPIRoutines.cpp](src/MPIRoutines.cpp).
2. Reuse host buffers for GPU staging in [src/cuda/CalculateRegularAcceleration.cpp](src/cuda/CalculateRegularAcceleration.cpp).
3. Use pinned host buffers and more async transfers in [src/cuda/cuda_my_acceleration.cu](src/cuda/cuda_my_acceleration.cu).
4. Release build flags in [src/Makefile](src/Makefile) without fast-math.

### Phase 2: CPU/GPU parallelism (medium confidence)
1. OpenMP parallel loops in [src/IrregularRoutines.cpp](src/IrregularRoutines.cpp) and [src/RegularRoutines.cpp](src/RegularRoutines.cpp).
2. Kernel parameter tuning and unroll in [src/cuda/cuda_kernels.cu](src/cuda/cuda_kernels.cu).
3. Parallel reduction kernel in [src/cuda/cuda_kernels.cu](src/cuda/cuda_kernels.cu).

### Phase 3: Structural refactors (lower confidence / higher effort)
1. SoA data layout for Jparticles in [src/cuda/cuda_kernels.cu](src/cuda/cuda_kernels.cu) and related staging in [src/cuda/CalculateRegularAcceleration.cpp](src/cuda/CalculateRegularAcceleration.cpp).
2. Task batching in queue distribution (if present in scheduler paths).

## Suggested verification
- Build: make clean && make USE_CUDA=1
- Regression: existing tests under tests/ (target a baseline run such as tests/test1)
- Performance: compare wall-clock and GPU kernel timings before/after each phase
