# Research Summary: AoS → SoA Conversion

## Key Findings

### Stack Recommendation

**Manual SoA container with accessor functions** is the recommended approach for ABYSS:
- C++11 compatible (no modern library dependencies)
- Full control over memory layout
- Explicit GPU transfers
- Each array gets its own MPI shared memory window

Reference implementations: [mini-nbody](https://github.com/harrism/mini-nbody), [NVIDIA GPU Gems 3](https://developer.nvidia.com/gpugems/gpugems3/part-v-physics-simulation/chapter-31-fast-n-body-simulation-cuda)

### Conversion Scope

| Priority | Arrays | Description |
|----------|--------|-------------|
| Critical | ~42 | Position, velocity, mass, accelerations — accessed in O(N²) loops |
| Important | ~12 | Timestep data — accessed per-timestep |
| Low | ~8 | Metadata — sparse access |
| Skip | ~10 | Pointers, SDAR/SEVN objects |

**Total: ~62 new arrays** for full conversion

### Architecture Highlights

1. **ParticleData class**: Central SoA container with accessor functions
2. **MPI**: Multiple `MPI_Win` objects (one per array) for shared memory
3. **GPU**: Separate device arrays with batched async transfers
4. **SDAR**: Proxy/adapter pattern to maintain compatibility

### Critical Pitfalls

| Pitfall | Impact | Prevention |
|---------|--------|------------|
| Runtime AoS↔SoA conversion | Defeats performance gains | Store in SoA permanently |
| Power-of-2 array sizes | Cache conflict misses | Pad allocations |
| SDAR breakage | Few-body integration fails | Use proxy pattern |
| Energy conservation | Silent correctness bug | Validate after each phase |

## Recommended Phase Structure

1. **Core SoA container** — ParticleData class with all arrays
2. **Accessor functions** — Complete coverage for all fields
3. **MPI integration** — Shared memory windows
4. **GPU integration** — Device arrays, kernel updates
5. **CPU routines** — Force calculations, integration
6. **SDAR compatibility** — Proxy pattern
7. **I/O updates** — HDF5 output
8. **Validation** — Full test suite

## Performance Expectations

Based on research:
- **GPU kernels**: 2-4x improvement from memory coalescing
- **CPU loops**: 1.5-2x improvement from cache efficiency
- **Irregular patterns**: May not benefit (SDAR, sparse access)

Measure baseline before changes to quantify actual gains.

## Sources

- [AoS vs SoA in practice: particle simulation](https://isocpp.org/blog/2025/05/aos-vs-soa-in-practice-particle-simulation-vittorio-romeo)
- [Annotation-guided AoS-to-SoA conversions](https://onlinelibrary.wiley.com/doi/full/10.1002/cpe.70199)
- [AoS and SoA - Algorithmica](https://en.algorithmica.org/hpc/cpu-cache/aos-soa/)
- [NVIDIA GPU Gems 3: Fast N-Body Simulation](https://developer.nvidia.com/gpugems/gpugems3/part-v-physics-simulation/chapter-31-fast-n-body-simulation-cuda)
