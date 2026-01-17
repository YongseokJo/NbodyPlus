# SoA Implementation Patterns for C++/CUDA

## Recommended Approach

### Container Design Pattern

For ABYSS with C++11 and CUDA, use a **manual SoA container with accessor functions**:

```cpp
class ParticleData {
private:
    size_t capacity_;
    size_t count_;

    // Core arrays (heap allocated)
    double* pos_x_;
    double* pos_y_;
    double* pos_z_;
    double* vel_x_;
    double* vel_y_;
    double* vel_z_;
    double* mass_;
    // ... all other fields

public:
    // Accessor functions
    void set_position(size_t i, double x, double y, double z);
    void get_position(size_t i, double& x, double& y, double& z) const;
    double* pos_x() { return pos_x_; }  // For GPU memcpy
    // etc.
};
```

### Why This Pattern

| Approach | Pros | Cons | Verdict |
|----------|------|------|---------|
| Manual SoA + accessors | Full control, C++11 compatible, explicit GPU transfers | More boilerplate | **Recommended** |
| soagen library | Less boilerplate, nice syntax | Requires C++17 | Not compatible |
| Ikra-Cpp/DynaSOAr | OOP syntax, SOA layout | Complex DSL, heavy dependency | Overkill |
| Proxy objects | AoS-like syntax | Performance overhead, complex | Avoid |

### GPU Memory Layout

Follow [NVIDIA GPU Gems 3](https://developer.nvidia.com/gpugems/gpugems3/part-v-physics-simulation/chapter-31-fast-n-body-simulation-cuda) guidance:

- Use `float4` where possible for coalesced access
- Store mass in position's w component for force calculations
- Separate arrays for each component (not `float3` arrays)

```cpp
// GPU-side SoA
struct ParticleGPU {
    float* pos_x;  // [N]
    float* pos_y;  // [N]
    float* pos_z;  // [N]
    float* mass;   // [N] (or combine with pos as float4)
    float* vel_x;  // [N]
    float* vel_y;  // [N]
    float* vel_z;  // [N]
};
```

### MPI Shared Memory Pattern

Each array gets its own `MPI_Win`:

```cpp
// Allocate multiple shared memory windows
MPI_Win_allocate_shared(N * sizeof(double), sizeof(double),
                        MPI_INFO_NULL, shared_comm, &pos_x_, &win_pos_x);
MPI_Win_allocate_shared(N * sizeof(double), sizeof(double),
                        MPI_INFO_NULL, shared_comm, &pos_y_, &win_pos_y);
// etc.
```

## What NOT To Do

1. **Don't use array-of-float3**: `float3 positions[N]` is still AoS-like, loses coalescing benefits
2. **Don't use proxy objects for hot paths**: Runtime overhead defeats purpose
3. **Don't convert AoS↔SoA at runtime**: Conversion overhead can exceed benefits for cheap kernels
4. **Don't use power-of-2 array sizes**: Cache associativity issues (pad to avoid)

## Libraries Reference

- [mini-nbody](https://github.com/harrism/mini-nbody) — Reference CUDA N-body with SoA
- [soagen](https://marzer.github.io/soagen/) — C++17 SoA library (for reference, not use)
- [DynaSOAr](https://github.com/prg-titech/dynasoar) — CUDA SoA framework (reference)

## Confidence Level

**High confidence** — These patterns are well-established in GPU computing literature and align with NVIDIA best practices.
