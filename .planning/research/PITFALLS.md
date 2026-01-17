# SoA Conversion Pitfalls

## Critical Pitfalls

### 1. Conversion Overhead for Cheap Kernels

**Problem**: If you convert AoS↔SoA at runtime for computationally cheap operations, the conversion overhead exceeds the performance benefit.

**Warning Signs**:
- Temporary SoA buffers created per timestep
- Converting data just before simple loops

**Prevention**:
- Store data in SoA permanently, no runtime conversion
- ABYSS approach: Full SoA with accessor functions

**Phase**: Address in Phase 1 (architecture design)

---

### 2. Cache Associativity with Power-of-2 Sizes

**Problem**: When array sizes are large powers of 2, cache associativity causes severe performance degradation due to conflict misses.

**Warning Signs**:
- Particle count is exactly 65536, 131072, etc.
- Unexplained slowdowns at specific N values

**Prevention**:
- Add padding to avoid power-of-2 array sizes
- Use `capacity = (requested + 15) & ~15` or similar

**Phase**: Address in Phase 1 (memory allocation)

```cpp
size_t padded_capacity(size_t n) {
    // Avoid exact powers of 2
    if ((n & (n - 1)) == 0 && n > 64) {
        return n + 64;  // Add padding
    }
    return n;
}
```

---

### 3. Incomplete Accessor Coverage

**Problem**: Missing accessor functions forces developers to bypass encapsulation, leading to mixed access patterns.

**Warning Signs**:
- Code that directly accesses `particle_data.pos_x_[i]`
- Inconsistent accessor usage across files

**Prevention**:
- Generate complete accessor set upfront
- Make internal arrays private
- Compile with warnings for unused functions

**Phase**: Address in Phase 2 (accessor implementation)

---

### 4. GPU Memory Transfer Fragmentation

**Problem**: Many small `cudaMemcpy` calls for individual arrays is slower than fewer large transfers.

**Warning Signs**:
- One `cudaMemcpy` per field
- Transfer overhead dominates kernel time

**Prevention**:
- Batch transfers where possible
- Use CUDA streams for overlapping
- Consider unified memory for small particle counts

**Phase**: Address in Phase 4 (GPU integration)

```cpp
// Bad: Many small transfers
cudaMemcpy(d_pos_x, h_pos_x, N*sizeof(float), ...);
cudaMemcpy(d_pos_y, h_pos_y, N*sizeof(float), ...);
cudaMemcpy(d_pos_z, h_pos_z, N*sizeof(float), ...);

// Better: Use streams for overlap
cudaMemcpyAsync(d_pos_x, h_pos_x, N*sizeof(float), ..., stream1);
cudaMemcpyAsync(d_pos_y, h_pos_y, N*sizeof(float), ..., stream2);
```

---

### 5. SDAR Integration Breakage

**Problem**: SDAR library expects `Particle` objects with specific field names. Converting to pure SoA breaks the interface.

**Warning Signs**:
- Compile errors in SDAR headers
- Runtime crashes in few-body integration

**Prevention**:
- Keep minimal `Particle` struct for SDAR compatibility
- Use proxy/adapter pattern for SDAR calls
- Sync data before/after SDAR operations

**Phase**: Address in Phase 6 (SDAR compatibility)

---

### 6. MPI Window Explosion

**Problem**: Creating 60+ MPI_Win objects (one per array) can stress MPI implementation limits or cause initialization overhead.

**Warning Signs**:
- Slow MPI initialization
- MPI errors about resource limits

**Prevention**:
- Test with target MPI implementation
- Consider grouping related arrays into larger allocations
- Monitor MPI_Win creation time in profiling

**Phase**: Address in Phase 3 (MPI integration)

---

### 7. Energy Conservation Regression

**Problem**: Subtle bugs in data access cause numerical differences that break energy conservation.

**Warning Signs**:
- dE/E0 suddenly increases after conversion
- Different results at high precision

**Prevention**:
- Run baseline test BEFORE any changes
- Compare bit-for-bit after each phase
- Use `analyze_energy.py` as regression test

**Phase**: Validate after EVERY phase

---

### 8. Mixed Float/Double Precision

**Problem**: GPU uses `cuda_real_t` (float) while CPU uses double. Conversion between precision can accumulate errors.

**Warning Signs**:
- Energy conservation differs between CPU-only and GPU runs
- Small timestep errors accumulate

**Prevention**:
- Be explicit about precision at each interface
- Minimize CPU↔GPU round-trips
- Consider double precision on GPU for critical paths

**Phase**: Address in Phase 4 (GPU integration)

---

## Testing Checklist

After each phase, verify:

- [ ] Code compiles with `-Wall -Wextra` without new warnings
- [ ] `tests/test1` runs to completion
- [ ] Energy conservation within baseline tolerance
- [ ] No segfaults or memory errors (run with `valgrind` or `cuda-memcheck`)
- [ ] Performance benchmark (don't regress)

## Phase-Pitfall Mapping

| Phase | Pitfalls to Address |
|-------|---------------------|
| 1. Container design | #1 Conversion overhead, #2 Cache associativity |
| 2. Accessors | #3 Incomplete coverage |
| 3. MPI | #6 Window explosion |
| 4. GPU | #4 Transfer fragmentation, #8 Precision |
| 5. CPU routines | #7 Energy conservation |
| 6. SDAR | #5 Integration breakage |
| 7. Testing | All validation |
