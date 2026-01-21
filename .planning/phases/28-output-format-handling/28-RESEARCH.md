# Phase 28: Output Format Handling - Research

**Researched:** 2026-01-20
**Domain:** McLuster output format transformation, unit conversion, ABYSS IC file compatibility
**Confidence:** HIGH

## Summary

This phase addresses a **CRITICAL unit mismatch** between McLuster output and ABYSS IC file expectations. The Phase 27 implementation already handles column reordering (mass first -> mass last), but the unit conversion is incorrect for the current transformation approach.

The core issue: ABYSS's `normalize_particle()` function expects input in specific units that differ from McLuster's astrophysical output units. The current implementation uses `-u 1` (astrophysical units) but the transformation does not account for the unit differences ABYSS expects.

**Key findings:**
1. ABYSS IC files expect: position in **kpc**, velocity in **km/s**, mass in **10^-9 Msun** units
2. McLuster `-u 1` outputs: position in **pc**, velocity in **km/s**, mass in **Msun**
3. McLuster `-u 0` outputs: **N-body units** (mass sums to 1, specific scaling)

**Primary recommendation:** Fix the unit conversion in `transformMclusterOutput()` to convert McLuster's astrophysical units (pc, km/s, Msun) to ABYSS expected units (kpc, km/s, 10^-9 Msun). This is a straightforward arithmetic fix requiring position/1000 and mass/1e9.

## Standard Stack

No additional libraries needed - this phase uses existing file I/O utilities.

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `<fstream>` | Standard C++ | File reading/writing | Already used in mcluster_runner.cpp |
| `<sstream>` | Standard C++ | Line parsing | Already used in transformation |
| `<iomanip>` | Standard C++ | Numeric formatting | Already used for scientific notation |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `<cmath>` | Standard C++ | Numeric operations if needed | Edge case handling |

## Architecture Patterns

### Current Implementation Structure
```
src/
├── mcluster_runner.cpp     # Contains transformMclusterOutput()
├── mcluster_runner.h       # Declares transformation function
├── particle.h              # Contains normalize_particle()
└── read_write.cpp          # Contains readData()
```

### Pattern 1: Unit Conversion in Transformation
**What:** Apply unit conversions during file transformation, before normalize_particle runs
**When to use:** When source and target have predictable unit systems
**Example:**
```cpp
// Source: Required fix for mcluster_runner.cpp
bool transformMclusterOutput(const std::string& mcluster_file,
                              const std::string& abyss_file) {
    // ... file open code ...

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        double mass, x, y, z, vx, vy, vz;
        if (!(iss >> mass >> x >> y >> z >> vx >> vy >> vz)) continue;

        // McLuster outputs: mass (Msun), position (pc), velocity (km/s)
        // ABYSS expects: mass (1e-9 Msun), position (kpc), velocity (km/s)
        double mass_abyss = mass / 1e9;    // Msun -> 1e-9 Msun units
        double x_abyss = x / 1000.0;       // pc -> kpc
        double y_abyss = y / 1000.0;
        double z_abyss = z / 1000.0;
        // Velocity stays the same (km/s -> km/s)

        out << std::scientific << std::setprecision(8)
            << x_abyss << " " << y_abyss << " " << z_abyss << " "
            << vx << " " << vy << " " << vz << " "
            << mass_abyss << "\n";

        ++lines_written;
    }
    // ... cleanup code ...
}
```

### Pattern 2: Alternative - Use N-body Units from McLuster
**What:** Configure McLuster to output N-body units (-u 0) and skip unit conversion entirely
**When to use:** When N-body normalized output is acceptable
**Trade-off:** Total mass sums to 1 in N-body units, requires understanding McLuster's internal scaling

### Anti-Patterns to Avoid
- **Double conversion:** Do NOT convert units in transformation AND expect normalize_particle to convert again
- **Ignoring units entirely:** The current implementation silently outputs wrong units
- **Hardcoded magic numbers without comments:** Always document what units are being converted

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Precision loss | Custom formatting | `std::setprecision(15)` with double | Scientific notation preserves precision |
| Line parsing | Character-by-character | `std::istringstream` | Already used, handles whitespace |
| File validation | Multiple stat calls | Existing validateMclusterOutput() | Already implemented |

**Key insight:** The core fix is simple arithmetic (divide by 1000, divide by 1e9). The complexity is in understanding the unit chain correctly.

## Common Pitfalls

### Pitfall 1: Unit Chain Confusion
**What goes wrong:** Applying wrong conversions because unit expectations are unclear
**Why it happens:** Three different unit systems involved (McLuster astrophysical, McLuster N-body, ABYSS internal)
**How to avoid:** Document the full unit chain in comments:
```
McLuster (-u 1): Msun, pc, km/s
      |
      v (transformation)
ABYSS input: 1e-9 Msun, kpc, km/s
      |
      v (normalize_particle)
ABYSS internal: code units (MASS_UNIT, POSITION_UNIT, VELOCITY_UNIT)
```
**Warning signs:** Particles clustered at origin (positions too small), extreme velocities, unrealistic masses

### Pitfall 2: Precision Loss in Output
**What goes wrong:** Positions near zero truncated, losing significant digits
**Why it happens:** Default double formatting doesn't preserve enough digits
**How to avoid:** Use `std::setprecision(15)` or at least `std::setprecision(8)` with scientific notation
**Warning signs:** Particles at exact same position, simulations exploding immediately

### Pitfall 3: Silent File Path Issues
**What goes wrong:** Generated IC file not found or placed in wrong directory
**Why it happens:** McLuster runs in working directory, ABYSS expects file at specific path
**How to avoid:** Use absolute paths or ensure consistent working directory
**Warning signs:** "File not found" errors, ABYSS loading old IC file

### Pitfall 4: Header Line Processing
**What goes wrong:** McLuster header line parsed as data, causing malformed particle
**Why it happens:** Header starts with `#` but code might not skip it
**How to avoid:** Explicit check for comment lines (current code does this correctly)
**Warning signs:** First particle has garbage values, NaN propagation

## Code Examples

Verified patterns from official sources:

### Current Transformation (BUGGY - for reference)
```cpp
// Source: mcluster_runner.cpp (current implementation)
// BUG: No unit conversion - outputs McLuster units directly
out << std::scientific << std::setprecision(8)
    << x << " " << y << " " << z << " "
    << vx << " " << vy << " " << vz << " "
    << mass << "\n";
```

### Fixed Transformation (REQUIRED)
```cpp
// Source: Proposed fix for mcluster_runner.cpp
// Convert McLuster (Msun, pc, km/s) to ABYSS input (1e-9 Msun, kpc, km/s)
out << std::scientific << std::setprecision(15)
    << (x / 1000.0) << " "        // pc -> kpc
    << (y / 1000.0) << " "
    << (z / 1000.0) << " "
    << vx << " " << vy << " " << vz << " "  // km/s unchanged
    << (mass / 1e9) << "\n";      // Msun -> 1e-9 Msun units
```

### ABYSS normalize_particle Reference
```cpp
// Source: particle.h lines 266-275
// This runs AFTER reading the IC file - expects specific input units
void normalize_particle() {
    this->mass *= 1e9;                        // Input in 1e-9 Msun -> Msun
    this->mass /= MASS_UNIT;                  // Msun -> code units
    for (int dim = 0; dim < DIM; dim++) {
        this->position[dim] *= 1000;          // Input in kpc -> pc
        this->position[dim] /= POSITION_UNIT; // pc -> code units
        this->velocity[dim] *= 1e5 * YR_TO_SEC / PC_TO_CM;  // km/s -> pc/yr
        this->velocity[dim] /= VELOCITY_UNIT; // pc/yr -> code units
    }
}
```

### McLuster Output Format Reference
```
#Mass_[Msun] x_[pc] y_[pc] z_[pc] vx_[km/s] vy_[km/s] vz_[km/s] ...
0.593310     0.29   -0.41  0.62   0.05      0.11      -0.07    ...
```

### ABYSS Expected IC File Format
```
x(kpc) y(kpc) z(kpc) vx(km/s) vy(km/s) vz(km/s) mass(1e-9 Msun)
0.00029 -0.00041 0.00062 0.05 0.11 -0.07 5.9331e-10
```

## Unit Conversion Reference

### McLuster Astrophysical Units (-u 1)
| Quantity | McLuster Unit | Symbol |
|----------|---------------|--------|
| Mass | Solar mass | Msun |
| Position | Parsec | pc |
| Velocity | km/s | km/s |

### ABYSS IC File Expected Units
| Quantity | ABYSS Input Unit | Symbol |
|----------|------------------|--------|
| Mass | 10^-9 Solar mass | 1e-9 Msun |
| Position | Kiloparsec | kpc |
| Velocity | km/s | km/s |

### Conversion Factors
| From | To | Multiply by |
|------|-----|-------------|
| Msun | 1e-9 Msun units | 1/1e9 = 1e-9 |
| pc | kpc | 1/1000 = 0.001 |
| km/s | km/s | 1 (unchanged) |

### Full Pipeline
```
McLuster Output (Msun, pc, km/s)
        |
        | transformMclusterOutput()
        | mass /= 1e9
        | pos /= 1000
        v
ABYSS IC File (1e-9 Msun, kpc, km/s)
        |
        | readData() -> normalize_particle()
        | mass *= 1e9 -> Msun -> /MASS_UNIT
        | pos *= 1000 -> pc -> /POSITION_UNIT
        | vel convert to pc/yr -> /VELOCITY_UNIT
        v
ABYSS Code Units
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Manual IC file creation | McLuster integration | Phase 27 | Automated workflow |
| No unit conversion | Astrophysical units | Phase 27 | **BUG: Wrong units** |

**Required change:**
- transformMclusterOutput() must apply position/1000 and mass/1e9 conversions

## Open Questions

Things that couldn't be fully resolved:

1. **Output file placement**
   - What we know: McLuster outputs to working directory, transformation writes to `mcluster_abyss.dat`
   - What's unclear: Should this be configurable? Should it match the config file Filename?
   - Recommendation: Keep current behavior for now, could add config option later

2. **Precision requirements**
   - What we know: Scientific notation with 8 digits used currently
   - What's unclear: Is 8 digits sufficient for all use cases? Very small positions?
   - Recommendation: Increase to 15 digits to ensure no precision loss for edge cases

3. **McLuster N-body units (-u 0) alternative**
   - What we know: N-body units have total mass = 1, specific scaling
   - What's unclear: Exact scaling factors McLuster uses internally
   - Recommendation: Stick with astrophysical units for clarity, apply conversion

## Sources

### Primary (HIGH confidence)
- `/gpfs/home/vjl4366/pkg/ABYSS/src/particle.h` lines 266-275 - normalize_particle() implementation
- `/gpfs/home/vjl4366/pkg/ABYSS/src/mcluster_runner.cpp` lines 254-311 - current transformation
- `/gpfs/home/vjl4366/pkg/ABYSS/src/def.h` lines 50-76 - unit constant definitions
- `/gpfs/home/vjl4366/pkg/ABYSS/tests/test1/nbody.dat` - actual IC file format verification

### Secondary (MEDIUM confidence)
- `/gpfs/home/vjl4366/pkg/ABYSS/mcluster/README` - McLuster parameter documentation

## Metadata

**Confidence breakdown:**
- Unit mismatch identification: HIGH - Verified by code inspection and test file analysis
- Conversion factors: HIGH - Derived from normalize_particle() source code
- Fix approach: HIGH - Simple arithmetic, well-understood pipeline
- Edge cases: MEDIUM - May need precision adjustments for very small clusters

**Research date:** 2026-01-20
**Valid until:** Until normalize_particle() or readData() changes

---

## Requirements Mapping

| Requirement | Finding | Confidence |
|-------------|---------|------------|
| OUTPUT-01: McLuster output format matches ABYSS nbody.dat expectations | **CRITICAL FIX NEEDED**: Unit conversion missing | HIGH |
| OUTPUT-02: Units conversion if needed (N-body vs astrophysical) | **CRITICAL FIX NEEDED**: Position and mass need conversion | HIGH |
| OUTPUT-03: Generated IC file placed in correct location for ABYSS | Already working - fname updated to point to generated file | HIGH |

## Implementation Summary

The fix is straightforward:

1. In `transformMclusterOutput()`, after parsing McLuster values:
   - Divide position components (x, y, z) by 1000.0 (pc -> kpc)
   - Divide mass by 1e9 (Msun -> 1e-9 Msun units)
   - Leave velocity unchanged (km/s -> km/s)

2. Increase precision to `std::setprecision(15)` for safety

3. Add unit documentation comments to make the pipeline clear

This is a minimal, surgical fix that corrects the unit mismatch without changing any other behavior.
