# Phase 26: Config Parser Extension - Research

**Researched:** 2026-01-20
**Domain:** TOML configuration parsing, parameter validation, McLuster integration
**Confidence:** HIGH

## Summary

This phase extends the existing ABYSS TOML parser (`read_parameter_file.cpp` with `toml.hpp`) to recognize and validate a new `[mcluster]` configuration section. The existing parser infrastructure is well-designed with a `Config` class providing typed accessors, validation helpers, and nested table support - exactly what's needed for the new section.

McLuster accepts 9 parameters relevant to ABYSS integration: N (star count), M (total mass), P (density profile), R (half-mass radius), f (IMF selection), Z (metallicity), b (binary fraction), e (stellar evolution epoch), plus a new ABYSS-specific `generate_only` flag. All parameters have well-defined types, ranges, and defaults from McLuster's source code.

**Primary recommendation:** Extend `read_parameter_file.cpp` with a new `McusterConfig` struct and parsing function following the exact patterns already established for `[numerics]`, `[output]`, and `[restart]` sections.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `toml.hpp` | Custom (ABYSS) | TOML parsing | Already integrated, lightweight header-only |
| Standard C++ | C++11+ | String processing, data structures | No external dependencies |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `<algorithm>` | Standard | std::min for Levenshtein distance | Typo detection |
| `<vector>` | Standard | Storing valid parameter names | Parameter suggestion |
| `<sstream>` | Standard | Error message formatting | Validation messages |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Custom toml.hpp | toml11 | More features but adds dependency; custom is sufficient |
| Hand-rolled Levenshtein | External library | Overkill for 10 parameter names; simple implementation is <30 lines |

**Installation:**
```bash
# No additional installation needed - all components already in codebase
```

## Architecture Patterns

### Recommended Project Structure
```
src/
├── read_parameter_file.cpp  # Extend with McusterConfig parsing
├── global.h                 # Add mcluster config extern declarations
├── toml.hpp                 # No changes needed - already supports nested tables
└── mcluster_config.h        # NEW: MclusterConfig struct definition
```

### Pattern 1: Config Struct with Validation
**What:** Dedicated struct holding all mcluster parameters with sensible defaults
**When to use:** Complex configuration with many related parameters
**Example:**
```cpp
// Source: Pattern from existing read_parameter_file.cpp Config class
struct MclusterConfig {
    // Required: user must specify N or M (not both)
    int N = 0;              // Number of stars (0 = use M instead)
    double M = 0.0;         // Total mass in Msun (0 = use N instead)

    // Optional with McLuster defaults
    int P = 0;              // Density profile: 0=Plummer, 1=King, 2=Subr, 3=EFF/Nuker, -1=none
    double R = 0.8;         // Half-mass radius [pc]
    int f = 1;              // IMF: 0=single mass, 1=Kroupa, 2=user-defined
    double Z = 0.02;        // Metallicity [0.0001-0.03, solar=0.02]
    double b = 0.0;         // Binary fraction [0.0-1.0]
    double e = 0.0;         // Stellar evolution epoch [Myr]
    bool generate_only = false;  // ABYSS-specific: exit after IC generation

    bool has_mcluster_section = false;  // Track if section was present
};
```

### Pattern 2: Two-Phase Validation
**What:** Basic type/range checks at parse time, semantic checks at runtime
**When to use:** When some validation requires runtime context
**Example:**
```cpp
// Source: Existing pattern in read_parameter_file.cpp
// Phase 1: Parse time - check types and basic ranges
void parseMclusterSection(const Config& config, MclusterConfig& mcluster) {
    mcluster.N = config.getNestedOr<int>("mcluster", "N", 0);
    if (mcluster.N < 0) {
        throw std::runtime_error("Parameter 'N' must be non-negative, got: " +
            std::to_string(mcluster.N));
    }
    // ... more parsing with immediate range checks
}

// Phase 2: Runtime - check semantic constraints
void validateMclusterConfig(const MclusterConfig& mcluster) {
    if (mcluster.N == 0 && mcluster.M == 0.0) {
        throw std::runtime_error("Must specify either N (star count) or M (total mass)");
    }
    if (mcluster.N > 0 && mcluster.M > 0.0) {
        // Warn but don't error - M takes precedence per CONTEXT.md
        if (my_rank == ROOT) {
            std::cerr << "Warning: Both N and M specified; M takes precedence, ignoring N\n";
        }
    }
}
```

### Pattern 3: Unknown Parameter Detection with Typo Suggestion
**What:** Iterate table keys, flag unknowns, suggest similar valid names
**When to use:** User-facing config with potential for typos
**Example:**
```cpp
// Source: Standard Levenshtein distance algorithm
int levenshteinDistance(const std::string& s1, const std::string& s2) {
    std::vector<std::vector<int>> dp(s1.size() + 1, std::vector<int>(s2.size() + 1));
    for (size_t i = 0; i <= s1.size(); ++i) dp[i][0] = i;
    for (size_t j = 0; j <= s2.size(); ++j) dp[0][j] = j;
    for (size_t i = 1; i <= s1.size(); ++i) {
        for (size_t j = 1; j <= s2.size(); ++j) {
            int cost = (s1[i-1] == s2[j-1]) ? 0 : 1;
            dp[i][j] = std::min({dp[i-1][j] + 1, dp[i][j-1] + 1, dp[i-1][j-1] + cost});
        }
    }
    return dp[s1.size()][s2.size()];
}

std::string suggestSimilar(const std::string& unknown, const std::vector<std::string>& valid) {
    int minDist = INT_MAX;
    std::string suggestion;
    for (const auto& v : valid) {
        int dist = levenshteinDistance(unknown, v);
        if (dist < minDist && dist <= 2) {  // Threshold: max 2 edits
            minDist = dist;
            suggestion = v;
        }
    }
    return suggestion;
}
```

### Anti-Patterns to Avoid
- **Silent default fallback:** Don't silently use defaults for unknown parameters - fail fast with helpful message
- **Stringly-typed config:** Don't pass config values as strings - use typed struct
- **Global state pollution:** Don't add 10 new global variables - use single config struct
- **Late validation:** Don't wait until McLuster subprocess fails to catch bad parameters

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| TOML parsing | Custom parser | Existing `toml.hpp` | Already handles edge cases (escapes, nested tables) |
| Nested table access | Manual string parsing | `Config::getNestedOr` | Already implemented with error handling |
| Type conversion | Manual atoi/atof | `toml::find<T>` templates | Handles int/double/bool with proper exceptions |
| Validation errors | ad-hoc strings | `std::runtime_error` with consistent format | Matches existing error pattern |

**Key insight:** The existing `read_parameter_file.cpp` already solves all the hard parsing problems. This phase is pure extension, not innovation.

## Common Pitfalls

### Pitfall 1: N vs M Mutual Exclusivity
**What goes wrong:** User specifies both N=1000 and M=5000, unclear which is used
**Why it happens:** Natural to want to specify both for "documentation"
**How to avoid:** Check for both at validation time, warn and prefer M per CONTEXT.md decision
**Warning signs:** Unexpected star counts in generated ICs

### Pitfall 2: Silent Unknown Parameters
**What goes wrong:** User types `n=1000` (lowercase) instead of `N=1000`, gets default
**Why it happens:** TOML is case-sensitive, McLuster uses uppercase single letters
**How to avoid:** Iterate all keys in `[mcluster]` table, error on any not in known list
**Warning signs:** "I set N but it didn't work"

### Pitfall 3: Integer vs Float Confusion
**What goes wrong:** User writes `N = 1000.0` (float), parser fails or truncates
**Why it happens:** TOML distinguishes integers from floats strictly
**How to avoid:** The existing `toml::find<int>` already handles this (truncates whole numbers)
**Warning signs:** Type mismatch exceptions

### Pitfall 4: Metallicity Range
**What goes wrong:** User specifies Z=0.05 (above solar), McLuster rejects silently or produces garbage
**Why it happens:** Valid range [0.0001-0.03] is not obvious
**How to avoid:** Validate at parse time with helpful message including valid range
**Warning signs:** McLuster crashes or produces unrealistic stellar populations

### Pitfall 5: Missing Section Detection
**What goes wrong:** Old config files without `[mcluster]` section cause errors
**Why it happens:** Parser assumes section exists
**How to avoid:** Use `config.hasTable("mcluster")` to check presence first
**Warning signs:** Regression in existing workflows

## Code Examples

Verified patterns from official sources:

### Reading Nested Table (Existing Pattern)
```cpp
// Source: read_parameter_file.cpp lines 142-144
eta = config.getDoubleOr("eta", 0.01);
if (config.hasTable("numerics")) {
    eta = config.getNestedOr<double>("numerics", "eta", eta);
}
```

### Validation with Range and Message (Existing Pattern)
```cpp
// Source: read_parameter_file.cpp lines 118-123
void validateRange(int value, int min, int max, const std::string& name) {
    if (value < min || value > max) {
        throw std::runtime_error("Parameter '" + name + "' must be between " +
            std::to_string(min) + " and " + std::to_string(max) + ", got: " + std::to_string(value));
    }
}
```

### Checking Table Existence (Existing Pattern)
```cpp
// Source: read_parameter_file.cpp lines 98-105
bool hasTable(const std::string& table) const {
    try {
        toml::find(data_, table);
        return true;
    } catch (...) {
        return false;
    }
}
```

### Complete Mcluster Section Example Config
```toml
# Example config with [mcluster] section
Filename = "generated_ic.dat"
StopTime = 1.0e7
OutputDirectory = "output"

[mcluster]
N = 10000           # Number of stars (OR use M instead)
# M = 5000.0        # Total mass in Msun (alternative to N)
P = 0               # Density profile: 0=Plummer
R = 0.8             # Half-mass radius in pc
f = 1               # IMF: 1=Kroupa (2001)
Z = 0.02            # Metallicity (solar)
b = 0.0             # Binary fraction
e = 0.0             # Stellar evolution epoch [Myr]
generate_only = false  # If true, exit after IC generation

[numerics]
eta = 0.01
FixNumNeighbor = 100
```

## McLuster Parameter Reference

Parameters extracted from `mcluster/main.c`:

| Parameter | TOML Key | Type | McLuster Default | Valid Range | Notes |
|-----------|----------|------|------------------|-------------|-------|
| N | `N` | int | 0 (use M) | 3 to ~10^6 | Number of stars |
| M | `M` | double | 1000.0 | > 0 | Total mass [Msun] |
| P | `P` | int | 0 | -1, 0, 1, 2, 3 | Density profile |
| R | `R` | double | 0.8 | > 0 or -1 | Half-mass radius [pc] |
| f | `f` | int | 1 | 0, 1, 2 | IMF selection |
| Z | `Z` | double | 0.02 | 0.0001-0.03 | Metallicity |
| b | `b` | double | 0.0 | 0.0-1.0 | Binary fraction |
| e | `e` | double | 0.0 | >= 0 | Epoch [Myr] |
| generate_only | `generate_only` | bool | false | true/false | ABYSS-specific |

### Profile Values (P parameter)
- P = -1: No density gradient
- P = 0: Plummer model (default)
- P = 1: King model (requires W0 parameter, not exposed)
- P = 2: Subr et al. (2007) mass-segregated
- P = 3: 2D EFF/Nuker template

### IMF Values (f parameter)
- f = 0: Single-mass stars
- f = 1: Kroupa (2001) IMF (default)
- f = 2: User-defined multi power-law

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Global variables for each param | Config struct | v1.0+ | Better encapsulation |
| Silently ignore unknowns | Error with suggestions | Best practice | Catches typos |
| Exit on first validation error | Continue with meaningful message | Best practice | Better UX |

**Deprecated/outdated:**
- None relevant - existing patterns are current

## Range Validation Strategy

Per CONTEXT.md, "Claude decides which ranges are worth checking early":

**Check at parse time (immediate feedback):**
- N: Must be non-negative (N < 0 is always wrong)
- M: Must be non-negative (M < 0 is always wrong)
- P: Must be -1, 0, 1, 2, or 3 (small enumeration)
- Z: Must be in [0.0001, 0.03] (McLuster hard limit)
- b: Must be in [0.0, 1.0] (fraction)

**Defer to runtime/McLuster:**
- R: Complex rules (can be -1 for auto-calculate, or positive)
- f: Valid values depend on IMF definitions
- e: Any non-negative value valid
- N upper limit: Depends on available memory

**Rationale:** Early validation for values that are always wrong; defer complex validation where McLuster itself will provide better error messages.

## Open Questions

Things that couldn't be fully resolved:

1. **Line number reporting in errors**
   - What we know: CONTEXT.md requests line numbers in errors
   - What's unclear: Current toml.hpp doesn't track line numbers in parsed values
   - Recommendation: Enhance error messages with section name (e.g., "[mcluster] section: Parameter 'N' must be positive"). Line-level tracking would require toml.hpp modifications - defer to future enhancement if users request.

2. **Unknown parameter iteration**
   - What we know: Need to detect unknown keys like "nstar" vs "N"
   - What's unclear: Current toml.hpp doesn't expose key iteration API
   - Recommendation: Add `keys()` method to toml::value class (simple extension) or check for specific known unknown patterns

## Sources

### Primary (HIGH confidence)
- `/gpfs/home/vjl4366/pkg/ABYSS/src/read_parameter_file.cpp` - Existing TOML parsing patterns
- `/gpfs/home/vjl4366/pkg/ABYSS/src/toml.hpp` - Custom TOML parser implementation
- `/gpfs/home/vjl4366/pkg/ABYSS/mcluster/main.c` lines 61-147 - McLuster parameter definitions with defaults and ranges
- `/gpfs/home/vjl4366/pkg/ABYSS/mcluster/README` - McLuster parameter documentation

### Secondary (MEDIUM confidence)
- [Levenshtein Distance Wikipedia](https://en.wikipedia.org/wiki/Levenshtein_distance) - Algorithm for typo detection
- [TOML validation best practices](https://realpython.com/python-toml/) - Error message patterns

### Tertiary (LOW confidence)
- WebSearch results on config validation patterns - General guidance only

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Using existing codebase patterns
- Architecture: HIGH - Direct extension of existing Config class
- Parameter ranges: HIGH - Verified from McLuster source code
- Typo detection: MEDIUM - Standard algorithm, implementation straightforward

**Research date:** 2026-01-20
**Valid until:** Indefinite (stable domain, no external dependencies changing)
