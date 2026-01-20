# Phase 25: Build System Integration - Research

**Researched:** 2026-01-20
**Domain:** GNU Make build system, C/Fortran mixed compilation
**Confidence:** HIGH

## Summary

This phase integrates McLuster compilation into the ABYSS build system. The existing ABYSS build uses a straightforward GNU Makefile in `src/Makefile` with pattern-based conditional compilation (e.g., `ifdef USE_CUDA`). McLuster has its own `mcluster/Makefile` that compiles Fortran sources with gfortran and links them with a C main program using gcc.

The standard approach is to create a top-level Makefile at the project root that orchestrates both builds. This Makefile will detect gfortran availability using `$(shell which gfortran)`, skip McLuster if unavailable (with warning), and provide explicit targets for McLuster operations. The existing `src/Makefile` remains unchanged.

**Primary recommendation:** Create a root-level `Makefile` that uses recursive make to build both ABYSS (via `src/Makefile`) and McLuster (via `mcluster/Makefile`), with shell-based gfortran detection and graceful degradation when Fortran compiler is unavailable.

## Standard Stack

The established tools for this domain:

### Core
| Tool | Version | Purpose | Why Standard |
|------|---------|---------|--------------|
| GNU Make | 3.81+ | Build orchestration | Already used by ABYSS, universal on HPC systems |
| gfortran | 4.8+ | Fortran compiler | Part of GCC suite, required for McLuster SSE/BSE |
| gcc | 4.8+ | C compiler | Part of GCC suite, required for McLuster main.c |

### Supporting
| Tool | Version | Purpose | When to Use |
|------|---------|---------|-------------|
| OpenMP | 2.0+ | Parallelization | Already enabled in McLuster's gcc flags (-fopenmp) |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Recursive make | Single Makefile with includes | More complex, would require refactoring existing Makefiles |
| Make | CMake | Overkill for this use case, ABYSS already uses Make |
| Shell detection | Hardcoded paths | Shell detection is more portable across HPC systems |

**No new installation required** - all tools are already present on the target HPC system.

## Architecture Patterns

### Recommended Project Structure
```
ABYSS/
|-- Makefile              # NEW: Top-level orchestration
|-- src/
|   |-- Makefile          # UNCHANGED: Existing ABYSS build
|   `-- abyss.exe         # Output binary
|-- mcluster/
|   |-- Makefile          # UNCHANGED: Existing McLuster build
|   `-- mcluster_sse      # Output binary (renamed to mcluster via symlink)
`-- bin/                  # OPTIONAL: Could symlink binaries here
```

### Pattern 1: Shell-Based Compiler Detection
**What:** Use `$(shell which gfortran 2>/dev/null)` to detect Fortran compiler availability
**When to use:** At Makefile parse time, before any build rules execute
**Example:**
```makefile
# Source: Standard GNU Make pattern
GFORTRAN := $(shell which gfortran 2>/dev/null)
```

### Pattern 2: Conditional Warning with Skip
**What:** Use `$(warning ...)` to alert user, then conditionally exclude targets
**When to use:** When optional component cannot be built
**Example:**
```makefile
# Source: GNU Make manual + existing ABYSS pattern
ifndef GFORTRAN
  ifndef DISABLE_MCLUSTER
    $(warning gfortran not found, skipping mcluster)
  endif
  MCLUSTER_AVAILABLE := 0
else
  ifdef DISABLE_MCLUSTER
    MCLUSTER_AVAILABLE := 0
  else
    MCLUSTER_AVAILABLE := 1
  endif
endif
```

### Pattern 3: Recursive Make with $(MAKE)
**What:** Use `$(MAKE)` variable (not bare `make`) to invoke sub-makes
**When to use:** Always when calling make recursively - preserves flags like -j for parallel builds
**Example:**
```makefile
# Source: GNU Make manual section 5.7
mcluster:
	$(MAKE) -C mcluster mcluster_sse
```

### Pattern 4: Quiet Mode with @ Prefix
**What:** Use `@` prefix to suppress command echo, controlled by variable
**When to use:** When user requests quiet output (QUIET=1)
**Example:**
```makefile
# Source: Standard GNU Make pattern
ifdef QUIET
  Q := @
else
  Q :=
endif

mcluster:
	$(Q)$(MAKE) -C mcluster mcluster_sse
```

### Pattern 5: Symlink Creation
**What:** Create symlink in target directory pointing to actual binary
**When to use:** After successful build to provide simplified binary name
**Example:**
```makefile
# After mcluster build succeeds
ln -sf ../mcluster/mcluster_sse src/mcluster
```

### Anti-Patterns to Avoid
- **Modifying src/Makefile:** Keep existing Makefile unchanged to avoid breaking current workflows
- **Using bare `make` in recipes:** Always use `$(MAKE)` to preserve parallel build flags
- **Hardcoding compiler paths:** Use `which` detection for portability
- **Silent failures:** Always warn user when skipping optional components
- **Modifying mcluster/Makefile:** Keep original intact as it's from external project

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Compiler detection | Custom path search | `$(shell which compiler)` | Standard pattern, handles PATH correctly |
| Recursive directory build | Manual cd commands | `$(MAKE) -C dir target` | Preserves make flags, handles errors |
| Conditional targets | Complex if/else in recipes | ifeq/ifdef at Makefile level | Cleaner, evaluated at parse time |
| Symlink handling | cp command | ln -sf | Preserves single source of truth |

**Key insight:** GNU Make has built-in features for all these use cases. Using shell functions and make conditionals is more portable and maintainable than custom scripts.

## Common Pitfalls

### Pitfall 1: Forgetting $(MAKE) in Recursive Calls
**What goes wrong:** Parallel build flags (-j) not propagated, causing serial builds
**Why it happens:** Using bare `make` instead of `$(MAKE)`
**How to avoid:** Always use `$(MAKE) -C subdir target`
**Warning signs:** Sub-builds run serially even with -j flag

### Pitfall 2: Using `=` Instead of `:=` for Shell Commands
**What goes wrong:** Shell command runs multiple times (once per use of variable)
**Why it happens:** `=` is lazy evaluation, `:=` is immediate
**How to avoid:** Use `:=` for shell results: `GFORTRAN := $(shell which gfortran)`
**Warning signs:** Slow Makefile parsing, repeated shell output

### Pitfall 3: Not Handling Missing Directory Gracefully
**What goes wrong:** Build fails with confusing error if mcluster/ doesn't exist
**Why it happens:** git submodule not initialized, or mcluster not yet added
**How to avoid:** Check for directory existence before attempting build
**Warning signs:** "No rule to make target" errors

### Pitfall 4: Symlink Pointing to Wrong Location
**What goes wrong:** Symlink breaks when running from different directory
**Why it happens:** Using absolute path that differs between systems
**How to avoid:** Use relative path from symlink location: `ln -sf ../mcluster/mcluster_sse src/mcluster`
**Warning signs:** "No such file or directory" when running mcluster

### Pitfall 5: Clean Target Not Cleaning Symlinks
**What goes wrong:** Old symlink remains after clean, pointing to nothing
**Why it happens:** `rm -f *.o` doesn't match symlinks in other directories
**How to avoid:** Explicitly remove symlinks in clean target
**Warning signs:** Broken symlinks after `make clean`

### Pitfall 6: Build Failure Not Stopping Top-Level Make
**What goes wrong:** Top-level make reports success even when mcluster fails
**Why it happens:** Using `-` prefix before recipe command (ignores errors)
**How to avoid:** Don't use `-` prefix; let errors propagate naturally
**Warning signs:** "McLuster built successfully" message even when build failed

## Code Examples

Verified patterns from research and existing ABYSS codebase:

### Root-Level Makefile Structure
```makefile
# Source: Pattern derived from GNU Make manual + existing ABYSS src/Makefile patterns
# =============================================================================
# ABYSS Top-Level Makefile
# =============================================================================
# Builds ABYSS and optionally McLuster
#
# Environment variables:
#   DISABLE_MCLUSTER - Set to 1 to skip mcluster even if gfortran available
#   QUIET           - Set to 1 to suppress build output
# =============================================================================

# Detect gfortran availability (immediate evaluation)
GFORTRAN := $(shell which gfortran 2>/dev/null)

# Determine if mcluster can/should be built
ifndef GFORTRAN
  ifndef DISABLE_MCLUSTER
    $(warning gfortran not found, skipping mcluster)
  endif
  BUILD_MCLUSTER := 0
else
  ifdef DISABLE_MCLUSTER
    BUILD_MCLUSTER := 0
  else
    BUILD_MCLUSTER := 1
  endif
endif

# Quiet mode
ifdef QUIET
  Q := @
  MAKE_QUIET := --no-print-directory
else
  Q :=
  MAKE_QUIET :=
endif

# Default target builds everything available
ifeq ($(BUILD_MCLUSTER),1)
all: abyss mcluster
else
all: abyss
endif

# Build ABYSS
abyss:
	$(Q)$(MAKE) $(MAKE_QUIET) -C src

# Build McLuster (only if available)
ifeq ($(BUILD_MCLUSTER),1)
mcluster:
	$(Q)$(MAKE) $(MAKE_QUIET) -C mcluster mcluster_sse
	$(Q)ln -sf ../mcluster/mcluster_sse src/mcluster
	$(Q)echo "McLuster built successfully"
else
mcluster:
	$(Q)echo "mcluster: gfortran not available or DISABLE_MCLUSTER=1"
	$(Q)exit 1
endif

# McLuster-specific targets
mcluster-clean:
	$(Q)$(MAKE) $(MAKE_QUIET) -C mcluster clean
	$(Q)rm -f src/mcluster

mcluster-rebuild: mcluster-clean mcluster

# Clean everything
clean: mcluster-clean
	$(Q)$(MAKE) $(MAKE_QUIET) -C src clean

.PHONY: all abyss mcluster mcluster-clean mcluster-rebuild clean
```

### Compiler Detection Pattern
```makefile
# Source: GNU Make manual, verified working on target HPC system
# Check for gfortran - immediate evaluation prevents repeated shell calls
GFORTRAN := $(shell which gfortran 2>/dev/null)
GCC := $(shell which gcc 2>/dev/null)

# Conditional based on detection
ifdef GFORTRAN
  $(info Fortran compiler: $(GFORTRAN))
else
  $(warning No Fortran compiler found)
endif
```

### Existing mcluster/Makefile Key Lines
```makefile
# Source: mcluster/Makefile (existing, verified working 2026-01-20)
FC = gfortran -O2
CC = gcc -O2 -fopenmp -Wall
CFLAGS = -L/usr/lib/ -lgfortran

# Fortran sources for SSE/BSE
SOURCE = deltat.f evolv1.f hrdiag.f kick.f mlwind.f mrenv.f \
         ran3.f star.f zcnsts.f zfuncs.f comenv.f corerd.f \
         dgcore.f evolv2.f gntage.f instar.f mix.f rl.f

OBJECTS = $(SOURCE:.f=.o)

mcluster_sse: $(OBJECTS) $(LFLAGS)
	$(CC) -c main.c -D SSE -lm
	$(CC) $(OBJECTS) main.o -o mcluster_sse -lm $(CFLAGS)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Manual builds | Integrated make | This phase | Single `make` builds everything |
| Separate binaries | Symlinked binaries | This phase | Simplified invocation |

**Deprecated/outdated:**
- Building mcluster manually in its directory: Will still work, but `make` from root is preferred

## Open Questions

Things that couldn't be fully resolved:

1. **gfortran library path (-L/usr/lib/ -lgfortran)**
   - What we know: Current mcluster Makefile hardcodes `/usr/lib/`
   - What's unclear: Whether this path works on all target systems
   - Recommendation: Keep existing path for now; if issues arise, use `$(shell gfortran -print-file-name=libgfortran.so | xargs dirname)` to detect dynamically

2. **OpenMP support detection**
   - What we know: mcluster uses -fopenmp flag
   - What's unclear: Whether OpenMP is always available with gfortran
   - Recommendation: Leave as-is; gfortran includes OpenMP support by default

## Sources

### Primary (HIGH confidence)
- `/gpfs/home/vjl4366/pkg/ABYSS/src/Makefile` - Existing ABYSS build patterns
- `/gpfs/home/vjl4366/pkg/ABYSS/mcluster/Makefile` - Existing McLuster build (verified working)
- `/gpfs/home/vjl4366/pkg/ABYSS/mcluster/README` - McLuster compilation instructions
- Verified build test on HPC system (2026-01-20)

### Secondary (MEDIUM confidence)
- [GNU Make Manual - Recursion](https://www.gnu.org/software/make/manual/html_node/Recursion.html) - Recursive make patterns
- [Makefile Tutorial by Example](https://makefiletutorial.com/) - Conditional and shell patterns
- [Writing Makefiles for Modern Fortran](https://aoterodelaroza.github.io/devnotes/modern-fortran-makefiles/) - Fortran build patterns
- [GNU Mixed-Language Programming](https://gcc.gnu.org/onlinedocs/gfortran/Mixed-Language-Programming.html) - C/Fortran linking

### Tertiary (LOW confidence)
- Various Stack Overflow and forum posts on compiler detection - patterns verified against official docs

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Using existing tools already on system
- Architecture: HIGH - Patterns derived from existing ABYSS Makefile and verified McLuster build
- Pitfalls: HIGH - Based on direct testing and GNU Make documentation

**Research date:** 2026-01-20
**Valid until:** 90 days (stable technology, no expected changes)
