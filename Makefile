# =============================================================================
# ABYSS Top-Level Makefile
# =============================================================================
# Builds ABYSS and optionally McLuster (if gfortran available)
#
# Targets:
#   all              - Build ABYSS + McLuster (default)
#   abyss            - Build only ABYSS
#   mcluster         - Build only McLuster
#   mcluster-clean   - Clean McLuster artifacts
#   mcluster-rebuild - Clean and rebuild McLuster
#   clean            - Clean all build artifacts
#
# Environment variables:
#   DISABLE_MCLUSTER - Set to 1 to skip mcluster even if gfortran available
#   QUIET            - Set to 1 to suppress build output
# =============================================================================

# Detect gfortran availability (immediate evaluation with :=)
GFORTRAN := $(shell which gfortran 2>/dev/null)

# Determine if mcluster can/should be built
ifndef GFORTRAN
  ifndef DISABLE_MCLUSTER
    $(warning gfortran not found, skipping mcluster build)
  endif
  BUILD_MCLUSTER := 0
else
  ifdef DISABLE_MCLUSTER
    BUILD_MCLUSTER := 0
  else
    BUILD_MCLUSTER := 1
  endif
endif

# Quiet mode configuration
ifdef QUIET
  Q := @
  MAKE_QUIET := --no-print-directory
else
  Q :=
  MAKE_QUIET :=
endif

# =============================================================================
# Default target - builds everything available
# =============================================================================
ifeq ($(BUILD_MCLUSTER),1)
all: abyss mcluster
else
all: abyss
endif

# =============================================================================
# ABYSS build
# =============================================================================
abyss:
	$(Q)$(MAKE) $(MAKE_QUIET) -C src

# =============================================================================
# McLuster build (conditional implementation)
# =============================================================================
ifeq ($(BUILD_MCLUSTER),1)
mcluster:
	$(Q)$(MAKE) $(MAKE_QUIET) -C mcluster mcluster_sse
	$(Q)ln -sf ../mcluster/mcluster_sse src/mcluster
	@echo "McLuster built successfully"
else
mcluster:
	@echo "mcluster: gfortran not available or DISABLE_MCLUSTER=1"
	@exit 1
endif

# =============================================================================
# McLuster maintenance targets
# =============================================================================
mcluster-clean:
	$(Q)$(MAKE) $(MAKE_QUIET) -C mcluster clean 2>/dev/null || true
	$(Q)rm -f src/mcluster

mcluster-rebuild: mcluster-clean mcluster

# =============================================================================
# Clean all artifacts
# =============================================================================
clean: mcluster-clean
	$(Q)$(MAKE) $(MAKE_QUIET) -C src clean

# =============================================================================
# Phony targets
# =============================================================================
.PHONY: all abyss mcluster mcluster-clean mcluster-rebuild clean
