# McLuster Integration Tests
#
# Quick tests use pre-generated fixtures (no gfortran needed)
# Full tests run McLuster live (requires gfortran)
#
# Usage:
#   make test        # Run quick tests (default)
#   make test-quick  # Run quick tests (fixture-based)
#   make test-full   # Run full tests (live McLuster)
#   make FAIL_FAST=1 test  # Stop on first failure

.PHONY: test test-quick test-full clean

SHELL := /bin/bash

# Default: quick tests
test: test-quick

# Quick tests with pre-generated fixtures
test-quick:
	@echo "=== Running McLuster Integration Tests (fixture mode) ==="
	@FIXTURE_MODE=1 ./run_tests.sh

# Full tests with live McLuster execution
test-full:
	@echo "=== Running McLuster Integration Tests (live mode) ==="
	@if ! which gfortran >/dev/null 2>&1; then \
		echo "ERROR: gfortran required for full tests"; \
		exit 1; \
	fi
	@FIXTURE_MODE=0 ./run_tests.sh

# Individual test targets for debugging
test-generate-only:
	@./test_generate_only.sh

test-run-only:
	@./test_run_only.sh

test-generate-run:
	@./test_generate_run.sh

test-energy:
	@./test_energy.sh

# Clean test output
clean:
	rm -rf output/
	rm -rf /tmp/abyss_test_*

# Create output directory
output:
	mkdir -p output
