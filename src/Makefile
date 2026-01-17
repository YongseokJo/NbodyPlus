# =============================================================================
# ABYSS Makefile
# =============================================================================
# All paths and options can be configured externally via environment variables.
# See workflow/config.sh for cluster-specific settings.
#
# Environment variables:
#   CXX         - C++ compiler (default: mpicxx)
#   CXXFLAGS    - Additional compiler flags (appended to defaults)
#   LDFLAGS     - Additional linker flags (appended to defaults)
#   HDF5_DIR    - HDF5 installation prefix (or HDF5_ROOT)
#   CUDA_HOME   - CUDA installation prefix (or CUDA_ROOT)
#   USE_CUDA    - Set to 1 to enable CUDA support
#   USE_SEVN    - Set to 1 to enable SEVN support
#   SEVN_DIR    - SEVN installation prefix
# =============================================================================

# Compiler (override via environment)
CXX ?= mpicxx

# Base compiler flags
CXXFLAGS_BASE = -O2 -march=native -ftree-vectorize -std=c++11 -Wall -Wextra -Wuninitialized -g
CXXFLAGS_ALL  = $(CXXFLAGS_BASE) $(CXXFLAGS)

# HDF5 (required)
HDF5_PATH := $(or $(HDF5_ROOT),$(HDF5_DIR))
ifdef HDF5_PATH
  CXXFLAGS_ALL += -I$(HDF5_PATH)/include
  LDFLAGS_HDF5  = -L$(HDF5_PATH)/lib -L$(HDF5_PATH)/lib64 -lhdf5_cpp -lhdf5
else
  $(warning HDF5_DIR or HDF5_ROOT not set; HDF5 headers/libs may not be found)
  LDFLAGS_HDF5 = -lhdf5_cpp -lhdf5
endif

# Base linker flags
LDFLAGS_BASE = -lm -lmpi $(LDFLAGS_HDF5)
LDFLAGS_ALL  = $(LDFLAGS_BASE) $(LDFLAGS)

# Source files
CXX_SRCS = main.cpp default_global.cpp mpi_routines.cpp read_parameter_file.cpp worker_routines.cpp root_routines.cpp initialization_routines.cpp irregular_routines.cpp regular_routines.cpp timestep_routines.cpp parser.cpp read_write.cpp particle_data.cpp particle_data_mpi.cpp Particle/initialize.cpp Particle/update_particle.cpp Particle/compute_acceleration.cpp

# SDAR (always included)
CXXFLAGS_ALL += -I./../SDAR/src -DAR_TTL -DAR_SLOWDOWN_TREE -DAR_SLOWDOWN_TIMESCALE -DFEWBODY
CXX_SRCS += FewBody/fb_check.cpp FewBody/fb_initialization.cpp FewBody/fb_integration.cpp FewBody/fb_termination.cpp FewBody/group_acceleration.cpp

# CUDA support (set USE_CUDA=1 to enable)
ifdef USE_CUDA
  CUDA_PATH := $(or $(CUDA_HOME),$(CUDA_ROOT))
  ifndef CUDA_PATH
    $(error USE_CUDA=1 but CUDA_HOME or CUDA_ROOT not set)
  endif
  CXXFLAGS_ALL += -I./cuda -I$(CUDA_PATH)/include -DCUDA -DNSIGHT
  NVCC          = $(CUDA_PATH)/bin/nvcc
  NVCCFLAGS     = -O2 -std=c++11 -Xcompiler -std=c++11 -Xcompiler -rdynamic -g -lineinfo -Xptxas=-v -I./../SDAR/src -I$(CUDA_PATH)/include
  NVCC_FLAGS_DEBUG = -g -G --device-debug -lineinfo
  LDFLAGS_ALL  += -L$(CUDA_PATH)/lib64 -lcudart -lcublas -lnvToolsExt -lcuda
  CXX_SRCS     += cuda/calculate_regular_acceleration.cpp
  CU_SRCS       = cuda/cuda_acceleration.cu cuda/cuda_kernels.cu
endif

# SEVN support (set USE_SEVN=1 and SEVN_DIR to enable)
ifdef USE_SEVN
  ifndef SEVN_DIR
    $(error USE_SEVN=1 but SEVN_DIR not set)
  endif
  CXXFLAGS_ALL += -DSEVN -DSEVN_BINARY -I$(SEVN_DIR)/include/sevn
  LDFLAGS_ALL  += -L$(SEVN_DIR)/lib64/sevn -lsevn_lib_shared
  CXX_SRCS     += stellar_evolution.cpp
endif

# Object files
OBJS = $(CXX_SRCS:.cpp=.o)
ifdef USE_CUDA
  OBJS += $(CU_SRCS:.cu=.o)
endif

# Executable name
EXEC = abyss.exe

# =============================================================================
# Build rules
# =============================================================================

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS_ALL) $(OBJS) -o $@ $(LDFLAGS_ALL)

%.o: %.cpp
	$(CXX) $(CXXFLAGS_ALL) -c $< -o $@

ifdef USE_CUDA
%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@
endif

clean:
	find . -type f -name "*.o" -delete
	rm -f $(EXEC)

.PHONY: all clean
