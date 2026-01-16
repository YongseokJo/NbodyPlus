# Compiler
#CXX = icc
#CXX = mpiicpc
CXX = mpicxx
#CXX = g++

# Compiler flags
CXXFLAGS  = -O2 -march=native -ftree-vectorize -std=c++11 -Wall -Wextra -Wuninitialized -g -DPERFORMANCETRACE #-fopenmp #-DDEBUG #-DMULTIMAP #-DDEBUG #-fopenmp  -fsanitize=thread
#CXXFLAGS  = -O0 -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g -qopenmp  #-fopenmp  -fsanitize=thread
#CXXFLAGS  = -O0 -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g -fopenmp  -fsanitize=thread
#CXXFLAGS  = -O2 -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g -fopenmp  #-fsanitize=thread

# TOML parsing uses the in-tree src/toml.hpp (no external include needed).

# HDF5 paths (adjust as needed for your system)
# Set HDF5_DIR environment variable or load HDF5 module before compiling
# Example: module load hdf5/1.14.1-2-openmpi-gcc-12.3.0
ifdef HDF5_ROOT
	HDF5_PATH = $(HDF5_ROOT)
else ifdef HDF5_DIR
	HDF5_PATH = $(HDF5_DIR)
else
	HDF5_PATH = /usr
endif
CXXFLAGS += -I$(HDF5_PATH)/include
LDFLAGS_HDF5 = -L$(HDF5_PATH)/lib -L$(HDF5_PATH)/lib64 -lhdf5_cpp -lhdf5

#NVCC_FLAGS_PROFILE = -g -lineinfo

# Libraries
#LDFLAGS = -L/usr/local/cuda-11.8/lib64/ -lcudart -lcublas -lnvToolsExt -lm -lcuda -lmpi #node14
LDFLAGS = -lm -lmpi $(LDFLAGS_HDF5)
# LDFLAGS = -L/usr/local/cuda/lib64/ -lm -lmpi #node15
#add -lnvToolsExt for <nvToolsExt.h>
#LDFLAGS = -L/mnt/sw/nix/store/zi2wc26znf75csf5hhz77p0d2bbz53ih-cuda-11.8.0/lib64/ -lcudart

# Source files
CXX_SRCS = main.cpp DefaultGlobal.cpp MPIRoutines.cpp readParameterFile.cpp WorkerRoutines.cpp RootRoutines.cpp InitializationRoutines.cpp IrregularRoutines.cpp RegularRoutines.cpp TimeStepRoutines.cpp parser.cpp ReadWrite.cpp Particle/Initialize.cpp Particle/UpdateParticle.cpp Particle/ComputeAcceleration.cpp
#cuda/cuda_routines.cpp


# cuda/cuda_routines.cu


# SDAR
CXXFLAGS += -I./../SDAR/src -DAR_TTL -DAR_SLOWDOWN_TREE -DAR_SLOWDOWN_TIMESCALE -DFEWBODY
CXX_SRCS += FewBody/FBCheck.cpp FewBody/FBInitialization.cpp FewBody/FBIntegration.cpp FewBody/FBTermination.cpp FewBody/GroupAccelerationRoutine.cpp

ifdef USE_CUDA
	# CUDA paths - use CUDA_HOME if set, otherwise try common locations
	CUDA_PATH ?= $(if $(CUDA_HOME),$(CUDA_HOME),$(if $(CUDA_ROOT),$(CUDA_ROOT),/usr/local/cuda))
	CXXFLAGS += -I./cuda -I$(CUDA_PATH)/include -DCUDA -DNSIGHT
	NVCC = $(CUDA_PATH)/bin/nvcc
	NVCCFLAGS = -O2 -std=c++11 -Xcompiler -std=c++11 -Xcompiler -rdynamic -g -lineinfo -Xptxas=-v -I./../SDAR/src -I$(CUDA_PATH)/include
	NVCC_FLAGS_DEBUG = -g -G --device-debug -lineinfo
	LDFLAGS += -L$(CUDA_PATH)/lib64 -lcudart -lcublas -lnvToolsExt -lcuda
	CXX_SRCS += cuda/CalculateRegularAcceleration.cpp
	CU_SRCS = cuda/cuda_my_acceleration.cu cuda/cuda_kernels.cu
endif

ifdef USE_SEVN
    CXXFLAGS += -DSEVN -DSEVN_BINARY -I/home/vinicius/install/sevn_custom/build_n15/include/sevn
    LDFLAGS += -L/home/vinicius/install/sevn_custom/build_n15/lib64/sevn -lsevn_lib_shared
    CXX_SRCS += StellarEvolution.cpp
endif

#
OBJS = $(CXX_SRCS:.cpp=.o)
ifdef USE_CUDA
	OBJS += $(CU_SRCS:.cu=.o) 
endif

# Executable name
EXEC = abyss.exe

# Main target
all: $(EXEC)


# Link object files into the executable
$(EXEC): $(OBJS)
	#$(NVCC) $(NVCCFLAGS) $(OBJS) -o $@ $(LDFLAGS)
	#$(NVCC) $(CXXFLAGS) $(NVCCFLAGS) $(OBJS) -o $(EXEC) -v
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ -v $(LDFLAGS)
	#$(NVCC) $(CXXFLAGS) $(NVCCFLAGS) $(OBJS) -o $@ -v

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

ifdef USE_CUDA
%.o: %.cu
	$(NVCC) $(NVCCFLAGS) $(NVCC_FLAGS_DEBUG) -c $< -o $@
endif


# Clean the project directory
clean:
	find . -type f -name "*.o" -delete
	rm -f $(EXEC)
