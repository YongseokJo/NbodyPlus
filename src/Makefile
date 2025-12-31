# Compiler
#CXX = icc
#CXX = mpiicpx
CXX = mpicxx
#CXX = g++

# Compiler flags
CXXFLAGS  = -O2 -march=native -ftree-vectorize -std=c++11 -Wall -Wextra -Wuninitialized -g -DPERFORMANCETRACE #-fopenmp #-DDEBUG #-DMULTIMAP #-DDEBUG #-fopenmp  -fsanitize=thread
#CXXFLAGS  = -O0 -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g -qopenmp  #-fopenmp  -fsanitize=thread
#CXXFLAGS  = -O0 -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g -fopenmp  -fsanitize=thread
#CXXFLAGS  = -O2 -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g -fopenmp  #-fsanitize=thread

#NVCC_FLAGS_PROFILE = -g -lineinfo

# Libraries
#LDFLAGS = -L/usr/local/cuda-11.8/lib64/ -lcudart -lcublas -lnvToolsExt -lm -lcuda -lmpi #node14
# NOTE: Do not add "-lmpi" manually; it can accidentally link against a different MPI (e.g. OpenMPI) than the compiler wrapper.
LDFLAGS = -lm
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

# CUDA (overridable; typically set by `module add cuda`)
CUDA_HOME ?= /usr/local/cuda

ifdef USE_CUDA
	CXXFLAGS += -I./cuda -I$(CUDA_HOME)/include -DCUDA -DNSIGHT
	NVCC = nvcc
	NVCCFLAGS = -O2 --std=c++11 -G -Xcompiler -rdynamic -g -lineinfo -Xptxas=-v -lcuda -lcudart -lcompute-sanitizer -I./../SDAR/src -I$(CUDA_HOME)/include
	NVCC_FLAGS_DEBUG = -g -G --device-debug -lineinfo
	LDFLAGS += -L$(CUDA_HOME)/lib64 -lcudart -lcublas -lnvToolsExt -lcuda
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
