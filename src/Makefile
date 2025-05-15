# Compiler
#CXX = icc
CXX = mpic++
#CXX = g++

# Compiler flags
CXXFLAGS  = -O2 -march=native -ftree-vectorize -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g -DPERFORMANCETRACE -DMultiNode #-DDEBUG #-DMULTIMAP #-DDEBUG #-fopenmp  -fsanitize=thread
#CXXFLAGS  = -O0 -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g -qopenmp  #-fopenmp  -fsanitize=thread
#CXXFLAGS  = -O0 -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g -fopenmp  -fsanitize=thread
#CXXFLAGS  = -O2 -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./ -g -fopenmp  #-fsanitize=thread

#NVCC_FLAGS_PROFILE = -g -lineinfo

# Libraries
#LDFLAGS = -L/usr/local/cuda-11.8/lib64/ -lcudart -lcublas -lnvToolsExt -lm -lcuda -lmpi #node14
LDFLAGS = -lm -lmpi
# LDFLAGS = -L/usr/local/cuda/lib64/ -lm -lmpi #node15
#add -lnvToolsExt for <nvToolsExt.h>
#LDFLAGS = -L/mnt/sw/nix/store/zi2wc26znf75csf5hhz77p0d2bbz53ih-cuda-11.8.0/lib64/ -lcudart

# Source files
CXX_SRCS = main.cpp DefaultGlobal.cpp MPIRoutines.cpp readParameterFile.cpp WorkerRoutines.cpp RootRoutines.cpp TimeStepRoutines.cpp parser.cpp ReadWrite.cpp Particle/Initialize.cpp Particle/UpdateParticle.cpp Particle/ComputeAcceleration.cpp
#cuda/cuda_routines.cpp


# cuda/cuda_routines.cu


# SDAR
CXXFLAGS += -I./../SDAR/src -DAR_TTL -DAR_SLOWDOWN_TREE -DAR_SLOWDOWN_TIMESCALE -DFEWBODY
CXX_SRCS += FewBody/FBCheck.cpp FewBody/FBInitialization.cpp FewBody/FBIntegration.cpp FewBody/FBTermination.cpp FewBody/GroupAccelerationRoutine.cpp

ifdef USE_CUDA
	CXXFLAGS += -I./cuda -DCUDA -DNSIGHT
	NVCC = nvcc
	NVCCFLAGS = -O2 --std=c++11 -G -Xcompiler -rdynamic -g -lineinfo -Xptxas=-v -lcuda -lcudart -lcompute-sanitizer
	NVCC_FLAGS_DEBUG = -g -G --device-debug -lineinfo
	LDFLAGS += -L/usr/local/cuda/lib64/ -lcudart -lcublas -lnvToolsExt -lcuda
	CXX_SRCS += cuda/CalculateRegularAcceleration.cpp
	CU_SRCS = cuda/cuda_my_acceleration.cu cuda/cuda_kernels.cu
endif

ifdef USE_SEVN
    CXXFLAGS += -DSEVN -I/data/vinicius/sevn/build/include/sevn
    LDFLAGS += -L/data/vinicius/sevn/build/lib64/sevn -lsevn_lib_shared
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
