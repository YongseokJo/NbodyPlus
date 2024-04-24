# Compiler
#CXX = icc
CXX = icc
NVCC = nvcc

# Compiler flags
CXXFLAGS  = -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./
NVCCFLAGS = -O2 -G -Xcompiler -rdynamic -lineinfo

# Libraries
LDFLAGS = -L/usr/local/cuda-11.8/lib64/ -lcudart
#LDFLAGS = -L/mnt/sw/nix/store/zi2wc26znf75csf5hhz77p0d2bbz53ih-cuda-11.8.0/lib64/ -lcudart

# Source files
CXX_SRCS = main.cpp ReadWrite.cpp Particle/ParticleRoutines.cpp ComputationChainRoutine.cpp SortComputationChain.cpp InitializeParticles.cpp TimeStepRoutines.cpp Particle/CalculateIrrForce.cpp IrregularAccelerationRoutine.cpp Particle/CalculateRegForce.cpp RegularAccelerationRoutine.cpp Particle/UpdateParticle.cpp Particle/CalculateTimeStep.cpp InitializeTimeStep.cpp Evolve.cpp Particle/CheckNeighborForEvolution.cpp Particle/UpdateEvolveParticle.cpp parser.cpp cuda/CalculateRegularAcceleration.cpp cuda/CalculateAccelerationForAll.cpp Particle/EvolveStar.cpp Binary/BinaryRoutines.cpp Binary/KSRegularlizationInitialization.cpp Binary/KSRegularlizationIntegration.cpp Binary/KSRegularlizationTermination.cpp BinaryAccelerationRoutine.cpp



CU_SRCS = cuda/cuda_my_acceleration.cu

#
OBJS = $(CU_SRCS:.cu=.o) $(CXX_SRCS:.cpp=.o)

# Executable name
EXEC = nbodyplus.exe

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

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@


# Clean the project directory
clean:
	rm -f $(OBJS) $(EXEC)


