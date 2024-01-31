# Compiler
#CXX = icc
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11 -Wall -Wextra -Wuninitialized -I./Particle -I./

# Source files
SRCS = main.cpp ReadWrite.cpp Particle/ParticleRoutines.cpp ComputationChainRoutine.cpp SortComputationChain.cpp InitializeParticles.cpp TimeStepRoutines.cpp Particle/CalculateIrrForce.cpp IrregularAccelerationRoutine.cpp Particle/CalculateRegForce.cpp RegularAccelerationRoutine.cpp Particle/UpdateParticle.cpp Particle/CalculateTimeStep.cpp InitializeTimeStep.cpp Evolve.cpp Particle/CheckNeighborForEvolution.cpp Particle/UpdateEvolveParticle.cpp parser.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
EXEC = nbodyplus.exe

# Main target
all: $(EXEC)

# Link object files into the executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(EXEC) -v

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean the project directory
clean:
	rm -f $(OBJS) $(EXEC)


