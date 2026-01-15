#ifdef SEVN
#include "sevn.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include "def.h"
#include "particle.h"
#include "GlobalVariable.h"
#include "global.h"
#include <mpi.h>
#ifdef CUDA
#include <cuda_runtime.h>
#include "cuda/cuda_functions.h"
#endif


void broadcastFromRoot(int &data);
void DefaultGlobal();
void initializeMPI(int argc, char *argv[]);
void WorkerRoutines();
void RootRoutines();
void Parser(int argc, char *argv[]);
bool readData();
void readParameterFile();


int main(int argc, char *argv[]) {

	/* Initialize global variables */
	DefaultGlobal();

	binout = fopen("binary_output.txt", "w");
	fprintf(binout, "Starting ABYSS - Binary OUTPUT\n");
	fflush(binout);
	mergerout = fopen("merger_output.txt", "w");
	fprintf(mergerout, "Starting ABYSS - Merger OUTPUT\n");
	fflush(mergerout);
#ifdef SEVN
	SEVNout = fopen("SEVN_output.txt", "w");
	fprintf(SEVNout, "Starting ABYSS - SEVN OUTPUT\n");
	fflush(SEVNout);
#endif

	/* MPI Initialization */
	initializeMPI(argc, argv);

#ifdef CUDA
	if (MyRank == ROOT) {
		OpenDevice();
		cudaDeviceSynchronize();
	}
#endif

	/* Input options */
	Parser(argc, argv);
	readParameterFile();

	// Write Particles
	// If shared windows are not actually shared (SharedCommSize==1), each rank must
	// load the initial condition data to avoid workers reading uninitialized memory.
	if (SharedCommSize == 1) {
		if (!readData()) {
			fprintf(stderr, "Read Data Failed!\n");
		}
		MPI_Barrier(MPI_COMM_WORLD);
	} else {
		if (MyRank == ROOT && !readData())
			fprintf(stderr, "Read Data Failed!\n");
		MPI_Barrier(MPI_COMM_WORLD);
	}
	

	if (MyRank == ROOT) {
		RootRoutines();
	} else {
		std::string filename = "worker_output_" + std::to_string(MyRank) + ".txt";
		workerout = fopen(filename.c_str(), "w");
		fprintf(workerout, "Starting ABYSS - WORKER OUTPUT\n");
		fflush(workerout);
		
		WorkerRoutines();
	}

	// Finalize the window and MPI environment
	MPI_Win_free(&win);
	MPI_Win_free(&win2);
	MPI_Win_free(&win3);
	MPI_Win_free(&win4);

	MPI_Comm_free(&shared_comm);

	MPI_Type_free(&QueueType);
    MPI_Type_free(&IparticleType);
    MPI_Type_free(&JparticleType);

	MPI_Finalize();

	return 0;
}

