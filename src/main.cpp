#ifdef SEVN
#include "sevn.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include "def.h"
#include "particle.h"
#include "global_state.h"
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

	bin_output_file = fopen("binary_output.txt", "w");
	fprintf(bin_output_file, "Starting ABYSS - Binary OUTPUT\n");
	fflush(bin_output_file);
	merger_output_file = fopen("merger_output.txt", "w");
	fprintf(merger_output_file, "Starting ABYSS - Merger OUTPUT\n");
	fflush(merger_output_file);
#ifdef SEVN
	sevn_output_file = fopen("SEVN_output.txt", "w");
	fprintf(sevn_output_file, "Starting ABYSS - SEVN OUTPUT\n");
	fflush(sevn_output_file);
#endif

	/* MPI Initialization */
	initializeMPI(argc, argv);

#ifdef CUDA
	if (my_rank == ROOT) {
		OpenDevice();
		cudaDeviceSynchronize();
	}
#endif

	/* Input options */
	Parser(argc, argv);
	readParameterFile();

	// Write Particles
	if (my_rank == ROOT && !readData())
		fprintf(stderr, "Read Data Failed!\n");
	

	if (my_rank == ROOT) {
		RootRoutines();
	} else {
		std::string filename = "worker_output_" + std::to_string(my_rank) + ".txt";
		worker_output_file = fopen(filename.c_str(), "w");
		fprintf(worker_output_file, "Starting ABYSS - WORKER OUTPUT\n");
		fflush(worker_output_file);
		
		WorkerRoutines();
	}

	// Finalize the window and MPI environment
	MPI_Win_free(&win);
	MPI_Win_free(&win2);
	MPI_Win_free(&win3);
	MPI_Win_free(&win4);

	MPI_Comm_free(&shared_comm);

	MPI_Type_free(&queue_type_mpi);
    MPI_Type_free(&iparticle_type_mpi);
    MPI_Type_free(&jparticle_type_mpi);

	MPI_Finalize();

	return 0;
}

