#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <sys/stat.h>
#include <iomanip>
#include "global.h"
#include "Queue.h"

int getLineNumber();
//void write_out(std::ofstream& outputFile, const Particle* ptcl);
void write_out(std::ofstream& outputFile, const Particle* ptcl, const double *pos, const double *vel);
void write_out_group(std::ofstream& outputFile, const Particle* ptcl, const Particle* members, const double *pos, const double *vel);
void write_neighbor(std::ofstream& outputFile, const Particle* ptcl);
#ifdef SEVN
void initializeStellarEvolution();
#endif
const int NUM_COLUMNS = 7; // Define the number of columns
const int width = 18;

bool readData() {

	fprintf(stdout, "Opening %s ...\n", fname);
	std::ifstream inputFile(fname);

	if (!inputFile) {
		std::cerr << "Error: Could not open the file." << std::endl;
		return false;
	}

	NumberOfParticle = getLineNumber();
	assert(NumberOfParticle < MaxNumberOfParticle);
	NewPID = NumberOfParticle;
	LastParticleIndex = NumberOfParticle - 1;
	global_variable->LastParticleIndex = LastParticleIndex;

	
	double** data = new double*[NumberOfParticle];
	for (int i = 0; i < NumberOfParticle; ++i) {
		data[i] = new double[NUM_COLUMNS];
	}

	// Initialization
	for (int i = 0; i < NumberOfParticle; ++i) {
		for (int j = 0; j < NUM_COLUMNS; ++j) {
			data[i][j] = 0;
		}
	}

	int row = 0;
	std::string line;
	while (std::getline(inputFile, line) && row < NumberOfParticle) { // Read lines from the file
		std::istringstream iss(line); // Create a stringstream for each line

		double value;
		int col = 0;
		while (iss >> value && col < NUM_COLUMNS) { // Read values from the stringstream
			data[row][col] = value;
			++col;
		}
		particles_original[row].initialize(data[row],row);
		++row;
	}

	/*
	for (int col = 0; col < NUM_COLUMNS; ++col) {
		std::cout << "Column " << col + 1 << " values: ";
		for (int r = 0; r < row; ++r) {
			std::cout << data[col][r] << ' ';
		}
		std::cout << std::endl;
	}
	*/


	// Normalize particles
	std::cout << "Particle normalizing." << std::endl;
	for (int i=0; i<NumberOfParticle; i++) {
		particles[i].normalizeParticle();
	}
	inputFile.close();

#ifdef SEVN
	initializeStellarEvolution();
#endif

	// Deallocate memory
	for (int i = 0; i < NumberOfParticle; ++i) {
		delete[] data[i];
	}
	delete[] data;

	return true;
}


int getLineNumber() {
	std::ifstream inputFile(fname); // Open the file

	if (!inputFile) {
		std::cerr << "Error: Could not open the file." << std::endl;
		return 1;
	}

	int lineCount = 0;
	std::string line;
	while (std::getline(inputFile, line)) { // Read lines from the file
		lineCount++;
	}

	std::cout << "Number of lines in the file: " << lineCount << std::endl;

	inputFile.close(); // Close the file

	return lineCount;
}


bool WriteData() {
	return true;
}



// Function to create a directory

bool createDirectory(const std::string& path) {
	// Create a folder with permissions 0777 (full access for user, group, others)
	int status = mkdir(path.c_str(), 0777);

	if (status == 0) {
		std::cout << "Folder created successfully." << std::endl;
	} else {
		std::cerr << "Error creating folder." << std::endl;
		// You can use perror to print the error message for more details
		perror("mkdir");
	}
	return true;
}



int writeParticle(double current_time, int outputNum) {

    std::cout << "Data is being written..." << std::endl;
    std::string directoryPath = "output";

    // Create the directory or check if it already exists
    if (!createDirectory(directoryPath)) {
        // Handle the error if necessary
        return 1;
    }

    // Now let's save the outputs in a new directory

    // Construct the filename with the timestamp
    std::string filename = directoryPath + "/" + foutput + "_" + std::to_string(outputNum) + ".txt";
    //std::string nn_fname = directoryPath + "/neighbor/nn_" + std::to_string(outputNum) + ".txt";

    // Open a file for writing
    std::ofstream outputFile(filename);
    //std::ofstream output_nn(nn_fname);


    // Check if the file is opened successfully
    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file!" << std::endl;
        return 1;
    }

	outputFile << "Time = " << current_time*EnzoTimeStep*1e10/1e6 << " Myr\n"; //

	Queue queue = {GetTotalEnergy, -1, -1.0};
	for (int i=0; i< NumberOfWorker; i++)
		MPI_Send(&queue, 1, QueueType, i+1, QUEUE_TAG, MPI_COMM_WORLD);

	MPI_Reduce(MPI_IN_PLACE, &E_binary,		1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
	MPI_Reduce(MPI_IN_PLACE, &E_binary_SD,	1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
	MPI_Reduce(MPI_IN_PLACE, &E_merger,		1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
	MPI_Reduce(MPI_IN_PLACE, &E_PN,			1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);

	outputFile << "NumberOfParticle = " << NumberOfParticle << "\n";
	double unit_energy = mass_unit * (velocity_unit/yr*pc/1e5) * (velocity_unit/yr*pc/1e5);
	outputFile << "E_binary = " << E_binary * unit_energy << "\n";
	outputFile << "E_binary_SD = " << E_binary_SD * unit_energy << "\n";
	outputFile << "E_merger = " << E_merger * unit_energy << "\n";
	outputFile << "E_PN = " << E_PN * unit_energy << "\n";

	E_binary = 0.0;
	E_binary_SD = 0.0;
	E_merger = 0.0;
	E_PN = 0.0;

    outputFile << std::left 
			<< std::setw(width) << "PID"
			<< std::setw(width) << "Mass (Msun)"
			<< std::setw(width) << "X (pc)"
			<< std::setw(width) << "Y (pc)"
			<< std::setw(width) << "Z (pc)"
			<< std::setw(width) << "Vx (km/s)"
		 	<< std::setw(width) << "Vy (km/s)" 
#ifdef SEVN
			<< std::setw(width) << "Vz (km/s)"
			<< std::setw(width) << "Type" << "\n";
#else
			<< std::setw(width) << "Vz (km/s)" << "\n";
#endif 


    // Write particle data to the file
	Particle *ptcl;
	double pos[Dim], vel[Dim];

	// for performance test by EW 2025.3.13
	int Index_minIrr, Index_minReg;
	double minTimeStepIrr = 1.;
	double minTimeStepReg = 1.;

	std::unordered_set<int> CMPtclsSet;

	for (int i=0; i<=LastParticleIndex; i++) {
		ptcl = &particles[i];

		if (!ptcl->isActive) continue;

		if (ptcl->TimeStepIrr < minTimeStepIrr) {
			Index_minIrr = ptcl->ParticleIndex;
			minTimeStepIrr = ptcl->TimeStepIrr;
		}
		if (ptcl->TimeStepReg < minTimeStepReg) {
			Index_minReg = ptcl->ParticleIndex;
			minTimeStepReg = ptcl->TimeStepReg;
		}

		if (ptcl->isCMptcl)
			CMPtclsSet.insert(i);
		
		ptcl->predictParticleSecondOrder(current_time - ptcl->CurrentTimeIrr, pos, vel);
		write_out(outputFile, ptcl, pos, vel);
		
// write_neighbor(output_nn, ptcl);
	}

	for (int i: CMPtclsSet) {
		ptcl = &particles[i];

		ptcl->predictParticleSecondOrder(current_time - ptcl->CurrentTimeIrr, pos, vel);

		Particle* members;
		for (int j=0; j < ptcl->NumberOfMember; j++) {
			members = &particles[ptcl->Members[j]];
			write_out_group(outputFile, ptcl, members, pos, vel);
		}
	}

	// Close the file
	outputFile.close();
	// output_nn.close();

	std::cout << "Data written to output.txt successfully!" << std::endl;

#ifdef PERFORMANCETRACE
	if (outputNum != 0) {
		std::cout << "--------------Performance-Summary--------------" << std::endl;
		fprintf(stdout, "Simulation Time: %f Myr\n", current_time*EnzoTimeStep*1e10/1e6);

		Particle* members = &particles[Index_minIrr];
		fprintf(stdout, "Particle Info with minimum TimeStepIrr...\n");
		members->printParticleInfo(stdout);

		members = &particles[Index_minReg];
		fprintf(stdout, "Particle Info with minimum TimeStepReg...\n");
		members->printParticleInfo(stdout);

		std::cout << std::fixed << std::setprecision(2);

		std::cout << "Elapsed time from the last output: " << performance.WholeRoutine*1e-9 << " s" << std::endl;

		std::cout << "Irregular Force: " << performance.IrregularForce*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.IrregularForce / performance.WholeRoutine << " %)" << std::endl;
		performance.IrregularForce = 0;
		std::cout << "Irregular Update: " << performance.IrregularUpdate*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.IrregularUpdate / performance.WholeRoutine << " %)" << std::endl;
		performance.IrregularUpdate = 0;

		std::cout << "FewBody Termination: " << performance.FewBodyTermination*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.FewBodyTermination / performance.WholeRoutine << " %)" << std::endl;
		performance.FewBodyTermination = 0;
		std::cout << "FewBody Search: " << performance.FewBodySearch*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.FewBodySearch / performance.WholeRoutine << " %)" << std::endl;
		performance.FewBodySearch = 0;
		std::cout << "FewBody Initialization: " << performance.FewBodyInitialization*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.FewBodyInitialization / performance.WholeRoutine << " %)" << std::endl;
		performance.FewBodyInitialization = 0;
#ifdef CUDA
		std::cout << "Regular SendToGPU: " << performance.RegularSendAllParticlesToGPU*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.RegularSendAllParticlesToGPU / performance.WholeRoutine << " %)" << std::endl;
		performance.RegularSendAllParticlesToGPU = 0;
		std::cout << "Regular GPU: " << performance.RegularGPU*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.RegularGPU / performance.WholeRoutine << " %)" << std::endl;
		performance.RegularGPU = 0;
		std::cout << "Regular Adjust: " << performance.RegularAdjust*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.RegularAdjust / performance.WholeRoutine << " %)" << std::endl;
		performance.RegularAdjust = 0;
		std::cout << "Regular Update: " << performance.RegularUpdate*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.RegularUpdate / performance.WholeRoutine << " %)" << std::endl;
		performance.RegularUpdate = 0;
#else
		std::cout << "Regular Force: " << performance.RegularForce*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.RegularForce / performance.WholeRoutine << " %)" << std::endl;
		performance.RegularForce = 0;
		std::cout << "Regular Update: " << performance.RegularUpdate*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.RegularUpdate / performance.WholeRoutine << " %)" << std::endl;
		performance.RegularUpdate = 0;
#endif
		std::cout << "SkipList Create: " << performance.SkipListCreate*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.SkipListCreate / performance.WholeRoutine << " %)" << std::endl;
		performance.SkipListCreate = 0;
		std::cout << "SkipList Update: " << performance.SkipListUpdate*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.SkipListUpdate / performance.WholeRoutine << " %)" << std::endl;
		performance.SkipListUpdate = 0;
#ifdef MULTIMAP
		std::cout << "RegularMap: " << performance.RegularMap*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.RegularMap / performance.WholeRoutine << " %)" << std::endl;
		performance.RegularMap = 0;
#else
		std::cout << "UpdateNextRegTime: " << performance.UpdateNextRegTime*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.UpdateNextRegTime / performance.WholeRoutine << " %)" << std::endl;
		performance.UpdateNextRegTime = 0;
#endif

#ifdef SEVN
		std::cout << "Stellar Evolution: " << performance.StellarEvolution*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.StellarEvolution / performance.WholeRoutine << " %)" << std::endl;
		performance.StellarEvolution = 0;
#ifdef SEVN_BINARY
		std::cout << "Binary Stellar Evolution: " << performance.BinaryStellarEvolution*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.BinaryStellarEvolution / performance.WholeRoutine << " %)" << std::endl;
		performance.BinaryStellarEvolution = 0;
#endif
#endif
		performance.WholeRoutine = 0;
		std::cout << "-----------------------------------------------" << std::endl;

		std::cout.unsetf(std::ios::fixed | std::ios::scientific);
	}
#endif

	return 0;

}


void write_out(std::ofstream& outputFile, const Particle* ptcl, const double *pos, const double *vel) {
        outputFile  << std::left << std::fixed << std::setprecision(8) // Eunwoo test
					<< std::setw(width) << ptcl->PID
					<< std::setw(width) << ptcl->Mass*mass_unit
                    << std::setw(width) << pos[0]*position_unit
                    << std::setw(width) << pos[1]*position_unit
                    << std::setw(width) << pos[2]*position_unit
                    << std::setw(width) << vel[0]*velocity_unit/yr*pc/1e5
                    << std::setw(width) << vel[1]*velocity_unit/yr*pc/1e5;
#ifdef SEVN
		outputFile	<< std::setw(width) << vel[2]*velocity_unit/yr*pc/1e5;
		if (ptcl->StellarEvolution == nullptr)
			outputFile << std::setw(width) << "1" << '\n'; // Main-sequence star;
		else if (!ptcl->StellarEvolution->amiremnant())
			outputFile << std::setw(width) << int(ptcl->StellarEvolution->getp(Phase::ID)) << '\n';
		else
			outputFile << std::setw(width) << 8+int(ptcl->StellarEvolution->getp(RemnantType::ID)) << '\n';
#else
		outputFile	<< std::setw(width) << vel[2]*velocity_unit/yr*pc/1e5 << '\n';
#endif
}

// This function is for group members cause group members have pos, vel in original frame, not predicted values.
void write_out_group(std::ofstream& outputFile, const Particle* ptclCM, const Particle* ptcl, const double *pos, const double *vel) {
        // outputFile  << std::left
		outputFile  << std::left << std::fixed << std::setprecision(8) // Eunwoo test
					<< std::setw(width) << ptcl->PID
					<< std::setw(width) << ptcl->Mass*mass_unit
                    << std::setw(width) << (pos[0] - ptclCM->Position[0] + ptcl->Position[0])*position_unit
                    << std::setw(width) << (pos[1] - ptclCM->Position[1] + ptcl->Position[1])*position_unit
                    << std::setw(width) << (pos[2] - ptclCM->Position[2] + ptcl->Position[2])*position_unit
                    << std::setw(width) << (vel[0] - ptclCM->Velocity[0] + ptcl->Velocity[0])*velocity_unit/yr*pc/1e5
                    << std::setw(width) << (vel[1] - ptclCM->Velocity[1] + ptcl->Velocity[1])*velocity_unit/yr*pc/1e5;
#ifdef SEVN
		outputFile	<< std::setw(width) << (vel[2] - ptclCM->Velocity[2] + ptcl->Velocity[2])*velocity_unit/yr*pc/1e5;
		if (ptcl->StellarEvolution == nullptr)
			outputFile << std::setw(width) << "1" << '\n';
		else if (!ptcl->StellarEvolution->amiremnant())
			outputFile << std::setw(width) << int(ptcl->StellarEvolution->getp(Phase::ID)) << '\n';
		else
			outputFile << std::setw(width) << 8+int(ptcl->StellarEvolution->getp(RemnantType::ID)) << '\n';
#else
		outputFile << std::setw(width) << (vel[2] - ptclCM->Velocity[2] + ptcl->Velocity[2])*velocity_unit/yr*pc/1e5 << '\n';
#endif
}

/*
void write_neighbor(std::ofstream& outputFile, const Particle* ptcl) {
	outputFile  << std::left\
			<< std::setw(width) << ptcl->PID << " = [ " ;
	for (Particle* nn:ptcl->ACList) {
			outputFile << nn->PID << "  ";
	}
	outputFile << "]\n";

}
*/


#ifdef time_trace
void output_time_trace() {

}
#endif
