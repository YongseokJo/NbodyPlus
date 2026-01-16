#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <sys/stat.h>
#include <iomanip>
#include "global.h"
#include "Queue.h"
#include "H5Cpp.h"

int getLineNumber();
void write_out(std::ofstream& outputFile, const Particle* ptcl, const double *pos, const double *vel);
void write_out_group(std::ofstream& outputFile, const Particle* ptcl, const Particle* members, const double *pos, const double *vel);
void write_neighbor(std::ofstream& outputFile, const Particle* ptcl);
int writeCheckpoint(double current_time, int checkpointNum);
bool readCheckpoint(const std::string& filename);
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
	assert(NumberOfParticle < MaxNumParticle);
	NewCMPID = NumberOfParticle;
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
		particles[row].initialize(data[row],row);
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

    // Construct the HDF5 filename
    std::string filename = directoryPath + "/" + foutput + ".h5";

	Queue queue = {GetTotalEnergy, -1, -1.0};
	for (int i=0; i< NumberOfWorker; i++)
		MPI_Send(&queue, 1, QueueType, i+1, QUEUE_TAG, MPI_COMM_WORLD);

	MPI_Reduce(MPI_IN_PLACE, &E_binary,		1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
	MPI_Reduce(MPI_IN_PLACE, &E_binary_SD,	1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
	MPI_Reduce(MPI_IN_PLACE, &E_merger,		1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
	MPI_Reduce(MPI_IN_PLACE, &E_PN,			1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);

	double unit_energy = mass_unit * (velocity_unit/yr*pc/1e5) * (velocity_unit/yr*pc/1e5);

    // Write particle data to HDF5
	Particle *ptcl;
	double pos[Dim], vel[Dim];

	// for performance test by EW 2025.3.13
	int Index_minIrr, Index_minReg;
	double minTimeStepIrr = 1.;
	double minTimeStepReg = 1.;

	std::unordered_set<int> CMPtclsSet;

	// Count active particles and collect data
	int activeCount = 0;
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

		activeCount++;
	}

	// Count group members
	for (int i: CMPtclsSet) {
		ptcl = &particles[i];
		activeCount += ptcl->NumberOfMember;
	}

	// Allocate arrays for particle data
	std::vector<int> pid_data(activeCount);
	std::vector<double> mass_data(activeCount);
	std::vector<double> x_data(activeCount);
	std::vector<double> y_data(activeCount);
	std::vector<double> z_data(activeCount);
	std::vector<double> vx_data(activeCount);
	std::vector<double> vy_data(activeCount);
	std::vector<double> vz_data(activeCount);
#ifdef SEVN
	std::vector<int> type_data(activeCount);
#endif

	// Collect particle data
	int idx = 0;
	for (int i=0; i<=LastParticleIndex; i++) {
		ptcl = &particles[i];
		if (!ptcl->isActive) continue;

		ptcl->predictParticleSecondOrder(current_time - ptcl->CurrentTimeIrr, pos, vel);

		pid_data[idx] = ptcl->PID;
		mass_data[idx] = ptcl->Mass * mass_unit;
		x_data[idx] = pos[0] * position_unit;
		y_data[idx] = pos[1] * position_unit;
		z_data[idx] = pos[2] * position_unit;
		vx_data[idx] = vel[0] * velocity_unit/yr*pc/1e5;
		vy_data[idx] = vel[1] * velocity_unit/yr*pc/1e5;
		vz_data[idx] = vel[2] * velocity_unit/yr*pc/1e5;
#ifdef SEVN
		if (ptcl->StellarEvolution == nullptr)
			type_data[idx] = 1;
		else if (!ptcl->StellarEvolution->amiremnant())
			type_data[idx] = int(ptcl->StellarEvolution->getp(Phase::ID));
		else
			type_data[idx] = 8+int(ptcl->StellarEvolution->getp(RemnantType::ID));
#endif
		idx++;
	}

	// Collect group member data
	for (int i: CMPtclsSet) {
		ptcl = &particles[i];
		ptcl->predictParticleSecondOrder(current_time - ptcl->CurrentTimeIrr, pos, vel);

		Particle* members;
		for (int j=0; j < ptcl->NumberOfMember; j++) {
			members = &particles[ptcl->Members[j]];

			pid_data[idx] = members->PID;
			mass_data[idx] = members->Mass * mass_unit;
			x_data[idx] = (pos[0] - ptcl->Position[0] + members->Position[0]) * position_unit;
			y_data[idx] = (pos[1] - ptcl->Position[1] + members->Position[1]) * position_unit;
			z_data[idx] = (pos[2] - ptcl->Position[2] + members->Position[2]) * position_unit;
			vx_data[idx] = (vel[0] - ptcl->Velocity[0] + members->Velocity[0]) * velocity_unit/yr*pc/1e5;
			vy_data[idx] = (vel[1] - ptcl->Velocity[1] + members->Velocity[1]) * velocity_unit/yr*pc/1e5;
			vz_data[idx] = (vel[2] - ptcl->Velocity[2] + members->Velocity[2]) * velocity_unit/yr*pc/1e5;
#ifdef SEVN
			if (members->StellarEvolution == nullptr)
				type_data[idx] = 1;
			else if (!members->StellarEvolution->amiremnant())
				type_data[idx] = int(members->StellarEvolution->getp(Phase::ID));
			else
				type_data[idx] = 8+int(members->StellarEvolution->getp(RemnantType::ID));
#endif
			idx++;
		}
	}

	// Write to HDF5 file
	try {
		H5::H5File file;

		// Open or create file
		if (outputNum == 0) {
			file = H5::H5File(filename, H5F_ACC_TRUNC);

			// Write root-level metadata (only on first output)
			hsize_t scalar_dim = 1;
			H5::DataSpace scalar_space(1, &scalar_dim);

			// Version info
			H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
			std::string version = "ABYSS 1.0";
			H5::Attribute version_attr = file.createAttribute("ABYSS_Version", str_type, scalar_space);
			version_attr.write(str_type, version);

			// Input filename
			std::string input_file = fname;
			H5::Attribute input_attr = file.createAttribute("InputFile", str_type, scalar_space);
			input_attr.write(str_type, input_file);

			// Units documentation as JSON string
			std::string units_json = "{"
				"\"Position\": \"pc\", "
				"\"Velocity\": \"km/s\", "
				"\"Mass\": \"Msun\", "
				"\"Time\": \"Myr\", "
				"\"Energy\": \"Msun*(km/s)^2\", "
				"\"CodeUnits\": {"
					"\"position_unit\": 4.0, "
					"\"time_unit\": 1e10, "
					"\"mass_unit\": 0.0001424198, "
					"\"velocity_unit\": 4e-10"
				"}"
			"}";
			H5::Attribute units_attr = file.createAttribute("Units", str_type, scalar_space);
			units_attr.write(str_type, units_json);

			// Simulation parameters
			H5::Attribute eta_attr = file.createAttribute("eta", H5::PredType::NATIVE_DOUBLE, scalar_space);
			eta_attr.write(H5::PredType::NATIVE_DOUBLE, &eta);

			H5::Attribute nneighbor_attr = file.createAttribute("FixNumNeighbor", H5::PredType::NATIVE_INT, scalar_space);
			nneighbor_attr.write(H5::PredType::NATIVE_INT, &FixNumNeighbor);

			double init_radius_pc = InitialNeighborRadius * position_unit;
			H5::Attribute radius_attr = file.createAttribute("InitialRadius_pc", H5::PredType::NATIVE_DOUBLE, scalar_space);
			radius_attr.write(H5::PredType::NATIVE_DOUBLE, &init_radius_pc);

			double rsearch_pc = RSearch * position_unit;
			H5::Attribute rsearch_attr = file.createAttribute("RSearch_pc", H5::PredType::NATIVE_DOUBLE, scalar_space);
			rsearch_attr.write(H5::PredType::NATIVE_DOUBLE, &rsearch_pc);

			double end_time_myr = endTime / 1e6;  // yr -> Myr
			H5::Attribute endtime_attr = file.createAttribute("EndTime_Myr", H5::PredType::NATIVE_DOUBLE, scalar_space);
			endtime_attr.write(H5::PredType::NATIVE_DOUBLE, &end_time_myr);

			H5::Attribute nproc_attr = file.createAttribute("NumberOfProcessors", H5::PredType::NATIVE_INT, scalar_space);
			nproc_attr.write(H5::PredType::NATIVE_INT, &NumberOfProcessor);

			int init_npart = NumberOfParticle;
			H5::Attribute initnpart_attr = file.createAttribute("InitialParticleCount", H5::PredType::NATIVE_INT, scalar_space);
			initnpart_attr.write(H5::PredType::NATIVE_INT, &init_npart);

			// Compression settings
			H5::Attribute compress_attr = file.createAttribute("CompressionEnabled", H5::PredType::NATIVE_HBOOL, scalar_space);
			hbool_t compress_val = UseCompression ? 1 : 0;
			compress_attr.write(H5::PredType::NATIVE_HBOOL, &compress_val);

			if (UseCompression) {
				H5::Attribute level_attr = file.createAttribute("CompressionLevel", H5::PredType::NATIVE_INT, scalar_space);
				level_attr.write(H5::PredType::NATIVE_INT, &CompressionLevel);
			}

		} else {
			file = H5::H5File(filename, H5F_ACC_RDWR);
		}

		// Create group for this timestep
		std::string groupName = "/Step_" + std::to_string(outputNum);
		H5::Group group = file.createGroup(groupName);

		// Write metadata
		double time_value = current_time*EnzoTimeStep*1e10/1e6;
		hsize_t scalar_dim = 1;
		H5::DataSpace scalar_space(1, &scalar_dim);

		H5::Attribute time_attr = group.createAttribute("Time_Myr", H5::PredType::NATIVE_DOUBLE, scalar_space);
		time_attr.write(H5::PredType::NATIVE_DOUBLE, &time_value);

		H5::Attribute npart_attr = group.createAttribute("NumberOfParticle", H5::PredType::NATIVE_INT, scalar_space);
		npart_attr.write(H5::PredType::NATIVE_INT, &activeCount);

		double e_bin = E_binary * unit_energy;
		double e_bin_sd = E_binary_SD * unit_energy;
		double e_mrg = E_merger * unit_energy;
		double e_pn = E_PN * unit_energy;

		H5::Attribute ebinary_attr = group.createAttribute("E_binary", H5::PredType::NATIVE_DOUBLE, scalar_space);
		ebinary_attr.write(H5::PredType::NATIVE_DOUBLE, &e_bin);

		H5::Attribute ebinary_sd_attr = group.createAttribute("E_binary_SD", H5::PredType::NATIVE_DOUBLE, scalar_space);
		ebinary_sd_attr.write(H5::PredType::NATIVE_DOUBLE, &e_bin_sd);

		H5::Attribute emerger_attr = group.createAttribute("E_merger", H5::PredType::NATIVE_DOUBLE, scalar_space);
		emerger_attr.write(H5::PredType::NATIVE_DOUBLE, &e_mrg);

		H5::Attribute epn_attr = group.createAttribute("E_PN", H5::PredType::NATIVE_DOUBLE, scalar_space);
		epn_attr.write(H5::PredType::NATIVE_DOUBLE, &e_pn);

		// Write datasets with optional compression and chunking
		hsize_t dims[1] = {(hsize_t)activeCount};
		H5::DataSpace dataspace(1, dims);

		// Create dataset creation property list with chunking and compression
		H5::DSetCreatPropList plist;
		if (UseCompression && activeCount > 0) {
			// Set chunk size (minimum of activeCount or 10000 for efficiency)
			hsize_t chunk_dims[1] = {std::min((hsize_t)activeCount, (hsize_t)10000)};
			plist.setChunk(1, chunk_dims);
			plist.setDeflate(CompressionLevel);  // GZIP compression
		}

		H5::DataSet pid_dataset = group.createDataSet("PID", H5::PredType::NATIVE_INT, dataspace, plist);
		pid_dataset.write(pid_data.data(), H5::PredType::NATIVE_INT);

		H5::DataSet mass_dataset = group.createDataSet("Mass_Msun", H5::PredType::NATIVE_DOUBLE, dataspace, plist);
		mass_dataset.write(mass_data.data(), H5::PredType::NATIVE_DOUBLE);

		H5::DataSet x_dataset = group.createDataSet("X_pc", H5::PredType::NATIVE_DOUBLE, dataspace, plist);
		x_dataset.write(x_data.data(), H5::PredType::NATIVE_DOUBLE);

		H5::DataSet y_dataset = group.createDataSet("Y_pc", H5::PredType::NATIVE_DOUBLE, dataspace, plist);
		y_dataset.write(y_data.data(), H5::PredType::NATIVE_DOUBLE);

		H5::DataSet z_dataset = group.createDataSet("Z_pc", H5::PredType::NATIVE_DOUBLE, dataspace, plist);
		z_dataset.write(z_data.data(), H5::PredType::NATIVE_DOUBLE);

		H5::DataSet vx_dataset = group.createDataSet("Vx_km_s", H5::PredType::NATIVE_DOUBLE, dataspace, plist);
		vx_dataset.write(vx_data.data(), H5::PredType::NATIVE_DOUBLE);

		H5::DataSet vy_dataset = group.createDataSet("Vy_km_s", H5::PredType::NATIVE_DOUBLE, dataspace, plist);
		vy_dataset.write(vy_data.data(), H5::PredType::NATIVE_DOUBLE);

		H5::DataSet vz_dataset = group.createDataSet("Vz_km_s", H5::PredType::NATIVE_DOUBLE, dataspace, plist);
		vz_dataset.write(vz_data.data(), H5::PredType::NATIVE_DOUBLE);

#ifdef SEVN
		H5::DataSet type_dataset = group.createDataSet("Type", H5::PredType::NATIVE_INT, dataspace, plist);
		type_dataset.write(type_data.data(), H5::PredType::NATIVE_INT);
#endif

		file.close();

	} catch (H5::Exception& e) {
		std::cerr << "HDF5 error: " << e.getDetailMsg() << std::endl;
		return 1;
	}

	E_binary = 0.0;
	E_binary_SD = 0.0;
	E_merger = 0.0;
	E_PN = 0.0;

	std::cout << "Data written to HDF5 file successfully!" << std::endl;

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
#ifdef unused
		std::cout << "FewBody Search: " << performance.FewBodySearch*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.FewBodySearch / performance.WholeRoutine << " %)" << std::endl;
		performance.FewBodySearch = 0;
#endif
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

		// New enhanced profiler output
		double sim_time_myr = current_time * EnzoTimeStep * 1e10 / 1e6;
		profiler().printIntervalSummary(std::cout, sim_time_myr);

		// Write CSV for analysis (appends to file)
		std::string csv_path = std::string(foutput) + "/profiling.csv";
		profiler().writeCSV(csv_path, sim_time_myr, outputNum);

		// Write JSON snapshot for this output
		std::string json_path = std::string(foutput) + "/profiling_" + std::to_string(outputNum) + ".json";
		profiler().writeJSON(json_path, sim_time_myr, outputNum);

		// Reset interval stats for next output period
		profiler().resetIntervalStats();
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


// ============================================================================
// Checkpoint/Restart Functions
// ============================================================================

int writeCheckpoint(double current_time, int checkpointNum) {
	if (!RestartEnabled) return 0;

	std::string directoryPath = foutput;
	createDirectory(directoryPath);

	std::string filename = directoryPath + "/checkpoint_" + std::to_string(checkpointNum) + ".h5";

	std::cout << "Writing checkpoint to " << filename << "..." << std::endl;

	try {
		H5::H5File file(filename, H5F_ACC_TRUNC);

		// Write global state as attributes
		hsize_t scalar_dim = 1;
		H5::DataSpace scalar_space(1, &scalar_dim);

		// Time variables
		H5::Attribute gt_attr = file.createAttribute("global_time", H5::PredType::NATIVE_DOUBLE, scalar_space);
		gt_attr.write(H5::PredType::NATIVE_DOUBLE, &global_time);

		H5::Attribute gti_attr = file.createAttribute("global_time_irr", H5::PredType::NATIVE_DOUBLE, scalar_space);
		gti_attr.write(H5::PredType::NATIVE_DOUBLE, &global_time_irr);

		H5::Attribute nrt_attr = file.createAttribute("NextRegTimeBlock", H5::PredType::NATIVE_ULLONG, scalar_space);
		nrt_attr.write(H5::PredType::NATIVE_ULLONG, &NextRegTimeBlock);

		H5::Attribute tb_attr = file.createAttribute("time_block", H5::PredType::NATIVE_INT, scalar_space);
		tb_attr.write(H5::PredType::NATIVE_INT, &time_block);

		H5::Attribute bm_attr = file.createAttribute("block_max", H5::PredType::NATIVE_ULLONG, scalar_space);
		bm_attr.write(H5::PredType::NATIVE_ULLONG, &block_max);

		// Particle counts
		H5::Attribute lpi_attr = file.createAttribute("LastParticleIndex", H5::PredType::NATIVE_INT, scalar_space);
		lpi_attr.write(H5::PredType::NATIVE_INT, &LastParticleIndex);

		H5::Attribute nop_attr = file.createAttribute("NumberOfParticle", H5::PredType::NATIVE_INT, scalar_space);
		nop_attr.write(H5::PredType::NATIVE_INT, &NumberOfParticle);

		H5::Attribute ncm_attr = file.createAttribute("NewCMPID", H5::PredType::NATIVE_INT, scalar_space);
		ncm_attr.write(H5::PredType::NATIVE_INT, &NewCMPID);

		// Output tracking
		H5::Attribute ot_attr = file.createAttribute("outputTime", H5::PredType::NATIVE_DOUBLE, scalar_space);
		ot_attr.write(H5::PredType::NATIVE_DOUBLE, &outputTime);

		H5::Attribute on_attr = file.createAttribute("outNum", H5::PredType::NATIVE_INT, scalar_space);
		on_attr.write(H5::PredType::NATIVE_INT, &outNum);

		// Count particles to save (active ones plus their data)
		int saveCount = 0;
		for (int i = 0; i <= LastParticleIndex; i++) {
			if (particles[i].PID >= 0) saveCount++;
		}

		// Create particle group
		H5::Group ptcl_group = file.createGroup("/Particles");

		hsize_t dims[1] = {(hsize_t)saveCount};
		H5::DataSpace dataspace(1, dims);

		// Allocate arrays
		std::vector<int> pid(saveCount), pindex(saveCount), ptype(saveCount);
		std::vector<int> nneighbor(saveCount), noffset(saveCount);
		std::vector<int> isactive(saveCount), iscm(saveCount), cmindex(saveCount);
		std::vector<int> nmember(saveCount);
		std::vector<double> mass(saveCount), radius_nb(saveCount);
		std::vector<double> px(saveCount), py(saveCount), pz(saveCount);
		std::vector<double> vx(saveCount), vy(saveCount), vz(saveCount);
		std::vector<double> ct_irr(saveCount), ct_reg(saveCount);
		std::vector<double> ts_irr(saveCount), ts_reg(saveCount);
		std::vector<ULL> cb_irr(saveCount), cb_reg(saveCount), nb_irr(saveCount);
		std::vector<ULL> tbl_irr(saveCount), tbl_reg(saveCount);
		std::vector<int> tl_irr(saveCount), tl_reg(saveCount);

		// Acceleration arrays (flattened)
		std::vector<double> a_tot_flat(saveCount * Dim * HERMITE_ORDER);
		std::vector<double> a_reg_flat(saveCount * Dim * HERMITE_ORDER);
		std::vector<double> a_irr_flat(saveCount * Dim * HERMITE_ORDER);

		// Collect data
		int idx = 0;
		for (int i = 0; i <= LastParticleIndex; i++) {
			Particle* ptcl = &particles[i];
			if (ptcl->PID < 0) continue;

			pid[idx] = ptcl->PID;
			pindex[idx] = ptcl->ParticleIndex;
			ptype[idx] = ptcl->ParticleType;
			mass[idx] = ptcl->Mass;
			px[idx] = ptcl->Position[0];
			py[idx] = ptcl->Position[1];
			pz[idx] = ptcl->Position[2];
			vx[idx] = ptcl->Velocity[0];
			vy[idx] = ptcl->Velocity[1];
			vz[idx] = ptcl->Velocity[2];
			nneighbor[idx] = ptcl->NumberOfNeighbor;
			noffset[idx] = ptcl->NeighborsOffset;
			radius_nb[idx] = ptcl->RadiusOfNeighbor;
			ct_irr[idx] = ptcl->CurrentTimeIrr;
			ct_reg[idx] = ptcl->CurrentTimeReg;
			cb_irr[idx] = ptcl->CurrentBlockIrr;
			cb_reg[idx] = ptcl->CurrentBlockReg;
			nb_irr[idx] = ptcl->NextBlockIrr;
			ts_irr[idx] = ptcl->TimeStepIrr;
			ts_reg[idx] = ptcl->TimeStepReg;
			tbl_irr[idx] = ptcl->TimeBlockIrr;
			tbl_reg[idx] = ptcl->TimeBlockReg;
			tl_irr[idx] = ptcl->TimeLevelIrr;
			tl_reg[idx] = ptcl->TimeLevelReg;
			isactive[idx] = ptcl->isActive ? 1 : 0;
			iscm[idx] = ptcl->isCMptcl ? 1 : 0;
			cmindex[idx] = ptcl->CMPtclIndex;
			nmember[idx] = ptcl->NumberOfMember;

			// Copy accelerations (flattened)
			for (int d = 0; d < Dim; d++) {
				for (int o = 0; o < HERMITE_ORDER; o++) {
					int flat_idx = idx * Dim * HERMITE_ORDER + d * HERMITE_ORDER + o;
					a_tot_flat[flat_idx] = ptcl->a_tot[d][o];
					a_reg_flat[flat_idx] = ptcl->a_reg[d][o];
					a_irr_flat[flat_idx] = ptcl->a_irr[d][o];
				}
			}
			idx++;
		}

		// Write datasets
		auto writeDataset = [&](const std::string& name, const H5::DataType& type, const void* data) {
			H5::DataSet ds = ptcl_group.createDataSet(name, type, dataspace);
			ds.write(data, type);
		};

		writeDataset("PID", H5::PredType::NATIVE_INT, pid.data());
		writeDataset("ParticleIndex", H5::PredType::NATIVE_INT, pindex.data());
		writeDataset("ParticleType", H5::PredType::NATIVE_INT, ptype.data());
		writeDataset("Mass", H5::PredType::NATIVE_DOUBLE, mass.data());
		writeDataset("Position_X", H5::PredType::NATIVE_DOUBLE, px.data());
		writeDataset("Position_Y", H5::PredType::NATIVE_DOUBLE, py.data());
		writeDataset("Position_Z", H5::PredType::NATIVE_DOUBLE, pz.data());
		writeDataset("Velocity_X", H5::PredType::NATIVE_DOUBLE, vx.data());
		writeDataset("Velocity_Y", H5::PredType::NATIVE_DOUBLE, vy.data());
		writeDataset("Velocity_Z", H5::PredType::NATIVE_DOUBLE, vz.data());
		writeDataset("NumberOfNeighbor", H5::PredType::NATIVE_INT, nneighbor.data());
		writeDataset("NeighborsOffset", H5::PredType::NATIVE_INT, noffset.data());
		writeDataset("RadiusOfNeighbor", H5::PredType::NATIVE_DOUBLE, radius_nb.data());
		writeDataset("CurrentTimeIrr", H5::PredType::NATIVE_DOUBLE, ct_irr.data());
		writeDataset("CurrentTimeReg", H5::PredType::NATIVE_DOUBLE, ct_reg.data());
		writeDataset("CurrentBlockIrr", H5::PredType::NATIVE_ULLONG, cb_irr.data());
		writeDataset("CurrentBlockReg", H5::PredType::NATIVE_ULLONG, cb_reg.data());
		writeDataset("NextBlockIrr", H5::PredType::NATIVE_ULLONG, nb_irr.data());
		writeDataset("TimeStepIrr", H5::PredType::NATIVE_DOUBLE, ts_irr.data());
		writeDataset("TimeStepReg", H5::PredType::NATIVE_DOUBLE, ts_reg.data());
		writeDataset("TimeBlockIrr", H5::PredType::NATIVE_ULLONG, tbl_irr.data());
		writeDataset("TimeBlockReg", H5::PredType::NATIVE_ULLONG, tbl_reg.data());
		writeDataset("TimeLevelIrr", H5::PredType::NATIVE_INT, tl_irr.data());
		writeDataset("TimeLevelReg", H5::PredType::NATIVE_INT, tl_reg.data());
		writeDataset("isActive", H5::PredType::NATIVE_INT, isactive.data());
		writeDataset("isCMptcl", H5::PredType::NATIVE_INT, iscm.data());
		writeDataset("CMPtclIndex", H5::PredType::NATIVE_INT, cmindex.data());
		writeDataset("NumberOfMember", H5::PredType::NATIVE_INT, nmember.data());

		// Write acceleration arrays (2D: saveCount x (Dim*HERMITE_ORDER))
		hsize_t acc_dims[2] = {(hsize_t)saveCount, (hsize_t)(Dim * HERMITE_ORDER)};
		H5::DataSpace acc_space(2, acc_dims);
		H5::DataSet atot_ds = ptcl_group.createDataSet("a_tot", H5::PredType::NATIVE_DOUBLE, acc_space);
		atot_ds.write(a_tot_flat.data(), H5::PredType::NATIVE_DOUBLE);
		H5::DataSet areg_ds = ptcl_group.createDataSet("a_reg", H5::PredType::NATIVE_DOUBLE, acc_space);
		areg_ds.write(a_reg_flat.data(), H5::PredType::NATIVE_DOUBLE);
		H5::DataSet airr_ds = ptcl_group.createDataSet("a_irr", H5::PredType::NATIVE_DOUBLE, acc_space);
		airr_ds.write(a_irr_flat.data(), H5::PredType::NATIVE_DOUBLE);

		// Save neighbor lists for particles that have neighbors
		H5::Group neighbor_group = file.createGroup("/Neighbors");
		for (int i = 0; i <= LastParticleIndex; i++) {
			Particle* ptcl = &particles[i];
			if (ptcl->PID < 0 || ptcl->NumberOfNeighbor <= 0) continue;

			std::string ds_name = std::to_string(ptcl->ParticleIndex);
			hsize_t nb_dims[1] = {(hsize_t)ptcl->NumberOfNeighbor};
			H5::DataSpace nb_space(1, nb_dims);
			H5::DataSet nb_ds = neighbor_group.createDataSet(ds_name, H5::PredType::NATIVE_INT, nb_space);
			nb_ds.write(&Neighbors[ptcl->NeighborsOffset], H5::PredType::NATIVE_INT);
		}

		file.close();
		std::cout << "Checkpoint written successfully: " << saveCount << " particles saved." << std::endl;

	} catch (H5::Exception& e) {
		std::cerr << "HDF5 error writing checkpoint: " << e.getDetailMsg() << std::endl;
		return 1;
	}

	return 0;
}


bool readCheckpoint(const std::string& filename) {
	std::cout << "Reading checkpoint from " << filename << "..." << std::endl;

	try {
		H5::H5File file(filename, H5F_ACC_RDONLY);

		// Read global state
		H5::Attribute gt_attr = file.openAttribute("global_time");
		gt_attr.read(H5::PredType::NATIVE_DOUBLE, &global_time);

		H5::Attribute gti_attr = file.openAttribute("global_time_irr");
		gti_attr.read(H5::PredType::NATIVE_DOUBLE, &global_time_irr);

		H5::Attribute nrt_attr = file.openAttribute("NextRegTimeBlock");
		nrt_attr.read(H5::PredType::NATIVE_ULLONG, &NextRegTimeBlock);

		H5::Attribute tb_attr = file.openAttribute("time_block");
		tb_attr.read(H5::PredType::NATIVE_INT, &time_block);

		H5::Attribute bm_attr = file.openAttribute("block_max");
		bm_attr.read(H5::PredType::NATIVE_ULLONG, &block_max);

		H5::Attribute lpi_attr = file.openAttribute("LastParticleIndex");
		lpi_attr.read(H5::PredType::NATIVE_INT, &LastParticleIndex);

		H5::Attribute nop_attr = file.openAttribute("NumberOfParticle");
		nop_attr.read(H5::PredType::NATIVE_INT, &NumberOfParticle);

		H5::Attribute ncm_attr = file.openAttribute("NewCMPID");
		ncm_attr.read(H5::PredType::NATIVE_INT, &NewCMPID);

		H5::Attribute ot_attr = file.openAttribute("outputTime");
		ot_attr.read(H5::PredType::NATIVE_DOUBLE, &outputTime);

		H5::Attribute on_attr = file.openAttribute("outNum");
		on_attr.read(H5::PredType::NATIVE_INT, &outNum);

		// Update global_variable
		global_variable->LastParticleIndex = LastParticleIndex;
		global_variable->NextRegTimeBlock = NextRegTimeBlock;
		time_step = std::pow(2, time_block);

		// Open particle group
		H5::Group ptcl_group = file.openGroup("/Particles");

		// Get dataset dimensions
		H5::DataSet pid_ds = ptcl_group.openDataSet("PID");
		H5::DataSpace space = pid_ds.getSpace();
		hsize_t dims[1];
		space.getSimpleExtentDims(dims);
		int saveCount = dims[0];

		// Allocate arrays
		std::vector<int> pid(saveCount), pindex(saveCount), ptype(saveCount);
		std::vector<int> nneighbor(saveCount), noffset(saveCount);
		std::vector<int> isactive(saveCount), iscm(saveCount), cmindex(saveCount);
		std::vector<int> nmember(saveCount);
		std::vector<double> mass(saveCount), radius_nb(saveCount);
		std::vector<double> px(saveCount), py(saveCount), pz(saveCount);
		std::vector<double> vx(saveCount), vy(saveCount), vz(saveCount);
		std::vector<double> ct_irr(saveCount), ct_reg(saveCount);
		std::vector<double> ts_irr(saveCount), ts_reg(saveCount);
		std::vector<ULL> cb_irr(saveCount), cb_reg(saveCount), nb_irr(saveCount);
		std::vector<ULL> tbl_irr(saveCount), tbl_reg(saveCount);
		std::vector<int> tl_irr(saveCount), tl_reg(saveCount);
		std::vector<double> a_tot_flat(saveCount * Dim * HERMITE_ORDER);
		std::vector<double> a_reg_flat(saveCount * Dim * HERMITE_ORDER);
		std::vector<double> a_irr_flat(saveCount * Dim * HERMITE_ORDER);

		// Read datasets
		auto readDataset = [&](const std::string& name, const H5::DataType& type, void* data) {
			H5::DataSet ds = ptcl_group.openDataSet(name);
			ds.read(data, type);
		};

		readDataset("PID", H5::PredType::NATIVE_INT, pid.data());
		readDataset("ParticleIndex", H5::PredType::NATIVE_INT, pindex.data());
		readDataset("ParticleType", H5::PredType::NATIVE_INT, ptype.data());
		readDataset("Mass", H5::PredType::NATIVE_DOUBLE, mass.data());
		readDataset("Position_X", H5::PredType::NATIVE_DOUBLE, px.data());
		readDataset("Position_Y", H5::PredType::NATIVE_DOUBLE, py.data());
		readDataset("Position_Z", H5::PredType::NATIVE_DOUBLE, pz.data());
		readDataset("Velocity_X", H5::PredType::NATIVE_DOUBLE, vx.data());
		readDataset("Velocity_Y", H5::PredType::NATIVE_DOUBLE, vy.data());
		readDataset("Velocity_Z", H5::PredType::NATIVE_DOUBLE, vz.data());
		readDataset("NumberOfNeighbor", H5::PredType::NATIVE_INT, nneighbor.data());
		readDataset("NeighborsOffset", H5::PredType::NATIVE_INT, noffset.data());
		readDataset("RadiusOfNeighbor", H5::PredType::NATIVE_DOUBLE, radius_nb.data());
		readDataset("CurrentTimeIrr", H5::PredType::NATIVE_DOUBLE, ct_irr.data());
		readDataset("CurrentTimeReg", H5::PredType::NATIVE_DOUBLE, ct_reg.data());
		readDataset("CurrentBlockIrr", H5::PredType::NATIVE_ULLONG, cb_irr.data());
		readDataset("CurrentBlockReg", H5::PredType::NATIVE_ULLONG, cb_reg.data());
		readDataset("NextBlockIrr", H5::PredType::NATIVE_ULLONG, nb_irr.data());
		readDataset("TimeStepIrr", H5::PredType::NATIVE_DOUBLE, ts_irr.data());
		readDataset("TimeStepReg", H5::PredType::NATIVE_DOUBLE, ts_reg.data());
		readDataset("TimeBlockIrr", H5::PredType::NATIVE_ULLONG, tbl_irr.data());
		readDataset("TimeBlockReg", H5::PredType::NATIVE_ULLONG, tbl_reg.data());
		readDataset("TimeLevelIrr", H5::PredType::NATIVE_INT, tl_irr.data());
		readDataset("TimeLevelReg", H5::PredType::NATIVE_INT, tl_reg.data());
		readDataset("isActive", H5::PredType::NATIVE_INT, isactive.data());
		readDataset("isCMptcl", H5::PredType::NATIVE_INT, iscm.data());
		readDataset("CMPtclIndex", H5::PredType::NATIVE_INT, cmindex.data());
		readDataset("NumberOfMember", H5::PredType::NATIVE_INT, nmember.data());

		// Read acceleration arrays
		H5::DataSet atot_ds = ptcl_group.openDataSet("a_tot");
		atot_ds.read(a_tot_flat.data(), H5::PredType::NATIVE_DOUBLE);
		H5::DataSet areg_ds = ptcl_group.openDataSet("a_reg");
		areg_ds.read(a_reg_flat.data(), H5::PredType::NATIVE_DOUBLE);
		H5::DataSet airr_ds = ptcl_group.openDataSet("a_irr");
		airr_ds.read(a_irr_flat.data(), H5::PredType::NATIVE_DOUBLE);

		// Restore particles
		for (int idx = 0; idx < saveCount; idx++) {
			int pi = pindex[idx];
			Particle* ptcl = &particles[pi];

			ptcl->PID = pid[idx];
			ptcl->ParticleIndex = pindex[idx];
			ptcl->ParticleType = ptype[idx];
			ptcl->Mass = mass[idx];
			ptcl->Position[0] = px[idx];
			ptcl->Position[1] = py[idx];
			ptcl->Position[2] = pz[idx];
			ptcl->Velocity[0] = vx[idx];
			ptcl->Velocity[1] = vy[idx];
			ptcl->Velocity[2] = vz[idx];
			ptcl->NumberOfNeighbor = nneighbor[idx];
			ptcl->NeighborsOffset = noffset[idx];
			ptcl->RadiusOfNeighbor = radius_nb[idx];
			ptcl->CurrentTimeIrr = ct_irr[idx];
			ptcl->CurrentTimeReg = ct_reg[idx];
			ptcl->CurrentBlockIrr = cb_irr[idx];
			ptcl->CurrentBlockReg = cb_reg[idx];
			ptcl->NextBlockIrr = nb_irr[idx];
			ptcl->TimeStepIrr = ts_irr[idx];
			ptcl->TimeStepReg = ts_reg[idx];
			ptcl->TimeBlockIrr = tbl_irr[idx];
			ptcl->TimeBlockReg = tbl_reg[idx];
			ptcl->TimeLevelIrr = tl_irr[idx];
			ptcl->TimeLevelReg = tl_reg[idx];
			ptcl->isActive = isactive[idx] != 0;
			ptcl->isCMptcl = iscm[idx] != 0;
			ptcl->CMPtclIndex = cmindex[idx];
			ptcl->NumberOfMember = nmember[idx];

			// Restore accelerations
			for (int d = 0; d < Dim; d++) {
				for (int o = 0; o < HERMITE_ORDER; o++) {
					int flat_idx = idx * Dim * HERMITE_ORDER + d * HERMITE_ORDER + o;
					ptcl->a_tot[d][o] = a_tot_flat[flat_idx];
					ptcl->a_reg[d][o] = a_reg_flat[flat_idx];
					ptcl->a_irr[d][o] = a_irr_flat[flat_idx];
				}
			}
		}

		// Read neighbor lists
		H5::Group neighbor_group = file.openGroup("/Neighbors");
		for (int idx = 0; idx < saveCount; idx++) {
			int pi = pindex[idx];
			Particle* ptcl = &particles[pi];
			if (ptcl->NumberOfNeighbor <= 0) continue;

			std::string ds_name = std::to_string(pi);
			try {
				H5::DataSet nb_ds = neighbor_group.openDataSet(ds_name);
				nb_ds.read(&Neighbors[ptcl->NeighborsOffset], H5::PredType::NATIVE_INT);
			} catch (...) {
				// Neighbor dataset may not exist for this particle
			}
		}

		file.close();
		std::cout << "Checkpoint loaded successfully: " << saveCount << " particles restored." << std::endl;
		std::cout << "Resuming from time: " << global_time * EnzoTimeStep * 1e10 / 1e6 << " Myr" << std::endl;

	} catch (H5::Exception& e) {
		std::cerr << "HDF5 error reading checkpoint: " << e.getDetailMsg() << std::endl;
		return false;
	}

	return true;
}