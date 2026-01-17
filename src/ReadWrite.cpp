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
		std::cerr << "TASK_ERROR: Could not open the file." << std::endl;
		return false;
	}

	num_particles = getLineNumber();
	assert(num_particles < MAX_NUM_PARTICLE);
	new_cm_pid = num_particles;
	last_particle_index = num_particles - 1;
	g_state->last_particle_index = last_particle_index;

	
	double** data = new double*[num_particles];
	for (int i = 0; i < num_particles; ++i) {
		data[i] = new double[NUM_COLUMNS];
	}

	// Initialization
	for (int i = 0; i < num_particles; ++i) {
		for (int j = 0; j < NUM_COLUMNS; ++j) {
			data[i][j] = 0;
		}
	}

	int row = 0;
	std::string line;
	while (std::getline(inputFile, line) && row < num_particles) { // Read lines from the file
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
	for (int i=0; i<num_particles; i++) {
		particles[i].normalize_particle();
	}
	inputFile.close();

#ifdef SEVN
	initializeStellarEvolution();
#endif

	// Deallocate memory
	for (int i = 0; i < num_particles; ++i) {
		delete[] data[i];
	}
	delete[] data;

	return true;
}


int getLineNumber() {
	std::ifstream inputFile(fname); // Open the file

	if (!inputFile) {
		std::cerr << "TASK_ERROR: Could not open the file." << std::endl;
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
		std::cerr << "TASK_ERROR creating folder." << std::endl;
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

	Queue queue = {TASK_GET_TOTAL_ENERGY, -1, -1.0};
	for (int i=0; i< num_workers; i++)
		MPI_Send(&queue, 1, queue_type_mpi, i+1, QUEUE_TAG, MPI_COMM_WORLD);

	MPI_Reduce(MPI_IN_PLACE, &energy_binary,		1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
	MPI_Reduce(MPI_IN_PLACE, &energy_binary_sd,	1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
	MPI_Reduce(MPI_IN_PLACE, &energy_merger,		1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
	MPI_Reduce(MPI_IN_PLACE, &energy_pn,			1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);

	double unit_energy = mass_unit * (velocity_unit/yr*pc/1e5) * (velocity_unit/yr*pc/1e5);

    // Write particle data to HDF5
	Particle *ptcl;
	double pos[DIM], vel[DIM];

	// for performance test by EW 2025.3.13
	int Index_minIrr, Index_minReg;
	double minTimeStepIrr = 1.;
	double minTimeStepReg = 1.;

	std::unordered_set<int> CMPtclsSet;

	// Count active particles and collect data
	int activeCount = 0;
	for (int i=0; i<=last_particle_index; i++) {
		ptcl = &particles[i];
		if (!ptcl->is_active) continue;

		if (ptcl->time_step_irr < minTimeStepIrr) {
			Index_minIrr = ptcl->particle_index;
			minTimeStepIrr = ptcl->time_step_irr;
		}
		if (ptcl->time_step_reg < minTimeStepReg) {
			Index_minReg = ptcl->particle_index;
			minTimeStepReg = ptcl->time_step_reg;
		}

		if (ptcl->is_cm_particle)
			CMPtclsSet.insert(i);

		activeCount++;
	}

	// Count group members
	for (int i: CMPtclsSet) {
		ptcl = &particles[i];
		activeCount += ptcl->num_members;
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
	for (int i=0; i<=last_particle_index; i++) {
		ptcl = &particles[i];
		if (!ptcl->is_active) continue;

		ptcl->predict_particle_second_order(current_time - ptcl->current_time_irr, pos, vel);

		pid_data[idx] = ptcl->pid;
		mass_data[idx] = ptcl->mass * mass_unit;
		x_data[idx] = pos[0] * position_unit;
		y_data[idx] = pos[1] * position_unit;
		z_data[idx] = pos[2] * position_unit;
		vx_data[idx] = vel[0] * velocity_unit/yr*pc/1e5;
		vy_data[idx] = vel[1] * velocity_unit/yr*pc/1e5;
		vz_data[idx] = vel[2] * velocity_unit/yr*pc/1e5;
#ifdef SEVN
		if (ptcl->stellar_evolution == nullptr)
			type_data[idx] = 1;
		else if (!ptcl->stellar_evolution->amiremnant())
			type_data[idx] = int(ptcl->stellar_evolution->getp(Phase::ID));
		else
			type_data[idx] = 8+int(ptcl->stellar_evolution->getp(RemnantType::ID));
#endif
		idx++;
	}

	// Collect group member data
	for (int i: CMPtclsSet) {
		ptcl = &particles[i];
		ptcl->predict_particle_second_order(current_time - ptcl->current_time_irr, pos, vel);

		Particle* members;
		for (int j=0; j < ptcl->num_members; j++) {
			members = &particles[ptcl->members[j]];

			pid_data[idx] = members->pid;
			mass_data[idx] = members->mass * mass_unit;
			x_data[idx] = (pos[0] - ptcl->position[0] + members->position[0]) * position_unit;
			y_data[idx] = (pos[1] - ptcl->position[1] + members->position[1]) * position_unit;
			z_data[idx] = (pos[2] - ptcl->position[2] + members->position[2]) * position_unit;
			vx_data[idx] = (vel[0] - ptcl->velocity[0] + members->velocity[0]) * velocity_unit/yr*pc/1e5;
			vy_data[idx] = (vel[1] - ptcl->velocity[1] + members->velocity[1]) * velocity_unit/yr*pc/1e5;
			vz_data[idx] = (vel[2] - ptcl->velocity[2] + members->velocity[2]) * velocity_unit/yr*pc/1e5;
#ifdef SEVN
			if (members->stellar_evolution == nullptr)
				type_data[idx] = 1;
			else if (!members->stellar_evolution->amiremnant())
				type_data[idx] = int(members->stellar_evolution->getp(Phase::ID));
			else
				type_data[idx] = 8+int(members->stellar_evolution->getp(RemnantType::ID));
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

			H5::Attribute nneighbor_attr = file.createAttribute("fixed_num_neighbors", H5::PredType::NATIVE_INT, scalar_space);
			nneighbor_attr.write(H5::PredType::NATIVE_INT, &fixed_num_neighbors);

			double init_radius_pc = initial_neighbor_radius * position_unit;
			H5::Attribute radius_attr = file.createAttribute("InitialRadius_pc", H5::PredType::NATIVE_DOUBLE, scalar_space);
			radius_attr.write(H5::PredType::NATIVE_DOUBLE, &init_radius_pc);

			double rsearch_pc = r_search * position_unit;
			H5::Attribute rsearch_attr = file.createAttribute("RSearch_pc", H5::PredType::NATIVE_DOUBLE, scalar_space);
			rsearch_attr.write(H5::PredType::NATIVE_DOUBLE, &rsearch_pc);

			double end_time_myr = end_time / 1e6;  // yr -> Myr
			H5::Attribute endtime_attr = file.createAttribute("EndTime_Myr", H5::PredType::NATIVE_DOUBLE, scalar_space);
			endtime_attr.write(H5::PredType::NATIVE_DOUBLE, &end_time_myr);

			H5::Attribute nproc_attr = file.createAttribute("NumberOfProcessors", H5::PredType::NATIVE_INT, scalar_space);
			nproc_attr.write(H5::PredType::NATIVE_INT, &num_processors);

			int init_npart = num_particles;
			H5::Attribute initnpart_attr = file.createAttribute("InitialParticleCount", H5::PredType::NATIVE_INT, scalar_space);
			initnpart_attr.write(H5::PredType::NATIVE_INT, &init_npart);

			// Compression settings
			H5::Attribute compress_attr = file.createAttribute("CompressionEnabled", H5::PredType::NATIVE_HBOOL, scalar_space);
			hbool_t compress_val = use_compression ? 1 : 0;
			compress_attr.write(H5::PredType::NATIVE_HBOOL, &compress_val);

			if (use_compression) {
				H5::Attribute level_attr = file.createAttribute("compression_level", H5::PredType::NATIVE_INT, scalar_space);
				level_attr.write(H5::PredType::NATIVE_INT, &compression_level);
			}

		} else {
			file = H5::H5File(filename, H5F_ACC_RDWR);
		}

		// Create group for this timestep
		std::string groupName = "/Step_" + std::to_string(outputNum);
		H5::Group group = file.createGroup(groupName);

		// Write metadata
		double time_value = current_time*enzo_time_step*1e10/1e6;
		hsize_t scalar_dim = 1;
		H5::DataSpace scalar_space(1, &scalar_dim);

		H5::Attribute time_attr = group.createAttribute("Time_Myr", H5::PredType::NATIVE_DOUBLE, scalar_space);
		time_attr.write(H5::PredType::NATIVE_DOUBLE, &time_value);

		H5::Attribute npart_attr = group.createAttribute("num_particles", H5::PredType::NATIVE_INT, scalar_space);
		npart_attr.write(H5::PredType::NATIVE_INT, &activeCount);

		double e_bin = energy_binary * unit_energy;
		double e_bin_sd = energy_binary_sd * unit_energy;
		double e_mrg = energy_merger * unit_energy;
		double e_pn = energy_pn * unit_energy;

		H5::Attribute ebinary_attr = group.createAttribute("energy_binary", H5::PredType::NATIVE_DOUBLE, scalar_space);
		ebinary_attr.write(H5::PredType::NATIVE_DOUBLE, &e_bin);

		H5::Attribute ebinary_sd_attr = group.createAttribute("energy_binary_sd", H5::PredType::NATIVE_DOUBLE, scalar_space);
		ebinary_sd_attr.write(H5::PredType::NATIVE_DOUBLE, &e_bin_sd);

		H5::Attribute emerger_attr = group.createAttribute("energy_merger", H5::PredType::NATIVE_DOUBLE, scalar_space);
		emerger_attr.write(H5::PredType::NATIVE_DOUBLE, &e_mrg);

		H5::Attribute epn_attr = group.createAttribute("energy_pn", H5::PredType::NATIVE_DOUBLE, scalar_space);
		epn_attr.write(H5::PredType::NATIVE_DOUBLE, &e_pn);

		// Write datasets with optional compression and chunking
		hsize_t dims[1] = {(hsize_t)activeCount};
		H5::DataSpace dataspace(1, dims);

		// Create dataset creation property list with chunking and compression
		H5::DSetCreatPropList plist;
		if (use_compression && activeCount > 0) {
			// Set chunk size (minimum of activeCount or 10000 for efficiency)
			hsize_t chunk_dims[1] = {std::min((hsize_t)activeCount, (hsize_t)10000)};
			plist.setChunk(1, chunk_dims);
			plist.setDeflate(compression_level);  // GZIP compression
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

	energy_binary = 0.0;
	energy_binary_sd = 0.0;
	energy_merger = 0.0;
	energy_pn = 0.0;

	std::cout << "Data written to HDF5 file successfully!" << std::endl;

#ifdef PERFORMANCETRACE
	if (outputNum != 0) {
		std::cout << "--------------Performance-Summary--------------" << std::endl;
		fprintf(stdout, "Simulation Time: %f Myr\n", current_time*enzo_time_step*1e10/1e6);

		Particle* members = &particles[Index_minIrr];
		fprintf(stdout, "Particle Info with minimum TimeStepIrr...\n");
		members->print_particle_info(stdout);

		members = &particles[Index_minReg];
		fprintf(stdout, "Particle Info with minimum TimeStepReg...\n");
		members->print_particle_info(stdout);

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
		std::cout << "Stellar Evolution: " << performance.stellar_evolution*1e-9 << " s";
		std::cout << " (" << 100.0 * performance.stellar_evolution / performance.WholeRoutine << " %)" << std::endl;
		performance.stellar_evolution = 0;
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
		double sim_time_myr = current_time * enzo_time_step * 1e10 / 1e6;
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
					<< std::setw(width) << ptcl->pid
					<< std::setw(width) << ptcl->mass*mass_unit
                    << std::setw(width) << pos[0]*position_unit
                    << std::setw(width) << pos[1]*position_unit
                    << std::setw(width) << pos[2]*position_unit
                    << std::setw(width) << vel[0]*velocity_unit/yr*pc/1e5
                    << std::setw(width) << vel[1]*velocity_unit/yr*pc/1e5;
#ifdef SEVN
		outputFile	<< std::setw(width) << vel[2]*velocity_unit/yr*pc/1e5;
		if (ptcl->stellar_evolution == nullptr)
			outputFile << std::setw(width) << "1" << '\n'; // Main-sequence star;
		else if (!ptcl->stellar_evolution->amiremnant())
			outputFile << std::setw(width) << int(ptcl->stellar_evolution->getp(Phase::ID)) << '\n';
		else
			outputFile << std::setw(width) << 8+int(ptcl->stellar_evolution->getp(RemnantType::ID)) << '\n';
#else
		outputFile	<< std::setw(width) << vel[2]*velocity_unit/yr*pc/1e5 << '\n';
#endif
}

// This function is for group members cause group members have pos, vel in original frame, not predicted values.
void write_out_group(std::ofstream& outputFile, const Particle* ptclCM, const Particle* ptcl, const double *pos, const double *vel) {
        // outputFile  << std::left
		outputFile  << std::left << std::fixed << std::setprecision(8) // Eunwoo test
					<< std::setw(width) << ptcl->pid
					<< std::setw(width) << ptcl->mass*mass_unit
                    << std::setw(width) << (pos[0] - ptclCM->position[0] + ptcl->position[0])*position_unit
                    << std::setw(width) << (pos[1] - ptclCM->position[1] + ptcl->position[1])*position_unit
                    << std::setw(width) << (pos[2] - ptclCM->position[2] + ptcl->position[2])*position_unit
                    << std::setw(width) << (vel[0] - ptclCM->velocity[0] + ptcl->velocity[0])*velocity_unit/yr*pc/1e5
                    << std::setw(width) << (vel[1] - ptclCM->velocity[1] + ptcl->velocity[1])*velocity_unit/yr*pc/1e5;
#ifdef SEVN
		outputFile	<< std::setw(width) << (vel[2] - ptclCM->velocity[2] + ptcl->velocity[2])*velocity_unit/yr*pc/1e5;
		if (ptcl->stellar_evolution == nullptr)
			outputFile << std::setw(width) << "1" << '\n';
		else if (!ptcl->stellar_evolution->amiremnant())
			outputFile << std::setw(width) << int(ptcl->stellar_evolution->getp(Phase::ID)) << '\n';
		else
			outputFile << std::setw(width) << 8+int(ptcl->stellar_evolution->getp(RemnantType::ID)) << '\n';
#else
		outputFile << std::setw(width) << (vel[2] - ptclCM->velocity[2] + ptcl->velocity[2])*velocity_unit/yr*pc/1e5 << '\n';
#endif
}


// ============================================================================
// Checkpoint/Restart Functions
// ============================================================================

int writeCheckpoint(double current_time, int checkpointNum) {
	if (!restart_enabled) return 0;

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

		H5::Attribute nrt_attr = file.createAttribute("next_reg_time_block", H5::PredType::NATIVE_ULLONG, scalar_space);
		nrt_attr.write(H5::PredType::NATIVE_ULLONG, &next_reg_time_block);

		H5::Attribute tb_attr = file.createAttribute("time_block", H5::PredType::NATIVE_INT, scalar_space);
		tb_attr.write(H5::PredType::NATIVE_INT, &time_block);

		H5::Attribute bm_attr = file.createAttribute("block_max", H5::PredType::NATIVE_ULLONG, scalar_space);
		bm_attr.write(H5::PredType::NATIVE_ULLONG, &block_max);

		// Particle counts
		H5::Attribute lpi_attr = file.createAttribute("last_particle_index", H5::PredType::NATIVE_INT, scalar_space);
		lpi_attr.write(H5::PredType::NATIVE_INT, &last_particle_index);

		H5::Attribute nop_attr = file.createAttribute("num_particles", H5::PredType::NATIVE_INT, scalar_space);
		nop_attr.write(H5::PredType::NATIVE_INT, &num_particles);

		H5::Attribute ncm_attr = file.createAttribute("new_cm_pid", H5::PredType::NATIVE_INT, scalar_space);
		ncm_attr.write(H5::PredType::NATIVE_INT, &new_cm_pid);

		// Output tracking
		H5::Attribute ot_attr = file.createAttribute("output_time", H5::PredType::NATIVE_DOUBLE, scalar_space);
		ot_attr.write(H5::PredType::NATIVE_DOUBLE, &output_time);

		H5::Attribute on_attr = file.createAttribute("output_num", H5::PredType::NATIVE_INT, scalar_space);
		on_attr.write(H5::PredType::NATIVE_INT, &output_num);

		// Count particles to save (active ones plus their data)
		int saveCount = 0;
		for (int i = 0; i <= last_particle_index; i++) {
			if (particles[i].pid >= 0) saveCount++;
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
		std::vector<ull_t> cb_irr(saveCount), cb_reg(saveCount), nb_irr(saveCount);
		std::vector<ull_t> tbl_irr(saveCount), tbl_reg(saveCount);
		std::vector<int> tl_irr(saveCount), tl_reg(saveCount);

		// Acceleration arrays (flattened)
		std::vector<double> a_tot_flat(saveCount * DIM * HERMITE_ORDER);
		std::vector<double> a_reg_flat(saveCount * DIM * HERMITE_ORDER);
		std::vector<double> a_irr_flat(saveCount * DIM * HERMITE_ORDER);

		// Collect data
		int idx = 0;
		for (int i = 0; i <= last_particle_index; i++) {
			Particle* ptcl = &particles[i];
			if (ptcl->pid < 0) continue;

			pid[idx] = ptcl->pid;
			pindex[idx] = ptcl->particle_index;
			ptype[idx] = ptcl->particle_type;
			mass[idx] = ptcl->mass;
			px[idx] = ptcl->position[0];
			py[idx] = ptcl->position[1];
			pz[idx] = ptcl->position[2];
			vx[idx] = ptcl->velocity[0];
			vy[idx] = ptcl->velocity[1];
			vz[idx] = ptcl->velocity[2];
			nneighbor[idx] = ptcl->num_neighbors;
			noffset[idx] = ptcl->neighbors_offset;
			radius_nb[idx] = ptcl->neighbor_radius_sq;
			ct_irr[idx] = ptcl->current_time_irr;
			ct_reg[idx] = ptcl->current_time_reg;
			cb_irr[idx] = ptcl->current_block_irr;
			cb_reg[idx] = ptcl->current_block_reg;
			nb_irr[idx] = ptcl->next_block_irr;
			ts_irr[idx] = ptcl->time_step_irr;
			ts_reg[idx] = ptcl->time_step_reg;
			tbl_irr[idx] = ptcl->time_block_irr;
			tbl_reg[idx] = ptcl->time_block_reg;
			tl_irr[idx] = ptcl->time_level_irr;
			tl_reg[idx] = ptcl->time_level_reg;
			isactive[idx] = ptcl->is_active ? 1 : 0;
			iscm[idx] = ptcl->is_cm_particle ? 1 : 0;
			cmindex[idx] = ptcl->cm_particle_index;
			nmember[idx] = ptcl->num_members;

			// Copy accelerations (flattened)
			for (int d = 0; d < DIM; d++) {
				for (int o = 0; o < HERMITE_ORDER; o++) {
					int flat_idx = idx * DIM * HERMITE_ORDER + d * HERMITE_ORDER + o;
					a_tot_flat[flat_idx] = ptcl->acc_total[d][o];
					a_reg_flat[flat_idx] = ptcl->acc_regular[d][o];
					a_irr_flat[flat_idx] = ptcl->acc_irregular[d][o];
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

		// Write acceleration arrays (2D: saveCount x (DIM*HERMITE_ORDER))
		hsize_t acc_dims[2] = {(hsize_t)saveCount, (hsize_t)(DIM * HERMITE_ORDER)};
		H5::DataSpace acc_space(2, acc_dims);
		H5::DataSet atot_ds = ptcl_group.createDataSet("a_tot", H5::PredType::NATIVE_DOUBLE, acc_space);
		atot_ds.write(a_tot_flat.data(), H5::PredType::NATIVE_DOUBLE);
		H5::DataSet areg_ds = ptcl_group.createDataSet("a_reg", H5::PredType::NATIVE_DOUBLE, acc_space);
		areg_ds.write(a_reg_flat.data(), H5::PredType::NATIVE_DOUBLE);
		H5::DataSet airr_ds = ptcl_group.createDataSet("a_irr", H5::PredType::NATIVE_DOUBLE, acc_space);
		airr_ds.write(a_irr_flat.data(), H5::PredType::NATIVE_DOUBLE);

		// Save neighbor lists for particles that have neighbors
		H5::Group neighbor_group = file.createGroup("/neighbors");
		for (int i = 0; i <= last_particle_index; i++) {
			Particle* ptcl = &particles[i];
			if (ptcl->pid < 0 || ptcl->num_neighbors <= 0) continue;

			std::string ds_name = std::to_string(ptcl->particle_index);
			hsize_t nb_dims[1] = {(hsize_t)ptcl->num_neighbors};
			H5::DataSpace nb_space(1, nb_dims);
			H5::DataSet nb_ds = neighbor_group.createDataSet(ds_name, H5::PredType::NATIVE_INT, nb_space);
			nb_ds.write(&neighbors[ptcl->neighbors_offset], H5::PredType::NATIVE_INT);
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

		H5::Attribute nrt_attr = file.openAttribute("next_reg_time_block");
		nrt_attr.read(H5::PredType::NATIVE_ULLONG, &next_reg_time_block);

		H5::Attribute tb_attr = file.openAttribute("time_block");
		tb_attr.read(H5::PredType::NATIVE_INT, &time_block);

		H5::Attribute bm_attr = file.openAttribute("block_max");
		bm_attr.read(H5::PredType::NATIVE_ULLONG, &block_max);

		H5::Attribute lpi_attr = file.openAttribute("last_particle_index");
		lpi_attr.read(H5::PredType::NATIVE_INT, &last_particle_index);

		H5::Attribute nop_attr = file.openAttribute("num_particles");
		nop_attr.read(H5::PredType::NATIVE_INT, &num_particles);

		H5::Attribute ncm_attr = file.openAttribute("new_cm_pid");
		ncm_attr.read(H5::PredType::NATIVE_INT, &new_cm_pid);

		H5::Attribute ot_attr = file.openAttribute("output_time");
		ot_attr.read(H5::PredType::NATIVE_DOUBLE, &output_time);

		H5::Attribute on_attr = file.openAttribute("output_num");
		on_attr.read(H5::PredType::NATIVE_INT, &output_num);

		// Update g_state
		g_state->last_particle_index = last_particle_index;
		g_state->next_reg_time_block = next_reg_time_block;
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
		std::vector<ull_t> cb_irr(saveCount), cb_reg(saveCount), nb_irr(saveCount);
		std::vector<ull_t> tbl_irr(saveCount), tbl_reg(saveCount);
		std::vector<int> tl_irr(saveCount), tl_reg(saveCount);
		std::vector<double> a_tot_flat(saveCount * DIM * HERMITE_ORDER);
		std::vector<double> a_reg_flat(saveCount * DIM * HERMITE_ORDER);
		std::vector<double> a_irr_flat(saveCount * DIM * HERMITE_ORDER);

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

			ptcl->pid = pid[idx];
			ptcl->particle_index = pindex[idx];
			ptcl->particle_type = ptype[idx];
			ptcl->mass = mass[idx];
			ptcl->position[0] = px[idx];
			ptcl->position[1] = py[idx];
			ptcl->position[2] = pz[idx];
			ptcl->velocity[0] = vx[idx];
			ptcl->velocity[1] = vy[idx];
			ptcl->velocity[2] = vz[idx];
			ptcl->num_neighbors = nneighbor[idx];
			ptcl->neighbors_offset = noffset[idx];
			ptcl->neighbor_radius_sq = radius_nb[idx];
			ptcl->current_time_irr = ct_irr[idx];
			ptcl->current_time_reg = ct_reg[idx];
			ptcl->current_block_irr = cb_irr[idx];
			ptcl->current_block_reg = cb_reg[idx];
			ptcl->next_block_irr = nb_irr[idx];
			ptcl->time_step_irr = ts_irr[idx];
			ptcl->time_step_reg = ts_reg[idx];
			ptcl->time_block_irr = tbl_irr[idx];
			ptcl->time_block_reg = tbl_reg[idx];
			ptcl->time_level_irr = tl_irr[idx];
			ptcl->time_level_reg = tl_reg[idx];
			ptcl->is_active = isactive[idx] != 0;
			ptcl->is_cm_particle = iscm[idx] != 0;
			ptcl->cm_particle_index = cmindex[idx];
			ptcl->num_members = nmember[idx];

			// Restore accelerations
			for (int d = 0; d < DIM; d++) {
				for (int o = 0; o < HERMITE_ORDER; o++) {
					int flat_idx = idx * DIM * HERMITE_ORDER + d * HERMITE_ORDER + o;
					ptcl->acc_total[d][o] = a_tot_flat[flat_idx];
					ptcl->acc_regular[d][o] = a_reg_flat[flat_idx];
					ptcl->acc_irregular[d][o] = a_irr_flat[flat_idx];
				}
			}
		}

		// Read neighbor lists
		H5::Group neighbor_group = file.openGroup("/neighbors");
		for (int idx = 0; idx < saveCount; idx++) {
			int pi = pindex[idx];
			Particle* ptcl = &particles[pi];
			if (ptcl->num_neighbors <= 0) continue;

			std::string ds_name = std::to_string(pi);
			try {
				H5::DataSet nb_ds = neighbor_group.openDataSet(ds_name);
				nb_ds.read(&neighbors[ptcl->neighbors_offset], H5::PredType::NATIVE_INT);
			} catch (...) {
				// Neighbor dataset may not exist for this particle
			}
		}

		file.close();
		std::cout << "Checkpoint loaded successfully: " << saveCount << " particles restored." << std::endl;
		std::cout << "Resuming from time: " << global_time * enzo_time_step * 1e10 / 1e6 << " Myr" << std::endl;

	} catch (H5::Exception& e) {
		std::cerr << "HDF5 error reading checkpoint: " << e.getDetailMsg() << std::endl;
		return false;
	}

	return true;
}