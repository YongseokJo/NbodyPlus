#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <stdexcept>
#include <cstring>
#include "toml.hpp"
#include "global.h"

// Static storage for C-string pointers (to maintain compatibility with existing code)
static std::string fname_storage;
static std::string foutput_storage;

class Config {
public:
	void read(const std::string& filename) {
		try {
			data_ = toml::parse(filename);
		} catch (const std::exception& e) {
			throw std::runtime_error("Unable to parse TOML configuration file: " + std::string(e.what()));
		}
	}

	// Required parameters (throw if not found)
	std::string getString(const std::string& key) const {
		try {
			return toml::find<std::string>(data_, key);
		} catch (const std::exception& e) {
			throw std::runtime_error("Required parameter not found: " + key);
		}
	}

	int getInt(const std::string& key) const {
		try {
			return toml::find<int>(data_, key);
		} catch (const std::exception& e) {
			throw std::runtime_error("Required parameter not found or invalid type: " + key);
		}
	}

	double getDouble(const std::string& key) const {
		try {
			return toml::find<double>(data_, key);
		} catch (const std::exception& e) {
			throw std::runtime_error("Required parameter not found or invalid type: " + key);
		}
	}

	// Optional parameters with defaults
	std::string getStringOr(const std::string& key, const std::string& defaultValue) const {
		try {
			return toml::find<std::string>(data_, key);
		} catch (...) {
			return defaultValue;
		}
	}

	int getIntOr(const std::string& key, int defaultValue) const {
		try {
			return toml::find<int>(data_, key);
		} catch (...) {
			return defaultValue;
		}
	}

	double getDoubleOr(const std::string& key, double defaultValue) const {
		try {
			return toml::find<double>(data_, key);
		} catch (...) {
			return defaultValue;
		}
	}

	bool getBoolOr(const std::string& key, bool defaultValue) const {
		try {
			return toml::find<bool>(data_, key);
		} catch (...) {
			return defaultValue;
		}
	}

	// Nested table access
	template<typename T>
	T getNestedOr(const std::string& table, const std::string& key, T defaultValue) const {
		try {
			const auto& tbl = toml::find(data_, table);
			return toml::find<T>(tbl, key);
		} catch (...) {
			return defaultValue;
		}
	}

	bool hasKey(const std::string& key) const {
		return data_.contains(key);
	}

	bool hasTable(const std::string& table) const {
		try {
			toml::find(data_, table);
			return true;
		} catch (...) {
			return false;
		}
	}

private:
	toml::value data_;
};

// Validation helpers
void validatePositive(double value, const std::string& name) {
	if (value <= 0) {
		throw std::runtime_error("Parameter '" + name + "' must be positive, got: " + std::to_string(value));
	}
}

void validateRange(int value, int min, int max, const std::string& name) {
	if (value < min || value > max) {
		throw std::runtime_error("Parameter '" + name + "' must be between " +
			std::to_string(min) + " and " + std::to_string(max) + ", got: " + std::to_string(value));
	}
}

void readParameterFile() {
	try {
		Config config;
		config.read(config_file);

		// === Required parameters ===
		fname_storage = config.getString("Filename");
		fname = const_cast<char*>(fname_storage.c_str());

		end_time = config.getDouble("StopTime");  // in years
		validatePositive(end_time, "StopTime");

		foutput_storage = config.getStringOr("OutputDirectory", "output");
		foutput = const_cast<char*>(foutput_storage.c_str());

		// === Numerics parameters (with defaults) ===
		eta = config.getDoubleOr("eta", 0.01);
		if (config.hasTable("numerics")) {
			eta = config.getNestedOr<double>("numerics", "eta", eta);
		}
		validatePositive(eta, "eta");

		fixed_num_neighbors = config.getIntOr("fixed_num_neighbors", 100);
		if (config.hasTable("numerics")) {
			fixed_num_neighbors = config.getNestedOr<int>("numerics", "fixed_num_neighbors", fixed_num_neighbors);
		}
		validateRange(fixed_num_neighbors, 10, MAX_NUM_NEIGHBOR, "fixed_num_neighbors");

		double initialRadius_pc = config.getDoubleOr("InitialRadius", 0.2);
		if (config.hasTable("numerics")) {
			initialRadius_pc = config.getNestedOr<double>("numerics", "InitialRadius", initialRadius_pc);
		}
		validatePositive(initialRadius_pc, "InitialRadius");
		initial_neighbor_radius = initialRadius_pc / position_unit;  // convert pc -> code units

		// Few-body search parameters
		double rsearch_pc = config.getDoubleOr("r_search", 2.5e-4);
		if (config.hasTable("numerics")) {
			rsearch_pc = config.getNestedOr<double>("numerics", "r_search", rsearch_pc);
		}
		validatePositive(rsearch_pc, "r_search");
		r_search = rsearch_pc / position_unit;  // convert pc -> code units

		double tsearch_myr = config.getDoubleOr("t_search", 1e-6);
		if (config.hasTable("numerics")) {
			tsearch_myr = config.getNestedOr<double>("numerics", "t_search", tsearch_myr);
		}
		validatePositive(tsearch_myr, "t_search");
		// t_search will be set properly after enzo_time_step is computed

		// === Output parameters ===
		output_time_step = config.getDoubleOr("dtOutput", end_time / 10.0);
		if (config.hasTable("output")) {
			output_time_step = config.getNestedOr<double>("output", "dtOutput", output_time_step);
		}
		validatePositive(output_time_step, "dtOutput");

		use_compression = config.getBoolOr("Compression", true);
		if (config.hasTable("output")) {
			use_compression = config.getNestedOr<bool>("output", "Compression", use_compression);
		}

		compression_level = config.getIntOr("compression_level", 6);
		if (config.hasTable("output")) {
			compression_level = config.getNestedOr<int>("output", "compression_level", compression_level);
		}
		validateRange(compression_level, 1, 9, "compression_level");

		// === Restart parameters ===
		restart_enabled = config.getBoolOr("restart_enabled", false);
		if (config.hasTable("restart")) {
			restart_enabled = config.getNestedOr<bool>("restart", "Enabled", restart_enabled);
		}

		checkpoint_file = config.getStringOr("checkpoint_file", "");
		if (config.hasTable("restart")) {
			checkpoint_file = config.getNestedOr<std::string>("restart", "checkpoint_file", checkpoint_file);
		}

		// === Compute derived quantities ===
		enzo_time_step = end_time / 1e10;  // end_time should be yr
		output_time_step = output_time_step / end_time;  // normalize to simulation time

		// Now set t_search properly with enzo_time_step
		t_search = tsearch_myr / (enzo_time_step * 1e4);  // Myr -> code units

		// === Print configuration summary ===
		if (my_rank == ROOT) {
			std::cout << "\n========== ABYSS Configuration ==========\n";
			std::cout << "Input file:        " << fname << std::endl;
			std::cout << "Output directory:  " << foutput << std::endl;
			std::cout << "\n--- Time ---\n";
			std::cout << "End time:          " << end_time / 1e6 << " Myr\n";
			std::cout << "Output interval:   " << output_time_step * enzo_time_step * 1e4 << " Myr\n";
			std::cout << "\n--- Numerics ---\n";
			std::cout << "eta:               " << eta << std::endl;
			std::cout << "fixed_num_neighbors:    " << fixed_num_neighbors << std::endl;
			std::cout << "InitialRadius:     " << initialRadius_pc << " pc\n";
			std::cout << "r_search:           " << rsearch_pc << " pc\n";
			std::cout << "t_search:           " << tsearch_myr << " Myr\n";
			std::cout << "\n--- Output ---\n";
			std::cout << "Compression:       " << (use_compression ? "enabled" : "disabled") << std::endl;
			if (use_compression) {
				std::cout << "Compression level: " << compression_level << std::endl;
			}
			std::cout << "\n--- Restart ---\n";
			std::cout << "Restart enabled:   " << (restart_enabled ? "yes" : "no") << std::endl;
			if (restart_enabled && !checkpoint_file.empty()) {
				std::cout << "Checkpoint file:   " << checkpoint_file << std::endl;
			}
			std::cout << "==========================================\n\n";
		}

	} catch (const std::exception& e) {
		std::cerr << "Configuration TASK_ERROR: " << e.what() << std::endl;
		std::exit(EXIT_FAILURE);
	}
}
