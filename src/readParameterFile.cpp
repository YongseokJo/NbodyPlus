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

		endTime = config.getDouble("StopTime");  // in years
		validatePositive(endTime, "StopTime");

		foutput_storage = config.getStringOr("OutputDirectory", "output");
		foutput = const_cast<char*>(foutput_storage.c_str());

		// === Numerics parameters (with defaults) ===
		eta = config.getDoubleOr("eta", 0.01);
		if (config.hasTable("numerics")) {
			eta = config.getNestedOr<double>("numerics", "eta", eta);
		}
		validatePositive(eta, "eta");

		FixNumNeighbor = config.getIntOr("FixNumNeighbor", 100);
		if (config.hasTable("numerics")) {
			FixNumNeighbor = config.getNestedOr<int>("numerics", "FixNumNeighbor", FixNumNeighbor);
		}
		validateRange(FixNumNeighbor, 10, MaxNumNeighbor, "FixNumNeighbor");

		double initialRadius_pc = config.getDoubleOr("InitialRadius", 0.2);
		if (config.hasTable("numerics")) {
			initialRadius_pc = config.getNestedOr<double>("numerics", "InitialRadius", initialRadius_pc);
		}
		validatePositive(initialRadius_pc, "InitialRadius");
		InitialNeighborRadius = initialRadius_pc / position_unit;  // convert pc -> code units

		// Few-body search parameters
		double rsearch_pc = config.getDoubleOr("RSearch", 2.5e-4);
		if (config.hasTable("numerics")) {
			rsearch_pc = config.getNestedOr<double>("numerics", "RSearch", rsearch_pc);
		}
		validatePositive(rsearch_pc, "RSearch");
		RSearch = rsearch_pc / position_unit;  // convert pc -> code units

		double tsearch_myr = config.getDoubleOr("TSearch", 1e-6);
		if (config.hasTable("numerics")) {
			tsearch_myr = config.getNestedOr<double>("numerics", "TSearch", tsearch_myr);
		}
		validatePositive(tsearch_myr, "TSearch");
		// TSearch will be set properly after EnzoTimeStep is computed

		// === Output parameters ===
		outputTimeStep = config.getDoubleOr("dtOutput", endTime / 10.0);
		if (config.hasTable("output")) {
			outputTimeStep = config.getNestedOr<double>("output", "dtOutput", outputTimeStep);
		}
		validatePositive(outputTimeStep, "dtOutput");

		UseCompression = config.getBoolOr("Compression", true);
		if (config.hasTable("output")) {
			UseCompression = config.getNestedOr<bool>("output", "Compression", UseCompression);
		}

		CompressionLevel = config.getIntOr("CompressionLevel", 6);
		if (config.hasTable("output")) {
			CompressionLevel = config.getNestedOr<int>("output", "CompressionLevel", CompressionLevel);
		}
		validateRange(CompressionLevel, 1, 9, "CompressionLevel");

		// === Restart parameters ===
		RestartEnabled = config.getBoolOr("RestartEnabled", false);
		if (config.hasTable("restart")) {
			RestartEnabled = config.getNestedOr<bool>("restart", "Enabled", RestartEnabled);
		}

		CheckpointFile = config.getStringOr("CheckpointFile", "");
		if (config.hasTable("restart")) {
			CheckpointFile = config.getNestedOr<std::string>("restart", "CheckpointFile", CheckpointFile);
		}

		// === Compute derived quantities ===
		EnzoTimeStep = endTime / 1e10;  // endTime should be yr
		outputTimeStep = outputTimeStep / endTime;  // normalize to simulation time

		// Now set TSearch properly with EnzoTimeStep
		TSearch = tsearch_myr / (EnzoTimeStep * 1e4);  // Myr -> code units

		// === Print configuration summary ===
		if (MyRank == ROOT) {
			std::cout << "\n========== ABYSS Configuration ==========\n";
			std::cout << "Input file:        " << fname << std::endl;
			std::cout << "Output directory:  " << foutput << std::endl;
			std::cout << "\n--- Time ---\n";
			std::cout << "End time:          " << endTime / 1e6 << " Myr\n";
			std::cout << "Output interval:   " << outputTimeStep * EnzoTimeStep * 1e4 << " Myr\n";
			std::cout << "\n--- Numerics ---\n";
			std::cout << "eta:               " << eta << std::endl;
			std::cout << "FixNumNeighbor:    " << FixNumNeighbor << std::endl;
			std::cout << "InitialRadius:     " << initialRadius_pc << " pc\n";
			std::cout << "RSearch:           " << rsearch_pc << " pc\n";
			std::cout << "TSearch:           " << tsearch_myr << " Myr\n";
			std::cout << "\n--- Output ---\n";
			std::cout << "Compression:       " << (UseCompression ? "enabled" : "disabled") << std::endl;
			if (UseCompression) {
				std::cout << "Compression level: " << CompressionLevel << std::endl;
			}
			std::cout << "\n--- Restart ---\n";
			std::cout << "Restart enabled:   " << (RestartEnabled ? "yes" : "no") << std::endl;
			if (RestartEnabled && !CheckpointFile.empty()) {
				std::cout << "Checkpoint file:   " << CheckpointFile << std::endl;
			}
			std::cout << "==========================================\n\n";
		}

	} catch (const std::exception& e) {
		std::cerr << "Configuration Error: " << e.what() << std::endl;
		std::exit(EXIT_FAILURE);
	}
}
