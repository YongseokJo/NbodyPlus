#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <stdexcept>
#include <cstring>
#include <algorithm>  // for std::min
#include <climits>    // for INT_MAX
#include "toml.hpp"
#include "global.h"
#include "mcluster_config.h"

// Static storage for C-string pointers (to maintain compatibility with existing code)
static std::string fname_storage;
static std::string foutput_storage;

// Global mcluster configuration
MclusterConfig mcluster_config;

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

	const toml::value& getData() const {
		return data_;
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

// ============================================================================
// McLuster config parsing and validation
// ============================================================================

// Levenshtein distance for typo detection in parameter names
int levenshteinDistance(const std::string& s1, const std::string& s2) {
	std::vector<std::vector<int>> dp(s1.size() + 1, std::vector<int>(s2.size() + 1));
	for (size_t i = 0; i <= s1.size(); ++i) dp[i][0] = static_cast<int>(i);
	for (size_t j = 0; j <= s2.size(); ++j) dp[0][j] = static_cast<int>(j);
	for (size_t i = 1; i <= s1.size(); ++i) {
		for (size_t j = 1; j <= s2.size(); ++j) {
			int cost = (s1[i-1] == s2[j-1]) ? 0 : 1;
			dp[i][j] = std::min({dp[i-1][j] + 1, dp[i][j-1] + 1, dp[i-1][j-1] + cost});
		}
	}
	return dp[s1.size()][s2.size()];
}

// Suggest similar parameter name for typo correction
std::string suggestSimilarParam(const std::string& unknown) {
	int minDist = INT_MAX;
	std::string suggestion;
	for (const auto& valid : MCLUSTER_VALID_PARAMS) {
		int dist = levenshteinDistance(unknown, valid);
		if (dist < minDist && dist <= 2) {
			minDist = dist;
			suggestion = valid;
		}
	}
	return suggestion;
}

// Range validation for mcluster double parameters
void validateMclusterRange(double value, double min, double max, const std::string& name) {
	if (value < min || value > max) {
		throw std::runtime_error("[mcluster] section: Parameter '" + name + "' must be between " +
			std::to_string(min) + " and " + std::to_string(max) + ", got: " + std::to_string(value));
	}
}

// Parse [mcluster] section from TOML config
void parseMclusterSection(const Config& config, MclusterConfig& mc) {
	if (!config.hasTable("mcluster")) {
		mc.has_mcluster_section = false;
		return;
	}
	mc.has_mcluster_section = true;

	// Get the mcluster table for key validation
	const auto& data = config.getData();
	const auto& mcluster_table = toml::find(data, "mcluster");
	std::vector<std::string> keys = mcluster_table.keys();

	// Check for unknown parameters
	for (const auto& key : keys) {
		bool found = false;
		for (const auto& valid : MCLUSTER_VALID_PARAMS) {
			if (key == valid) {
				found = true;
				break;
			}
		}
		if (!found) {
			std::string suggestion = suggestSimilarParam(key);
			std::string msg = "[mcluster] section: Unknown parameter '" + key + "'";
			if (!suggestion.empty()) {
				msg += ", did you mean '" + suggestion + "'?";
			}
			throw std::runtime_error(msg);
		}
	}

	// Parse N (number of stars)
	mc.N = config.getNestedOr<int>("mcluster", "N", 0);
	if (mc.N < 0) {
		throw std::runtime_error("[mcluster] section: Parameter 'N' must be non-negative, got: " +
			std::to_string(mc.N));
	}

	// Parse M (total mass in Msun)
	mc.M = config.getNestedOr<double>("mcluster", "M", 0.0);
	if (mc.M < 0) {
		throw std::runtime_error("[mcluster] section: Parameter 'M' must be non-negative, got: " +
			std::to_string(mc.M));
	}

	// Parse P (density profile)
	mc.P = config.getNestedOr<int>("mcluster", "P", 0);
	if (mc.P < -1 || mc.P > 3) {
		throw std::runtime_error("[mcluster] section: Parameter 'P' must be -1, 0, 1, 2, or 3, got: " +
			std::to_string(mc.P) + "\n  -1 = no density gradient\n  0 = Plummer (default)\n  " +
			"1 = King\n  2 = Subr et al.\n  3 = EFF/Nuker");
	}

	// Parse R (half-mass radius in pc)
	mc.R = config.getNestedOr<double>("mcluster", "R", 0.8);

	// Parse f (IMF selection)
	mc.f = config.getNestedOr<int>("mcluster", "f", 1);
	if (mc.f < 0 || mc.f > 2) {
		throw std::runtime_error("[mcluster] section: Parameter 'f' must be 0, 1, or 2, got: " +
			std::to_string(mc.f) + "\n  0 = single mass\n  1 = Kroupa (default)\n  2 = user-defined");
	}

	// Parse Z (metallicity)
	mc.Z = config.getNestedOr<double>("mcluster", "Z", 0.02);
	validateMclusterRange(mc.Z, 0.0001, 0.03, "Z");

	// Parse b (binary fraction)
	mc.b = config.getNestedOr<double>("mcluster", "b", 0.0);
	validateMclusterRange(mc.b, 0.0, 1.0, "b");

	// Parse e (stellar evolution epoch in Myr)
	mc.e = config.getNestedOr<double>("mcluster", "e", 0.0);
	if (mc.e < 0) {
		throw std::runtime_error("[mcluster] section: Parameter 'e' must be non-negative, got: " +
			std::to_string(mc.e));
	}

	// Parse generate_only (ABYSS-specific)
	mc.generate_only = config.getNestedOr<bool>("mcluster", "generate_only", false);
}

// Validate parsed mcluster configuration (N/M mutual exclusivity)
void validateMclusterConfig(const MclusterConfig& mc) {
	if (!mc.has_mcluster_section) {
		return;  // No validation needed if section not present
	}

	// N and M mutual exclusivity check
	if (mc.N == 0 && mc.M == 0.0) {
		throw std::runtime_error("[mcluster] section: Must specify either 'N' (star count) or 'M' (total mass in Msun)");
	}

	if (mc.N > 0 && mc.M > 0.0) {
		// Per CONTEXT.md: M takes precedence, warn user
		if (my_rank == ROOT) {
			std::cerr << "Warning: [mcluster] Both N and M specified; M takes precedence, ignoring N="
			          << mc.N << std::endl;
		}
	}

	// Minimum star count if using N
	if (mc.N > 0 && mc.N < 3) {
		throw std::runtime_error("[mcluster] section: Parameter 'N' must be at least 3 for N-body simulation, got: " +
			std::to_string(mc.N));
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

		// === McLuster parameters ===
		parseMclusterSection(config, mcluster_config);
		validateMclusterConfig(mcluster_config);

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
			if (mcluster_config.has_mcluster_section) {
				std::cout << "\n--- McLuster IC Generation ---\n";
				if (mcluster_config.M > 0) {
					std::cout << "M (total mass):    " << mcluster_config.M << " Msun\n";
				} else {
					std::cout << "N (star count):    " << mcluster_config.N << std::endl;
				}
				std::cout << "P (profile):       " << mcluster_config.P;
				switch (mcluster_config.P) {
					case -1: std::cout << " (none)"; break;
					case 0: std::cout << " (Plummer)"; break;
					case 1: std::cout << " (King)"; break;
					case 2: std::cout << " (Subr)"; break;
					case 3: std::cout << " (EFF/Nuker)"; break;
				}
				std::cout << std::endl;
				std::cout << "R (half-mass):     " << mcluster_config.R << " pc\n";
				std::cout << "f (IMF):           " << mcluster_config.f;
				switch (mcluster_config.f) {
					case 0: std::cout << " (single mass)"; break;
					case 1: std::cout << " (Kroupa)"; break;
					case 2: std::cout << " (user-defined)"; break;
				}
				std::cout << std::endl;
				std::cout << "Z (metallicity):   " << mcluster_config.Z << std::endl;
				std::cout << "b (binary frac):   " << mcluster_config.b << std::endl;
				std::cout << "e (epoch):         " << mcluster_config.e << " Myr\n";
				std::cout << "generate_only:     " << (mcluster_config.generate_only ? "yes" : "no") << std::endl;
			}
			std::cout << "==========================================\n\n";
		}

	} catch (const std::exception& e) {
		std::cerr << "Configuration TASK_ERROR: " << e.what() << std::endl;
		std::exit(EXIT_FAILURE);
	}
}
