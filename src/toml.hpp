#ifndef ABYSS_TOML_HPP
#define ABYSS_TOML_HPP

#include <cerrno>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace toml {

class value {
public:
	enum class kind {
		kTable,
		kString,
		kInteger,
		kFloat,
		kBoolean
	};

	value() : kind_(kind::kTable), integer_(0), floating_(0.0), boolean_(false) {}

	static value make_table() {
		return value();
	}

	static value make_string(std::string v) {
		value out;
		out.kind_ = kind::kString;
		out.string_ = std::move(v);
		return out;
	}

	static value make_integer(long long v) {
		value out;
		out.kind_ = kind::kInteger;
		out.integer_ = v;
		return out;
	}

	static value make_float(double v) {
		value out;
		out.kind_ = kind::kFloat;
		out.floating_ = v;
		return out;
	}

	static value make_boolean(bool v) {
		value out;
		out.kind_ = kind::kBoolean;
		out.boolean_ = v;
		return out;
	}

	kind type() const {
		return kind_;
	}

	bool is_table() const {
		return kind_ == kind::kTable;
	}

	bool contains(const std::string& key) const {
		if (!is_table()) {
			return false;
		}
		return table_.find(key) != table_.end();
	}

	const value& at(const std::string& key) const {
		if (!is_table()) {
			throw std::runtime_error("TOML value is not a table");
		}
		auto it = table_.find(key);
		if (it == table_.end()) {
			throw std::runtime_error("TOML key not found: " + key);
		}
		return it->second;
	}

	value& at(const std::string& key) {
		if (!is_table()) {
			throw std::runtime_error("TOML value is not a table");
		}
		auto it = table_.find(key);
		if (it == table_.end()) {
			throw std::runtime_error("TOML key not found: " + key);
		}
		return it->second;
	}

	value& ensure_table(const std::string& key) {
		if (!is_table()) {
			throw std::runtime_error("TOML value is not a table");
		}
		auto it = table_.find(key);
		if (it == table_.end()) {
			auto res = table_.emplace(key, value::make_table());
			return res.first->second;
		}
		if (!it->second.is_table()) {
			throw std::runtime_error("TOML key is not a table: " + key);
		}
		return it->second;
	}

	void set(const std::string& key, value v) {
		if (!is_table()) {
			throw std::runtime_error("TOML value is not a table");
		}
		table_[key] = std::move(v);
	}

	const std::string& as_string() const {
		if (kind_ != kind::kString) {
			throw std::runtime_error("TOML value is not a string");
		}
		return string_;
	}

	long long as_integer() const {
		if (kind_ != kind::kInteger) {
			throw std::runtime_error("TOML value is not an integer");
		}
		return integer_;
	}

	double as_float() const {
		if (kind_ != kind::kFloat) {
			throw std::runtime_error("TOML value is not a float");
		}
		return floating_;
	}

	bool as_boolean() const {
		if (kind_ != kind::kBoolean) {
			throw std::runtime_error("TOML value is not a boolean");
		}
		return boolean_;
	}

private:
	kind kind_;
	std::string string_;
	long long integer_;
	double floating_;
	bool boolean_;
	std::unordered_map<std::string, value> table_;
};

namespace detail {

inline std::string trim(const std::string& input) {
	size_t start = 0;
	while (start < input.size() &&
	       std::isspace(static_cast<unsigned char>(input[start]))) {
		++start;
	}
	if (start == input.size()) {
		return std::string();
	}
	size_t end = input.size();
	while (end > start &&
	       std::isspace(static_cast<unsigned char>(input[end - 1]))) {
		--end;
	}
	return input.substr(start, end - start);
}

inline std::string strip_comments(const std::string& input) {
	bool in_string = false;
	bool escaped = false;
	for (size_t i = 0; i < input.size(); ++i) {
		const char c = input[i];
		if (escaped) {
			escaped = false;
			continue;
		}
		if (in_string && c == '\\') {
			escaped = true;
			continue;
		}
		if (c == '"') {
			in_string = !in_string;
			continue;
		}
		if (c == '#' && !in_string) {
			return input.substr(0, i);
		}
	}
	return input;
}

inline size_t find_unquoted(const std::string& input, char needle) {
	bool in_string = false;
	bool escaped = false;
	for (size_t i = 0; i < input.size(); ++i) {
		const char c = input[i];
		if (escaped) {
			escaped = false;
			continue;
		}
		if (in_string && c == '\\') {
			escaped = true;
			continue;
		}
		if (c == '"') {
			in_string = !in_string;
			continue;
		}
		if (!in_string && c == needle) {
			return i;
		}
	}
	return std::string::npos;
}

inline std::string remove_underscores(const std::string& input) {
	std::string out;
	out.reserve(input.size());
	for (char c : input) {
		if (c != '_') {
			out.push_back(c);
		}
	}
	return out;
}

inline std::string parse_basic_string(const std::string& input, size_t line_no) {
	if (input.size() < 2 || input.front() != '"' || input.back() != '"') {
		throw std::runtime_error("TOML parse error at line " +
			std::to_string(line_no) + ": invalid string literal");
	}
	std::string out;
	out.reserve(input.size() - 2);
	bool escaped = false;
	for (size_t i = 1; i + 1 < input.size(); ++i) {
		const char c = input[i];
		if (escaped) {
			switch (c) {
			case '"':
				out.push_back('"');
				break;
			case '\\':
				out.push_back('\\');
				break;
			case 'n':
				out.push_back('\n');
				break;
			case 't':
				out.push_back('\t');
				break;
			case 'r':
				out.push_back('\r');
				break;
			default:
				throw std::runtime_error("TOML parse error at line " +
					std::to_string(line_no) + ": unsupported escape sequence");
			}
			escaped = false;
			continue;
		}
		if (c == '\\') {
			escaped = true;
			continue;
		}
		out.push_back(c);
	}
	if (escaped) {
		throw std::runtime_error("TOML parse error at line " +
			std::to_string(line_no) + ": unterminated escape sequence");
	}
	return out;
}

inline value parse_value(const std::string& input, size_t line_no) {
	if (input.empty()) {
		throw std::runtime_error("TOML parse error at line " +
			std::to_string(line_no) + ": missing value");
	}
	if (input.front() == '"') {
		return value::make_string(parse_basic_string(input, line_no));
	}
	if (input == "true") {
		return value::make_boolean(true);
	}
	if (input == "false") {
		return value::make_boolean(false);
	}

	std::string cleaned = remove_underscores(input);
	const char first = cleaned.empty() ? '\0' : cleaned[0];
	if (std::isdigit(static_cast<unsigned char>(first)) || first == '+' || first == '-') {
		const bool has_dot = cleaned.find('.') != std::string::npos;
		const bool has_exp = cleaned.find_first_of("eE") != std::string::npos;
		char* end = nullptr;
		errno = 0;
		if (has_dot || has_exp) {
			const double val = std::strtod(cleaned.c_str(), &end);
			if (!end || *end != '\0' || errno == ERANGE) {
				throw std::runtime_error("TOML parse error at line " +
					std::to_string(line_no) + ": invalid float");
			}
			return value::make_float(val);
		}
		const long long val = std::strtoll(cleaned.c_str(), &end, 10);
		if (!end || *end != '\0' || errno == ERANGE) {
			throw std::runtime_error("TOML parse error at line " +
				std::to_string(line_no) + ": invalid integer");
		}
		return value::make_integer(val);
	}

	throw std::runtime_error("TOML parse error at line " +
		std::to_string(line_no) + ": unsupported value");
}

}  // namespace detail

inline value parse(const std::string& filename) {
	std::ifstream input(filename.c_str());
	if (!input) {
		throw std::runtime_error("Unable to open TOML file: " + filename);
	}

	value root = value::make_table();
	value* current = &root;
	std::string line;
	size_t line_no = 0;

	while (std::getline(input, line)) {
		++line_no;
		std::string cleaned = detail::trim(detail::strip_comments(line));
		if (cleaned.empty()) {
			continue;
		}

		if (cleaned.front() == '[') {
			if (cleaned.back() != ']') {
				throw std::runtime_error("TOML parse error at line " +
					std::to_string(line_no) + ": malformed table header");
			}
			std::string header = detail::trim(cleaned.substr(1, cleaned.size() - 2));
			if (header.empty()) {
				throw std::runtime_error("TOML parse error at line " +
					std::to_string(line_no) + ": empty table name");
			}
			current = &root;
			std::stringstream ss(header);
			std::string segment;
			while (std::getline(ss, segment, '.')) {
				std::string name = detail::trim(segment);
				if (name.empty()) {
					throw std::runtime_error("TOML parse error at line " +
						std::to_string(line_no) + ": empty table segment");
				}
				current = &current->ensure_table(name);
			}
			continue;
		}

		const size_t eq = detail::find_unquoted(cleaned, '=');
		if (eq == std::string::npos) {
			throw std::runtime_error("TOML parse error at line " +
				std::to_string(line_no) + ": expected '='");
		}
		std::string key = detail::trim(cleaned.substr(0, eq));
		std::string val = detail::trim(cleaned.substr(eq + 1));
		if (key.empty()) {
			throw std::runtime_error("TOML parse error at line " +
				std::to_string(line_no) + ": empty key");
		}
		if (val.empty()) {
			throw std::runtime_error("TOML parse error at line " +
				std::to_string(line_no) + ": empty value");
		}
		current->set(key, detail::parse_value(val, line_no));
	}

	return root;
}

inline const value& find(const value& data, const std::string& key) {
	return data.at(key);
}

template <typename T>
T find(const value& data, const std::string& key);

template <>
inline std::string find<std::string>(const value& data, const std::string& key) {
	return data.at(key).as_string();
}

template <>
inline int find<int>(const value& data, const std::string& key) {
	const value& v = data.at(key);
	if (v.type() == value::kind::kFloat) {
		const double val = v.as_float();
		if (std::floor(val) != val) {
			throw std::runtime_error("TOML value is not an integer");
		}
		return static_cast<int>(val);
	}
	return static_cast<int>(v.as_integer());
}

template <>
inline double find<double>(const value& data, const std::string& key) {
	const value& v = data.at(key);
	if (v.type() == value::kind::kInteger) {
		return static_cast<double>(v.as_integer());
	}
	return v.as_float();
}

template <>
inline bool find<bool>(const value& data, const std::string& key) {
	return data.at(key).as_boolean();
}

}  // namespace toml

#endif  // ABYSS_TOML_HPP
