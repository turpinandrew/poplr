#ifndef READ_UTILS_HPP
#define READ_UTILS_HPP

#include <vector>
#include <string>
#include "poplr.hpp"

std::vector<Eye> read_json(const std::string& filename);

std::vector<std::vector<long double>> read_csv(const std::string filename, int nrows, int ncols, bool row_sum_check);

#endif // READ_UTILS_HPP