#ifndef POPLR_HPP
#define POPLR_HPP

#include <iostream>
#include <vector>
#include <set>
#include <string>

#define MIN_VISITS 6

using Location = std::vector<double>;   // number of visits long
using Series = std::pair<std::vector<double>, std::vector<Location>>;   // <x[n], y[locs][n]>
using Eye = std::vector<Series>;

double PoPLR0(const Series series, const int max_perm_count, const float p_threshold, int verbosity_level, double aging_slope);
double PoPLR4(const Series series, const int max_perm_count, const float p_threshold, int verbosity_level, double aging_slope);
double PoPLR5(const Series series, const int max_perm_count, const float p_threshold, int verbosity_level, double aging_slope);

    // Utility print helpers for debugging
template<typename T>
void print(const std::string& name, const std::vector<T>& v) {
    std::cout << name << ": ";
    for (T x : v)
        std::cout << x << " ";
    std::cout << std::endl;
}
template<typename T>
void print(const std::string& name, const std::vector<std::vector<T>>& vv) {
    std::cout << name << ": ";
    for (auto v : vv) {
        std::cout << "\n\t";
        for (T x : v)
            std::cout << x << " ";
    }
    std::cout << std::endl;
}
template<typename T>
void print(const std::string& name, const std::set<T>& vv) {
    std::cout << name << ": ";
    for (auto v : vv) {
            std::cout << v << " ";
    }
    std::cout << std::endl;
}

template<typename T, typename Y>
void print(const std::string& name, const std::vector<std::pair<T, Y>>& v) {
    std::cout << name << ": ";
    for (auto p : v)
        std::cout << p.first << "," << p.second << " ";
    std::cout << std::endl;
}

#endif // POPLR_HPP