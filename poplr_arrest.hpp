#ifndef POPLR_ARREST_HPP
#define POPLR_ARREST_HPP

#include <vector>
#include <set>

using Location = std::vector<double>;   // number of visits long
using Series = std::vector<Location>;   // number of locations long
using Eye = std::vector<Series>;

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

    // encoding ARREST values into dB values
double R() { return( -999); }
double Y(double db) { return(-100 - db); }
double O(double db) { return( -200 - db); }

    // decoding ARREST values into dB values
static bool is_green(double x) { return(x > -60); }
static bool is_yellow_orange(double x) { return(x < -60 && x > -999); }
static bool is_orange(double x) { return(x < -160 && x > -999); }
static bool is_red(double x) { return(x == -999); }
#define GET_DB_FROM_YELLOW(x) ((x) < -160 ? -(x) - 200 : -(x) - 100)

static double get_db(double x) { return (is_green(x) ? x : (is_red(x) ? -1 : GET_DB_FROM_YELLOW(x))); }


#endif // POPLR_ARREST_HPP