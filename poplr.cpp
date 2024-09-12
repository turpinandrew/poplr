/*
** PoPLR implemntation in c.

# Based on.
# Visual Field Progression in Glaucoma: Estimating the Overall Significance of Deterioration 
# with Permutation Analyses of Pointwise Linear Regression (PoPLR)
#   Neil O'Leary; Balwantray C. Chauhan; Paul H. Artes
#   Investigative Ophthalmology & Visual Science*
#   October 2012, Vol.53, 6776-6784. https://doi.org/10.1167/iovs.12-10049


 Only use permutations where NAs are at the start
 Only use locations with at least 6 visits.
 If arrest = TRUE, exclude any locations with <= 16 at penultimate visit.
 ASSUMES NAs are at the start.

 @param series A matrix of dB values with `num_locations` rows and `num_visits` columns.
               (ordered in time so that col1 is earliest
               and the final column the most recent visit)
               Use NAN to signal missing data
 @param perm_count Number of permutations of series to get p-value for S

 @return Probability series is stable  (-1 if data has too few rows or visits)

*/
#include <vector>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <random>
#include <fstream>
#include <boost/math/distributions/students_t.hpp>
#include <boost/container_hash/hash.hpp>
#include "permutation.hpp"
//#include <omp.h>

typedef std::unordered_set<std::vector<int>, boost::hash<std::vector<int>>> PermSet;

typedef std::vector<double> Location;   // number of visits long
typedef std::vector<Location> Series;   // number of locations long
typedef std::vector<Series> Eye;

    // Utility print helpers
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

// Function to get perm_count permutations of 1..n
// Do include "first" perm 1..n
PermSet get_perms(int n, int perm_count) {
    if (perm_count > std::tgamma(n + 1))
        perm_count = std::tgamma(n + 1);

    std::unordered_set<std::vector<int>, boost::hash<std::vector<int>>> perms;

    std::random_device rd;
    std::mt19937 g{ rd() };

    std::vector<int> base_perm(n);
    std::iota(base_perm.begin(), base_perm.end(), 1);

    perms.insert(base_perm);
    while (perms.size() < perm_count) {
        std::shuffle(base_perm.begin(), base_perm.end(), g);
        perms.insert(base_perm);
    }

    return perms;
}

// Function to get the boundaries of NA blocks
std::vector<int> get_boundaries(const std::vector<std::vector<double>>& series) {
    int num_visits = series[0].size();
    std::vector<int> num_NAs(num_visits, 0);

    // Count NAs for each visit
    for (int visit = 0; visit < num_visits; ++visit) {
        for (const auto& loc : series) {
            if (std::isnan(loc[visit])) {
                num_NAs[visit]++;
            }
        }
    }

    if (std::all_of(num_NAs.begin(), num_NAs.end(), [](int n) { return n == 0; })) {
        return {};
    }

    int curr_num = num_NAs[0];
    std::vector<int> block_boundaries = {1};
    for (int i = 1; i < num_visits; i++) {
        if (num_NAs[i] == 0) {
            break;
        }

        int n = block_boundaries.size();
        if (num_NAs[i] == curr_num) {
            block_boundaries[n - 1] = i + 1;
        } else {
            block_boundaries.push_back(i + 1);
            curr_num = num_NAs[i];
        }
    }
    return block_boundaries;
}

/*
 @param start Index of visit to start perm - 1
 @param end   Index of visit to end perm
 @return Perms of (start, end] = 1:(end - start) + start. Resample to get exactly perm_count perms.
*/
std::vector<std::vector<int>> gp(int start, int end, int perm_count) {
    PermSet p = get_perms(end - start, perm_count);

    std::vector<std::vector<int>> pp(p.begin(), p.end());
    for (int i = 0; i < pp.size(); i++)
        for (int j = 0; j < pp[i].size(); j++)
            pp[i][j] += start;

    while (pp.size() < perm_count) {
        pp.push_back(pp[rand() % pp.size()]);
    }

    return pp;
};

// Function to get permutations based on NA block boundaries
// Perm 1...n should be the first element in returned vector
std::vector<std::vector<int>> get_perms4(const std::vector<std::vector<double>>& series, int perm_count) {
    auto n = series[0].size();

    std::vector<std::vector<int>> res;

        // Look for block boundaries defined by NaNs
    std::vector<int> bb = get_boundaries(series);

    if (bb.empty()) {
        res = gp(0, n, perm_count);
    } else {
//print("bb", bb);
            // get permutations for each block
        std::vector<std::vector<std::vector<int>>> perm_blocks(bb.size());
        for (size_t i = 0; i < bb.size(); i++)
            perm_blocks[i] = gp(i == 0 ? 0 : bb[i - 1], bb[i], perm_count);

        if (bb.back() != n)
            perm_blocks.push_back(gp(bb.back(), n, perm_count));
//for (size_t i = 0; i < perm_blocks.size(); i++)
//print("pb" + std::to_string(i), perm_blocks[i]);

        res = perm_blocks[0];
        for (size_t i = 1; i < perm_blocks.size(); i++) {
            for (int j = 0; j < perm_count; ++j) {
                res[j].insert(res[j].end(), perm_blocks[i][j].begin(), perm_blocks[i][j].end());
            }
        }
    }

        // Find the 'first' perm 1..n and move it to the front
    auto v = std::vector<int>(n);
    std::iota(v.begin(), v.end(), 1);
    auto it = std::find(res.begin(), res.end(), v);
    assert (it != res.end()); 
    iter_swap(res.begin(), it);

    return res;
}

// @return -1 For no valid data (too few visits or locations)
double PoPLR4(Series& series, int num_locations, int num_visits, int perm_count) {
    if (perm_count > std::tgamma(num_visits + 1))
        perm_count = std::tgamma(num_visits + 1);

        // Throw out locations with 5 or less visits
    for (int loc = 0; loc < num_locations;) {
        int na_count = 0;
        for (int visit = 0; visit < num_visits; visit++)
            if (std::isnan(series[loc][visit]))
                na_count++;
        if (num_visits - na_count < 6) {
                // swap last location up to here and shrink num_locations
            for (int visit = 0; visit < num_visits; visit++)
                series[loc][visit] = series[num_locations - 1][visit];
            num_locations--;
        } else
            loc++;
    }

    if (num_locations <= 1) return -1;

              // Prune each location to num_visits entries
    for (int loc = 0; loc < num_locations;)
        series[loc].resize(num_visits);

        // skip over any columns at the start that are all NA
    int min_visit_index = 0;
    for (int visit = 0; visit < num_visits; visit++) {
        int loc;
        for (loc = 0; loc < num_locations; loc++)
            if (!std::isnan(series[loc][visit]))
                break;
        if (loc == num_locations)
            min_visit_index = visit + 1;
    }

    if (min_visit_index >= num_visits) return -1;

        // get permutation vectors for largest n
    std::vector<std::vector<int>> perms = get_perms4(series, perm_count - 1);
    perm_count = perms.size() - 1;  // does not include first 1:n

        // work out the x and ybar values for each location
    std::vector<int> n(num_locations);
    std::vector<double> xbar(num_locations);
    std::vector<double> sumXdSq(num_locations);
    std::vector<double> sqrtSumXdSq(num_locations);
    std::vector<int> y_skip(num_locations);    // number of leading NaNs
    std::vector<double> ybar(num_locations);
    std::vector<std::vector<double>> yds(num_locations);  // for all including leading NaNs

    for (int loc = 0; loc < num_locations; loc++) {
        n[loc] = 0;
        y_skip[loc] = 0;
        ybar[loc]= 0.0;
        for (double y : series[loc])
            if (std::isnan(y))
                y_skip[loc]++;
            else {
                ybar[loc] += y;
                n[loc]++;
            }

        ybar[loc] /= (double)n[loc];

        for (double y : series[loc])
            yds[loc].push_back(y - ybar[loc]);

        xbar[loc] = 0;
        for (int x = 1 ; x <= n[loc]; x++)
            xbar[loc] += x;
        xbar[loc] /= (double)n[loc];

        for (int x = 1 ; x <= n[loc]; x++)
            sumXdSq[loc] += (x - xbar[loc]) * (x - xbar[loc]);

        sqrtSumXdSq[loc] = std::sqrt(sumXdSq[loc]);
    }
//print("n", n);
//print("xbar", xbar);
//print("ybar", ybar);
//print("sumXdSq", sumXdSq);
//print("sqrtSumXdSq", sqrtSumXdSq);

    // Now compute all the p-values for each loc in each appropriate permutation
    std::vector<double> S(perm_count);
    std::vector<double> p_vals(num_locations);
    for (int perm = 0; perm <= perm_count; perm++) {
//print("perm:", perms[perm]);

//std::cout << "num_locations " << num_locations << std::endl;
        for (int loc = 0; loc < num_locations; loc++) {
            double beta = 0;
            for (int x = 1 ; x <= num_visits ; x++) {
                int ip = perms[perm][x - 1];
                if (ip - 1 < y_skip[loc]) continue; // skip the NaN values
                beta += (x - xbar[loc]) * yds[loc][ip - 1];
//std::cout << "x= " << x << " ip= " << ip << " yd= " << yds[loc][ip - 1] << std::endl;
            }
            beta /= sumXdSq[loc];

            double alpha = ybar[loc] - beta * xbar[loc];

            double sse = 0, se, t;
            for (int x = 1 ; x <= num_visits ; x++) {
                int ip = perms[perm][x - 1];
                if (ip - 1 < y_skip[loc]) continue; // skip the NaN values
                sse += (series[loc][ip - 1] - alpha - beta * x) * 
                       (series[loc][ip - 1] - alpha - beta * x);
            }
//std::cout << "beta " << beta << " alpha " << alpha << " sse " << sse;

            if (sse < 1e-10) {
                if (beta >= 0) {
                    p_vals[loc] = 1.0; // t -> Inf
                } else {
                    p_vals[loc] = 0.0; // t -> -Inf
                }
            } else {
                se = std::sqrt(sse / (n[loc] - 2)) / sqrtSumXdSq[loc];
                t = beta / se;
                boost::math::students_t dist(n[loc] - 2);
                p_vals[loc] = cdf(dist, t);
            }
//std::cout << " se " << se << " t " << t;
//std::cout << " p " << p_vals[loc] << std::endl;
        }

//std::cout << "ps: "; for (double p : p_vals) std::cout << p << " "; std::cout << std::endl;
            // Now compute the S value for this permutation
        S[perm] = 0.0;
        for (double p : p_vals)
            S[perm] += -std::log(p);
    }

//print("S", S);
    int less_than = 0, equal = 0;
    for(int i = 1; i < S.size(); i++) {
        if (S[i] < S[0])
            less_than++;
        else if (S[i] == S[0])
            equal++;
    }

    return 1 - (less_than + equal / 2.0) / perm_count;
}

// Function to read CSV file (generated by co-pilot) 11 Sep 2024
std::vector<std::vector<double>> readCSV(const std::string& filename) {
    std::vector<std::vector<double>> data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return data;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::stringstream ss(line);
        std::string value;

        while (std::getline(ss, value, ',')) {
            try {
                row.push_back(std::stod(value));
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid number: " << value << std::endl;
            } catch (const std::out_of_range& e) {
                std::cerr << "Number out of range: " << value << std::endl;
            }
        }

        data.push_back(row);
    }

    file.close();
    return data;
}
    
std::vector<Eye> read_json(const std::string& filename) {
    FILE *file = fopen(filename.c_str(), "r");

    std::vector<Eye> eyes;

    #define ADD_NUM {\
        eyes[eye][ser][loc].push_back(num); \
        num = 0;\
        tens = 1; \
    }

    int ch;
    int state = 0;  // 0 = null ; 1 = whole data; 2 = in eye ; 3 = in rep ; 4 = in row
    int eye, ser, loc;
    double num = 0;
    double tens = 1;
    while ((ch = getc(file)) != EOF) {
        if (ch == '[') {
            switch (state) { 
            case 0: { state++; eye = -1; break; }
            case 1: { state++; eye++; eyes.push_back(Eye()); ser = -1; break; }
            case 2: { state++; ser++; eyes[eye].push_back(Series()); loc = -1; break; }
            case 3: { state++; loc++; eyes[eye][ser].push_back(Location()); break; }
            default: {
                std::cout << "parse ERROR [" << std::endl;
                return(std::vector<Eye>());
            }};
        } else if (ch == ']') {
            switch (state) { 
            case 0: { std::cout << "parse ERROR [" << std::endl; break;}
            case 1: { state--; break; }
            case 2: { state--; break; }
            case 3: { state--; break; }
            case 4: { state--; ADD_NUM; break; }
            };
        } else if (ch == ',' && state == 4) {
            ADD_NUM;
        } else if (ch == ',') {
            continue;
        } else if (ch == ' ' || ch == '\n') {
            continue;
        } else if (ch == '.') {
            tens = 0.1;
        } else {
            if (tens < 1) {
                num += ((double)(ch - '0')) * tens;
                tens /= 10.0;
            } else
                num = num * 10.0 + (double)(ch - '0');
        }
    }

    return eyes;
}


int main(int argc, char *argv[]) {
    
    PermutationIterator pi(4, 4);
    while (pi.hasNext()) {
        std::cout << "p ";
        for (int x : pi.next())
            std::cout << x << " ";
        std::cout << std::endl;
    }
    
    // Tests
    /*
    std::vector<std::vector<std::vector<double>>> tests = {
       {{  1,   2, 3},
        {  1,   2, 3},
        {  1,   2, 3},
        {  1,   2, 3}},

       {{  1,   2,    3, 4},
        {  1,   2,    3, 4},
        {NAN, NAN,    3, 4},
        {NAN, NAN,    3, 4},
        {NAN, NAN,  NAN, 4},
        {NAN, NAN,  NAN, 4}},

       {{  8,   7,    6, 5, 4, 3, 2, 1},
        {  8,   7,    7, 7, 7, 7, 7, 7},
        {NAN, NAN,    3, 4, 4, 4, 4, 4},
        {NAN, NAN,  NAN, 4, 4, 4, 4, 4}},

    };

    for (std::vector<std::vector<double>> series : tests) {
        std::vector<std::vector<int>> perms = get_perms4(series, 10);

        // Print permutations
        std::cout << "************" << std::endl;
        for (const auto& perm : perms) {
            for (int val : perm) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }

    std::cout << PoPLR4(tests[2], tests[2].size(), tests[2][1].size(), 5000) << std::endl;
    */


    return(-1);
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename.json>" << std::endl;
        return 1;
    }

    std::vector<Eye> eyes = read_json(argv[1]);

    std::vector<std::vector<std::vector<double>>> results(126, 
         std::vector<std::vector<double>>(13, std::vector<double>(100)));
 
     #pragma omp parallel for
     //for (int i_eye = 0 ; i_eye < eyes.size(); i_eye++) {
     //for (int i_rep = 0 ; i_rep < eyes[i_eye].size(); i_rep++) {
     //for (int visit = eyes[i_eye][i_rep][0].size(); visit >= 6; visit--) {
     for (int i_eye = 0 ; i_eye < 1; i_eye++) {
     for (int i_rep = 0 ; i_rep < 1; i_rep++) {
     for (int visit = eyes[i_eye][i_rep][0].size(); visit >= 12; visit--) {
         results[i_eye][visit][i_rep] = PoPLR4(eyes[i_eye][i_rep], eyes[i_eye][i_rep].size(), visit, 5000);
     }}}

     for (int i_eye = 0 ; i_eye < eyes.size(); i_eye++) {
     for (int i_rep = 0 ; i_rep < eyes[i_eye].size(); i_rep++) {
     for (int visit = eyes[i_eye][i_rep][0].size(); visit >= 6; visit--) {
         std::cout << i_eye << "," << i_rep << "," << visit << ",";
         std::cout << results[i_eye][visit][i_rep];
         std::cout << std::endl;
     }}}
    
    return 0;
}
