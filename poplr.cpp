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

typedef std::vector<double> Location;   // number of visits long
typedef std::vector<Location> Series;   // number of locations long
typedef std::vector<Series> Eye;

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

// @param series A matrix of dB values with `num_locations` rows and `num_visits` columns.
// @return The right boundaries of NaN blocks. Must end in num_visits.
std::vector<int> get_boundaries(const Series& series) {
    int num_visits = series[0].size();

    // Count NAs for each visit
    std::vector<int> num_NAs(num_visits, 0);
    for (int visit = 0; visit < num_visits; ++visit) {
        for (const auto& loc : series) {
            if (std::isnan(loc[visit])) {
                num_NAs[visit]++;
            }
        }
    }

        // invariant: block_boundaries[-1] has curr_num NAs to the left
    std::vector<int> block_boundaries;
    if (num_NAs[0] > 0) {
        int curr_num = num_NAs[0];
        block_boundaries.push_back(1);

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
    }

    block_boundaries.push_back(num_visits);

    return block_boundaries;
}

#define MIN_VISITS 6
/*
    Resize matrix to have num_locations rows and num_visits columns.
    Ignore any data before a column of NANs. (presumably at the start)
    Remove locations with less than MIN_VISITS in their row.

    @param series A matrix of dB values with at least `num_locations` rows and `num_visits` columns.
    @return Updated series matrix or empty Series if cannot meet requirements.
*/
bool is_nan (double d) { return std::isnan(d); }
Series preProcess(Series series, int num_locations, int num_visits) {
    if (num_visits < MIN_VISITS || num_locations < 2) 
        return Series();

    if (series.size() < num_locations || series[0].size() < num_visits)
        return Series();

    int min_visit_index = 0;

        // Skip over any columns at the start that are all NA to determine min_visit_index
    for (int visit = 0; visit < num_visits; visit++) {
        int loc;
        for (loc = 0; loc < num_locations; loc++)
            if (!std::isnan(series[loc][visit]))
                break;
        if (loc == num_locations)  // all locations have this visit as NAN
            min_visit_index = visit + 1;
    }

    if (num_visits - min_visit_index + 1 < MIN_VISITS) 
        return Series();

        // Prune each location to num_visits from min_visit_index
    for (int loc = 0; loc < num_locations; loc++) {
        if (min_visit_index > 0)
            series[loc].erase(series[loc].begin(), series[loc].begin() + min_visit_index);
        series[loc].resize(num_visits);
    }

        // Throw out locations with less than MIN_VISITS non-nan
    for (int loc = 0; loc < num_locations;) {
        int nan_count = count_if(series[loc].begin(), series[loc].end(), [](double x){return(isnan(x));});
        if (num_visits - nan_count < MIN_VISITS) {
                // swap last location up to here and shrink num_locations
            for (int visit = 0; visit < num_visits; visit++)
                series[loc][visit] = series[num_locations - 1][visit];
            num_locations--;
        } else
            loc++;
    }
    series.erase(series.begin() + num_locations, series.end());

    if (num_locations < 2) 
        return Series();

    return series;
}

/*
    @param series A matrix of dB values; rows are locations, columns are visits
    @return -1 For no valid data (too few visits or locations) else estimated probability of stability
*/
double PoPLR4(Series series, int perm_count) {
    int num_locations = series.size();
    int num_visits = series[0].size();

    if (perm_count > std::tgamma(num_visits + 1))
        perm_count = std::tgamma(num_visits + 1);

        // work out the x*, ybar, yds, y_skip values for each location
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
    std::vector<double> S;
    std::vector<double> p_vals(num_locations);
    std::vector<int> boundaries = get_boundaries(series);
//print("boundaries", boundaries);
    PermutationIterator pi(perm_count, boundaries);
    perm_count = 0;
    while (pi.hasNext()) {
        std::vector<int> perm = pi.next();
        perm_count++;
//print("perm:", perm);

//std::cout << "num_locations " << num_locations << std::endl;
        for (int loc = 0; loc < num_locations; loc++) {
            double beta = 0;
            for (int x = 1 ; x <= num_visits ; x++) {
                int ip = perm[x - 1];
                if (ip - 1 < y_skip[loc]) continue; // skip the NaN values
                beta += (x - xbar[loc]) * yds[loc][ip - 1];
//std::cout << "x= " << x << " ip= " << ip << " yd= " << yds[loc][ip - 1] << std::endl;
            }
            beta /= sumXdSq[loc];

            double alpha = ybar[loc] - beta * xbar[loc];

            double sse = 0, se, t;
            for (int x = 1 ; x <= num_visits ; x++) {
                int ip = perm[x - 1];
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
        double temp = 0.0;
        for (double p : p_vals)
            temp += -std::log(p);
        S.push_back(temp);
    }

//print("S", S);
    int less_than = 0, equal = 0;
    for(int i = 1; i < S.size(); i++) {
        if (S[i] < S[0])
            less_than++;
        else if (S[i] == S[0])
            equal++;
    }
//cout << "less_than " << less_than << " equal " << equal << " perm_count " << perm_count << endl;

    return 1 - (less_than + equal / 2.0) / (perm_count - 1);
}
    
/*
Read in JSON file that has series of VFs. Format is
[                   # whole lot
  [                   # begin eye 1
    [                   # begin repeat 1
      [20,19,18],         # location 1  (length "num_visits")
      [20,19,18]          # location 2
    ],
    [
        [21,19,18],     # begin repeat 2
        [21,19,18]        # location 1
    ]                     # location 2
  ],
  [ eye 2...]
]

Hand crafted parser with simple state.
Ignores ' ' and '\n'.
Converts any ".*" to nan.
Assumes no negative numbers.
*/
std::vector<Eye> read_json(const std::string& filename) {
    FILE *file = fopen(filename.c_str(), "r");

    std::vector<Eye> eyes;

    #define ADD_NUM {\
        eyes[eye][ser][loc].push_back(sign * num); \
        num = 0;\
        tens = 1;\
        sign = +1;\
    }

    int ch;
    int state = 0;  // 0 = null ; 1 = whole data; 2 = in eye ; 3 = in rep ; 4 = in row; 5 = in "..."
    int eye, ser, loc;
    double num = 0;
    double tens = 1;
    double sign = +1;
    while ((ch = getc(file)) != EOF) {
//std::cout << (char)ch << " " << "st " << state << " e " << eye << " s " << ser << " l " << loc << std::endl;
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
        } else if (ch == '-' && state == 4) {
            sign = -1;
        } else if (ch == ',' && state == 4) {
            ADD_NUM;
        } else if (ch == '"' && state == 4) {
            state++;
            num = std::numeric_limits<double>::quiet_NaN();
        } else if (ch == '"' && state == 5) {
            state--;  // ADD_NUM will be called at , or ]
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
    // Tests

    /*
    PermutationIterator pi(15, 3);
    while (pi.hasNext()) {
        std::cout << "p ";
        for (int x : pi.next())
            std::cout << x << " ";
        std::cout << std::endl;
    }
    */

    /*
    std::vector<std::vector<std::vector<double>>> tests = {
       {{  1,   2, 3},
        {  1,   2, 3},
        {  1,   2, 3},
        {  1,   2, 3}},

       {{  1,    2,     3, 4},
        {  1,    2,     3, 4},
        {"NA", "NA",    3, 4},
        {"NA", "NA",    3, 4},
        {"NA", "NA",  "NA", 4},
        {"NA", "NA",  "NA", 4}},

       {{  8,    7,    6,  5, 4, 3, 2, 1},
        {  8,    7,    7,  7, 7, 7, 7, 7},
        {"NA", "NA",   3,  4, 4, 4, 4, 4},
        {"NA", "NA", "NA", 4, 4, 4, 4, 4}},

    };


    for (std::vector<std::vector<double>> series : tests) {
            // Print input
        std::cout << "************" << std::endl;
        for (std::vector<double> loc : series) {
            for (double val : loc) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }

            // Print permutations
        std::cout << "-------------" << std::endl;
        auto bs = get_boundaries(series);
        print("bs", bs);
        PermutationIterator pi(10, bs);
        while (pi.hasNext()) {
            vector<int> perm = pi.next();
            print("", perm);
        }
    }

    std::cout << PoPLR4(tests[2], 5000) << std::endl;
    */


    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename.json>" << std::endl;
        return 1;
    }

    std::vector<Eye> eyes = read_json(argv[1]);

    std::cout << "eye,rep,visit,p" << std::endl;
    //#pragma omp parallel for
    for (int i_eye = 0 ; i_eye < eyes.size(); i_eye++) {
        for (int i_rep = 0 ; i_rep < eyes[i_eye].size(); i_rep++) {
            for (int visit = eyes[i_eye][i_rep][0].size(); visit >= 6; visit--) {
                Series series = preProcess(eyes[i_eye][i_rep], eyes[i_eye][i_rep].size(), visit);
                #pragma omp critical
                {
                std::cout << i_eye << "," << i_rep << "," << visit << ",";
                std::cout << PoPLR4(series, 5000) << std::endl;
                std::cout << endl;
                //print("Series:", series);
                }
            }
        }
    }
    
    return 0;
}
