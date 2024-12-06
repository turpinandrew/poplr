/*
** PoPLR implementation for ARREST

# Based on.
# Visual Field Progression in Glaucoma: Estimating the Overall Significance of Deterioration 
# with Permutation Analyses of Pointwise Linear Regression (PoPLR)
#   Neil O'Leary; Balwantray C. Chauhan; Paul H. Artes
#   Investigative Ophthalmology & Visual Science*
#   October 2012, Vol.53, 6776-6784. https://doi.org/10.1167/iovs.12-10049


 Only use permutations where NAs are at the start
 Only use locations with at least 6 visits.
 For locations that are:
    G.....G Apply PoPLR 
    G..Y..R Apply pr_yellow_red_sequence
    Y..Y..R Apply pr_yellow_red_sequence
    Y...Y|O Apply PoPLR using last ZEST value for Y/O
    O..Y/O  Apply PoPLR using last ZEST value for Y/O
    G...Y|O Apply pr_green_yellow

 @param series A matrix of dB values with `num_locations` rows and `num_visits` columns.
               (ordered in time so that col1 is earliest
               and the final column the most recent visit)
               Use NAN to signal missing data
 @param perm_count Number of permutations of series to get p-value for S

 @return Probability series is stable  (-1 if data has too few rows or visits)


/opt/homebrew/Cellar/gcc/14.2.0/bin/g++-14 -Wc++11-extensions -std= c++14 -I /opt/homebrew/include/ poplr_arrest.cpp
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
#include "poplr_arrest.hpp"
#include "arrestSeries.hpp"
#include "read_json.hpp"

#include <_stdio.h>
#include <functional>
//#include <omp.h>

using namespace std;

// @param series A matrix of dB values with `num_locations` rows and `num_visits` columns.
// @return The right boundaries of NaN blocks. Must end in num_visits.
vector<int> get_boundaries(const Series& series) {
    int num_visits = series[0].size();

    // Count NAs for each visit
    vector<int> num_NAs(num_visits, 0);
    for (int visit = 0; visit < num_visits; ++visit) {
        for (const auto& loc : series) {
            if (isnan(loc[visit])) {
                num_NAs[visit]++;
            }
        }
    }

        // invariant: block_boundaries[-1] has curr_num NAs to the left
    vector<int> block_boundaries;
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
    Resize matrix to have at least 2 rows and min(num_visits, max loc seq length) columns.
    Remove locations with less than MIN_VISITS non-NaN in their row.

    @param series A matrix of dB values with at least `num_locations` rows and up to `num_visits` columns.
    @param num_locations Number of rows in series
    @param num_visits Max number of columns in series

    @return Updated series matrix or empty Series if cannot meet requirements.
*/
bool is_nan (double d) { return isnan(d); }
Series preProcess(Series series, int num_visits) {
    int n = series.size();

    if (num_visits < MIN_VISITS || n < 2)
        return Series();

        // Cap each location at num_visits in size (including NaNs)
    for (int loc = 0; loc < n; loc++)
        if (series[loc].size() > num_visits)
            series[loc].erase(series[loc].begin() + num_visits, series[loc].end());

        // Throw out locations with less than MIN_VISITS non-nan
    for (int loc = 0; loc < n; ) {
        int nan_count = count_if(series[loc].begin(), series[loc].end(), [](double x){return(isnan(x));});
        if (series[loc].size() - nan_count < MIN_VISITS) {
                // swap last location up to here and shrink num_locations
            for (int visit = 0; visit < num_visits; visit++)
                series[loc][visit] = series[n - 1][visit];
            n--;
        } else
            loc++;
    }
    if (n < 2) 
        return Series();
    else 
        series.erase(series.begin() + n, series.end()); // throw away junk swapped to the end 

    return series;
}

/*
    This version partitions columns based on NaNs at the start of rows.
    ASSUMES: locations have at least 6 non-NaN values in their row.

    @param series A matrix of dB values; rows are locations, columns are visits
    @param perm_count Number of permutations of series to get p-value for S

    @return Estimated probability of stability of the series
*/
double PoPLR4(Series series, int perm_count) {
    int num_locations = series.size();
    if (num_locations < 2) 
        return 1.0;

    int num_visits = series[0].size();

    if (perm_count > tgamma(num_visits + 1))
        perm_count = tgamma(num_visits + 1);

        // work out the x*, ybar, yds, y_skip values for each location
    vector<int> n(num_locations);
    vector<long double> xbar(num_locations);
    vector<long double> sumXdSq(num_locations);
    vector<long double> sqrtSumXdSq(num_locations);
    vector<int> y_skip(num_locations);    // number of leading NaNs
    vector<long double> ybar(num_locations);
    vector<vector<long double>> yds(num_locations);  // for all including leading NaNs

    for (int loc = 0; loc < num_locations; loc++) {
        n[loc] = 0;
        y_skip[loc] = 0;
        ybar[loc]= 0.0;
        for (double y : series[loc])
            if (isnan(y))
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

        sqrtSumXdSq[loc] = sqrt((long double)sumXdSq[loc]);
    }
//print("n", n);
//print("xbar", xbar);
//print("ybar", ybar);
//print("sumXdSq", sumXdSq);
//print("sqrtSumXdSq", sqrtSumXdSq);

    // Now compute all the p-values for each loc in each appropriate permutation
    vector<long double> S;
    vector<long double> p_vals(num_locations);
    vector<int> boundaries = get_boundaries(series);
//print("boundaries", boundaries);
    PermutationIterator pi(perm_count, boundaries);
    perm_count = 0;
    while (pi.hasNext()) {
        vector<int> perm = pi.next();
        perm_count++;
//print("perm:", perm);

//cout << "num_locations " << num_locations << endl;
        for (int loc = 0; loc < num_locations; loc++) {
            long double beta = 0;
            for (int x = 1 ; x <= num_visits ; x++) {
                int ip = perm[x - 1];
                if (ip - 1 < y_skip[loc]) continue; // skip the NaN values
                beta += (x - xbar[loc]) * yds[loc][ip - 1];
//cout << "x= " << x << " ip= " << ip << " yd= " << yds[loc][ip - 1] << endl;
            }
            beta /= sumXdSq[loc];

            long double alpha = ybar[loc] - beta * xbar[loc];

            long double sse = 0, se, t;
            for (int x = 1 ; x <= num_visits ; x++) {
                int ip = perm[x - 1];
                if (ip - 1 < y_skip[loc]) continue; // skip the NaN values
                sse += (series[loc][ip - 1] - alpha - beta * x) * 
                       (series[loc][ip - 1] - alpha - beta * x);
            }
//cout << "beta " << beta << " alpha " << alpha << " sse " << sse;

            if (sse < 1e-10) {
                if (beta >= 0) {
                    p_vals[loc] = 1.0; // t -> Inf
                } else {
                    p_vals[loc] = 0.0; // t -> -Inf
                }
            } else {
                se = sqrt(sse / (n[loc] - 2)) / sqrtSumXdSq[loc];
                t = beta / se;
                boost::math::students_t dist(n[loc] - 2);
                p_vals[loc] = cdf(dist, t);
            }
//cout << " se " << se << " t " << t;
//cout << " p " << p_vals[loc] << endl;
        }

//cout << "ps: "; for (double p : p_vals) cout << p << " "; cout << endl;
            // Now compute the S value for this permutation
        long double temp = 0.0;
        for (long double p : p_vals)
            temp += -log(p);
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
    This version breaks series into collections of rows that have the 
    same number of leading NaNs. The end result is the min of PoPLR
    on each collection.

    @param series A matrix of dB values; rows are locations, columns are visits
    @param perm_count Number of permutations of series to get p-value for S

    @return Estimated probability of stability of the series
*/
double PoPLR5(Series series, int perm_count) {
    int num_locations = series.size();

        // Work out number of leading NaNs (skips)
    vector<int> nan_count(num_locations);    // number of leading NaNs
    set<int> nan_count_set = set<int>();

    for (int loc = 0; loc < num_locations; loc++) {
        nan_count[loc] = 0;
        for (double y : series[loc])
            if (isnan(y))
                nan_count[loc]++;
        nan_count_set.insert(nan_count[loc]);
    }

    double minP = 1.0; // default stable
    for (int count : nan_count_set) {
        vector<vector<double>> sub_series;
        for (int loc = 0; loc < num_locations; loc++) {
            if (nan_count[loc] == count)
                sub_series.push_back(series[loc]);
        }
        int num_visits = sub_series[0].size() - count;

        if (num_visits < MIN_VISITS)
            continue;

        if (perm_count > tgamma(num_visits + 1))
            perm_count = tgamma(num_visits + 1);

        if (sub_series.size() > 0) {
            double p = PoPLR4(sub_series, perm_count);
            if (p < minP) 
                minP = p;
        }
    }
    return minP;
}
    
    // print usage and exit
void usage() {
    cerr << "Usage: poplr_arrest part_style db_treatment filename.json" << endl;
    cerr << "       where" << endl;
    cerr << "           part_style   is H[orizontal] | V[ertical] grouping" << endl;
    cerr << "           db_treatment is 0 = Not arrest (no converting values, use every line)" << endl;
    cerr << "                           1 = ARREST with probs (convert values, apply PLR to Green lines, special probs to others)" << endl;
    cerr << "                           2 = ARREST no probs (convert values, only all-Green locations)" << endl;
    cerr << "        Horizontal partitioning groups rows/locations with the same number of visits." << endl;
    cerr << "        Vertical partitioning groups columns/visits so that locations all have the same number of NAs." << endl;
    cerr << "        Permuting happens within groups and group order remains the same." << endl;
    exit(-1);
}

int main(int argc, char *argv[]) {
    // Tests

    /*
    PermutationIterator pi(15, 3);
    PermutationIterator pi(0, 5000);
    while (pi.hasNext()) {
        cout << "p ";
        for (int x : pi.next())
            cout << x << " ";
        cout << endl;
    }
    return 0;
    */

    /*
    vector<vector<vector<double>>> tests = {
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


    for (vector<vector<double>> series : tests) {
            // Print input
        cout << "************" << endl;
        for (vector<double> loc : series) {
            for (double val : loc) {
                cout << val << " ";
            }
            cout << endl;
        }

            // Print permutations
        cout << "-------------" << endl;
        auto bs = get_boundaries(series);
        print("bs", bs);
        PermutationIterator pi(10, bs);
        while (pi.hasNext()) {
            vector<int> perm = pi.next();
            print("", perm);
        }
    }

    cout << PoPLR4(tests[2], 5000) << endl;
    */

   /*
   Series a7_stable_21_e101_r1 = {
    {22.0682,17.0677,21.2230,11.2733,24.7591,22.0360,28.7629,23.6009,17.0347,18.6527,26.6499,16.1782},
    {25.7968,26.8750,27.2849,24.4909,25.7968,28.1244,25.2500,22.0360,25.9485,24.7572,26.1599,26.5921},
    {22.9376,16.1782,30.1011,30.1011,20.8235,27.2849,33.1565,25.9485,18.2311,24.4909,23.2762,18.2311},
    {22.7690,25.9485,24.4909,30.1398,24.6636,25.9485,26.5921,24.6636,25.8902,25.9485,17.0347,18.2311},
    {25.2500,25.9485,24.4909,28.1244,33.1565,24.4909,26.0985,30.1011,28.7629,25.7018,26.1599,30.1011},
    {22.8387,28.1244,22.7690,26.1599,27.6529,18.2311,26.5921,28.1244,18.2311,18.2311,18.2311,28.7629},
    {23.6009,27.2849,15.8378,22.8387,23.9373,23.2351,30.9964,25.2500,30.9964,23.2351,22.7690,30.9964},
    {27.2849,12.3794,24.9259,15.9625,13.1613,16.7070,12.6871,20.5077,24.4909,17.0677,30.9964,25.1069},
    {28.7629,17.1913,30.9964,17.0677,23.6009,20.7474,23.6009,24.6636,17.0347,15.7726,10.7447,17.0677},
    {17.0878, 7.7520,24.4909,22.9436,18.2311,11.8247,17.5191,18.2311,24.0696,24.6636, 9.0223,19.6088},
    {22.8387,15.7726,17.0677,16.3946,17.0677,22.1201,14.8890,26.6499,24.0696,18.6527,25.9485,22.8387},
    {24.0696,22.8387,22.8387,24.0696,20.6300,19.6088,24.4909,14.5603,25.2500,25.9485,23.6866,22.8387},
    {15.4363,18.9521,16.1782,15.6878,18.2311,16.0244,24.6636,26.9530,26.5921,26.0395,16.1782,26.5921},
    {16.3946,28.7629,24.7572,21.2230,16.1782,21.2230,24.4909,24.4909,23.8218,25.1069,23.6009,22.0360},
    {25.2500,25.2500,20.7474,22.7017,27.6529,22.8219,21.2230,19.5588,23.6009,23.8218,22.0360,19.6088},
    {19.5588,16.3872,19.5588,14.8890,28.7629,16.5376,24.6636,25.9485,23.2351,21.6479,23.6009,20.0873},
    {17.0347,28.7629,26.5921,30.1011,30.1011,32.4865,20.5077,14.3871,17.1865,18.0365,27.6529,18.2311},
    {18.2311,26.0985,25.7018,24.4909,28.1301,16.1782,20.2171,18.3376,14.3871,27.6529,16.3946,20.2164},
    {32.4865,22.8387,19.6088,26.5921,24.6636,27.2849,25.9485,25.2500,26.6499,19.6163,27.2849,24.9259},
    {24.7572,27.2849,25.1069,26.5921,24.7572,23.6009,25.2500,26.5921,28.1244,27.6529,22.1028,24.4909},
    {25.9485,17.5191,16.1782,18.3376,18.8111,26.1599,23.6866,20.3476,28.7629,23.2351,24.4909,15.2355},
    {18.2311,17.0677,16.1782,17.0442,24.4909,27.6529,18.2311,15.5183,20.9824,23.9373,11.8658,27.6529},
    {15.2355,17.0677,16.3946,26.6499,18.2311,18.3376,24.0696,20.8786,30.5894,21.2230,14.3871,19.6088},
    {27.6529,28.1244,24.0696,16.1782,28.7629,20.0873,30.1011,12.8327,19.6579,17.0372,19.6088,27.2849},
    {16.3387,22.1028,23.2762,33.1565,23.2351,15.2355,26.9530,24.4909,19.6088,22.8387,25.8902,19.6088},
    {30.5894,30.1011,25.9485,26.9530,30.1011,24.9259,24.4909,25.9485,28.7629,28.7629,25.9485,25.2500},
    {24.5945,25.7968,27.6529,27.2849,22.8387,21.2230,26.5921,24.1527,28.1244,23.6009,25.9485,28.1244},
    {23.2351,30.1011,24.4909,27.6529,25.2500,17.0347,26.1599,17.0677,23.8218,22.9436,23.6009,19.6088},
    {15.8378,15.2355,25.7968,23.6009,16.7194,22.8387,25.2500,16.1782,18.2311,21.9394,16.3171,21.2230}};

    ArrestSeries *testAS = new ArrestSeries(a7_stable_21_e101_r1, true);
    Series series = preProcess(testAS->get_series(), testAS->get_series().size(), 12);
    cout << "poplr: " << PoPLR4(series, 5000);
    cout << "other_ps: " << testAS->get_min_p();
    cout << endl;
    print("as", testAS->get_series());

    return -1;
    */
    

    if (argc != 4)
        usage();

    function<double(Series, int)> PoPLR;
    if (argv[1][0] == 'H' || argv[1][0] == 'h')
        PoPLR = PoPLR5;
    else if (argv[1][0] == 'V' || argv[1][0] == 'v')
        PoPLR = PoPLR4;
    else {
        cerr << "Invalid part_type." << endl;
        usage();
    }

    ArrestSeries::ArrestProcessType arrest;
    switch (stoi(argv[2])) {
        case 0: arrest = ArrestSeries::ArrestProcessType::NOT_ARREST; break;
        case 1: arrest = ArrestSeries::ArrestProcessType::ARREST_WITH_PROBS; break;
        case 2: arrest = ArrestSeries::ArrestProcessType::ARREST_NO_PROBS; break;
        default: {
            cerr << "Invalid db_type." << endl;
            usage();
        }
    }

    vector<Eye> eyes = read_json(argv[3]);

    cout << "eye,rep,visit,p" << endl;
    //#pragma omp parallel for
    for (int i_eye = 0 ; i_eye < eyes.size(); i_eye++) {
        for (int i_rep = 0 ; i_rep < eyes[i_eye].size(); i_rep++) {
            for (int visit = eyes[i_eye][i_rep][0].size(); visit >= 6; visit--) {
                Series series = preProcess(eyes[i_eye][i_rep], visit);
                ArrestSeries *as = new ArrestSeries(series, arrest);
                double poplr_p = PoPLR(as->get_series(), 5000);
                double other_ps = as->get_min_p();
                #pragma omp critical
                {
                cout << i_eye << "," << i_rep << "," << visit << ",";
                cout << (poplr_p < other_ps ? poplr_p :other_ps);
                cout << endl;
                //print("Series:", series);
                }
            }
        }
    }
    
    return 0;
}
