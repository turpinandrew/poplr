/* PoPLR implementation for ARREST

Based on.
 Visual Field Progression in Glaucoma: Estimating the Overall Significance of Deterioration 
 with Permutation Analyses of Pointwise Linear Regression (PoPLR)
   Neil O'Leary; Balwantray C. Chauhan; Paul H. Artes
   Investigative Ophthalmology & Visual Science*
   October 2012, Vol.53, 6776-6784. https://doi.org/10.1167/iovs.12-10049

*/
#ifndef POPLR_ARREST_HPP
#define POPLR_ARREST_HPP

#include <vector>
#include <set>
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
#include "types.hpp"
#include "arrestSeries.hpp"

#include <_stdio.h>
#include <functional>
//#include <omp.h>

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

    @param series A matrix of dB values; rows are locations, columns are visits. First row is x values.
    @param perm_count Number of permutations of series to get p-value for S
    @param slope_limit Only include p-values in S statistic when PLR slope is <= slope_limit

    @return Estimated probability of stability of the series
*/
double PoPLR4(Series series, int perm_count, double slope_limit) {
    int num_locations = series.size() - 1; // leading xs
    if (num_locations < 2) 
        return 1.0;

    vector<double> xs = get_xs(&series);

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
        xbar[loc] = 0;
        for (int i = 0 ; i < num_visits; i++) {
            double y = series[loc][i];
            if (isnan(y))
                y_skip[loc]++;
            else {
                ybar[loc] += y;
                n[loc]++;
                xbar[loc] += xs[i];
            }
        }

        ybar[loc] /= (double)n[loc];
        xbar[loc] /= (double)n[loc];

        for (double y : series[loc])
            yds[loc].push_back(y - ybar[loc]);

        for (int i = y_skip[loc] ; i < num_visits; i++)
            sumXdSq[loc] += (xs[i] - xbar[loc]) * (xs[i] - xbar[loc]);

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
            for (int i = 1 ; i <= num_visits ; i++) {
                int ip = perm[i - 1];
                if (ip - 1 < y_skip[loc]) continue; // skip the NaN values
                beta += (xs[i - 1] - xbar[loc]) * yds[loc][ip - 1];
//cout << "x= " << i << " ip= " << ip << " yd= " << yds[loc][ip - 1] << endl;
            }
            beta /= sumXdSq[loc];

            long double alpha = ybar[loc] - beta * xbar[loc];

            long double sse = 0, se, t;
            for (int i = 1 ; i <= num_visits ; i++) {
                int ip = perm[i - 1];
                if (ip - 1 < y_skip[loc]) continue; // skip the NaN values
                sse += (series[loc][ip - 1] - alpha - beta * xs[i - 1]) * 
                       (series[loc][ip - 1] - alpha - beta * xs[i - 1]);
            }
//cout << "beta " << beta << " alpha " << alpha << " sse " << sse;

            if (beta >= slope_limit) {
                p_vals[loc] = 1.0; // t -> Inf, no contribution to S
            } else if (sse < 1e-10) {
                p_vals[loc] = 0.0; // t -> -Inf, massive contribution to S
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

    if (equal == S.size() - 1)
        return 1.0;
    else
        return 1.0 - (less_than + equal / 2.0) / (perm_count - 1.0);
}

/*
    This version breaks series into collections of rows that have the 
    same number of leading NaNs. The end result is the min of PoPLR
    on each collection.

    @param series A matrix of dB values; rows are locations, columns are visits
    @param perm_count Number of permutations of series to get p-value for S
    @param slope_limit Only include p-values in S statistic when PLR slope is <= slope_limit

    @return Estimated probability of stability of the series
*/
double PoPLR5(Series series, int perm_count, double slope_limit) {
    int num_locations = series.size() - 1;  // leading xs

        // Work out number of leading NaNs (skips)
    vector<int> nan_count(num_locations);    // number of leading NaNs
    set<int> nan_count_set = set<int>();

    for (int loc = 0; loc < num_locations; loc++) {
        nan_count[loc] = 0;
        for (double y : series[loc + 1])  // skip xs
            if (isnan(y))
                nan_count[loc]++;
        nan_count_set.insert(nan_count[loc]);
    }

    double minP = 1.0; // default stable
    for (int count : nan_count_set) {
        Series sub_series;
        sub_series.push_back(series[0]);  // xs

        for (int loc = 0; loc < num_locations; loc++) {
            if (nan_count[loc] == count)
                sub_series.push_back(series[loc + 1]);  // skip xs
        }
        int num_visits = sub_series[0].size() - count;

        if (num_visits < MIN_VISITS)
            continue;

        if (perm_count > tgamma(num_visits + 1))
            perm_count = tgamma(num_visits + 1);

        if (sub_series.size() > 0) {
            double p = PoPLR4(sub_series, perm_count, slope_limit);
            if (p < minP) 
                minP = p;
        }
    }
    return minP;
}

/*
    This version does not partition based on NaNs, it allows missing values in the PLR.
    ASSUMES: locations have at least 6 non-NaN values in their row.

    @param series A matrix of dB values; rows are locations, columns are visits
    @param perm_count Number of permutations of series to get p-value for S
    @param slope_limit Only include p-values in S statistic when PLR slope is <= slope_limit

    @return Estimated probability of stability of the series
*/
double PoPLR6(Series series, int perm_count, double slope_limit) {
    int num_locations = series.size() - 1;  // xs are first
//cout << "num_locations " << num_locations << endl;
    if (num_locations < 2) 
        return 1.0;

    vector<double> xs = get_xs(&series);

    int num_visits = series[0].size();

    if (perm_count > tgamma(num_visits + 1))
        perm_count = tgamma(num_visits + 1);

        // Work out ybar, yds values for all locations (excluding NaNs)
        // Also record which locs have NaNs
    vector<long double> ybar(num_locations);
    vector<vector<long double>> yds(num_locations);
    vector<bool> has_NaNs(num_locations);
    for (int loc = 0; loc < num_locations; loc++) {
        has_NaNs[loc] = false;
        ybar[loc]= 0.0;
        int n = 0;
        for (double y : series[loc]) {
            if (!isnan(y)) {
                ybar[loc] += y;
                n++;
            } else
                has_NaNs[loc] = true;
        }
        ybar[loc] /= (double)n;

        for (double y : series[loc])
            yds[loc].push_back(y - ybar[loc]);
    }
//print("ybars", ybar);

        // Work out the x* for locations without NaNs
    double xbar_all = 0;
    for (double x : xs) 
        xbar_all += x;
    xbar_all /= (double)num_visits;

    double sumXdSq_all = 0;
    for (double x : xs)
        sumXdSq_all += (x - xbar_all) * (x - xbar_all);

    double sqrtSumXdSq_all = sqrt(sumXdSq_all);

    // Now compute all the p-values for each loc in each appropriate permutation
    vector<long double> S;
    vector<long double> p_vals(num_locations);
    PermutationIterator pi(perm_count, num_visits);
    perm_count = 0;
    while (pi.hasNext()) {
        vector<int> perm = pi.next();
        perm_count++;
//print("perm:", perm);
        for (int loc = 0; loc < num_locations; loc++) {
                // (1) get the n and X values for this location
            int n;
            double xbar = 0, sumXdSq = 0, sqrtSumXdSq;
            if (has_NaNs[loc]) {
                n = 0;
                for (int i_x : perm) {
                    if (!isnan(series[loc][i_x - 1])) {
                        xbar += xs[i_x - 1];
                        n++;
                    }
                }
                xbar /= (double)n;

                for (int i_x : perm)
                    if (!isnan(series[loc][i_x - 1]))
                        sumXdSq += (xs[i_x - 1] - xbar) * (xs[i_x - 1] - xbar);
                
                sqrtSumXdSq = sqrt(sumXdSq);
            } else {
                n = num_visits;
                xbar = xbar_all;
                sumXdSq = sumXdSq_all;
                sqrtSumXdSq = sqrtSumXdSq_all;
            }
//cout << "n " << n << endl;
//cout << "xbar " << xbar << endl;
//cout << "ybar " << ybar[loc] << endl;
//cout << "sumXdSq " << sumXdSq << endl;
//cout << "sqrtSumXdSq " << sqrtSumXdSq << endl;

                // (2) Compute beta, alpha, etc
            long double beta = 0;
            for (int ix = 1 ; ix <= num_visits ; ix++) {
                int ip = perm[ix - 1];
                if (!isnan(series[loc][ip - 1])) {
                    beta += (xs[ix - 1] - xbar) * yds[loc][ip - 1];
//cout << "ix= " << ix << " ip= " << ip << " x= " << xs[ix - 1] << " y=" << series[loc][ip - 1] << " ybar= " << ybar[loc] << " yd= " << yds[loc][ip - 1] << endl;
                }
            }
            beta /= sumXdSq;

            long double alpha = ybar[loc] - beta * xbar;

            long double sse = 0, se, t;
            for (int ix = 1 ; ix <= num_visits ; ix++) {
                int ip = perm[ix - 1];
                double y = series[loc][ip - 1];
                if (!isnan(y))
                    sse += (y - alpha - beta * xs[ix - 1]) * (y - alpha - beta * xs[ix - 1]);
            }
//cout << "beta " << beta << " alpha " << alpha << " sse " << sse << " n " << n;

            if (beta >= slope_limit) {
                p_vals[loc] = 1.0; // t -> Inf, no contribution to S
            } else if (sse < 1e-10) {
                p_vals[loc] = 0.0; // t -> -Inf, massive contribution to S
            } else {
                se = sqrt(sse / (n - 2)) / sqrtSumXdSq;
                t = beta / se;
                boost::math::students_t dist(n - 2);
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

//cout <<  "S[0] " << S[0] << endl;
//print("S", S);
    int less_than = 0, equal = 0;
    for(int i = 1; i < S.size(); i++) {
        if (S[i] < S[0])
            less_than++;
        else if (S[i] == S[0])
            equal++;
    }
//cout << "less_than " << less_than << " equal " << equal << " perm_count " << perm_count << endl;

    if (equal == S.size() - 1)
        return 1.0;
    else
        return 1.0 - (less_than + equal / 2.0) / (perm_count - 1.0);
}

#endif // POPLR_ARREST_HPP