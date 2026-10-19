#include <vector>
#include <numeric>
#include <cmath>
#include <boost/math/distributions/students_t.hpp>
#include <iostream>
#include "permutation.hpp"
#include "poplr.hpp"

// @param series A pair whose 2nd element is a matrix of dB values with `num_locations` rows and `num_visits` columns.
// @return The right boundaries of NaN blocks. Must end in num_visits.
vector<int> get_boundaries(const Series& series) {
    vector<Location> yvals = series.second;

    int num_visits = yvals[0].size();

    // Count NAs for each visit
    vector<int> num_NAs(num_visits, 0);
    for (int visit = 0; visit < num_visits; ++visit) {
        for (const auto& loc : yvals) {
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

/*
    Performs PoPLR as published on the series.
    Does not partition based on NaNs.
    NaNs are just skipped when computing means and sums.

    ASSUMES: locations have at least 6 non-NaN values in their row.

    @param series A pair <vector of x values, A matrix of dB values; rows are locations, columns are visits>
    @param perm_count Number of permutations of series to get p-value for S
    @param verbosity_level Verbosity level for debugging and output

    @return Estimated probability of stability of the series
*/
double PoPLR(Series series, PermutationIterator *perms, float p_threshold, int verbosity_level) {
    int num_locations = series.second.size();
    if (num_locations < 2) 
        return 1.0;

    int num_visits = series.first.size();

        // work out the x*, ybar, yds, y_skip values for each location
    vector<int> n(num_locations);
    vector<long double> xbar(num_locations);
    vector<long double> sumXdSq(num_locations);
    vector<long double> sqrtSumXdSq(num_locations);
    vector<long double> ybar(num_locations);
    vector<vector<long double>> yds(num_locations);  // for all including leading NaNs

    for (int loc = 0; loc < num_locations; loc++) {
        n[loc] = 0;
        xbar[loc]= 0;
        ybar[loc]= 0;
        sumXdSq[loc] = 0;
            // get non-nans
            // xs come from series.first (the vector of x values)
            // ys come from series.second[loc] (the row of dB values)
        for (int i = 0 ; i < num_visits ; i++) {
            double y = series.second[loc][i];
            if (!isnan(y)) {
                ybar[loc] += y;
                xbar[loc] += series.first[i];
                n[loc]++;
            }
        }

        xbar[loc]= xbar[loc] / n[loc];
        ybar[loc]= ybar[loc] / n[loc];

        for (int i = 0 ; i < num_visits ; i++) {
            double y = series.second[loc][i];
            yds[loc].push_back(y - ybar[loc]);
            if (!isnan(y)) {
                double x = series.first[i];
                sumXdSq[loc] += (x - xbar[loc]) * (x - xbar[loc]);
            }
        }

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

    int perm_count = 0;
    while (perms->hasNext()) {
        vector<int> perm = perms->next();
        perm_count++;
//print("perm:", perm);

//cout << "num_locations " << num_locations << endl;
        for (int loc = 0; loc < num_locations; loc++) {
            long double beta = 0;
            for (int ip : perm) {
                if (isnan(yds[loc][ip - 1])) continue;
                double x = series.first[ip - 1];
                beta += (x - xbar[loc]) * yds[loc][ip - 1];
//cout << "x= " << x << " ip= " << ip << " yd= " << yds[loc][ip - 1] << endl;
            }
            beta /= sumXdSq[loc];

            long double alpha = ybar[loc] - beta * xbar[loc];

            long double sse = 0, se, t;
            for (int ip : perm) {
                if (isnan(series.second[loc][ip - 1])) continue; // skip the NaN values
                double x = series.first[ip - 1];
                sse += (series.second[loc][ip - 1] - alpha - beta * x) * 
                       (series.second[loc][ip - 1] - alpha - beta * x);
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
            if (p <= p_threshold)
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
    This version does no partitioning of the series based on NaNs.
    ASSUMES: locations have at least 6 non-NaN values in their row.

    @param series A pair <vector of x values, A matrix of dB values; rows are locations, columns are visits>
    @param perm_count Number of permutations of series to get p-value for S
    @param verbosity_level Verbosity level for debugging and output

    @return Estimated probability of stability of the series
*/
double PoPLR0(const Series series, const int perm_count, const float p_threshold, int verbosity_level) {
    int num_locations = series.second.size();
    if (num_locations < 2) 
        return 1.0;

    int num_visits = series.first.size();

    int this_pc = perm_count;
    if (this_pc > tgamma(num_visits + 1))
        this_pc = tgamma(num_visits + 1);

    PermutationIterator pi(this_pc, num_visits);
    return PoPLR(series, &pi, p_threshold, verbosity_level);
}

/*
    This version vertically partitions series based on NaNs at the start of rows.
    ASSUMES: locations have at least 6 non-NaN values in their row.

    @param series A pair <vector of x values, A matrix of dB values; rows are locations, columns are visits>
    @param perm_count Number of permutations of series to get p-value for S
    @param verbosity_level Verbosity level for debugging and output

    @return Estimated probability of stability of the series
*/
double PoPLR4(const Series series, const int perm_count, const float p_threshold, int verbosity_level) {
    int num_locations = series.second.size();
    if (num_locations < 2) 
        return 1.0;

    int num_visits = series.first.size();

    int this_pc = perm_count;
    if (this_pc > tgamma(num_visits + 1))
        this_pc = tgamma(num_visits + 1);

    vector<int> boundaries = get_boundaries(series);
    PermutationIterator pi(this_pc, boundaries);
    return PoPLR(series, &pi, p_threshold, verbosity_level);
}

/*
    This version breaks series into collections of rows that have the 
    same number of leading NaNs. The end result is the min of PoPLR
    on each collection.

    @param series A pair <vector of x values, A matrix of dB values; rows are locations, columns are visits>
    @param perm_count Number of permutations of series to get p-value for S
    @param verbosity_level Verbosity level for debugging and output

    @return Estimated probability of stability of the series
*/
double PoPLR5(const Series series, const int perm_count, const float p_threshold, int verbosity_level) {
    int num_locations = series.second.size();

        // Work out number of leading NaNs (skips)
    vector<int> nan_count(num_locations);    // number of leading NaNs
    set<int> nan_count_set = set<int>();

    for (int loc = 0; loc < num_locations; loc++) {
        nan_count[loc] = 0;
        for (double y : series.second[loc])
            if (isnan(y))
                nan_count[loc]++;
        nan_count_set.insert(nan_count[loc]);
    }

        // For each distinct count, call PoPLR
    double minP = 1.0; // default stable
    for (int count : nan_count_set) {
        int num_visits = series.second[0].size() - count;

        if (num_visits < MIN_VISITS)
            continue;

        Series sub_series = Series();
        sub_series.first = series.first;

        for (int loc = 0; loc < num_locations; loc++)
            if (nan_count[loc] == count)
                sub_series.second.push_back(series.second[loc]);

        int this_perm_count = perm_count;
        if (this_perm_count > tgamma(num_visits + 1))
            this_perm_count = tgamma(num_visits + 1);

        if (sub_series.second.size() > 0) {
            PermutationIterator pi(this_perm_count, num_visits);
            double p = PoPLR(sub_series, &pi, p_threshold, verbosity_level);
            if (p < minP) 
                minP = p;
        }
    }
    return minP;
}
