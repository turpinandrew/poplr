#ifndef ARREST_SERIES_HPP
#define ARREST_SERIES_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <numeric>
#include <boost/math/distributions/students_t.hpp>
#include "types.hpp"

using namespace std;

/*
    Class for taking in some ARREST output (as a location * visits matrix)
    and altering it for use by PoPLR by
        1) Removing any locations with a red value at the end and computing the probability of that sequence.
        2) Removing any locations with a yellow or orange value at the end and computing the probability of that sequence.
        3) Converting all yellow/orange sequences to dB.
*/
class ArrestSeries {
public:

    enum class ArrestProcessType {
        NOT_ARREST = 1,          // Not expecting -999, -104, etc. All db >= -1
        ARREST_WITH_PROBS = 2,   // Treat seqs ending in non-Green with probs
        ARREST_NO_PROBS = 3      // Treat seqs ending in non-green as stable
    };

    /*
        @param series A matrix of dB or ARREST values.
        @param type If type == NOT_ARREST, do nothing, else call process
        @param pr_tt_given_mt_filename The filename of the CSV file containing P(tt|mt) values (one row per mt, col per tt).
        @param pr_mtlt17_given_tt_filename The filename of the CSV file containing P(mt < 17 | tt) values (one row, one column per tt).
    */
    ArrestSeries(Series series, ArrestProcessType type, string pr_tt_given_mt_filename, string pr_mtlt17_given_tt_filename) {
        this->series = series;

        if (type == ArrestProcessType::NOT_ARREST)
            return;

        this->pr_yrs = std::vector<std::pair<double, int>>();
        this->pr_gys = std::vector<std::pair<double, int>>();

        cache_pr_yellow_red = std::vector<std::vector<double>>(pr_yellow_red_cache_size, std::vector<double>(pr_yellow_red_cache_size, -1));

        read_csv(pr_tt_given_mt_filename, mt_domain_size, tt_domain_size, true, pr_tt_given_mt);
        read_csv(pr_mtlt17_given_tt_filename, 1, tt_domain_size, false, pr_mtlt17_given_tt);

        process(type);
    }
    //~ArrestSeries();  // TODO

    Series get_series() { return series; }

    /*
        Return the minimum p in pr_yrs and pr_gys. 1.0 if they are empty.
    */
    double get_min_p() {
        double p = 1.0;
        for (auto& pr : pr_gys) if (pr.first < p) p = pr.first;
        for (auto& pr : pr_yrs) if (pr.first < p) p = pr.first;
        return p;
    }

    vector<pair<double, int>> get_pr_gys() { return pr_gys; }
    vector<pair<double, int>> get_pr_yrs() { return pr_yrs; }

    /*
        Test out the probability functions.
    */
    void test() {
        for (int no = 1 ; no <= 5 ; no++) {
            for (int ny = 1 ; ny <= 5 ; ny++) {
                cout << "no = " << no << ", ny = " << ny << " : " << pr_yellow_red(no, ny, vector<long double>(0)) << endl;
            }
        }

        vector<double> y1 = {30, 29, 28, 27, 26, 25, Y(16), Y(16), Y(16), Y(16)};
        print("y1", y1);
        cout << pr_green_yellow_sequence(y1) << endl ;

        vector<double> y2 = {30, 31, 32, 33, 34, 35, Y(16), Y(16), Y(16), Y(16)};
        print("y2", y2);
        cout << pr_green_yellow_sequence(y2) << endl ;

        vector<double> y3 = {30, 28, 26, 24, 22, Y(14), Y(14), Y(14), Y(14), Y(14)};
        print("y3", y3);
        cout << pr_green_yellow_sequence(y3) << endl;

        double NaN = numeric_limits<double>::quiet_NaN();
        vector<double> y4 = {NaN, NaN, 26, 24, 22, Y(14), Y(14), Y(14), Y(14), Y(14)};
        print("y4", y4);
        cout << pr_green_yellow_sequence(y3) << endl;
    }

private:
    Series series;            // After locations removed and values converted to dB
    vector<pair<double, int>> pr_yrs; // Probabilities and locations for any pr_yellow_red_sequence
    vector<pair<double, int>> pr_gys; // Probabilities and locations for any pr_green_yellow_sequence

    /*
      For ArrestProcessType::ARREST_WITH_PROBS
        Find any 
            G..Y..R Apply pr_yellow_red_sequence
            Y..Y..R Apply pr_yellow_red_sequence
            G...Y|O Apply pr_green_yellow
        sequences and remove them from series, adding their prob to pr_yrs or pr_gys.
        Convert any remaining
          Y|O...Y|O sequences to dB
      For ArrestProcessType::ARREST_NO_PROBS
          Chuck out all sequences not ending in Green
    */
    void process(ArrestProcessType type) {
        int num_locations = series.size();
        vector<int> old_location_numbers(series.size()); // index into the original series locations/rows
        for (int loc = 0; loc < num_locations; loc++)
            old_location_numbers[loc] = loc;

        for (int loc = 0; loc < num_locations; ) {
            Location &loc_series = series[loc];
            bool keep = 1;
            if (type == ArrestProcessType::ARREST_NO_PROBS) {
                keep = is_green(loc_series.back());
            } else if (type == ArrestProcessType::ARREST_WITH_PROBS) {

                    // strip any leading NaNs
                int front;
                for ( front = 0; front < loc_series.size() && isnan(loc_series[front]) ; front++);

                if (loc_series.size() - front < 3) { // if series too short, chuck it, no judgement
                    keep = 0;
                } else if (get_db(loc_series[front]) < 6 && !is_green(loc_series.back())) {  // if the first measurement is a low yellow/orange, not prog
                    pr_yrs.push_back(std::make_pair(1.0, old_location_numbers[loc]));
                    keep = 0;
                } else if (is_red(loc_series.back())) { // if red at end, call pr_yrs
                    double p = pr_yellow_red_sequence(loc_series);
                    pr_yrs.push_back(std::make_pair(p, old_location_numbers[loc]));
                    keep = 0;
                } else if (is_yellow_orange(loc_series.back()) && is_green(loc_series[front])) { // G..Y|O
                    double p = pr_green_yellow_sequence(loc_series, front);
                    pr_gys.push_back(std::make_pair(p, old_location_numbers[loc]));
                    keep = 0;
                } else if (is_yellow_orange(loc_series.back()) && is_yellow_orange(loc_series[front])) { // Y|O..Y|O = stable
                    pr_gys.push_back(std::make_pair(1.0, old_location_numbers[loc]));
                    keep = 0;
                }
            }
            if (keep) {
                loc++;
            } else {
                    // swap last location up to here and shrink num_locations
                for (int visit = 0; visit < loc_series.size(); visit++)
                    series[loc][visit] = series[num_locations - 1][visit];
                num_locations--;

                old_location_numbers[loc] = num_locations; // swap location number up too
            }
        }
        series.erase(series.begin() + num_locations, series.end());  // erase any junk swapped to the end
    }

    /*
        Cummulative gaussian function (I hope!) Thanks co-pilot.
    */
    double pnorm(float x, float mean, float sd) {
        return 0.5 * (1 + std::erf((x - mean) / (sd * (float)std::sqrt(2))));
    }

        // mt's are -1:38
    static const int mt_domain_size = 37;
    static const int mt_domain_min = -1;
    static const int mt_domain_max = 35;

        // tt's are -10:40
    static const int tt_domain_size = 51;
    static const int tt_domain_min = -10;
    static const int tt_domain_max = 40;

        // Assumes not fp or fn, using Henson combined sd
        // Pr (not seeing  0 | tt in -10, -9, ..., 40, fpr=0%, fnr=0%)^2 = pr_not_see_0_twice[tt + 10]
    const long double pr_not_see_0_twice[tt_domain_size] {
        0.90670321322473, 0.870848799603662, 0.825897047152145,
        0.77145918922064, 0.707860981737141, 0.636280011816083,
        0.558767406423038, 0.478120335351116, 0.39760422347177,
        0.320564132309972, 0.25, 0.188196467088164, 0.136486903835298,
        0.0951954128030899, 0.0637524815168843, 0.0409367737433691,
        0.0251714896000551, 0.0148041983694028, 0.00831948660388055,
        0.00446320214137771, 0.00228391777035972, 0.00111399125855936,
        0.000517568503659564, 0.000228921136729339, 9.63406760883043e-05,
        3.85599434581467e-05, 1.46718152925694e-05, 5.30503487341715e-06,
        1.82222469579883e-06, 1.46433941203573e-07, 3.7537082074113e-09,
        3.74367563853235e-11, 1.14255734379317e-13, 7.91740029845402e-17,
        8.6019090265088e-21, 9.26058731040435e-26, 4.93038065763132e-32,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	// Pr (! seeing 16 | tt in -10, -9, ..., 40, fpr = 0%, fnr = 0%)^2 = pr_notsee_16_twice[tt + 10]
    const long double pr_not_see_16_twice[tt_domain_size] =
        {0.999985313206252, 0.999969091645071, 0.999936658519401,
        0.99987358553167, 0.999754282316219, 0.999534795958217,
        0.999142063426376, 0.998458624848598, 0.997302026161436,
        0.995398772771482, 0.992353910680113, 0.987619229291906,
        0.980465683418798, 0.969968641116258, 0.955017304607301,
        0.934360976088925, 0.90670321322473, 0.870848799603662,
        0.825897047152145, 0.77145918922064, 0.707860981737141,
        0.636280011816083, 0.558767406423038, 0.478120335351116,
        0.39760422347177, 0.320564132309972, 0.25, 0.188196467088164,
        0.136486903835298, 0.0885645429143801, 0.0489224924948571,
        0.0221640345757331, 0.00769441969238858, 0.0018707113417954,
        0.000282999985214457, 2.28328107278297e-05, 8.04744325390379e-07,
        9.58544998915119e-09, 2.77870011451845e-11, 1.29068888875652e-14,
        5.65418537516187e-19, 1.19565183446615e-24, 4.93038065763132e-32,
        0, 0, 0, 0, 0, 0, 0, 0};

    long double pr_tt_given_mt[mt_domain_size][tt_domain_size];
    long double pr_mtlt17_given_tt[1][tt_domain_size];

    /*
        Read csv file with nrows row, ncols columns

        @param filename Filename of a CSV file with nrows row and ncols columns
        @param nrows Number of rows to read
        @param ncols Number of cols to read
        @param row_sum_check Whether to check if rows sum to 1

        @return 0 on success, -1 on failure
    */
    int read_csv(const std::string filename, int nrows, int ncols, bool row_sum_check, long double data[][tt_domain_size]) {
        FILE *f = fopen(filename.c_str(), "r");
        if (!f) {
            perror(("Error opening file " + filename).c_str());
            return -1;
        }

            // read csv from f (loop via copilot)
        char line[2048];
        int row = 0;
        while (fgets(line, sizeof(line), f) && row < nrows) {
            int col = 0;
            char *token = strtok(line, ",");
            while (token && col < ncols) {
                if (strcmp(token, "NaN") == 0) {
                    data[row][col] = std::numeric_limits<double>::quiet_NaN();
                } else {
                    data[row][col] = strtold(token, nullptr);
                }
                token = strtok(nullptr, ",");
                col++;
            }
            if (col != ncols) {
                fprintf(stderr, "Row %d in %s has incorrect number of columns: %d expected %d\n", row, filename.c_str(), col, ncols);
                return(-1);
            }
            row++;
        }
        if (row != nrows) {
            fprintf(stderr, "File %s has incorrect number of rows: %d expected %d\n", filename.c_str(), row, nrows);
            return(-1);
        }

        fclose(f);
    
        // check rows sum to 1
        if (row_sum_check)
        for (int row = 0; row < nrows; row++) {
            long double row_sum = 0;
            for (int col = 0; col < ncols; col++)
                row_sum += data[row][col];
            if (abs(row_sum - 1.0) > 1e-12) {
                fprintf(stderr, "Row %d does not sum to 1 in %s: %Lf\n", row, filename.c_str(), row_sum);
                return(-1);
            }
        }

        return(0);
    }

    /*
        @param loc_series A location series (vector) of ARREST values.
        @param front Index of first non NaN in loc_series

        @return Probability of a G...G sequence at the start of loc_series 
                occurring if underlying true is not progressing. (stable or increasing)
                If loc_series is empty, returns uniform prior.

        @note Alters value of *front to the index of the first non-green value.
    */
    vector<long double> pr_tt_if_stable(const Location& loc_series, int *front) {
        vector<long double> pr_tt(tt_domain_size);
        for (int i = 0 ; i < pr_tt.size(); i++)
            pr_tt[i] = 1;

//print("loc_series", loc_series);
            // accumulate product of P(tt | green value at visit)
        int visit = *front;
        for ( ; visit < loc_series.size() && is_green(loc_series[visit]); visit++)
            for (int i = 0 ; i < pr_tt.size(); i++) {
                //cout << "pr_tt_given_mt " << (int)std::round((float)loc_series[visit]) + 1;
                //cout << " " << pr_tt_given_mt[(int)std::round((float)loc_series[visit]) + 1][i] << endl;
                pr_tt[i] *= pr_tt_given_mt[(int)std::round((float)loc_series[visit]) + 1][i];
            }

        *front = visit;
// print("pr_tt", pr_tt);
            // normalise pr_tt
        long double sum = 0;
        for (int i = 0 ; i < pr_tt.size(); i++)
            sum += pr_tt[i];
        for (int i = 0 ; i < pr_tt.size(); i++)
            pr_tt[i] /= sum;

        return pr_tt;
    }

    /*
        For a G..GY..Y|O sequence, 
            1) Find prob distribution of stable tt by p = \prod_G pr_tt_given_mt / Z
            2) Find E[pr_notsee_16_twice] = \sum p * pr_notsee_16_twice

        @param loc_series A location series (vector) of ARREST values. 
        @param front Index of first non NaN in loc_series
        
        ASSUMES loc_series starts with Green and terminates with a Yellow or Orange.

        @return Probability of a G...Y|O sequence occurring if underlying true is not progressing. (stable or increasing)
    */

    double pr_green_yellow_sequence(const Location& loc_series, int front = 0) {
        int *visit = &front;

        vector<long double> pr_tt = pr_tt_if_stable(loc_series, visit);

        if (*visit < 2)
           return 1;  // We need 2 green to get a good handle on the underlying pr_tt... (do we!? need more than 1, for sure.)
                      // 1 green only at start could be learning effect

            // Compute expectation of pr_notsee_16_twice using pr_tt as the probs.
        long double p = 0;
        for (int i = 0 ; i < pr_tt.size(); i++)
            p += pr_tt[i] * pr_not_see_16_twice[i];

        return (double)p;
    }


    /*
        Use any leading greens to set a prior on true threshold, tt.

        @param loc_series A location series (vector) of ARREST values.

        @return Probability of a sequence ending in red if it is stable.
    */
    double pr_yellow_red_sequence(const Location& loc_series) {
        int num_orange = count_if(loc_series.begin(), loc_series.end(), is_orange);
        int num_yellow = count_if(loc_series.begin(), loc_series.end(), is_yellow_orange) - num_orange;

        //return std::pow(0.15, num_orange) * std::pow(0.85, num_yellow);

        vector<long double> pr_tt = vector<long double>(0);
        if (is_green(loc_series[0])) {
            int visit = 0;
            pr_tt = pr_tt_if_stable(loc_series, &visit);

            if (visit < 2)           // need at least 2 greens to avoid "learning effect"
                pr_tt.resize(0);
        }

        //print("pr_tt: ", pr_tt);

        return pr_yellow_red(num_orange, num_yellow, pr_tt);
    }

    /*
        # Let p be prob of not seeing 0dB twice assuming stable at tt
        # Working forwards towards current Red
        #     Yellow -> Orange = p
        #     Orange -> Red = p
        #     Yellow -> Yellow = * (1 - p)
        #     Orange -> Yellow = * (1 - p)
        # For a sequence ending in Red+, every Yellow and Orange had to have a ->
        # Assuming there are `ny` Yellow and `no` Orange leading up to the Red
        # Each Orange had to have a Y->O coming in, so they have prob p^{no}
        # Each of the Yellows but the first had to have a (1-p) coming in, so they have prob (1-p)^{ny - 1}
        # Pr of getting to the first Yellow = (1 - pr_see(16))^2
        # Pr getting a ZEST result less than 17 (tt_prior)

        # Assuming all stable tt's are, equally likely this can be evaluated as...
        @param num_orange Number of orange values in the sequence.
        @param num_yellow Number of yellow values in the sequence.
        @param tt_prior Prior probability of each true threshold. If size() == 0, then assume uniform prior.

        @return Probability of a yellow or orange sequence ending in red.
    */
    double pr_yellow_red(int num_orange, int num_yellow, vector<long double> tt_prior) {
        if (num_orange < 1 || num_yellow < 0)
            return 0.0;

        if (num_yellow == 0)   // in this incantation of ARREST7, it is possible to get Orange on first visit, so spoof YO
            num_yellow = 1;

        double *p, pp = -1;
        p = &pp;

        if (tt_prior.empty()) {  // don't use tt_prior and cache result
            if (num_orange < pr_yellow_red_cache_size && num_yellow < pr_yellow_red_cache_size)
                p = &cache_pr_yellow_red[num_orange][num_yellow];

            if (*p == -1) {
                *p = 0;
                for (int i = 0 ; i < tt_domain_size ; i++) {
                    *p += 1.0 / tt_domain_size *
                        pr_mtlt17_given_tt[0][i] *
                        pr_not_see_16_twice[i] *
                        std::pow(pr_not_see_0_twice[i], num_orange) *
                        std::pow(1 - pr_not_see_0_twice[i], num_yellow - 1);
                }
            }
        } else { // no caching and use tt_prior
//print("tt_prior in p_yr", tt_prior);
            *p = 0;
            for (int i = 0 ; i < tt_domain_size ; i++) {
                *p += tt_prior[i] * 
                    pr_mtlt17_given_tt[0][i] *
                    pr_not_see_16_twice[i] *
                    std::pow(pr_not_see_0_twice[i], num_orange) *
                    std::pow(1 - pr_not_see_0_twice[i], num_yellow - 1);
            }
        }

        return *p;
    }

        // cache for pr_yellow_red when tt_prior is uniform
    const int pr_yellow_red_cache_size = 20;
    vector<vector<double>> cache_pr_yellow_red;
};

#endif // ARREST_SERIES_HPP