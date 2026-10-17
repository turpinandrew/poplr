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


/opt/homebrew/Cellar/gcc/14.2.0/bin/g++-14 -Wc++11-extensions -std=c++14 -I /opt/homebrew/include/ poplr_arrest.cpp -o poplr_arrest
*/
#include <vector>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <random>
#include <fstream>
#include <functional>
#include <boost/math/distributions/students_t.hpp>
#include <boost/container_hash/hash.hpp>
#include "permutation.hpp"
#include "read_utils.hpp"
#include "poplr.hpp"
#include "arrestSeries.hpp"

//#include <_stdio.h>
//#include <functional>
#include <omp.h>

using namespace std;

/*
    Resize matrix to have at least 2 rows and min(num_visits, max loc seq length) columns.
    Remove locations with less than MIN_VISITS non-NaN in their row.

    @param series A pair whose 1st element is list of num_visits x values (years) and 2nd element is a matrix of dB values 
                  with at least `num_locations` rows and up to `num_visits` columns.
    @param num_locations Number of rows in series
    @param num_visits Max number of columns in series

    @return Updated series pair or empty Series if cannot meet requirements.
*/
Series preProcess(Series series, int num_visits) {
    int n = series.second.size();

    if (num_visits < MIN_VISITS || n < 2)
        return Series();

        // Cap number of x values at num_visits in size (including NaNs)
    if (series.first.size() > num_visits)
        series.first.erase(series.first.begin() + num_visits, series.first.end());

        // Cap length of each location at num_visits in size (including NaNs)
    for (int loc = 0; loc < n; loc++)
        if (series.second[loc].size() > num_visits)
            series.second[loc].erase(series.second[loc].begin() + num_visits, series.second[loc].end());

        // Throw out locations with less than MIN_VISITS non-nan
    for (int loc = 0; loc < n; ) {
        int nan_count = count_if(series.second[loc].begin(), series.second[loc].end(), [](double x){return(isnan(x));});
        if (series.second[loc].size() - nan_count < MIN_VISITS) {
                // swap last location up to here and shrink num_locations
            for (int visit = 0; visit < num_visits; visit++)
                series.second[loc][visit] = series.second[n - 1][visit];
            n--;
        } else
            loc++;
    }
    if (n < 2) 
        return Series();
    else 
        series.second.erase(series.second.begin() + n, series.second.end()); // throw away junk swapped to the end 

    return series;
}
    
    // print usage and exit
void usage() {
    std::cerr << "Usage: poplr_arrest [-c n] [-t p] part_style db_treatment tt_file mt_file filename.json" << std::endl;
    std::cerr << "       where" << std::endl;
    std::cerr << "           part_style   is H[orizontal] | V[ertical | N[one]] grouping" << std::endl;
    std::cerr << "           db_treatment is 0 = Not arrest (no converting values, use every line)" << std::endl;
    std::cerr << "                           1 = ARREST with probs (convert values, apply PLR to Green lines, special probs to others)" << std::endl;
    std::cerr << "                           2 = ARREST no probs (convert values, only all-Green locations)" << std::endl;
    std::cerr << "           -c n         Restrict permutation count to n (default 5000 or visit!)." << std::endl;
    std::cerr << "           -t p         Only include p-values <= p in the computation of S (default 1.1)." << std::endl;
    std::cerr << "          tt_file is a csv file containing P(tt|mt): one row per mt (" << ArrestSeries::mt_domain_size << " rows ";
    std::cerr << ArrestSeries::mt_domain_min << ":" << ArrestSeries::mt_domain_max << "), one column per tt (";
    std::cerr << ArrestSeries::tt_domain_size << " cols, " << ArrestSeries::tt_domain_min << ":" << ArrestSeries::tt_domain_max << ")." << std::endl;
    std::cerr << "          mt_file is a csv file containing P(mt<17|tt): one row, one column per tt (";
    std::cerr << ArrestSeries::tt_domain_size << " cols, " << ArrestSeries::tt_domain_min << ":" << ArrestSeries::tt_domain_max << ")." << std::endl;
    std::cerr << "          filename.json is the JSON input file." << std::endl;
    std::cerr << std::endl;
    std::cerr << "        Horizontal partitioning groups rows/locations with the same number of visits." << std::endl;
    std::cerr << "        Vertical partitioning groups columns/visits so that locations all have the same number of NAs at the start." << std::endl;
    std::cerr << "        None partitioning does no grouping and just ignores NAs in regressions." << std::endl;
    std::cerr << "        Permuting happens within groups and group order remains the same." << std::endl;
    std::cerr << std::endl;
    std::cerr << "        JSON file is [[" << std::endl;
    std::cerr << "                        [ // rep 1" << std::endl;
    std::cerr << "                           [x1, x2, ..., x_n]," << std::endl;
    std::cerr << "                           [loc_1_visit1, loc_1_visit_2, ..., loc_1_visit_n]," << std::endl;
    std::cerr << "                           [loc_2_visit1, loc_2_visit_2, ..., loc_2_visit_n]," << std::endl;
    std::cerr << "                            ..." << std::endl;
    std::cerr << "                           [loc_n_visit1, loc_n_visit_2, ..., loc_n_visit_n]" << std::endl;
    std::cerr << "                        ]," << std::endl;
    std::cerr << "                        [ // rep 2" << std::endl;
    std::cerr << "                           [x1, x2, ..., x_n]," << std::endl;
    std::cerr << "                           [loc_1_visit1, loc_1_visit_2, ..., loc_1_visit_n]," << std::endl;
    std::cerr << "                           [loc_2_visit1, loc_2_visit_2, ..., loc_2_visit_n]," << std::endl;
    std::cerr << "                           ..." << std::endl;
    std::cerr << "                           [loc_n_visit1, loc_n_visit_2, ..., loc_n_visit_n]" << std::endl;
    std::cerr << "                        ], ..." << std::endl;
    std::cerr << "                     ]]" << std::endl;    std::cerr << std::endl;
    exit(-1);
}

int main(int argc, char *argv[]) {
    if (argc < 6) {
        std::cerr << "\nMissing a mandatory argument." << std::endl << std::endl;
        usage();
    }

    int max_perm_count = 5000;
    float p_threshold = 1.1;
    string filename;
    string pr_tt_given_mt_filename;
    string pr_mtlt17_given_tt_filename;
    int optional_arg_pos = argc - 7; // p_a _ _ _ _  HV 012 tt mt json
                                     //              -5 -4 -3  -2 argc-1

    while (optional_arg_pos > 0) {
        if (strcmp(argv[optional_arg_pos], "-c") == 0)
            max_perm_count = stoi(argv[optional_arg_pos + 1]);
        else if (strcmp(argv[optional_arg_pos], "-t") == 0)
            p_threshold = stof(argv[optional_arg_pos + 1]);
        else {
            std::cerr << "Unknown optional argument: " << argv[optional_arg_pos] << endl;
            usage();
        }

        optional_arg_pos -= 2;
    }

        // argc - 5 is style of partitioning
    function<double(Series, int, float)> PoPLR_to_use;
    if (argv[argc - 5][0] == 'H' || argv[argc - 5][0] == 'h')
        PoPLR_to_use = PoPLR5;
    else if (argv[argc - 5][0] == 'V' || argv[argc - 5][0] == 'v')
        PoPLR_to_use = PoPLR4;
    else if (argv[argc - 5][0] == 'N' || argv[argc - 5][0] == 'n')
        PoPLR_to_use = PoPLR0;
    else {
        std::cerr << "\nInvalid part_type." << std::endl << std::endl;
        usage();
    }

        // argc - 4 is type of data
    ArrestSeries::ArrestProcessType arrest;
    switch (stoi(argv[argc - 4])) {
        case 0: arrest = ArrestSeries::ArrestProcessType::NOT_ARREST; break;
        case 1: arrest = ArrestSeries::ArrestProcessType::ARREST_WITH_PROBS; break;
        case 2: arrest = ArrestSeries::ArrestProcessType::ARREST_NO_PROBS; break;
        default: {
            std::cerr << "\nInvalid db_type." << std::endl << std::endl;
            usage();
        }
    }
    
    pr_tt_given_mt_filename = argv[argc - 3];
    pr_mtlt17_given_tt_filename = argv[argc - 2];
    filename = argv[argc - 1];

    if (filename.empty()) {
        std::cerr << "\nNo P(tt | mt) filename given." << std::endl << std::endl;
        usage();
    }
    if (pr_tt_given_mt_filename.empty()) {
        std::cerr << "\nNo P(mt < 17 | tt) filename given." << std::endl << std::endl;
        usage();
    }
    if (pr_mtlt17_given_tt_filename.empty()) {
        std::cerr << "\nNo P(mt < 17 | tt) filename given." << std::endl << std::endl;
        usage();
    }

    vector<vector<long double>> pr_tt_given_mt = read_csv(pr_tt_given_mt_filename, ArrestSeries::mt_domain_size, ArrestSeries::tt_domain_size, true);
    vector<vector<long double>> pr_mtlt17_given_tt_2d = read_csv(pr_mtlt17_given_tt_filename, 1, ArrestSeries::tt_domain_size, false);
    vector<long double> pr_mtlt17_given_tt;

    if (pr_tt_given_mt.empty()) {
        std::cerr << "\nTrouble reading " << pr_tt_given_mt_filename << std::endl << std::endl;
        usage();
    }
    if (pr_mtlt17_given_tt_2d.empty()) {
        std::cerr << "\nTrouble reading " << pr_mtlt17_given_tt_filename << std::endl << std::endl;
        usage();
    } else {
        pr_mtlt17_given_tt = pr_mtlt17_given_tt_2d[0];
    }

    vector<Eye> eyes = read_json(argv[argc - 1]);

    std::cout << "eye,rep,visit,p" << std::endl;
    #pragma omp parallel for
    for (int i_eye = 0 ; i_eye < eyes.size(); i_eye++) {
        for (int i_rep = 0 ; i_rep < eyes[i_eye].size(); i_rep++) {
            for (int visit = eyes[i_eye][i_rep].second[0].size(); visit >= 6; visit--) {
                Series series = preProcess(eyes[i_eye][i_rep], visit);
                ArrestSeries *as = new ArrestSeries(series, arrest, pr_tt_given_mt, pr_mtlt17_given_tt);
                double poplr_p = PoPLR_to_use(as->get_series(), max_perm_count, p_threshold);
                double other_ps = as->get_min_p();
                #pragma omp critical
                {
                std::cout << i_eye << "," << i_rep << "," << visit << ",";
                std::cout << (poplr_p < other_ps ? poplr_p :other_ps);
                std::cout << std::endl;
                //print("Series:", series);
                }
            }
        }
    }
    
    return 0;
}