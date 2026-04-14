/*
Driver.
*/
#include "poplr_arrest.hpp"
#include "read_json.hpp"

using namespace std;

    // print usage and exit
void usage() {
    cerr << "Usage: poplr_arrest [-c n] [-{V | H | I}] [-d db_treatment] [-s slope_limit] p(tt|mt)_filename p(mt<17|tt)_filename filename.json" << endl;
    cerr << "       where" << endl;
    cerr << "           -c n restrict permutation count to n (default 5000 or visit!)." << endl;
    cerr << "           -V uses vertical grouping. (default uses vertical grouping)." << endl;
    cerr << "           -H uses horizontal grouping. (default uses vertical grouping)." << endl;
    cerr << "           -I uses independent grouping. (default uses vertical grouping)." << endl;
    cerr << "           -d db_treatment is 0 = Not arrest (no converting values, use every line)." << endl;
    cerr << "                           1 = ARREST with probs (convert values, apply PLR to," << endl;
    cerr << "                               Green lines, special probs to others)." << endl;
    cerr << "                           2 = ARREST no probs (convert values, only all-Green locations)." << endl;
    cerr << "                          (default = 0)." << endl;
    cerr << "           -s slope_limit Only include p-values in S statistic when PLR slope is < slope_limit." << endl;
    cerr << "                          (default = 0)." << endl;
    cerr << endl;
    cerr << "        Horizontal partitioning groups rows/locations with the same number of visits." << endl;
    cerr << "        Vertical partitioning groups columns/visits so that locations" << endl;
    cerr << "                              all have the same number of NAs at the start." << endl;
    cerr << "        Permuting happens within groups and group order remains the same." << endl;
    cerr << endl;
    cerr << "        P(tt|mt) file is a csv, one row per mt, one column per tt." << endl;
    cerr << endl;
    cerr << "        P(mt<17|tt) file is a csv, one row, one column per tt." << endl;
    cerr << endl;
    cerr << "        JSON file is [[" << endl;
    cerr << "                        [ // VF 1" << endl;
    cerr << "                           [loc_1_visit1, loc_1_visit_2, ..., loc_1_visit_n]," << endl;
    cerr << "                           [loc_2_visit1, loc_2_visit_2, ..., loc_2_visit_n]," << endl;
    cerr << "                            ..." << endl;
    cerr << "                           [loc_n_visit1, loc_n_visit_2, ..., loc_n_visit_n]" << endl;
    cerr << "                        ]," << endl;
    cerr << "                        [ // VF 2" << endl;
    cerr << "                           [loc_1_visit1, loc_1_visit_2, ..., loc_1_visit_n]," << endl;
    cerr << "                           [loc_2_visit1, loc_2_visit_2, ..., loc_2_visit_n]," << endl;
    cerr << "                           ..." << endl;
    cerr << "                           [loc_n_visit1, loc_n_visit_2, ..., loc_n_visit_n]]," << endl;
    cerr << "                        ], ..." << endl;
    cerr << "                     ]]" << endl;
    cerr << endl;
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

        // set default arguments
    int max_perm_count = 5000;                       // for -c
    function<double(Series, int, double)> PoPLR = PoPLR4;    // for -H  or -I
    ArrestSeries::ArrestProcessType arrest = ArrestSeries::ArrestProcessType::NOT_ARREST; // for -d
    double slope_upper_limit = 0;  // for -s
    string filename;
    string pr_tt_given_mt_filename;
    string pr_mtlt17_given_tt_filename;

        // process options
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-c") == 0) {
            max_perm_count = stoi(argv[++i]);
        } else if (strcmp(argv[i], "-V") == 0) {
            PoPLR = PoPLR4;
        } else if (strcmp(argv[i], "-H") == 0) {
            PoPLR = PoPLR5;
        } else if (strcmp(argv[i], "-I") == 0) {
            PoPLR = PoPLR6;
        } else if (strcmp(argv[i], "-d") == 0) {
            switch (stoi(argv[++i])) {
                case 0: arrest = ArrestSeries::ArrestProcessType::NOT_ARREST; break;
                case 1: arrest = ArrestSeries::ArrestProcessType::ARREST_WITH_PROBS; break;
                case 2: arrest = ArrestSeries::ArrestProcessType::ARREST_NO_PROBS; break;
                default: {
                    cerr << "Invalid db_type." << endl;
                    usage();
                }
            }
        } else if (strcmp(argv[i], "-s") == 0) {
            slope_upper_limit = stod(argv[++i]);
        } else if (i == argc - 3) {
            pr_tt_given_mt_filename = argv[i];
        } else if (i == argc - 2) {
            pr_mtlt17_given_tt_filename = argv[i];
        } else if (i == argc - 1) {
            filename = argv[i];
        }
    }
    if (filename.empty()) {
        cerr << "No json filename given." << endl;
        usage();
    }
    if (pr_tt_given_mt_filename.empty()) {
        cerr << "No P(mt|tt) filename given." << endl;
        usage();
    }

    vector<Eye> eyes = read_json(filename);

    cout << "eye,rep,visit,p" << endl;
    //#pragma omp parallel for
    for (int i_eye = 0 ; i_eye < eyes.size(); i_eye++) {
        for (int i_rep = 0 ; i_rep < eyes[i_eye].size(); i_rep++) {
            for (int visit = eyes[i_eye][i_rep][0].size(); visit >= 6; visit--) {
                Series series = preProcess(eyes[i_eye][i_rep], visit);
                ArrestSeries *as = new ArrestSeries(series, arrest, pr_tt_given_mt_filename, pr_mtlt17_given_tt_filename);
//cout << as->get_min_p() << endl;
//print("greens", as->get_series());
                double poplr_p = PoPLR(as->get_series(), max_perm_count, slope_upper_limit);
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
