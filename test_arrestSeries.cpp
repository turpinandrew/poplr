
#include <iostream>
#include <string>
#include <fstream>

#include "read_utils.hpp"
#include "arrestSeries.hpp"
#include "poplr.hpp"
#include "permutation.hpp"

using namespace std;

const string PTT="/Users/aturpin/src/hIPPOS/analysis/probs/arrest11_fp15_fn03_pr_tt_given_mt.csv";
const string F17="/Users/aturpin/src/hIPPOS/analysis/probs/arrest11_fp15_fn03_pr_mtlt17_given_tt.csv";

vector<vector<long double>> pr_tt_given_mt;
vector<long double> pr_mtlt17_given_tt;


Series s1 {{1,2,3,4,5,6,7,8,9,10}, {
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30}
}};

Series a1 {{}, {
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30}
}};

    // end in Red
Series s2 {{1,2,3,4,5,6,7,8,9,10}, {
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, Y(17), Y(17), Y(17), O(17), R(), R()}
}};

Series a2 {{}, {
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30}
}};

    // all Yellow
Series s3 {{1,2,3,4,5,6,7,8,9,10}, {
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10)}
}};

Series a3 {{}, {
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30}
}};

    // degenerate
Series s4 {{1,2,3,4,5,6,7,8,9,10}, {
    {30, Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10)},
    {30, Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10)},
    {30, Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10)},
    {30, Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10)},
    {30, Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10), Y(10)}
}};

Series a4 {{}, {}};

    // G G G Y Y Y 
Series s5 {{1,2,3,4,5,6,7,8,9,10}, {
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {16, 16, 16, 16, 16, Y(16), Y(16), Y(16), Y(16), Y(16)},
    {30, 28, 26, 24, 22, Y(14), Y(14), Y(14), Y(14), Y(14)},
    {30, 28, 26, 24, Y(14), Y(14), Y(14), Y(14), Y(14), Y(14)},
    {30, 28, 26, Y(14), Y(14), Y(14), Y(14), Y(14), Y(14), Y(14)},
    {30, 28, Y(14), Y(14), Y(14), Y(14), Y(14), Y(14), Y(14), Y(14)},
    {30, Y(14), Y(14), Y(14), Y(14), Y(14), Y(14), Y(14), Y(14), Y(14)}
}};

Series a5 {{1,2,3,4,5,6,7,8,9,10}, {
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
    {30, 30, 30, 30, 30, 30, 30, 30, 30, 30}
}};

Series s6 {{1,2,3,4,5,6}, {
    {30, 30, 30, Y(14), O(14), -999},
    {25, 25, 25, Y(14), O(14), -999},
    {20, 20, 20, Y(14), O(14), -999},
    {15, 15, 15, Y(14), O(14), -999},
    {10, 10, 10, Y(14), O(14), -999}
}};

double NaN = numeric_limits<double>::quiet_NaN();

Series a_ip006 {{0, 0.364134154688569, 0.747433264887064, 1.13073237508556, 1.34154688569473, 1.76317590691307, 2.03148528405202, 2.34360027378508, 2.57357973990418, 2.95687885010267}, {
{23,27,25,27,24,24,17,22,25,24},
{25,26,22,23,22,25,22,25,21,25},
{22,27,25,25,17,22,21,23,25,22},
{27,25,22,28,17,23,18,24,27,22},
{27,28,25,28,27,25,25,28,27,27},
{23,27,25,27,23,30,25,24,28,28},
{24,27,22,25,24,25,22,22,25,25},
{27,25,25,25,22,27,21,22,27,21},
{22,27,22,27,25,24,23,22,27,24},
{27,24,25,28,22,27,22,23,28,27},
{27,30,28,27,25,30,27,27,27,27},
{27,28,28,30,25,30,27,25,27,28},
{30,27,28,28,27,28,25,27,28,27},
{27,27,28,25,27,30,24,27,28,23},
{27,25,25,30,28,27,27,30,27,30},
{23,25,23,27,24,24,24,22,27,25},
{27,30,25,28,27,30,25,23,28,27},
{27,28,25,30,25,28,25,25,27,27},
{27,28,28,30,28,30,27,28,28,28},
{28,28,27,30,27,30,27,25,30,24},
{27,27,27,27,27,28,27,27,27,28},
{30,30,28,28,28,27,28,28,30,30},
{30,28,28,30,27,30,28,30,27,30},
{28,30,28,30,28,28,30,30,30,28},
{30,28,30,30,30,30,27,30,30,27},
{28,25,25,30,27,30,30,28,30,28},
{28,28,28,30,30,30,28,28,30,30},
{22,27,28,27,28,25,27,25,25,30},
{22,25,27,27,27,27,27,25,26,27},
{28,25,25,27,27,27,28,23,27,25},
{25,27,28,28,27,27,27,25,28,28},
{NaN,25,28,27,28,27,22,21,21,25},
{25,25,27,30,27,27,28,27,30,28},
{27,27,27,30,28,28,30,26,27,28},
{28,27,25,28,26,30,28,25,28,30},
{28,27,27,27,28,30,30,27,28,28},
{27,28,25,30,30,28,30,30,27,27},
{NaN,NaN,NaN,NaN,NaN,NaN,22,6,6,6},
{NaN,NaN,30,28,25,30,25,30,28,28},
{28,25,26,27,25,27,24,27,27,25},
{24,23,22,25,25,22,11,6,22,22},
{6,27,27,27,6,22,6,22,19,27},
{27,25,28,28,24,25,27,28,27,30},
{27,27,30,28,30,28,28,25,27,28},
{NaN,NaN,24,23,6,15,22,11,24,22},
{NaN,NaN,23,11,27,17,24,25,28,28},
{NaN,27,22,27,21,25,27,27,28,17},
{NaN,28,28,30,27,30,27,27,27,28},
{30,28,30,27,30,30,30,27,30,30},
{30,30,25,30,28,25,30,30,30,30},
{27,30,28,30,30,30,30,27,28,30},
{28,28,27,30,30,28,27,27,30,30},
{NaN,NaN,22,22,22,23,22,18,22,25}
}};

/*
    @param s Series to test
    @param a Expected series after processing
    @param expected_p Expected minimum p
    @param expected_size Expected size of resulting series (should be nrow(a))
    @param header Header for output
*/
void tester(Series s, Series a, double expected_p, int expected_size, string header) {
    ArrestSeries *as = new ArrestSeries(s, ArrestSeries::ArrestProcessType::ARREST_WITH_PROBS, pr_tt_given_mt, pr_mtlt17_given_tt);

    double poplr_p = PoPLR0(as->get_series(), 1000, 0.05, 1, -0.1);
    double other_ps = as->get_min_p();
    double p = poplr_p < other_ps ? poplr_p :other_ps;

    bool pass = abs(p - expected_p) < 0.00001;
    pass &= as->get_series().second.size() == expected_size;
    
    vector<double> first_fail = {-1, 1, -2, -2, -1, -1};
    if (pass)
        for(int i = 0 ; i < a.second.size() ; i++)
            for(int j = 0 ; j < a.second[i].size() ; j++) {
                if (isnan(a.second[i][j])) {
                    if (!isnan(as->get_series().second[i][j])) {
                        first_fail = {(double)i, (double)j, a.second[i][j], as->get_series().second[i][j], (double)isnan(a.second[i][j]), (double)isnan(as->get_series().second[i][j])};
                        pass = false;
                    }
                } else if (a.second[i][j] != round(as->get_series().second[i][j])) {
                    first_fail = {(double)i, (double)j, a.second[i][j], as->get_series().second[i][j], (double)isnan(a.second[i][j]), (double)isnan(as->get_series().second[i][j])};
                    pass = false;
                }
            }

    cout << "\n*********************\n" << header << (pass ?  " Pass" : " Fail") << endl;
    cout << "PopLR p: " << poplr_p << endl;
    cout << "Other p: " << as->get_min_p() << endl;
    if (!pass) {
        cout << "Size arrests: " << as->get_series().second.size() << endl;
        print("first fail ", first_fail);
        print("s ", as->get_series().second);
        print("a ", a.second);
    }
    print("pr gys:", as->get_pr_gys());
    print("pr yrs:", as->get_pr_yrs());
}


int main() {
    pr_tt_given_mt = read_csv(PTT, ArrestSeries::mt_domain_size, ArrestSeries::tt_domain_size, true);
    vector<vector<long double>> pr_mtlt17_given_tt_2d = read_csv(F17, 1, ArrestSeries::tt_domain_size, false);
    pr_mtlt17_given_tt = pr_mtlt17_given_tt_2d[0];

    tester(s1, a1, 1, 5, "s1 ");
    tester(s2, a2, 9.10979e-16, 4, "s2 ");
    tester(s3, a3, 1, 4, "s3 ");
    tester(s4, a4, 1, 0, "s4 ");
    tester(s5, a5, 5.04921e-05, 2, "s5 ");
    tester(s6, {}, 3.88467e-10, 0, "s6 ");

    for (string filename : vector<string> {"ip031.right_final_a.json"}) {
        vector<Eye> eyes = read_json(filename);
        //cout << "\nRead " << eyes.size() << " eyes from " << filename << endl;
        //cout << "\nEye 1: " << eyes[0].size() << " series" << endl;
        //print("series 1 xs:", eyes[0][0].first);
        //print("series 1 ys:", eyes[0][0].second);
        //print("series last ys:", eyes[0][0].second.back());
        tester(eyes[0][0], a_ip006, 4.51377e-05, 53, filename);
    }

    return 0;
}
