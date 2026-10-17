#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>

#include "arrestSeries.hpp"
#include "read_utils.hpp"

using namespace std;

/*
Read in JSON file that has series of VFs. Format is
[                   # whole lot
  [                   # begin eye 1
    [                   # begin repeat 1
      [0,0.4,1.2],        # x values in years (length "num_visits")
      [20,19,18],         # location 1
      [20,19,18]          # location 2
    ],
    [
        [0,0.4,1.2],        # x values in years (length "num_visits")
        [21,19,18],     # begin repeat 2
        [21,19,18]        # location 1
    ]                     # location 2
  ],
  [ eye 2...]
]

Hand crafted parser with simple state.
Ignores ' ' and '\n'.
Converts any ".*" to nan.
*/
vector<Eye> read_json(const string& filename) {
    FILE *file = fopen(filename.c_str(), "r");

    vector<Eye> eyes;

    #define ADD_NUM {\
        if (loc == 0) \
            eyes[eye][ser].first.push_back(sign * num); \
        else \
            eyes[eye][ser].second[loc].push_back(sign * num); \
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
//cout << (char)ch << " " << "st " << state << " e " << eye << " s " << ser << " l " << loc << endl;
        if (ch == '[') {
            switch (state) { 
            case 0: { state++; eye = -1; break; }
            case 1: { state++; eye++; eyes.push_back(Eye()); ser = -1; break; }
            case 2: { state++; ser++; eyes[eye].push_back(Series()); loc = -1; break; }
            case 3: { state++; loc++; eyes[eye][ser].second.push_back(Location()); break; }
            default: {
                cout << "parse ERROR [" << endl;
                return(vector<Eye>());
            }};
        } else if (ch == ']') {
            switch (state) { 
            case 0: { cout << "parse ERROR [" << endl; break;}
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
            num = numeric_limits<double>::quiet_NaN();
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

    // The first one of each series.second will be empty as we put that in series.first
    for (int eye = 0; eye < eyes.size(); eye++) {
        for (int ser = 0; ser < eyes[eye].size(); ser++) {
            if (!eyes[eye][ser].second[0].size() == 0) {
                perror(("Non-empty first location in series " + to_string(ser) + " of eye " + to_string(eye)).c_str());
                exit(-1);
            }
            eyes[eye][ser].second.erase(eyes[eye][ser].second.begin());
        }
    }
    return eyes;
}

/*
    Read csv file with nrows row, ncols columns

    @param filename Filename of a CSV file with nrows row and ncols columns
    @param nrows Number of rows to read
    @param ncols Number of cols to read
    @param row_sum_check Whether to check if rows sum to 1

    @return vector<vector<long double>>  of size nrows x ncols
*/
vector<vector<long double>> read_csv(const std::string filename, int nrows, int ncols, bool row_sum_check) {
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        perror(("Error opening file " + filename).c_str());
        return vector<vector<long double>>(); // return empty vector on error
    }

    vector<vector<long double>> data(nrows, vector<long double>(ncols));

    int row = 0;
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;

        int col = 0;
        while (std::getline(ss, cell, ',')) {
            if (cell ==  "NaN") {
                data[row][col] = std::numeric_limits<double>::quiet_NaN();
            } else {
                data[row][col] = strtold(cell.c_str(), nullptr);
            }
            col++;
        }
        if (col != ncols) {
            fprintf(stderr, "Row %d in %s has incorrect number of columns: %d expected %d\n", row, filename.c_str(), col, ncols);
            file.close();
            return vector<vector<long double>>(); // return empty vector on error
        }
        row++;
    }

    file.close();

    if (row != nrows) {
        fprintf(stderr, "File %s has incorrect number of rows: %d expected %d\n", filename.c_str(), row, nrows);
        return vector<vector<long double>>(); // return empty vector on error
    }

    // check rows sum to 1
    if (row_sum_check)
    for (int row = 0; row < nrows; row++) {
        long double row_sum = 0;
        for (int col = 0; col < ncols; col++)
            row_sum += data[row][col];
        if (abs(row_sum - 1.0) > 1e-12) {
            fprintf(stderr, "Row %d does not sum to 1 in %s: %Lf\n", row, filename.c_str(), row_sum);
            return vector<vector<long double>>(); // return empty vector on error
        }
    }

    return data;
}