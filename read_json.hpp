#ifndef READ_JSON_HPP
#define READ_JSON_HPP

#include <string>
#include <vector>
#include <iostream>
#include <limits>

#include "poplr_arrest.hpp"

using namespace std;

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
vector<Eye> read_json(const string& filename) {
    FILE *file = fopen(filename.c_str(), "r");

    vector<Eye> eyes;

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
//cout << (char)ch << " " << "st " << state << " e " << eye << " s " << ser << " l " << loc << endl;
        if (ch == '[') {
            switch (state) { 
            case 0: { state++; eye = -1; break; }
            case 1: { state++; eye++; eyes.push_back(Eye()); ser = -1; break; }
            case 2: { state++; ser++; eyes[eye].push_back(Series()); loc = -1; break; }
            case 3: { state++; loc++; eyes[eye][ser].push_back(Location()); break; }
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

    return eyes;
}

#endif