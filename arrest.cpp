
#include "arrest.hpp"

    // encoding ARREST values into dB values
double R() { return( -999); }
double Y(double db) { return(-100 - db); }
double O(double db) { return( -200 - db); }

    // decoding ARREST values into dB values
bool is_green(double x) { return(x > -60); }
bool is_yellow_orange(double x) { return(x < -60 && x > -999); }
bool is_orange(double x) { return(x < -160 && x > -999); }
bool is_red(double x) { return(x == -999); }

double get_db(double x) { return (is_green(x) ? x : (is_red(x) ? -1 : GET_DB_FROM_YELLOW(x))); }
