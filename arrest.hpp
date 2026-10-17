#ifndef ARREST_HPP
#define ARREST_HPP

    // encoding ARREST values into dB values
double R();
double Y(double db);
double O(double db);

    // decoding ARREST values into dB values
bool is_green(double x);
bool is_yellow_orange(double x);
bool is_orange(double x);
bool is_red(double x);

#define GET_DB_FROM_YELLOW(x) ((x) < -160 ? -(x) - 200 : -(x) - 100)

double get_db(double x);

#endif // ARREST_HPP
