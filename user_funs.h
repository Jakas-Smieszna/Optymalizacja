#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix ff1T(matrix, matrix = NAN, matrix = NAN);
matrix ff1R(matrix, matrix = NAN, matrix = NAN);
matrix df1(double, matrix, matrix = NAN, matrix = NAN);
matrix ff2T(matrix, matrix = NAN, matrix = NAN);
matrix ff2R(matrix, matrix = NAN, matrix = NAN);
matrix df2(double, matrix, matrix, matrix);
matrix ff3T(matrix, matrix = 0, matrix = 0);
typedef double(*boundF)(matrix, double);
double g3T1(matrix, double);
double g3T2(matrix, double);
double g3T3(matrix, double);
matrix df3(double, matrix, matrix = 0, matrix = 0);
matrix ff3R(matrix, matrix = 0, matrix = 0);
matrix df3(matrix, matrix = 0, matrix = 0);
matrix ff4R(matrix, matrix = 0, matrix = 0);
matrix gf4T(matrix x, matrix ud1 = 0, matrix ud2 = 0);
matrix Hf4T(matrix x, matrix ud1 = 0, matrix ud2 = 0);
matrix ff3T_zewn(matrix, matrix = 0, matrix = 0);
matrix ff3T_wewn(matrix, matrix = 0, matrix = 0);
matrix ff4T(matrix, matrix = 0, matrix = 0);
