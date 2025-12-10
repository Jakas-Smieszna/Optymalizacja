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
matrix ff3R(matrix, matrix = 0, matrix = 0);
matrix df3(matrix, matrix = 0, matrix = 0);
matrix ff4R(matrix, matrix = 0, matrix = 0);
