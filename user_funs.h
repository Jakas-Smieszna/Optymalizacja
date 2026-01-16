#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include"ode_solver.h"

#define A5 1// 1 lub 10 lub 100. Nie mam jak dac 4. argumentu do gg5T by dalo sie tego uzywac goldenem.

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
matrix ff4T(matrix, matrix = 0, matrix = 0);
matrix gf4T(matrix, matrix = 0, matrix = 0);
matrix Hf4T(matrix, matrix = 0, matrix = 0);
void load_data_lab4(matrix*& xD, matrix*& yD);
matrix ff4R(matrix, matrix = 0, matrix = 0);
matrix gf4R(matrix, matrix = 0, matrix = 0);
matrix zlotf4T(matrix a, matrix d, matrix x);
matrix zlotf4R(matrix a, matrix d, matrix x);
double poprawne4R(matrix theta);
matrix ff5T1(matrix, matrix, matrix);
matrix ff5T2(matrix, matrix, matrix);
matrix gg5T1(matrix, matrix, matrix);
matrix gg5T2(matrix, matrix, matrix);
matrix ff5TX(matrix, matrix, matrix);
matrix gg5TX(matrix, matrix, matrix);
matrix ff6T(matrix, matrix, matrix);
