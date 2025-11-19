//Ten plik nie powinien byc edytowany

#ifndef _H_ODE_SOLVER
#define _H_ODE_SOLVER

#include"matrix.h"
#include"user_funs.h"

matrix* solve_ode(matrix(*)(double, matrix, matrix, matrix), double, double, double, matrix, matrix = NAN, matrix = NAN); // throw (string);

#endif