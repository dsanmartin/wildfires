#ifndef PDE_H
#define PDE_H

#include <stdlib.h>
#include "utils.h"
#include "structures.h"
#include "functions.h"
#include "output.h"

void Phi(double t, double *R_old, double *R_new, Parameters *parameters);
void euler(double t_n, double *y_n, double *y_np1, double *F, double dt, int size, Parameters *parameters);
void create_y_0(double *u, double *v, double *w, double *T, double *Y, double *y_0, Parameters *parameters);
void solve_PDE(double *y_n, Parameters *parameters);

#endif