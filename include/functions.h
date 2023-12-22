#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <math.h>
#include "structures.h"
#include "utils.h"

double power_law(double z, double u_r, double z_r, double alpha_u);
void power_law_initial_condition(double *x, double *y, double *z, double *u, double *v, double *w, Parameters *parameters);

#endif