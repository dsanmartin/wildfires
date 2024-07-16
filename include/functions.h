#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <math.h>
#include "structures.h"
#include "utils.h"

double power_law(double z, double u_r, double z_r, double alpha_u);
double gaussian(double x, double y, double z, double x_0, double y_0, double z_0, double sx, double sy, double sz);
double K(double T, double A, double T_a);
// double H(double x);
double H(double x, double T_pc);
double f_damping(double z, double u_tau, double nu);
double source(double T, double Y, double H_R, double A, double T_a, double h, double a_v, double T_inf, double c_p, double rho, double T_pc);
void power_law_initial_condition(double *x, double *y, double *z, double *u, double *v, double *w, Parameters *parameters);
void gaussian_temperature_initial_condition(double *x, double *y, double *z, double *T, Parameters *parameters);
void fuel_initial_condition(double *x, double *y, double *z, double *Y, Parameters *parameters);
void initial_conditions(double *x, double *y, double *z, double *u, double *v, double *w, double *T, double *Y, double *p, Parameters *parameters);


#endif