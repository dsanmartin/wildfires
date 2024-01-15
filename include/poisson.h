#ifndef POISSON_H
#define POISSON_H

#include <complex.h>
#include "structures.h"

void solve_pressure(double *U, double *P_RHS, Parameters *parameters);
void thomas_algorithm(double complex *a, double complex *b, double complex *c, double complex *d, double complex *x, double complex *l, double complex *u, double complex *y, int N);
void fftfd(double complex *f, double complex *F, Parameters *parameters, int N);

#endif