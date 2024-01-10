#ifndef POISSON_H
#define POISSON_H

#include <complex.h>

void thomas_algorithm(double complex *a, double complex *b, double complex *c, double complex *d, double complex *x, double complex *l, double complex *u, double complex *y, int N);

#endif