#ifndef POISSON_H
#define POISSON_H

#include <complex.h>
#include <fftw3.h>
#include "structures.h"
#include "utils.h"

void fftfd(double complex *f, double complex *F, Parameters *parameters, int N);
void prepare_data_fft(fftw_complex *A, double *B, Parameters *parameters);
void restore_data_fft(fftw_complex *B, double *A, Parameters *parameters);
void solve_pressure(double *U, double *P_RHS, Parameters *parameters);
void thomas_algorithm(double complex *a, double complex *b, double complex *c, double complex *d, double complex *x, double complex *l, double complex *u, double complex *y, int N);


#endif