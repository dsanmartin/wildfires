#ifndef POISSON_H
#define POISSON_H

#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include "structures.h"
#include "utils.h"

void fftfd(double complex *f, double complex *F, Parameters *parameters, int N);
void prepare_data_fft(fftw_complex *A, double *B, Parameters *parameters);
void restore_data_fft(fftw_complex *B, double *A, Parameters *parameters);
// void solve_pressure_v1(double *U, double *p, Parameters *parameters);
// void solve_pressure(double *U, double *p, double complex *a, double complex *b, double complex *c, double complex *d, double complex *l, double complex *u, double complex *y, double complex *pk, fftw_complex *f_in, fftw_complex *f_out, fftw_complex *p_top_in, fftw_complex *p_top_out, fftw_complex *p_in, fftw_complex *p_out, Parameters *parameters);
// void solve_pressure(double *U, double *p, double complex *a, double complex *b, double complex *c, double complex *d, double complex *l, double complex *u, double complex *y, double complex *pk, Parameters *parameters);
void solve_pressure(double *U, double *p, double complex *a, double complex *b, double complex *c, double complex *d, double complex *l, double complex *u, double complex *y, double complex *pk, fftw_plan p_plan, fftw_plan f_plan, fftw_plan p_top_plan, fftw_complex *f_in, fftw_complex *f_out, fftw_complex *p_top_in, fftw_complex *p_top_out, fftw_complex *p_in, fftw_complex *p_out, Parameters *parameters);
void thomas_algorithm(double complex *a, double complex *b, double complex *c, double complex *d, double complex *x, double complex *l, double complex *u, double complex *y, int N);


#endif