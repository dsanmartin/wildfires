/**
 * @file pressure.h
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Functions for solving the Pressure equation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef PRESSURE_H
#define PRESSURE_H

#include <cufftw.h>
#include <math.h>
#include <stdlib.h>
#include "../c/structures.h"
#include "functions.cuh"
#include "logs.cuh"
#include "utils.cuh"

__global__
void compute_f(double *U, double *p, double *z, cufftDoubleComplex *f_in, cufftDoubleComplex *p_top_in, Parameters parameters);

__global__
void compute_f_density(double *y_np1, double *y_n, double *p, double *z, cufftDoubleComplex *f_in, cufftDoubleComplex *p_top_in, Parameters parameters);

__global__
void gammas_and_coefficients(double *kx, double *ky, double *gamma, double *a, double *b, double *c, double *z, Parameters parameters);

__global__
void post_fft(double *p, cufftDoubleComplex *p_out, Parameters parameters);

/**
 * @brief Implements the Thomas algorithm for solving a tridiagonal linear system.
 *
 * This file contains the implementation of the Thomas algorithm, also known as the tridiagonal matrix algorithm,
 * for solving a tridiagonal linear system of equations. The algorithm consists of three steps: decomposition,
 * forward substitution, and backward substitution.
 *
 * @param a The lower diagonal elements of the tridiagonal matrix.
 * @param b The main diagonal elements of the tridiagonal matrix.
 * @param c The upper diagonal elements of the tridiagonal matrix.
 * @param d The right-hand side vector.
 * @param x The solution vector.
 * @param l Temporary array to store the values of the lower diagonal elements during decomposition.
 * @param u Temporary array to store the values of the main diagonal elements during decomposition.
 * @param y Temporary array to store the values of the intermediate solution during forward substitution.
 * @param parameters The parameters of the mathematical model.
 */
 __global__
void thomas_algorithm(double *a, double *b, double *c, cufftDoubleComplex *d, cufftDoubleComplex *x, cufftDoubleComplex *l, cufftDoubleComplex *u, cufftDoubleComplex *y, Parameters parameters);

__global__ 
void update_coefficients(double *gamma, double *b, double *z, cufftDoubleComplex *d, cufftDoubleComplex *f_out, cufftDoubleComplex *p_top_out, Parameters parameters);

/**
 * @brief Solves the pressure equation.
 *
 * This function solves the pressure equation using FFT-FD method.
 * FFT for periodic domain and FD for non-periodic domain.
 *
 * @param U The input array with the velocity field.
 * @param p The output array with the pressure field.
 * @param gamma The array with FFT-FD coefficients.
 * @param a The array with the lower diagonal elements of the tridiagonal matrix.
 * @param b The array with the main diagonal elements of the tridiagonal matrix.
 * @param c The array with the upper diagonal elements of the tridiagonal matrix.
 * @param d The array with the right-hand side vector.
 * @param l The array for temporary storage of the lower diagonal elements, L in LU decomposition.
 * @param u The array for temporary storage of the main diagonal elements, U in LU decomposition.
 * @param y The array for temporary storage of the intermediate solution during forward substitution.
 * @param data_in The input array for FFT .
 * @param data_out The output array for IFFT.
 * @param p_top_in The input array for p_top FFT.
 * @param p_top_out The output array for p_top IFFT.
 * @param parameters The parameters of the mathematical model.
 */
void solve_pressure(double *y_np1, double *y_n, double *p, double *z, double *gamma, double *a, double *b, double *c, cufftDoubleComplex *d, cufftDoubleComplex *l, cufftDoubleComplex *u, cufftDoubleComplex *y, cufftDoubleComplex *data_in, cufftDoubleComplex *data_out, cufftDoubleComplex *p_top_in, cufftDoubleComplex *p_top_out, Parameters parameters);

void solve_pressure_iterative(double *y_np1, double *y_n, double *p, double *z, double *gamma, double *a, double *b, double *c, cufftDoubleComplex *d, cufftDoubleComplex *l, cufftDoubleComplex *u, cufftDoubleComplex *y, cufftDoubleComplex *data_in, cufftDoubleComplex *data_out, cufftDoubleComplex *p_top_in, cufftDoubleComplex *p_top_out, Parameters parameters, double *error, int *max_iter);

#endif