/**
 * @file poisson.h
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Functions for solving the Poisson equation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef POISSON_H
#define POISSON_H

#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include "structures.h"
#include "utils.h"

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
 * @param N The size of the linear system.
 */
void thomas_algorithm(fftw_complex *a, fftw_complex *b, fftw_complex *c, fftw_complex *d, fftw_complex *x, fftw_complex *l, fftw_complex *u, fftw_complex *y, int N);

/**
 * @brief Solves the pressure equation.
 *
 * This function solves the pressure equation using FFT-FD method.
 * FFT for periodic domain and FD for non-periodic domain.
 *
 * @param U The input array U with the velocity field.
 * @param p The output array p with the pressure field.
 * @param a The array a with the lower diagonal elements of the tridiagonal matrix.
 * @param b The array b with the main diagonal elements of the tridiagonal matrix.
 * @param c The array c with the upper diagonal elements of the tridiagonal matrix.
 * @param d The array d with the right-hand side vector.
 * @param l The array l for temporary storage of the lower diagonal elements, L in LU decomposition.
 * @param u The array u for temporary storage of the main diagonal elements, U in LU decomposition.
 * @param y The array y for temporary storage of the intermediate solution during forward substitution.
 * @param pk The array pk to store the solution of the Poisson equation.
 * @param p_plan The FFTW plan for p.
 * @param f_plan The FFTW plan for f.
 * @param p_top_plan The FFTW plan for p_top.
 * @param f_in The input array for f.
 * @param f_out The output array for f.
 * @param p_top_in The input array for p_top.
 * @param p_top_out The output array for p_top.
 * @param p_in The input array for p.
 * @param p_out The output array for p.
 * @param parameters The parameters for the solve_pressure function.
 */
void solve_pressure(double *U, double *p, fftw_complex *a, fftw_complex *b, fftw_complex *c, fftw_complex *d, fftw_complex *l, fftw_complex *u, fftw_complex *y, fftw_complex *pk, fftw_plan p_plan, fftw_plan f_plan, fftw_plan p_top_plan, fftw_complex *f_in, fftw_complex *f_out, fftw_complex *p_top_in, fftw_complex *p_top_out, fftw_complex *p_in, fftw_complex *p_out, Parameters *parameters);


#endif