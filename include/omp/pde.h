/**
 * @file pde.h
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Functions for solving the partial differential equations of the wildfire simulation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef PDE_H
#define PDE_H

#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <sys/time.h>
#include "../c/functions.h"
#include "../c/logs.h"
#include "../c/output.h"
#include "../c/poisson.h"
#include "../c/turbulence.h"
#include "../c/structures.h"
#include "../c/utils.h"

/**
 * @brief Calculates the right-hand side (RHS) of the partial differential equation (PDE).
 *
 * This function calculates the RHS of the PDE at a given time `t` using the provided
 * arrays `R_old`, `R_new`, and `R_turbulence`, as well as the `Parameters` struct.
 *
 * @param t The current time.
 * @param R_old Pointer to the array containing the old values of R.
 * @param R_new Pointer to the array where the new values of R will be stored.
 * @param R_turbulence Pointer to the array containing the turbulence values of R.
 * @param parameters Pointer to the Parameters struct containing additional parameters.
 */
void RHS(double t, double *R_old, double *R_new, double *R_turbulence, Parameters *parameters);

/**
 * @brief Applies boundary conditions to the given array.
 *
 * This function applies boundary conditions to the array `R_new` based on the provided `parameters`.
 *
 * @param R_new The array to which the boundary conditions will be applied.
 * @param parameters The parameters used to determine the boundary conditions.
 */
void boundary_conditions(double *R_new, Parameters *parameters);

/**
 * @brief Function evaluation for method of lines (MOL) time integration.
 *
 * This function compute the RHS, add the turbulence and add the boundary conditions.
 *
 * @param t The current time.
 * @param R_old Pointer to the array containing the old values of R.
 * @param R_new Pointer to the array where the new values of R will be stored.
 * @param U_turbulence Pointer to the array containing the turbulence values of U.
 * @param parameters Pointer to the Parameters struct containing additional parameters.
 */
void Phi(double t, double *R_old, double *R_new, double *U_turbulence, Parameters *parameters);

/**
 * @brief Applies velocity correction of Chorin's projection method to the given arrays.
 *
 * This function applies velocity correction using Chorin's projection method to the arrays. 
 * Also, it applies the boundary conditions to the velocity field and bounds the temperature and fuel values.
 *
 * @param R_new Pointer to the array representing the new values.
 * @param p Pointer to the array representing the previous values.
 * @param dt The time step used for the velocity correction.
 * @param parameters Pointer to the structure containing the parameters.
 */
void velocity_correction_boundaries_bounding(double *R_new, double *p, double dt, Parameters *parameters);

#endif