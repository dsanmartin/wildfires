/**
 * @file pde.cuh
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Functions to evaluate RHS of the PDE
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef PDE_CUH
#define PDE_CUH

#include "../c/structures.h"
#include "functions.cuh"
#include "turbulence.cuh"
#include "utils.cuh"

 /**
 * @brief Applies boundary conditions to the given array.
 *
 * This function applies boundary conditions to the array `R_new` based on the provided `parameters`.
 *
 * @param R_new The array to which the boundary conditions will be applied.
 * @param z The array containing the z values.
 * @param Nz_Y The number of nodes in the z direction for fuel fraction Y.
 * @param cut_nodes The array containing the cut nodes.
 * @param parameters Parameters struct containing additional parameters.
 */
__global__
void boundary_conditions(double *R_new, double *z, int *Nz_Y, int *cut_nodes, Parameters parameters);

/**
 * @brief Bounds the temperature and fuel arrays.
 *
 * This function modifies the given temperature and fuel arrays to enforce the
 * specified bounding values.
 *
 * @param R_new The array of temperature values and fuel.
 * @param parameters Parameters struct containing additional parameters.
 */
__global__
void bounds(double *R_new, Parameters parameters);

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
 * @param z Pointer to the array containing the z values.
 * @param z_ibm Pointer to the array containing the z values of the immersed boundary.
 * @param Nz_Y Pointer to the number of nodes in the z direction for fuel fraction Y.
 * @param parameters Parameters struct containing additional parameters.
 */
__global__
void RHS(double t, double *R_old, double *R_new, double *R_turbulence, double *z, double *z_ibm, int *Nz_Y, Parameters parameters);

 /**
  * @brief Applies velocity correction of Chorin's projection method to the given arrays.
  *
  * This function applies velocity correction using Chorin's projection method to the arrays. 
  * Also, it applies the boundary conditions to the velocity field and bounds the temperature and fuel values.
  *
  * @param R_new Pointer to the array representing the new values.
  * @param p Pointer to the array representing the previous values.
  * @param z Pointer to the array containing the z values.
  * @param fd_z Second-order finite difference option in the z direction. -1 backward, 0 centered, 1 forward.
  * @param parameters Pointer to the structure containing the parameters.
  */
  __global__
  void velocity_correction(double *R_new, double *p, double *z, int fd_z, Parameters parameters);
  
/**
 * @brief Function evaluation for method of lines (MOL) time integration.
 *
 * This function compute the RHS, add the turbulence and add the boundary conditions.
 *
 * @param t The current time.
 * @param R_old Pointer to the array containing the old values of R.
 * @param R_new Pointer to the array where the new values of R will be stored.
 * @param R_turbulence Pointer to the array containing the turbulence values of the approximation.
 * @param z Pointer to the array containing the z values.
 * @param z_ibm Pointer to the array containing the z values of the immersed boundary.
 * @param Nz_Y Pointer to the number of nodes in the z direction for fuel fraction Y.
 * @param cut_nodes Pointer to the array containing the cut nodes.
 * @param parameters Parameters struct containing additional parameters.
 */
 void Phi(double t, double *R_old, double *R_new, double *R_turbulence, double *z, double *z_ibm, int *Nz_Y, int *cut_nodes, Parameters parameters);

 
#endif