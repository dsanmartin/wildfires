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
 __global__
 void RHS(double t, double *R_old, double *R_new, double *R_turbulence, double *z, double *z_ibm, int *Nz_Y, Parameters parameters);
 
 /**
  * @brief Applies boundary conditions to the given array.
  *
  * This function applies boundary conditions to the array `R_new` based on the provided `parameters`.
  *
  * @param R_new The array to which the boundary conditions will be applied.
  * @param parameters The parameters used to determine the boundary conditions.
  */
 __global__
 void boundary_conditions(double *R_new, int *Nz_Y, int *cut_nodes, Parameters parameters);
 
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
 void Phi(double t, double *R_old, double *R_new, double *R_turbulence, double *z, double *z_ibm, int *Nz_Y, int *cut_nodes, Parameters parameters);
 
 /**
  * @brief Bounds the temperature and fuel arrays.
  *
  * This function modifies the given temperature and fuel arrays to enforce the
  * specified bounding values.
  *
  * @param R_new The array of temperature values and fuel.
  * @param parameters The parameters struct containing information about the
  *                   computational domain and the boundary conditions.
  */
 __global__
 void bounds(double *R_new, Parameters parameters);
 
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
 __global__
 void velocity_correction(double *R_new, double *p, double dt, Parameters parameters);
 
 __global__
 void velocity_correction_fw(double *R_new, double *p, double dt, Parameters parameters);
 
 __global__
 void velocity_correction_bw(double *R_new, double *p, double dt, Parameters parameters);
 
 #endif