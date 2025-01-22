/**
 * @file solver.cuh
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Functions for solving the partial differential equations of the wildfire simulation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef SOLVER_CUH
#define SOLVER_CUH

#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include "../c/structures.h"
#include "output.cuh"
#include "pressure.cuh"
#include "functions.cuh"
#include "logs.cuh"
#include "pde.cuh"
#include "utils.cuh"

/**
 * Performs a single Euler step for solving a partial differential equation (PDE).
 *
 * @param dt The time step size.
 * @param y_n The current state of the solution.
 * @param y_np1 The updated state of the solution after the Euler step.
 * @param F The function used in the PDE.
 * @param size The size of the solution vector.
 */
__global__
void euler_step(double dt, double *y_n, double *y_np1, double *F, int size);

/**
 * @brief Performs an Euler integration step for a partial differential equation (PDE).
 *
 * This function calculates the solution at time t_n+1 using the Euler method, given the solution at time t_n,
 * the function evaluation F, the turbulence vector U_turbulence, the time step dt, the size of the problem, and the parameters.
 *
 * @param t_n The current time.
 * @param y_n The solution vector at time t_n.
 * @param y_np1 The solution vector at time t_n+1 (output parameter).
 * @param F The function vector.
 * @param U_turbulence The turbulence vector.
 * @param dt The time step.
 * @param size The size of the problem.
 * @param parameters The parameters for the PDE.
 */
void euler(double t_n, double *y_n, double *y_np1, double *F, double *U_turbulence, double *z_ibm, int *Nz_Y, int *cut_nodes, double dt, int size, Parameters parameters);

/**
 * @brief Performs a single step of the second-order Runge-Kutta method (RK2).
 *
 * This function calculates the next state of the system using the RK2 method.
 *
 * @param dt The time step size.
 * @param y_n The current state of the system.
 * @param y_np1 The next state of the system (output).
 * @param k1 Temporary variable for intermediate calculations.
 * @param k2 Temporary variable for intermediate calculations.
 * @param size The size of the arrays.
 */
__global__
void RK2_step(double dt, double *y_n, double *y_np1, double *k1, double *k2, int size);

/**
 * @brief Performs a second-order Runge-Kutta (RK2) integration step.
 *
 * This function calculates the values of y_np1 at time t_n+1 using the values of y_n at time t_n,
 * along with other parameters such as k, F, U_turbulence, dt, size, and parameters.
 *
 * @param t_n The current time.
 * @param y_n The array of values at time t_n.
 * @param y_np1 The array to store the calculated values at time t_n+1.
 * @param k The array used for intermediate calculations.
 * @param F The array of function evaluation.
 * @param U_turbulence The array of turbulence values.
 * @param dt The time step size.
 * @param size The size of the arrays.
 * @param parameters The struct containing additional parameters.
 */
void RK2(double t_n, double *y_n, double *y_np1, double *k, double *F, double *U_turbulence, double *z_ibm, int *Nz_Y, int *cut_nodes, double dt, int size, Parameters parameters);

/**
 * @brief Performs a single step of the fourth-order Runge-Kutta method (RK4).
 *
 * This function calculates the next state of the system using the RK4 method.
 *
 * @param dt The time step size.
 * @param y_n The current state of the system.
 * @param y_np1 The next state of the system (output).
 * @param k1 Temporary variable for intermediate calculations.
 * @param k2 Temporary variable for intermediate calculations.
 * @param k3 Temporary variable for intermediate calculations.
 * @param k4 Temporary variable for intermediate calculations.
 * @param size The size of the arrays.
 */
__global__
void RK4_step(double dt, double *y_n, double *y_np1, double *k1, double *k2, double *k3, double *k4, int size);

/**
 * @brief Performs a fourth-order Runge-Kutta (RK4) integration step.
 *
 * This function calculates the values of y_np1 at time t_n+1 using the values of y_n at time t_n,
 * along with other parameters such as k, F, U_turbulence, dt, size, and parameters.
 *
 * @param t_n The current time.
 * @param y_n The array of values at time t_n.
 * @param y_np1 The array to store the calculated values at time t_n+1.
 * @param k The array used for intermediate calculations.
 * @param F The array of function evaluation.
 * @param U_turbulence The array of turbulence values.
 * @param dt The time step size.
 * @param size The size of the arrays.
 * @param parameters The struct containing additional parameters.
 */
void RK4(double t_n, double *y_n, double *y_np1, double *k, double *F, double *U_turbulence, double *z_ibm, int *Nz_Y, int *cut_nodes, double dt, int size, Parameters parameters);

/**
 * @brief Creates the initial condition for the variable y_0.
 *
 * This function takes arrays of variables u, v, w, T, Y, and computes the initial condition for the variable y_0.
 * The computed result is stored in the array y_0.
 *
 * @param u         Array of variables u.
 * @param v         Array of variables v.
 * @param w         Array of variables w.
 * @param T         Array of variables T.
 * @param Y         Array of variables Y.
 * @param y_0       Array to store the computed initial condition for y_0.
 * @param parameters Pointer to the Parameters struct containing additional parameters.
 */
void create_y_0(double *u, double *v, double *w, double *T, double *Y, double *y_0, Parameters parameters);

/**
 * @brief Solves the partial differential equation (PDE) using the method of lines (MOL).
 *
 * This function solves the PDE using the method of lines (MOL) with the given initial conditions and parameters.
 *
 * @param y_n       Array containing the initial condition for the PDE.
 * @param p_0       Array containing the initial condition for the pressure field.
 * @param parameters Pointer to the Parameters struct containing additional parameters.
 */
void solve_PDE(double *y_n, double *p_0, Parameters parameters);

#endif