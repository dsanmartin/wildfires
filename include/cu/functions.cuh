/**
 * @file functions.cuh
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Implementation of various functions used in the wildfire simulation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef FUNCTIONS_CUH
#define FUNCTIONS_CUH

#include <math.h>
#include "../c/structures.h"
#include "utils.cuh"

/**
 * @brief Calculates the power law wind profile.
 * 
 * @param z The height above ground level.
 * @param u_r The reference wind speed.
 * @param z_r The reference height.
 * @param alpha_u The power law exponent.
 * @return The wind speed at the given height.
 */
double power_law(double z, double u_r, double z_r, double alpha_u);

/**
 * @brief Calculates the Gaussian function.
 * 
 * @param x The x-coordinate.
 * @param y The y-coordinate.
 * @param z The z-coordinate.
 * @param x_0 The x-coordinate of the center.
 * @param y_0 The y-coordinate of the center.
 * @param z_0 The z-coordinate of the center.
 * @param sx The standard deviation in the x-direction.
 * @param sy The standard deviation in the y-direction.
 * @param sz The standard deviation in the z-direction.
 * @return The value of the Gaussian function at the given coordinates.
 */
double gaussian(double x, double y, double z, double x_0, double y_0, double z_0, double sx, double sy, double sz);

/**
 * @brief Calculates the temperature-dependent reaction rate constant.
 * 
 * @param T The temperature.
 * @param A The pre-exponential factor.
 * @param T_a The activation energy.
 * @return The reaction rate constant.
 */
double K(double T, double A, double T_a);

/**
 * @brief Calculates the Heaviside step function.
 * 
 * @param x The input value.
 * @param T_pc The threshold value.
 * @return 1 if x > T_pc, 0 otherwise.
 */
double H(double x, double T_pc);

/**
 * @brief Calculates the damping factor for the wind velocity.
 * 
 * @param z The height above ground level.
 * @param u_tau The friction velocity.
 * @param nu The kinematic viscosity.
 * @return The damping factor.
 */
double f_damping(double z, double u_tau, double nu);

/**
 * @brief Calculates the source term for the temperature equation.
 * 
 * @param T The temperature.
 * @param Y The fuel mass fraction.
 * @param H_R The heat of reaction.
 * @param A The pre-exponential factor.
 * @param T_a The activation energy.
 * @param h The heat transfer coefficient.
 * @param a_v The vegetation area index.
 * @param T_inf The ambient temperature.
 * @param c_p The specific heat capacity.
 * @param rho_inf The air density.
 * @param T_pc The threshold temperature for combustion.
 * @return The source term.
 */
double source(double T, double Y, double H_R, double A, double T_a, double h, double a_v, double T_inf, double c_p, double rho_inf, double T_pc);

/**
 * @brief Calculates the timestep reports based on the given parameters.
 *
 * @param y_n The array of values.
 * @param CFL The CFL (Courant-Friedrichs-Lewy) number.
 * @param Y_min The minimum value for Y (fuel).
 * @param Y_max The maximum value for Y (fuel).
 * @param T_min The minimum value for T (temperature).
 * @param T_max The maximum value for T (temperature).
 * @param parameters The pointer to the Parameters struct.
 */
void timestep_reports(double *y_n, double *CFL, double *Y_min, double *Y_max, double *T_min, double *T_max, Parameters parameters);

/**
 * @brief Sets the initial conditions for the simulation.
 * 
 * @param u The x-component of the velocity.
 * @param v The y-component of the velocity.
 * @param w The z-component of the velocity.
 * @param T The temperature.
 * @param Y The fuel mass fraction.
 * @param p The pressure.
 * @param parameters The simulation parameters.
 */
void initial_conditions(double *u, double *v, double *w, double *T, double *Y, double *p, Parameters parameters);

__global__
void temperature_source(double *x, double *y, double *z, double *y_n, Parameters paramenters);

__global__
void norm(double *v1, double *v2, double *result, double p, int size);

void equispaced_domain(double *z, int Nz, double z_min, double dz);

void non_equispaced_domain(double *z, int Nz, double z_min, double z_max, double k);

#endif