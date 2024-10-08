/**
 * @file turbulence.h
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Functions for adding turbulence to the wildfire simulation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef TURBULENCE_H
#define TURBULENCE_H

#include <math.h>
#include <stdlib.h>
#include "utils.h"
#include "functions.h"
#include "structures.h"

/**
 * @brief Calculates the turbulence of a given set of parameters.
 *
 * This function calculates the turbulence based on the given parameters and updates the value of R_turbulence.
 *
 * @param R_turbulence Pointer with partial derivatives used to calculate the turbulence.
 * @param R_new Pointer to the new value.
 * @param parameters Pointer to the Parameters struct containing the necessary parameters.
 */
void turbulence(double *R_turbulence, double *R_new, Parameters *parameters);

void turbulence_debug(double *R_turbulence, double *R_new, Parameters *parameters);

void turbulence_v1(double *R_turbulence, double *R_new, Parameters *parameters); 

void turbulence_v1_debug(double *R_turbulence, double *R_new, Parameters *parameters); 

#endif