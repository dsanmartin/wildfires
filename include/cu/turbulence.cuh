/**
 * @file turbulence.cuh
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Functions for adding turbulence to the wildfire simulation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef TURBULENCE_CUH
#define TURBULENCE_CUH

#include <math.h>
#include <stdlib.h>
#include "../c/structures.h"
#include "utils.cuh"
#include "functions.cuh"


/**
 * @brief Calculates the turbulence of a given set of parameters.
 *
 * This function calculates the turbulence based on the given parameters and updates the value of R_turbulence.
 *
 * @param R_turbulence Pointer with partial derivatives used to calculate the turbulence.
 * @param R_new Pointer to the new value.
 * @param parameters Pointer to the Parameters struct containing the necessary parameters.
 */
__global__
void turbulence(double *R_turbulence, double *R_new, double *z, double *z_ibm, Parameters parameters);

#endif