/**
 * @file ibm.h
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Implementation of the Immersed Boundary Method (IBM) for the wildfire simulation.
 * @version 0.1
 * @date 2024-08-20
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef IBM_CUH
#define IBM_CUH

#include "../c/structures.h"
#include "functions.cuh"
#include "utils.cuh"

void flat_terrain(Parameters *parameters);

void simple_hill(Parameters *parameters);

void cube(Parameters *parameters);

void ibm_parameters(Parameters *parameters);

#endif