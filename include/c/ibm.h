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

#ifndef IBM_H
#define IBM_H

#include "functions.h"
#include "structures.h"
#include "utils.h"

void flat_terrain(Parameters *parameters);

void simple_hill(Parameters *parameters);

void ibm_parameters(Parameters *parameters);

// void topography(Parameters parameters);

// void ibm_nodes(Parameters parameters);


#endif