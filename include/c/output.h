/**
 * @file output.h
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Functions for saving simulation data to files.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef OUTPUT_H
#define OUTPUT_H

#define FILENAME_SIZE 128

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structures.h"

/**
 * @brief Saves the simulation data to a file.
 * 
 * This function saves the simulation data to a file.
 * The data is saved in a binary format.
 * 
 * @param data       A pointer to the array containing the simulation data.
 * @param p          A pointer to the array containing the pressure data.
 * @param n          The number of spatial points.
 * @param parameters A pointer to the Parameters struct containing the simulation parameters.
 */
void save_data(double *data, double *p, int n, Parameters *parameters);

/**
 * @brief Saves the domain to a file.
 * 
 * This function saves the domain to a file in the parameter path.
 * The domain is saved in a binary format.
 * 
 * @param parameters A pointer to the Parameters struct containing the simulation parameters.
 */
void save_domain(Parameters *parameters);

void save_topography(Parameters *parameters);

#endif