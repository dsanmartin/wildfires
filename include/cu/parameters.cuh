/**
 * @file parameters.h
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Functions for reading the simulation parameters from a file.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#define MAX_LINE_LENGTH 512

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h> 
#include "../c/structures.h"
#include "utils.cuh"

/**
 * @brief Reads the simulation parameters from a file.
 * 
 * This function reads the simulation parameters from a file in the given path.
 * 
 * @param file_path The path to the file where the parameters are stored.
 * @return The Parameters struct containing the simulation parameters.
 */
Parameters read_parameters_file(const char *file_path); 

/**
 * @brief Free parameters memory.
 * 
 */
void free_parameters(Parameters parameters);

#endif