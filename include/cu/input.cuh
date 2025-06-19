/**
 * @file input.cuh
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Procedures for input data.
 * @version 0.1
 * @date 2025-06-19
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#ifndef INPUT_CUH
#define INPUT_CUH

#define FILENAME_SIZE 128

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "logs.cuh"
#include "../c/structures.h"

void initial_conditions_input(double *u, double *v, double *w, double *T, double *Y, double *p, Parameters parameters);

#endif