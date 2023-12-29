#ifndef OUTPUT_H
#define OUTPUT_H

#include <stdio.h>
#include "utils.h"
#include "structures.h"

void save_scalar(char *path, double *x, double *y, double *z, double *f, int Nx, int Ny, int Nz);
void save_data(char *path, double *data, int n, Parameters *parameters);

#endif