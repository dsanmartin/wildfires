#ifndef OUTPUT_H
#define OUTPUT_H

#include <stdio.h>
#include "utils.h"

void save_data(char *path, double *x, double *y, double *z, double *f, int Nx, int Ny, int Nz);
void save_data_periodic_xy(char *path, double *x, double *y, double *z, double *f, int Nx, int Ny, int Nz);

#endif