#ifndef OUTPUT_H
#define OUTPUT_H

#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "structures.h"

void save_scalar(char *path, double *x, double *y, double *z, double *f, int Nx, int Ny, int Nz);
void save_data_txt(char *path, double *data, double *p, int n, Parameters *parameters);
void save_data(char *path, double *data, double *p, int n, Parameters *parameters);
void save_domain(char *path, Parameters *parameters);
void save_time(char *save_path, double *t, int Nt, int NT);

#endif