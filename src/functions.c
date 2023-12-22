#include "../include/functions.h"
#include <stdio.h>

double power_law(double z, double u_r, double z_r, double alpha_u) {
    return u_r * pow(z / z_r, alpha_u);
}

void power_law_initial_condition(double *x, double *y, double *z, double *u, double *v, double *w, Parameters *parameters) {
    // Set u log wind, v and w as 0
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    double u_r = parameters->u_r;
    double z_r = parameters->z_r;
    double alpha_u = parameters->alpha_u;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                u[idx(i, j, k, Nx, Ny, Nz)] = power_law(z[k], u_r, z_r, alpha_u);
                v[idx(i, j, k, Nx, Ny, Nz)] = 0;
                w[idx(i, j, k, Nx, Ny, Nz)] = 0;
            }
        }
    }
}