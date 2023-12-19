#ifndef STRUCTURES_H
#define STRUCTURES_H

typedef struct _parameters {
    int Nx;
    int Ny;
    int Nz;
    int Nt;
    double x_min;
    double x_max;
    double y_min;
    double y_max;
    double z_min;
    double z_max;
    double t_min;
    double t_max;
    double dx;
    double dy;
    double dz;
    double dt;
} Parameters;

#endif