#include "../include/functions.h"

double power_law(double z, double u_r, double z_r, double alpha_u) {
    return u_r * pow(z / z_r, alpha_u);
}

double gaussian(double x, double y, double z, double x_0, double y_0, double z_0, double sx, double sy, double sz) {
    return exp(-pow((x - x_0) / sx, 2) - pow((y - y_0) / sy, 2) - pow((z - z_0) / sz, 2));
}

double K(double T, double A, double T_a) {
    return A * exp(-T_a / T);
}

double H(double x) {
    if (x > 0) {
        return 1;
    } else {
        return 0;
    }
}

double source(double T, double Y, double H_R, double A, double T_a, double h, double a_v, double T_inf, double c_p, double rho) {
    return H_R * Y * K(T, A, T_a) * H(T - T_a) / c_p - h * a_v * (T - T_inf) / (c_p * rho);
}

void power_law_initial_condition(double *x, double *y, double *z, double *u, double *v, double *w, Parameters *parameters) {
    // Set u log wind, v and w as 0
    int Nx = parameters->Nx;// - 1;
    int Ny = parameters->Ny;// - 1;
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

void gaussian_temperature_initial_condition(double *x, double *y, double *z, double *T, Parameters *parameters) {
    // Set T as a gaussian
    int Nx = parameters->Nx;// - 1;
    int Ny = parameters->Ny;// - 1;
    int Nz = parameters->Nz;
    double x_0 = parameters->T0_x_center;
    double y_0 = parameters->T0_y_center;
    double z_0 = parameters->T0_z_center;
    double sx = parameters->T0_length;
    double sy = parameters->T0_width;
    double sz = parameters->T0_height;
    double T_hot = parameters->T_hot;
    double T_inf = parameters->T_inf;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                T[idx(i, j, k, Nx, Ny, Nz)] = T_inf + (T_hot - T_inf) * gaussian(x[i], y[j], z[k], x_0, y_0, z_0, sx, sy, sz);
            }
        }
    }
}