/**
 * @file functions.c
 * @brief Implementation of various functions used in the wildfire simulation.
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 */

#include "../../include/cu/functions.cuh"

double power_law(double z, double u_r, double z_r, double alpha_u) {
    return u_r * pow(z / z_r, alpha_u);
}

double gaussian(double x, double y, double z, double x_0, double y_0, double z_0, double sx, double sy, double sz) {
    return exp(-pow((x - x_0) / sx, 2.0) - pow((y - y_0) / sy, 2.0) - pow((z - z_0) / sz, 2.0));
}

double K(double T, double A, double T_a) {
    return A * exp(-T_a / T);
}

double H(double x, double T_pc) {
    if (x > T_pc) {
        return 1.0;
    } else {
        return 0.0;
    }
}

double f_damping(double z, double u_tau, double nu) {
    return 1 - exp(-z * u_tau / 25 / nu);
}

double source(double T, double Y, double H_R, double A, double T_a, double h, double a_v, double T_inf, double c_p, double rho, double T_pc) {
    return H_R * Y * K(T, A, T_a) * H(T, T_pc) / c_p - h * a_v * (T - T_inf) / (c_p * rho);
}

void timestep_reports(double *y_n, double *CFL, double *Y_min, double *Y_max, double *T_min, double *T_max, Parameters parameters) {
    // printf("Timestep reports\n");
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int Nz_Y_max = parameters.Nz_Y_max;
    int u_index = parameters.field_indexes.u;
    int v_index = parameters.field_indexes.v;
    int w_index = parameters.field_indexes.w;
    int T_index = parameters.field_indexes.T;
    int Y_index = parameters.field_indexes.Y;
    // int size = Nx * Ny * Nz;
    double dx = parameters.dx;
    double dy = parameters.dy;
    double dz = parameters.dz;
    double dt = parameters.dt;
    double max_u = 0.0;
    double max_v = 0.0;
    double max_w = 0.0;
    double abs_u, abs_v, abs_w;
    // double CFL_tmp = 0.0;
    double Y_min_tmp = 0.0, Y_max_tmp = -1e9;
    double T_min_tmp = 1e9, T_max_tmp = -1e9;
    // int idx = threadIdx.x + blockIdx.x * blockDim.x;
    // int stride = gridDim.x * blockDim.x;
    // for (int ijk = idx; ijk < size; ijk += stride) {
        // int i = ijk / (Ny * Nz);
        // int j = (ijk % (Ny * Nz)) / Nz;
        // int k = ijk % Nz;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                abs_u = fabs(y_n[u_index + IDX(i, j, k, Nx, Ny, Nz)]);
                abs_v = fabs(y_n[v_index + IDX(i, j, k, Nx, Ny, Nz)]);
                abs_w = fabs(y_n[w_index + IDX(i, j, k, Nx, Ny, Nz)]);
                max_u = MAX(max_u, abs_u);
                max_v = MAX(max_v, abs_v);
                max_w = MAX(max_w, abs_w);
                if (k < Nz_Y_max) {
                    Y_min_tmp = MIN(Y_min_tmp, y_n[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y_max)]);
                    Y_max_tmp = MAX(Y_max_tmp, y_n[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y_max)]);
                }
                T_min_tmp = MIN(T_min_tmp, y_n[T_index + IDX(i, j, k, Nx, Ny, Nz)]);
                T_max_tmp = MAX(T_max_tmp, y_n[T_index + IDX(i, j, k, Nx, Ny, Nz)]);
            }
        }
    }
    *CFL = dt * (max_u / dx + max_v / dy + max_w / dz);
    *Y_min = Y_min_tmp;
    *Y_max = Y_max_tmp;
    *T_min = T_min_tmp;
    *T_max = T_max_tmp;
}

void initial_conditions(double *u, double *v, double *w, double *T, double *Y, double *p, Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int Nz_Y_max = parameters.Nz_Y_max;
    int *Nz_Y = parameters.Nz_Y;
    /* Velocity parameters */
    double u_r = parameters.u_r;
    double z_r = parameters.z_r;
    double alpha_u = parameters.alpha_u;
    /* Temperature parameters */
    double x_0 = parameters.T0_x_center;
    double y_0 = parameters.T0_y_center;
    double z_0 = parameters.T0_z_center;
    double sx = parameters.T0_length;
    double sy = parameters.T0_width;
    double sz = parameters.T0_height;
    double T_hot = parameters.T_hot;
    double T_inf = parameters.T_inf;
    /* Spatial domain */
    double *x = parameters.x;
    double *y = parameters.y;
    double *z = parameters.z;
    /* Pressure paramenters */
    double p_top = parameters.p_top;
    /* Fill arrays */
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                u[IDX(i, j, k, Nx, Ny, Nz)] = power_law(z[k], u_r, z_r, alpha_u);
                v[IDX(i, j, k, Nx, Ny, Nz)] = 0.0;
                w[IDX(i, j, k, Nx, Ny, Nz)] = 0.0;
                T[IDX(i, j, k, Nx, Ny, Nz)] = T_inf +  (T_hot - T_inf) * gaussian(x[i], y[j], z[k], x_0, y_0, z_0, sx, sy, sz);
                if (k == Nz - 1) {
                    p[IDX(i, j, k, Nx, Ny, Nz)] = p_top;
                } else {
                    p[IDX(i, j, k, Nx, Ny, Nz)] = 0.0;
                }
                if (k < Nz_Y[IDX(i, j, 0, Nx, Ny, 1)]) {
                    Y[IDX(i, j, k, Nx, Ny, Nz_Y_max)] = 1.0;
                }
            }
        }
    }
}