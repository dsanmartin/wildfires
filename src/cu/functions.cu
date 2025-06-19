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

double parallelepiped(double x, double y, double z, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max) {
    if (x >= x_min && x <= x_max && y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
        return 1.0;
    } else {
        return 0.0;
    }
}

double K(double T, double A, double T_act) {
    return A * exp(-T_act / T);
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

double source(double T, double Y, double H_R, double A, double T_a, double h, double a_v, double T_inf, double c_p, double rho_inf, double T_pc) {
    return H_R * Y * K(T, A, T_a) * H(T, T_pc) / c_p - h * a_v * (T - T_inf) / (c_p * rho_inf);
}

double non_uniform_z(double z, double z_min, double z_max, double k) {
    return (z_max - z_min) * (exp(k * (z - z_min) / (z_max - z_min)) - 1) / (exp(k) - 1) + z_min;
}

void equispaced_domain(double *z, int Nz, double z_min, double dz) {
    for (int k = 0; k < Nz; k++) 
        z[k] = z_min + k * dz;
}

void non_equispaced_domain(double *z, int Nz, double z_min, double z_max, double nu) {
    double dz = (z_max - z_min) / (Nz - 1);
    for (int k = 0; k < Nz; k++) {
        z[k] = (z_max - z_min) * (exp(nu * k * dz / (z_max - z_min)) - 1) / (exp(nu) - 1) + z_min;
    }
}

void transition_domain(double *z, int Nz, double z_min, double z_max, double z_t, int k_t) {
    double dz_1 = (z_t - z_min) / (k_t - 1);
    double dz_2 = (z_max - z_t) / (Nz - k_t - 1);
    for (int k = 0; k < Nz; k++) {
        if (k <= k_t) {
            z[k] = z_min + k * dz_1;
        } else {
            z[k] = z_t + (k - k_t) * dz_2;
        }
    }
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
    double u, v, w, T, Y;
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
                u = y_n[u_index + IDX(i, j, k, Nx, Ny, Nz)];
                v = y_n[v_index + IDX(i, j, k, Nx, Ny, Nz)];
                w = y_n[w_index + IDX(i, j, k, Nx, Ny, Nz)];
                T = y_n[T_index + IDX(i, j, k, Nx, Ny, Nz)];
                Y = 0;
                if (k < Nz_Y_max) {
                    Y = y_n[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y_max)];
                    Y_min_tmp = MIN(Y_min_tmp, Y);
                    Y_max_tmp = MAX(Y_max_tmp, Y);
                }
                // Check if any value is NaN
                if (isnan(u) || isnan(v) || isnan(w) || isnan(T) || isnan(Y)) {
                    log_message(parameters, "NaN value found. Exiting...");
                    exit(1);
                }
                abs_u = fabs(u);
                abs_v = fabs(v);
                abs_w = fabs(w);
                max_u = MAX(max_u, abs_u);
                max_v = MAX(max_v, abs_v);
                max_w = MAX(max_w, abs_w);
                T_min_tmp = MIN(T_min_tmp, T);
                T_max_tmp = MAX(T_max_tmp, T);
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
    double T_source = parameters.T_source;
    double T_inf = parameters.T_inf;
    double T0_x_start = parameters.T0_x_start;
    double T0_x_end = parameters.T0_x_end;
    double T0_y_start = parameters.T0_y_start;
    double T0_y_end = parameters.T0_y_end;
    double T0_z_start = parameters.T0_z_start;
    double T0_z_end = parameters.T0_z_end;
    double Y0_x_start = parameters.Y0_x_start;
    double Y0_x_end = parameters.Y0_x_end;
    double Y0_y_start = parameters.Y0_y_start;
    double Y0_y_end = parameters.Y0_y_end;
    // int fuel_relax = parameters.fuel_relax;
    double xa = parameters.Y0_xa;
    double ya = parameters.Y0_ya;
    double Y_h = parameters.Y_h;
    double z_bottom, z_top, z_left, z_right;
    /* Spatial domain */
    double *x = parameters.x;
    double *y = parameters.y;
    double *z = parameters.z;
    // double dx = parameters.dx;
    // double dy = parameters.dy;
    /* Pressure paramenters */
    double p_top = parameters.p_top;
    /* IBM parameters */
    int *cut_nodes = parameters.cut_nodes;
    double u_dead_nodes = parameters.u_dead_nodes;
    double v_dead_nodes = parameters.v_dead_nodes;
    double w_dead_nodes = parameters.w_dead_nodes;
    double T_dead_nodes = parameters.T_dead_nodes;
    double Y_dead_nodes = parameters.Y_dead_nodes;
    /* Fill arrays */
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                u[IDX(i, j, k, Nx, Ny, Nz)] = power_law(z[k], u_r, z_r, alpha_u);
                v[IDX(i, j, k, Nx, Ny, Nz)] = 0.0;
                w[IDX(i, j, k, Nx, Ny, Nz)] = 0.0;
                if (strcmp(parameters.T0_shape, "gaussian") == 0) {
                    T[IDX(i, j, k, Nx, Ny, Nz)] = T_inf +  (T_source - T_inf) * gaussian(x[i], y[j], z[k], x_0, y_0, z_0, sx, sy, sz);
                } else if (strcmp(parameters.T0_shape, "parallelepiped") == 0) {
                    T[IDX(i, j, k, Nx, Ny, Nz)] = T_inf +  (T_source - T_inf) * parallelepiped(x[i], y[j], z[k], T0_x_start, T0_x_end, T0_y_start, T0_y_end, T0_z_start, T0_z_end);
                }
                // T[IDX(i, j, k, Nx, Ny, Nz)] = T_inf +  (T_source - T_inf) * gaussian(x[i], y[j], z[k], x_0, y_0, z_0, sx, sy, sz);
                if (k == Nz - 1) {
                    p[IDX(i, j, k, Nx, Ny, Nz)] = p_top;
                } else {
                    p[IDX(i, j, k, Nx, Ny, Nz)] = 0.0;
                }
                // if (k < Nz_Y[IDX(i, j, 0, Nx, Ny, 1)] && x[i] > -1.0 && x[i] < 201 && y[j] > -1.0 && y[j] < 201) {
                //     Y[IDX(i, j, k, Nx, Ny, Nz_Y_max)] = 1.0;
                // }
                if (k < Nz_Y[IDX(i, j, 0, Nx, Ny, 1)]) {
                    z_bottom = Y_h * (y[j] - Y0_y_start + ya) / ya;
                    z_top = Y_h * (Y0_y_end + ya - y[j]) / ya;
                    z_left = Y_h * (x[i] - Y0_x_start + xa) / xa;
                    z_right = Y_h * (Y0_x_end + xa - x[i]) / xa;
                    if ((Y0_x_start <= x[i] && x[i] <= Y0_x_end) && (Y0_y_start <= y[j] && y[j] <= Y0_y_end) && (0 <= z[k] && z[k] <= Y_h))
                        Y[IDX(i, j, k, Nx, Ny, Nz_Y_max)] = 1.0;
                    else if ((Y0_y_start - ya <= y[j] && y[j] < Y0_y_start) && (Y0_x_start - xa <= x[i] && x[i] <= Y0_x_end + xa) && (0 <= z[k] && z[k] <= z_bottom))
                        Y[IDX(i, j, k, Nx, Ny, Nz_Y_max)] = 1.0; //z_bottom / Y_h;
                    else if ((Y0_y_end < y[j] && y[j] <= Y0_y_end + ya) && (Y0_x_start - xa <= x[i] && x[i] <= Y0_x_end + xa) && (0 <= z[k] && z[k] <= z_top))
                        Y[IDX(i, j, k, Nx, Ny, Nz_Y_max)] = 1.0; //z_top / Y_h;
                    else if ((Y0_x_start - xa <= x[i] && x[i] < Y0_x_start) && (Y0_y_start - ya <= y[j] && y[j] <= Y0_y_end + ya) && (0 <= z[k] && z[k] <= z_left))
                        Y[IDX(i, j, k, Nx, Ny, Nz_Y_max)] = 1.0; //z_left / Y_h;
                    else if ((Y0_x_end < x[i] && x[i] <= Y0_x_end + xa) && (Y0_y_start - ya <= y[j] && y[j] <= Y0_y_end + ya) && (0 <= z[k] && z[k] <= z_right))
                        Y[IDX(i, j, k, Nx, Ny, Nz_Y_max)] = 1.0; //z_right / Y_h;
                }
                // if (k < Nz_Y[IDX(i, j, 0, Nx, Ny, 1)] && x[i] >= Y0_x_start && x[i] <= Y0_x_end && y[j] >= Y0_y_start && y[j] <= Y0_y_end) {
                //     Y[IDX(i, j, k, Nx, Ny, Nz_Y_max)] = 1.0;
                // }
                // if (k < Nz_Y[IDX(i, j, 0, Nx, Ny, 1)] - 1 && 
                //     x[i] >= Y0_x_start - fuel_relax * dx && x[i] < Y0_x_start && x[i] > Y0_x_end && x[i] < Y0_x_end + fuel_relax * dx &&
                //     y[j] >= Y0_y_start - fuel_relax * dy && y[j] < Y0_y_start && y[j] > Y0_y_end && y[j] < Y0_y_end+ fuel_relax * dy) {
                //         Y[IDX(i, j, k, Nx, Ny, Nz_Y_max)] = 0.75;
                // }
                // if (k < Nz_Y[IDX(i, j, 0, Nx, Ny, 1)] - 2 && 
                //     x[i] >= Y0_x_start - 2 * fuel_relax * dx && x[i] < Y0_x_start - fuel_relax * dx && x[i] > Y0_x_end + fuel_relax * dx && x[i] < Y0_x_end + 2 * fuel_relax * dx &&
                //     y[j] >= Y0_y_start - 2 * fuel_relax * dy && y[j] < Y0_y_start - fuel_relax * dy && y[j] > Y0_y_end + fuel_relax * dy && y[j] < Y0_y_end + 2 * fuel_relax * dy) {
                //         Y[IDX(i, j, k, Nx, Ny, Nz_Y_max)] = 0.5;
                // }
                // if (k < Nz_Y[IDX(i, j, 0, Nx, Ny, 1)] - 3 && 
                //     x[i] >= Y0_x_start - 3 * fuel_relax * dx && x[i] < Y0_x_start - 2 * fuel_relax * dx && x[i] > Y0_x_end + 2 * fuel_relax * dx && x[i] < Y0_x_end + 3 * fuel_relax * dx &&
                //     y[j] >= Y0_y_start - 3 * fuel_relax * dy && y[j] < Y0_y_start - 2 * fuel_relax * dy && y[j] > Y0_y_end + 2 * fuel_relax * dy && y[j] < Y0_y_end + 3 * fuel_relax * dy) {
                //         Y[IDX(i, j, k, Nx, Ny, Nz_Y_max)] = 0.25;
                // }
                // } else {
                //     if (k < Nz_Y[IDX(i, j, 0, Nx, Ny, 1)] - 1 && 
                //         x[i] >= 0.0 - fuel_relax * dx && x[i] < 0.0 && x[i] > 200 && x[i] < 200 + fuel_relax * dx &&
                //         y[j] >= 0.0 - fuel_relax * dy && y[j] < 0.0 && y[j] > 200 && y[j] < 200 + fuel_relax * dy) {
                //             Y[IDX(i, j, k, Nx, Ny, Nz_Y_max)] = 0.75;
                //     } else {
                //         if (k < Nz_Y[IDX(i, j, 0, Nx, Ny, 1)] - 2 && 
                //             x[i] >= 0.0 - 2 * fuel_relax * dx && x[i] < - fuel_relax * dx && x[i] > 200 + fuel_relax * dx && x[i] < 200 + 2 * fuel_relax * dx &&
                //             y[j] >= 0.0 - 2 * fuel_relax * dy && y[j] < - fuel_relax * dy && y[j] > 200 + fuel_relax * dy && y[j] < 200 + 2 * fuel_relax * dy) {
                //                 Y[IDX(i, j, k, Nx, Ny, Nz_Y_max)] = 0.5;
                //         } else {
                //             if (k < Nz_Y[IDX(i, j, 0, Nx, Ny, 1)] - 3 && 
                //                 x[i] >= 0.0 - 3 * fuel_relax * dx && x[i] < -2 * fuel_relax * dx && x[i] > 200 + 2 * fuel_relax * dx && x[i] < 200 + 3 * fuel_relax * dx &&
                //                 y[j] >= 0.0 - 3 * fuel_relax * dy && y[j] < -2 * fuel_relax * dy && y[j] > 200 + 2 * fuel_relax * dy && y[j] < 200 + 3 * fuel_relax * dy) {
                //                     Y[IDX(i, j, k, Nx, Ny, Nz_Y_max)] = 0.25;
                //             }
                //         }
                //     }
                // }
                // Set dead nodes
                if (k < cut_nodes[IDX(i, j, 0, Nx, Ny, 1)]) {
                    u[IDX(i, j, k, Nx, Ny, Nz)] = u_dead_nodes;
                    v[IDX(i, j, k, Nx, Ny, Nz)] = v_dead_nodes;
                    w[IDX(i, j, k, Nx, Ny, Nz)] = w_dead_nodes;
                    T[IDX(i, j, k, Nx, Ny, Nz)] = T_dead_nodes;
                    if (k < Nz_Y[IDX(i, j, 0, Nx, Ny, 1)])
                        Y[IDX(i, j, k, Nx, Ny, Nz_Y_max)] = Y_dead_nodes;
                }
            }
        }
    }
}

__global__
void temperature_source(double *x, double *y, double *z, double *y_n, Parameters paramenters) {
    int Nx = paramenters.Nx;
    int Ny = paramenters.Ny;
    int Nz = paramenters.Nz;
    int T_index = paramenters.field_indexes.T;
    double T_source = paramenters.T_source;
    double T0_x_start = paramenters.T0_x_start;
    double T0_x_end = paramenters.T0_x_end;
    double T0_y_start = paramenters.T0_y_start;
    double T0_y_end = paramenters.T0_y_end;
    double T0_z_start = paramenters.T0_z_start;
    double T0_z_end = paramenters.T0_z_end;
    int size = Nx * Ny * Nz;
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;
    for (int ijk = idx; ijk < size; ijk += stride) {
        int i = ijk / (Ny * Nz);
        int j = (ijk % (Ny * Nz)) / Nz;
        int k = ijk % Nz;
        if (x[i] >= T0_x_start && x[i] <= T0_x_end && y[j] >= T0_y_start && y[j] <= T0_y_end && z[k] >= T0_z_start && z[k] <= T0_z_end) {
            y_n[T_index + IDX(i, j, k, Nx, Ny, Nz)] = T_source;
        }
    }
}

__global__
void norm(double *v1, double *v2, double *result, double p, int size) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;
    double res = 0;
    if (p == INFINITY) {
        double max = 0.0;
        for (int i = idx; i < size; i += stride) {
            max = MAX(max, fabs(v1[i] - v2[i]));
        }
        res = max;
    } else {
        double sum = 0.0;
        for (int i = idx; i < size; i += stride) {
            sum += pow(fabs(v1[i] - v2[i]), p);
        }
        res = pow(sum, 1.0 / p);
    }
    *result = res;
}