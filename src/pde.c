#include "../include/pde.h"

void Phi(float t, float *R_n, float *R_np1, Parameters *parameters) {
    int Nx = parameters->Nx - 1;
    int Ny = parameters->Ny - 1;
    int Nz = parameters->Nz;
    double dx = parameters->dx;
    double dy = parameters->dy;
    double dz = parameters->dz;
    double *x = parameters->x;
    double *y = parameters->y;
    double *z = parameters->z;
    double alpha = parameters->alpha;
    double Y_h = parameters->Y_h;
    double H_R = parameters->H_R;
    double A = parameters->A;
    double T_a = parameters->T_a;
    double h = parameters->h;
    double a_v = parameters->a_v;
    double T_inf = parameters->T_inf;
    double c_p = parameters->c_p;
    double rho = parameters->rho;
    int u_index = 0;
    int v_index = Nx * Ny * Nz;
    int w_index = 2 * Nx * Ny * Nz;
    int T_index = 3 * Nx * Ny * Nz;
    double u_ijk, u_ip1jk, u_im1jk, u_ijp1k, u_ijm1k, u_ijkp1, u_ijkm1;
    double v_ijk, v_ip1jk, v_im1jk, v_ijp1k, v_ijm1k, v_ijkp1, v_ijkm1;
    double w_ijk, w_ip1jk, w_im1jk, w_ijp1k, w_ijm1k, w_ijkp1, w_ijkm1;
    double T_ijk, T_ip1jk, T_im1jk, T_ijp1k, T_ijm1k, T_ijkp1, T_ijkm1;
    double Tx, Ty, Tz, Txx, Tyy, Tzz;
    double S = 0.0;    
    int ii, jj;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                ii = i % Nx;
                jj = j % Ny;
                T_ijk = R_n[T_index + idx(ii, jj, k, Nx, Ny, Nz)];
                T_ip1jk = R_n[T_index + idx(ii + 1, jj, k, Nx, Ny, Nz)];
                T_im1jk = R_n[T_index + idx(ii - 1, jj, k, Nx, Ny, Nz)];
                T_ijp1k = R_n[T_index + idx(ii, jj + 1, k, Nx, Ny, Nz)];
                T_ijm1k = R_n[T_index + idx(ii, jj - 1, k, Nx, Ny, Nz)];
                T_ijkp1 = R_n[T_index + idx(ii, jj, k + 1, Nx, Ny, Nz)];
                T_ijkm1 = R_n[T_index + idx(ii, jj, k - 1, Nx, Ny, Nz)];
                Tx = (T_ip1jk - T_im1jk) / (2 * dx);
                Ty = (T_ijp1k - T_ijm1k) / (2 * dy);
                Tz = (T_ijkp1 - T_ijkm1) / (2 * dz);
                Txx = (T_ip1jk - 2 * T_ijk + T_im1jk) / (dx * dx);
                Tyy = (T_ijp1k - 2 * T_ijk + T_ijm1k) / (dy * dy);
                Tzz = (T_ijkp1 - 2 * T_ijk + T_ijkm1) / (dz * dz);
                if (z[k] <= Y_h)
                    S = source(T_ijk, R_n[idx(i, j, k, Nx, Ny, Nz)], H_R, A, T_a, h, a_v, T_inf, c_p, rho);
                R_np1[T_index + idx(ii, jj, k, Nx, Ny, Nz)] = alpha * (Txx + Tyy + Tzz) - (u_ijk * Tx + v_ijk * Ty + w_ijk * Tz) + S;
            }
        }
    }
}

void euler(float t_n, float *y_n, float *y_np1, float *F, float dt, int size) {
    for (int i = 0; i < size; i++) {
        y_np1[i] = y_n[i] + dt * F[i];
    }
}

void rk4(float t_n, float *y_n, float *y_np1, float *k1, float *k2, float *k3, float *k4, float dt, int size) {
    for (int i = 0; i < size; i++) {
        y_np1[i] = y_n[i] + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
}