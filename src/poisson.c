#include "../include/poisson.h"

/**
 * @file poisson.c
 * @brief Implements the Thomas algorithm for solving a tridiagonal linear system.
 *
 * This file contains the implementation of the Thomas algorithm, also known as the tridiagonal matrix algorithm,
 * for solving a tridiagonal linear system of equations. The algorithm consists of three steps: decomposition,
 * forward substitution, and backward substitution.
 *
 * @param a The lower diagonal elements of the tridiagonal matrix.
 * @param b The main diagonal elements of the tridiagonal matrix.
 * @param c The upper diagonal elements of the tridiagonal matrix.
 * @param d The right-hand side vector.
 * @param x The solution vector.
 * @param l Temporary array to store the values of the lower diagonal elements during decomposition.
 * @param u Temporary array to store the values of the main diagonal elements during decomposition.
 * @param y Temporary array to store the values of the intermediate solution during forward substitution.
 * @param N The size of the linear system.
 */
void thomas_algorithm(double complex *a, double complex *b, double complex *c, double complex *d, double complex *x, double complex *l, double complex *u, double complex *y, int N) {
    // Create L and U
    l[0] = 0;
    u[0] = b[0];
    for (int i = 1; i < N; i++) {
        l[i] = a[i] / u[i - 1];
        u[i] = b[i] - l[i] * c[i - 1];
    }

    // Solve Ly = d
    y[0] = d[0];
    for (int i = 1; i < N; i++) {
        y[i] = d[i] - l[i] * y[i - 1];
    }

    // Solve Ux = y
    x[N - 1] = y[N - 1] / u[N - 1];
    for (int i = N - 2; i >= 0; i--) {
        x[i] = (y[i] - c[i] * x[i + 1]) / u[i];
    }
}

void solve_pressure(double *U, double *P_RHS, Parameters *parameters) {
    /* Compute div(U) */
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int size = Nx * Ny * Nz;
    int u_index = 0;
    int v_index = size;
    int w_index = 2 * size;
    double dx = parameters->dx;
    double dy = parameters->dy;
    double dz = parameters->dz;
    double dt = parameters->dt;
    double rho = parameters->rho;
    double u_ijk, u_ip1jk, u_im1jk, v_ijk, v_ijp1k, v_ijm1k, w_ijk, w_ijkp1, w_ijkm1;
    double ux, vy, wz;
    /* Compute rho / dt * div(U)*/
    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int k = 1; k < Nz - 1; k++) {
                u_ijk   = U[u_index + IDX(i, j, k, Nx, Ny, Nz)];
                u_ip1jk = U[u_index + IDX(i + 1, j, k, Nx, Ny, Nz)];
                u_im1jk = U[u_index + IDX(i - 1, j, k, Nx, Ny, Nz)];
                v_ijk   = U[v_index + IDX(i, j, k, Nx, Ny, Nz)];
                v_ijp1k = U[v_index + IDX(i, j + 1, k, Nx, Ny, Nz)];
                v_ijm1k = U[v_index + IDX(i, j - 1, k, Nx, Ny, Nz)];
                w_ijk   = U[w_index + IDX(i, j, k, Nx, Ny, Nz)];
                w_ijkp1 = U[w_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
                w_ijkm1 = U[w_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
                ux = (u_ip1jk - u_im1jk) / (2 * dx);
                vy = (v_ijp1k - v_ijm1k) / (2 * dy);
                wz = (w_ijkp1 - w_ijkm1) / (2 * dz);
                P_RHS[IDX(i, j, k, Nx, Ny, Nz)] = rho * (ux + vy + wz) / dt;
            }
        }
    }
}