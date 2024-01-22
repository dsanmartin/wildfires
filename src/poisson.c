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

/**
 * @brief Prepare data for Fast Fourier Transform (FFT).
 *
 * This function prepares the data for FFT by creating contiguous z slices.
 * It takes an array of complex numbers A, an array of double numbers B, and a Parameters struct.
 * The Parameters struct contains the dimensions of the data (Nx, Ny, Nz).
 *
 * @param A The array of complex numbers to store the prepared data.
 * @param B The array of double numbers containing the original data.
 * @param parameters The Parameters struct containing the dimensions of the data.
 */
void prepare_data_fft(fftw_complex *A, double *B, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    // Create contiguous_z_slices
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                A[j + Ny * i + Nx * Ny * k] = B[IDX(i, j, k, Nx, Ny, Nz)] + I * 0.0;
            }
        }
    }
}

/**
 * @brief Restores data from complex array to real array using FFT.
 *
 * This function restores data from a complex array to a real array using the Fast Fourier Transform (FFT).
 * The restored data is stored in the array A.
 *
 * @param B The complex array containing the data to be restored.
 * @param A The real array where the restored data will be stored.
 * @param parameters The parameters struct containing the dimensions of the arrays.
 */
void restore_data_fft(fftw_complex *B, double *A, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    // Restore data
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                A[IDX(i, j, k, Nx, Ny, Nz)] = creal(B[j + Ny * i + Nx * Ny * k]);
            }
        }
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
    double u_ip1jk, u_im1jk, v_ijp1k, v_ijm1k, w_ijk, w_ijkp1, w_ijkm1, w_ijkp2, w_ijkm2;
    double ux, vy, wz;
    fftw_complex *in = fftw_alloc_complex(size);
    fftw_complex *out = fftw_alloc_complex(size);
    int n[] = {Nx, Ny};
    int idist = Nx * Ny;
    int odist = Nx * Ny;
    int istride = 1;
    int ostride = 1;
    int howmany = Nz;
    int *inembed = NULL;
    int *onembed = NULL;
    fftw_plan p = fftw_plan_many_dft(
        2, n, howmany, 
        in, inembed, istride, idist,
        out, onembed, ostride, odist, 
        FFTW_FORWARD, FFTW_ESTIMATE);
    /* Compute rho / dt * div(U)*/
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                // Periodic boundary conditions on xy
                if (i == 0) { // Left boundary
                    u_im1jk = U[u_index + IDX(Nx - 2, j, k, Nx, Ny, Nz)];
                } else { // Interior
                    u_im1jk = U[u_index + IDX(i - 1, j, k, Nx, Ny, Nz)];
                }
                if (i == Nx - 1) { // Right boundary
                    u_ip1jk = U[u_index + IDX(1, j, k, Nx, Ny, Nz)];
                } else { // Interior
                    u_ip1jk = U[u_index + IDX(i + 1, j, k, Nx, Ny, Nz)];
                }
                if (j == 0) { // Back boundary
                    v_ijm1k = U[v_index + IDX(i, Ny - 2, k, Nx, Ny, Nz)];
                } else { // Interior
                    v_ijm1k = U[v_index + IDX(i, j - 1, k, Nx, Ny, Nz)];
                }
                if (j == Ny - 1) { // Front boundary
                    v_ijp1k = U[v_index + IDX(i, 1, k, Nx, Ny, Nz)];
                } else {
                    v_ijp1k = U[v_index + IDX(i, j + 1, k, Nx, Ny, Nz)];
                }
                // dw/dz 
                if (k == 0) { // Bottom boundary
                    w_ijk   = U[w_index + IDX(i, j, k, Nx, Ny, Nz)];
                    w_ijkp1 = U[w_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
                    w_ijkp2 = U[w_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
                    wz = (-3 * w_ijk + 4 * w_ijkp1 - w_ijkp2) / (2 * dz); // dw/dz at z = z_min
                } else if (k == Nz - 1) { // Top boundary
                    w_ijk   = U[w_index + IDX(i, j, k, Nx, Ny, Nz)];
                    w_ijkm1 = U[w_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
                    w_ijkm2 = U[w_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
                    wz = (3 * w_ijk - 4 * w_ijkm1 + w_ijkm2) / (2 * dz); // dw/dz at z = z_max
                } else { // Interior
                    w_ijkm1 = U[w_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
                    w_ijkp1 = U[w_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
                    wz = (w_ijkp1 - w_ijkm1) / (2 * dz); // dw/dz at z = z_k
                }
                ux = (u_ip1jk - u_im1jk) / (2 * dx); // du/dx
                vy = (v_ijp1k - v_ijm1k) / (2 * dy); // dv/dy
                P_RHS[IDX(i, j, k, Nx, Ny, Nz)] = rho * (ux + vy + wz) / dt;
            }
        }
    }
    // Prepare data for FFT
    prepare_data_fft(in, P_RHS, parameters);
    // Compute FFT
    fftw_execute(p);
    // Restore data
    restore_data_fft(out, P_RHS, parameters);
    // Destroy plan
    fftw_destroy_plan(p);
    // Free memory
    fftw_free(in);
    fftw_free(out);
}