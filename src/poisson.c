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
void thomas_algorithm(double complex *a, double complex *b, double complex *c, double complex *d, fftw_complex *x, double complex *l, double complex *u, double complex *y, int N) {
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

void solve_pressure(double *U, double *p, Parameters *parameters) {
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
    double gamma_rs;
    double *kx = parameters->kx;
    double *ky = parameters->ky;
    double complex *a = malloc((Nz - 2) * sizeof(double complex));
    double complex *b = malloc((Nz - 1) * sizeof(double complex));
    double complex *c = malloc((Nz - 2) * sizeof(double complex));
    double complex *d = malloc((Nz - 1) * sizeof(double complex));
    double complex *l = malloc((Nz - 1) * sizeof(double complex));
    double complex *u = malloc((Nz - 1) * sizeof(double complex));
    double complex *y = malloc((Nz - 1) * sizeof(double complex));
    double complex *pk = malloc((Nz - 1) * sizeof(double complex));
    // double *P_k = malloc((Nx - 1) * (Ny - 1) * (Nz - 1) * sizeof(double));
    fftw_complex *f_in = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    fftw_complex *f_out = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    fftw_complex *p_top_in = fftw_alloc_complex((Nx - 1) * (Ny - 1));
    fftw_complex *p_top_out = fftw_alloc_complex((Nx - 1) * (Ny - 1));
    fftw_complex *p_in = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    fftw_complex *p_out = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    int n[] = {Nx - 1, Ny - 1};
    int idist = (Nx - 1) * (Ny - 1);
    int odist = (Nx - 1) * (Ny - 1);
    int istride = 1;
    int ostride = 1;
    int howmany = Nz - 1;
    int *inembed = NULL;
    int *onembed = NULL;
    fftw_plan p_plan, f_plan, p_top_plan;
    
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
                // RHS of Poisson problem
                // f[IDX(i, j, k, Nx, Ny, Nz)] = rho * (ux + vy + wz) / dt;                
                // Fill p with zeros
                if (i < Nx - 1 && j < Ny - 1 && k < Nz - 1) {
                    // Compute rho / dt * div(U) and store it for many DFT (contiguous z slices)
                    f_in[j + Ny * i + Nx * Ny * k] = rho * (ux + vy + wz) / dt; //f[IDX(i, j, k, Nx, Ny, Nz)] + I * 0.0;
                    p[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)] = 0.0;
                }
                // Fill a and c
                if (i == 0 && j == 0 && k < Nz - 2) {
                    a[k] = 1 / (dz * dz);
                    c[k] = 1 / (dz * dz);
                }
                // Fill p_top
                if (k == Nz - 1) {
                    p_top_in[j + Ny * i] = p[IDX(i, j, k, Nx, Ny, Nz)];
                }
            }
        }
    }

    // Prepare data for FFT
    // Plan for FFT2(f) for each z slice
    f_plan = fftw_plan_many_dft(2, n, howmany, f_in, inembed, istride, idist, f_out, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
    // Plan for FFT2(p_top)
    p_top_plan = fftw_plan_dft_2d(Nx - 1, Ny - 1, p_top_in, p_top_out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Compute FFT
    fftw_execute(p_top_plan); // FFT2(p_top)
    fftw_execute(f_plan); // FFT2(f) for each z slice
    
    // Compute r,s systems of equations
    for (int r = 0; r < Nx - 1; r++) {
        for (int s = 0; s < Ny - 1; s++) {
            gamma_rs = - 2 - kx[r] * kx[r] - ky[s] * ky[s];
            // First equation
            // f[IDX(r, s, 0, Nx, Ny, Nz)] = 0.0 + 0.5 * dz * f[IDX(r, s, 1, Nx, Ny, Nz)];
            f_out[s + (Ny - 1) * r] = 0.0 + 0.5 * dz * f_out[s + (Ny - 1) * r + (Nx - 1) * (Ny - 1) * 1];
            // Last equation
            // f[IDX(r, s, Nz - 2, Nx, Ny, Nz)] -= p_top_out[IDX(r, s, 0, Nx, Ny, 1)] / (dz * dz);
            f_out[s + (Ny - 1) * r + (Nx - 1) * (Ny - 1) * (Nz - 2)] -= p_top_out[s + (Ny - 1) * r] / (dz * dz);
            // Fill diagonal elements of A and RHS, and temporal arrays
            for (int k = 0; k < Nz - 1; k++) {
                b[k] = gamma_rs / (dz * dz);
                d[k] = f_out[s + (Ny - 1) * r + (Nx - 1) * (Ny - 1) * k];
                l[k] = 0.0;
                u[k] = 0.0;
                y[k] = 0.0;
                pk[k] = 0.0; // To store the solution
            }
            // Fix first coefficient of b and c
            b[0] =  -1 / dz;
            c[0] = (2 + 0.5 * gamma_rs) /dz;
            // Solve tridiagonal system
            // thomas_algorithm(a, b, c, d, p_in + IDX(r, s, 0, Nx - 1, Ny - 1, Nz - 1), p_out + IDX(r, s, 0, Nx, Ny, Nz), a, b, c, Nz - 1);
            thomas_algorithm(a, b, c, d, pk, l, u, y, Nz - 1);
            // thomas_algorithm(a, b, c, f + IDX(r, s, 0, Nx, Ny, Nz), pk_in + IDX(r, s, 0, Nx, Ny, Nz), a, b, c, Nz - 1);

            // Fill p_in with solution
            for (int k = 0; k < Nz - 1; k++) {
                p_in[s + (Ny - 1) * r + (Nx - 1) * (Ny - 1) * k] = pk[k];
            }
        }
    }

    /*
    // Plan for IFFT2(p) for each z slice
    p_plan = fftw_plan_many_dft(2, n, howmany, p_in, inembed, istride, idist, p_out, onembed, ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    // Compute IFFT
    fftw_execute(p_plan);

    // Restore data shape, compute real part and fill boundary conditions
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                if (i < Nx - 1 && j < Ny - 1 && k < Nz - 1) {
                    p[IDX(i, j, k, Nx, Ny, Nz)] = creal(p_out[j + Ny * i + Nx * Ny * k]);
                }
                // Periodic boundary conditions on xy
                if (i == Nx - 1) { // Right boundary
                    p[IDX(i, j, k, Nx, Ny, Nz)] = p[IDX(0, j, k, Nx, Ny, Nz)];
                }
                if (j == Ny - 1) { // Front boundary
                    p[IDX(i, j, k, Nx, Ny, Nz)] = p[IDX(i, 0, k, Nx, Ny, Nz)];
                }
                // p_top on z = z_max
                if (k == Nz - 1) {
                    p[IDX(i, j, k, Nx, Ny, Nz)] = p[IDX(i, j, k, Nx, Ny, Nz)];
                }
            }
        }
    }
    */

    // Destroy plans
    fftw_destroy_plan(p_top_plan);
    fftw_destroy_plan(f_plan);
    fftw_destroy_plan(p_plan);
    // Free memory
    // fftw_free(p_top_in);
    // fftw_free(p_in);
    // fftw_free(f_in);
    // fftw_free(p_top_out);
    // fftw_free(f_out);
    // fftw_free(p_out);
    free(a);
    free(b);
    free(c);
    free(d);
    free(l);
    free(u);
    free(y);
}