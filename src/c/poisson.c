/**
 * @file poisson.c
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Functions for solving the Poisson equation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "../../include/c/poisson.h"

void thomas_algorithm(fftw_complex *a, fftw_complex *b, fftw_complex *c, fftw_complex *d, fftw_complex *x, fftw_complex *l, fftw_complex *u, fftw_complex *y, int N) {
    // Create L and U
    u[0] = b[0];
    for (int i = 1; i < N; i++) {
        l[i-1] = a[i-1] / u[i - 1];
        u[i] = b[i] - l[i-1] * c[i - 1];
    }
    // Solve Ly = d
    y[0] = d[0];
    for (int i = 1; i < N; i++) {
        y[i] = d[i] - l[i-1] * y[i - 1];
    }
    // Solve Ux = y
    x[N - 1] = y[N - 1] / u[N - 1];
    for (int i = N - 2; i >= 0; i--) {
        x[i] = (y[i] - c[i] * x[i+1]) / u[i];
    }
}

void solve_pressure(double *U, double *p, fftw_complex *a, fftw_complex *b, fftw_complex *c, fftw_complex *d, fftw_complex *l, fftw_complex *u, fftw_complex *y, fftw_complex *pk, fftw_plan p_plan, fftw_plan f_plan, fftw_plan p_top_plan, fftw_complex *f_in, fftw_complex *f_out, fftw_complex *p_top_in, fftw_complex *p_top_out, fftw_complex *p_in, fftw_complex *p_out, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int u_index = parameters->field_indexes.u;
    int v_index = parameters->field_indexes.v;
    int w_index = parameters->field_indexes.w;
    double dx = parameters->dx;
    double dy = parameters->dy;
    double dz = parameters->dz;
    double dt = parameters->dt;
    double rho = parameters->rho;
    double u_ijk, u_ip1jk, u_im1jk, u_iphjk, u_imhjk;
    double v_ijk, v_ijp1k, v_ijm1k, v_ijphk, v_ijmhk;
    double w_ijk, w_ijkp1, w_ijkm1, w_ijkp2, w_ijkm2, w_ijkph, w_ijkmh;
    double ux, vy, wz, f;
    double gamma_rs;
    double *kx = parameters->kx;
    double *ky = parameters->ky;
    int n[] = {Nx - 1, Ny - 1};
    int idist = (Nx - 1) * (Ny - 1);
    int odist = (Nx - 1) * (Ny - 1);
    int istride = 1;
    int ostride = 1;
    int howmany = Nz - 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int im1, ip1, jm1, jp1;
    // Loop over nodes to compute f
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                // Indexes for periodic boundary conditions
                im1 = (i - 1 + Nx - 1) % (Nx - 1);
                jm1 = (j - 1 + Ny - 1) % (Ny - 1);
                ip1 = (i + 1) % (Nx - 1);
                jp1 = (j + 1) % (Ny - 1);
                // Local nodes
                u_ijk = U[u_index + IDX(i, j, k, Nx, Ny, Nz)];
                v_ijk = U[v_index + IDX(i, j, k, Nx, Ny, Nz)];
                w_ijk = U[w_index + IDX(i, j, k, Nx, Ny, Nz)];
                // Periodic boundary conditions on xy
                u_im1jk = U[u_index + IDX(im1, j, k, Nx, Ny, Nz)];
                u_ip1jk = U[u_index + IDX(ip1, j, k, Nx, Ny, Nz)];
                v_ijm1k = U[v_index + IDX(i, jm1, k, Nx, Ny, Nz)];
                v_ijp1k = U[v_index + IDX(i, jp1, k, Nx, Ny, Nz)];
                // dw/dz 
                if (k == 0) { // Bottom boundary                    
                    w_ijkp1 = U[w_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
                    w_ijkp2 = U[w_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
                    wz = (-3 * w_ijk + 4 * w_ijkp1 - w_ijkp2) / (2 * dz); // dw/dz at z = z_min
                } else if (k == Nz - 1) { // Top boundary
                    w_ijkm1 = U[w_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
                    w_ijkm2 = U[w_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
                    wz = (3 * w_ijk - 4 * w_ijkm1 + w_ijkm2) / (2 * dz); // dw/dz at z = z_max
                } else { // Interior
                    w_ijkp1 = U[w_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
                    w_ijkm1 = U[w_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
                    w_ijkph = 0.5 * (w_ijk + w_ijkp1);
                    w_ijkmh = 0.5 * (w_ijk + w_ijkm1);
                    wz = (w_ijkph - w_ijkmh) / dz; // dw/dz at z = z_k
                }
                // Half derivatives
                u_iphjk = 0.5 * (u_ip1jk + u_ijk);
                u_imhjk = 0.5 * (u_ijk + u_im1jk);
                v_ijphk = 0.5 * (v_ijp1k + v_ijk);
                v_ijmhk = 0.5 * (v_ijk + v_ijm1k);
                ux = (u_iphjk - u_imhjk) / dx; // du/dx
                vy = (v_ijphk - v_ijmhk) / dy; // dv/dy
                if (i < Nx - 1 && j < Ny - 1 && k < Nz - 1) {
                    // Compute rho / dt * div(U) and store it for many DFT (contiguous z slices)
                    f = rho * (ux + vy + wz) / dt;
                    f_in[FFTWIDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)] = f;
                }
                // Fill a and c
                if (i == 0 && j == 0 && k < Nz - 2) {
                    a[k] = 1.0 / (dz * dz) + 0.0 * I;
                    c[k] = 1.0 / (dz * dz) + 0.0 * I;
                }
                // Fill p_top
                if (k == Nz - 1 && j < Ny - 1 && i < Nx - 1) {
                    p_top_in[FFTWIDX(i, j, 0, Nx - 1, Ny - 1, Nz - 1)] = p[IDX(i, j, k, Nx, Ny, Nz)];
                }
            }
        }
    }
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
            // First equation k=0
            // f_out[s + (Ny - 1) * r] = 0.5 * dz * f_out[s + (Ny - 1) * r + (Nx - 1) * (Ny - 1) * 1];
            f_out[FFTWIDX(r, s, 0, Nx - 1, Ny - 1, Nz - 1)] =  0.5 * dz * f_out[FFTWIDX(r, s, 1, Nx - 1, Ny - 1, Nz - 1)];
            // Last equation k=Nz-2
            // f_out[s + (Ny - 1) * r + (Nx - 1) * (Ny - 1) * (Nz - 2)] -= p_top_out[s + (Ny - 1) * r] / (dz * dz);
            f_out[FFTWIDX(r, s, Nz - 2, Nx - 1, Ny - 1, Nz - 1)] -= p_top_out[FFTWIDX(r, s, 0, Nx - 1, Ny - 1, Nz - 1)] / (dz * dz);
            // Fill diagonal elements of A and RHS, and temporal arrays
            for (int k = 0; k < Nz - 1; k++) {
                b[k] = gamma_rs / (dz * dz) + 0.0 * I;
                d[k] = f_out[FFTWIDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)];
            }
            // Fix first coefficient of b and c
            b[0] =  -1.0 / dz + 0.0 * I;
            c[0] = (2.0 + 0.5 * gamma_rs) / dz + 0.0 * I;
            // Solve tridiagonal system
            thomas_algorithm(a, b, c, d, pk, l, u, y, Nz - 1);
            // Fill p_in with solution 
            for (int k = 0; k < Nz - 1; k++) {
                p_in[FFTWIDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = pk[k];                
            }
        }
    }
    // Plan for IFFT2(p) for each z slice
    p_plan = fftw_plan_many_dft(2, n, howmany, p_in, inembed, istride, idist, p_out, onembed, ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);
    // Compute IFFT
    fftw_execute(p_plan);
    // Restore data shape, compute real part and fill boundary conditions
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                if (i < Nx - 1 && j < Ny - 1 && k < Nz - 1) {
                    // p[IDX(i, j, k, Nx, Ny, Nz)] = creal(p_out[j + (Ny - 1) * i + (Nx - 1) * (Ny - 1) * k] / ((Nx - 1) * (Ny - 1)));
                    p[IDX(i, j, k, Nx, Ny, Nz)] = creal(p_out[FFTWIDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)] / ((Nx - 1) * (Ny - 1)));   
                }
                // Periodic boundary conditions on xy
                if (i == Nx - 1) { // Right boundary
                    p[IDX(i, j, k, Nx, Ny, Nz)] = p[IDX(0, j, k, Nx, Ny, Nz)];
                }
                if (j == Ny - 1) { // Front boundary
                    p[IDX(i, j, k, Nx, Ny, Nz)] = p[IDX(i, 0, k, Nx, Ny, Nz)];
                }
            }
        }
    }
    // Destroy plans
    fftw_destroy_plan(p_top_plan);
    fftw_destroy_plan(f_plan);
    fftw_destroy_plan(p_plan);
}