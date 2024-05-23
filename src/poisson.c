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
void thomas_algorithm(double *a, double *b, double *c, double complex *d, double complex *x, double complex *l, double complex *u, double complex *y, int N) {
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
        x[i] = (y[i] - c[i+1] * x[i+1]) / u[i];
    }
}
/*
void thomas_algorithm(double complex *a, double complex *b, double complex *c, double complex *d, fftw_complex *x, double complex *l, double complex *u, double complex *y, int N) {
    // Create L and U
    // printf("Thomas algorithm\n");
    u[0] = b[0];
    // printf("u[%d]: %f + %fi\n", 0, creal(u[0]), cimag(u[0]));
    for (int i = 1; i < N; i++) {
        l[i-1] = a[i-1] / u[i - 1];
        u[i] = b[i] - l[i-1] * c[i - 1];
        // printf("l[%d]: %f + %fi\n", i-1, creal(l[i-1]), cimag(l[i-1]));
        // printf("u[%d]: %f + %fi\n", i, creal(u[i]), cimag(u[i]));
    }
    // Solve Ly = d
    y[0] = d[0];
    // printf("y[%d]: %f + %fi\n", 0, creal(y[0]), cimag(y[0]));
    for (int i = 1; i < N; i++) {
        y[i] = d[i] - l[i-1] * y[i - 1];
        // printf("y[%d]: %f + %fi\n", i, creal(y[i]), cimag(y[i]));
    }
    // Solve Ux = y
    x[N - 1] = y[N - 1] / u[N - 1];
    // printf("x[%d]: %f + %fi\n", N-1, creal(x[N-1]), cimag(x[N-1]));
    // for (int i = N - 2; i >= 0; i--) {
    //     x[i] = (y[i] - c[i] * x[i]) / u[i-1];
    //     // printf("x[%d]: %f + %fi\n", i, creal(x[i]), cimag(x[i]));
    // }
    // for (int i = N - 1; i > 0; i--) {
    //     x[i] = (y[i] - c[i] * x[i+1]) / u[i-1];
    //     // printf("x[%d]: %f + %fi\n", i, creal(x[i]), cimag(x[i]));
    // }
    for (int i = N - 2; i >= 0; i--) {
        x[i] = (y[i] - c[i+1] * x[i+1]) / u[i];
        // printf("x[%d]: %f + %fi\n", i, creal(x[i]), cimag(x[i]));
    }
    // printf("Thomas algorithm... OK!\n");
}
*/

void solve_pressure(double *U, double *p, double *a, double *b, double *c, double complex *d, double complex *l, double complex *u, double complex *y, double complex *pk, fftw_plan p_plan, fftw_plan f_plan, fftw_plan p_top_plan, fftw_complex *f_in, fftw_complex *f_out, fftw_complex *p_top_in, fftw_complex *p_top_out, fftw_complex *p_in, fftw_complex *p_out, Parameters *parameters) {
// void solve_pressure(double *U, double *p, double complex *a, double complex *b, double complex *c, double complex *d, double complex *l, double complex *u, double complex *y, double complex *pk, Parameters *parameters) {
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
    double u_ip1jk = 0.0, u_im1jk = 0.0, u_ikj = 0.0, u_iphjk = 0.0, u_imhjk = 0.0;
    double v_ijp1k = 0.0, v_ijm1k = 0.0, v_ijk = 0.0, v_ijphk = 0.0, v_ijmhk = 0.0;
    double w_ijk = 0.0, w_ijkp1 = 0.0, w_ijkm1 = 0.0, w_ijkp2 = 0.0, w_ijkm2 = 0.0;
    double ux = 0.0, vy = 0.0, wz = 0.0;
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
    // fftw_plan p_plan, f_plan, p_top_plan;
    // fftw_complex *f_in = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    // fftw_complex *f_out = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    // fftw_complex *p_top_in = fftw_alloc_complex((Nx - 1) * (Ny - 1));
    // fftw_complex *p_top_out = fftw_alloc_complex((Nx - 1) * (Ny - 1));
    // fftw_complex *p_in = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    // fftw_complex *p_out = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    double f;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                // Indexes for periodic boundary conditions
                im1 = (i - 1 + Nx - 1) % (Nx - 1);
                jm1 = (j - 1 + Ny - 1) % (Ny - 1);
                ip1 = (i + 1) % (Nx - 1);
                jp1 = (j + 1) % (Ny - 1);
                // Periodic boundary conditions on xy
                u_im1jk = U[u_index + IDX(im1, j, k, Nx, Ny, Nz)];
                u_ip1jk = U[u_index + IDX(ip1, j, k, Nx, Ny, Nz)];
                v_ijm1k = U[v_index + IDX(i, jm1, k, Nx, Ny, Nz)];
                v_ijp1k = U[v_index + IDX(i, jp1, k, Nx, Ny, Nz)];
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
                // Half derivatives
                u_ikj = U[u_index + IDX(i, j, k, Nx, Ny, Nz)];
                v_ijk = U[v_index + IDX(i, j, k, Nx, Ny, Nz)];
                u_iphjk = 0.5 * (u_ip1jk + u_ikj);
                u_imhjk = 0.5 * (u_ikj + u_im1jk);
                v_ijphk = 0.5 * (v_ijp1k + v_ijk);
                v_ijmhk = 0.5 * (v_ijk + v_ijm1k);
                ux = (u_iphjk - u_imhjk) / dx; // du/dx
                vy = (v_ijphk - v_ijmhk) / dy; // dv/dy
                // ux = (u_ip1jk - u_im1jk) / (2 * dx); // du/dx
                // vy = (v_ijp1k - v_ijm1k) / (2 * dy); // dv/dy
                // RHS of Poisson problem
                // f[IDX(i, j, k, Nx, Ny, Nz)] = rho * (ux + vy + wz) / dt;                
                // Fill p with zeros
                // if (ux > 10)
                //     printf("ux[%d,%d,%d]: %.14f\n", i, j, k, ux);
                // if (vy > 10)
                //     printf("vy[%d,%d,%d]: %.14f\n", i, j, k, vy);
                // if (wz > 10)
                //     printf("wz[%d,%d,%d]: %.14f\n", i, j, k, wz);
                if (i < Nx - 1 && j < Ny - 1 && k < Nz - 1) {
                    // Compute rho / dt * div(U) and store it for many DFT (contiguous z slices)
                    // if (i > 0 && j > 0 && k > 0) {
                        // printf("u[%d,%d,%d]: %.14f\n", i, j, k, U[u_index + IDX(i, j, k, Nx, Ny, Nz)]);
                        // printf("u_ip1jk = %.14f, u_im1jk = %.14f\n", u_ip1jk, u_im1jk);
                        // printf("ux[%d,%d,%d]: %.14f\n", i, j, k, ux);
                        // printf("v[%d,%d,%d]: %.14f\n", i, j, k, U[v_index + IDX(i, j, k, Nx, Ny, Nz)]);
                        // printf("v_ijp1k = %.14f, v_ijm1k = %.14f\n", v_ijp1k, v_ijm1k);
                        // printf("vy[%d,%d,%d]: %.14f\n", i, j, k, vy);
                        // printf("w[%d,%d,%d]: %.14f\n", i, j, k, U[w_index + IDX(i, j, k, Nx, Ny, Nz)]);
                        // printf("w_ijkp1 = %.14f, w_ijkm1 = %.14f\n", w_ijkp1, w_ijkm1);
                        // printf("wz[%d,%d,%d]: %.14f\n", i, j, k, wz);
                    // }
                    f = rho * (ux + vy + wz) / dt;
                    // if (f != 0)
                    //     printf("f[%d,%d,%d] = %.8f\n", i, j, k, f);
                    f_in[j + (Ny - 1) * i + (Nx - 1) * (Ny - 1) * k] = f; //f[IDX(i, j, k, Nx, Ny, Nz)] + I * 0.0;
                    f_out[j + (Ny - 1) * i + (Nx - 1) * (Ny - 1) * k] = 0.0;
                    p[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)] = 0.0;
                    // if (abs(f_in[j + (Ny - 1) * i + (Nx - 1) * (Ny - 1) * k]) > 1000)
                    // if (k == 0)
                    // printf("f_in[%d,%d,%d]: %.14f + %.14fi\n", i, j, k, creal(f_in[j + (Ny - 1) * i + (Nx - 1) * (Ny - 1) * k]), cimag(f_in[j + (Ny - 1) * i + (Nx - 1) * (Ny - 1) * k]));
                }
                // Fill a and c
                if (i == 0 && j == 0 && k < Nz - 2) {
                    a[k] = 1.0 / (dz * dz);
                    c[k] = 1.0 / (dz * dz);
                }
                // Fill p_top
                // if (k == Nz - 1) {
                if (k == Nz - 1 && j < Ny - 1 && i < Nx - 1) {
                    // p_top_in[j + Ny * i] = p[IDX(i, j, k, Nx, Ny, Nz)];
                    p_top_in[j + (Ny - 1) * i] = p[IDX(i, j, k, Nx, Ny, Nz)];
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
    // Destroy plans
    // fftw_destroy_plan(p_top_plan);
    // fftw_destroy_plan(f_plan);
    
    // Compute r,s systems of equations
    for (int r = 0; r < Nx - 1; r++) {
        for (int s = 0; s < Ny - 1; s++) {
            gamma_rs = - 2 - kx[r] * kx[r] - ky[s] * ky[s];
            // printf("r: %d, s: %d\n", r, s);
            // // printf("kx[r]: %f\n", kx[r]);
            // // printf("ky[s]: %f\n", ky[s]);
            // printf("gamma_rs: %f\n", gamma_rs);
            // First equation
            // f[IDX(r, s, 0, Nx, Ny, Nz)] = 0.0 + 0.5 * dz * f[IDX(r, s, 1, Nx, Ny, Nz)];
            f_out[s + (Ny - 1) * r] = 0.0 + 0.0 * I + 0.5 * dz * f_out[s + (Ny - 1) * r + (Nx - 1) * (Ny - 1) * 1];
            // Last equation
            // f[IDX(r, s, Nz - 2, Nx, Ny, Nz)] -= p_top_out[IDX(r, s, 0, Nx, Ny, 1)] / (dz * dz);
            f_out[s + (Ny - 1) * r + (Nx - 1) * (Ny - 1) * (Nz - 2)] -= p_top_out[s + (Ny - 1) * r] / (dz * dz);
            // printf("r: %d, s: %d, k: 0, f_out: %f + %fi\n", r, s, creal(f_out[s + (Ny - 1) * r]), cimag(f_out[s + (Ny - 1) * r]));
            // printf("r: %d, s: %d, k: Nz - 2, f_out: %f + %fi\n", r, s, creal(f_out[s + (Ny - 1) * r + (Nx - 1) * (Ny - 1) * (Nz - 2)]), cimag(f_out[s + (Ny - 1) * r + (Nx - 1) * (Ny - 1) * (Nz - 2)]));
            // Fill diagonal elements of A and RHS, and temporal arrays
            for (int k = 0; k < Nz - 1; k++) {
                b[k] = gamma_rs / (dz * dz);
                d[k] = f_out[s + (Ny - 1) * r + (Nx - 1) * (Ny - 1) * k];
                u[k] = 0.0;
                y[k] = 0.0;
                pk[k] = 0.0; // To store the solution
                if (k < Nz - 2) 
                    l[k] = 0.0;
                // printf("d[k] = %f + %fi\n", creal(d[k]), cimag(d[k]));
                // if (r > 0 && s > 0 && k > 0)
                //     printf("F_out[%d,%d,%d] = %.14f + %.14fi\n", r, s, k, creal(f_out[s + (Ny - 1) * r + (Nx - 1) * (Ny - 1) * k]), cimag(f_out[s + (Ny - 1) * r + (Nx - 1) * (Ny - 1) * k]));
            }
            // Fix first coefficient of b and c
            b[0] =  -1.0 / dz;
            c[0] = (2.0 + 0.5 * gamma_rs) / dz;

            // if (abs(f_out[s + (Ny - 1) * r]) > 100)
            //     printf("r: %d, s: %d, f_out: %f + %fi\n", r, s, creal(f_out[s + (Ny - 1) * r]), cimag(f_out[s + (Ny - 1) * r]));

            // Solve tridiagonal system
            thomas_algorithm(a, b, c, d, pk, l, u, y, Nz - 1);
            
            // printf("gamma_rs = %.14f\n", gamma_rs);
            // Fill p_in with solution
            for (int k = 0; k < Nz - 1; k++) {
                // Check if the solution is nan
                // if (isnan(creal(pk[k]))) {
                //     printf("NaN in: r = %d, s = %d, k = %d\n", r, s, k);
                // }
                p_in[s + (Ny - 1) * r + (Nx - 1) * (Ny - 1) * k] = pk[k];
                // if (cabs(pk[k]) > 100) {
                // printf("p[%d, %d, %d] = %e + %ei\n", r, s, k, creal(pk[k]), cimag(pk[k]));                    
                // if (k < Nz - 2)
                //     printf("a[%d] = %.14f\n", k, a[k]);
                // printf("b[%d] = %.14f\n", k, b[k]);
                // if (k < Nz - 2)
                //     printf("c[%d] = %.14f\n", k, c[k]);
                // // printf("d[%d] = %e + %ei\n", k, creal(d[k]), cimag(d[k]));
                // printf("d[%d] = %.14f\n", k, creal(d[k]));
                // printf("\n");
                // }
            }
        }
    }

    // Plan for IFFT2(p) for each z slice
    p_plan = fftw_plan_many_dft(2, n, howmany, p_in, inembed, istride, idist, p_out, onembed, ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);
    // Compute IFFT
    fftw_execute(p_plan);
    // Destroy plan
    // fftw_destroy_plan(p_plan);

    // Restore data shape, compute real part and fill boundary conditions
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                if (i < Nx - 1 && j < Ny - 1 && k < Nz - 1) {
                    p[IDX(i, j, k, Nx, Ny, Nz)] = creal(p_out[j + (Ny - 1) * i + (Nx - 1) * (Ny - 1) * k] / ((Nx - 1) * (Ny - 1)));
                }
                // Periodic boundary conditions on xy
                if (i == Nx - 1) { // Right boundary
                    p[IDX(i, j, k, Nx, Ny, Nz)] = p[IDX(0, j, k, Nx, Ny, Nz)];
                }
                if (j == Ny - 1) { // Front boundary
                    p[IDX(i, j, k, Nx, Ny, Nz)] = p[IDX(i, 0, k, Nx, Ny, Nz)];
                }
                // if (p[IDX(i, j, k, Nx, Ny, Nz)] > 200) {
                //     printf("p[%d,%d,%d]: %.14f\n", i, j, k, p[IDX(i, j, k, Nx, Ny, Nz)]);
                //     // int r = i;
                //     // int s = j;
                //     // printf("r: %d, s: %d, f_out: %f + %fi\n", r, s, creal(f_out[s + (Ny - 1) * r]), cimag(f_out[s + (Ny - 1) * r]));
                // }
            }
        }
    }

    // fftw_destroy_plan(p_top_plan);
    // fftw_destroy_plan(f_plan);
    // fftw_destroy_plan(p_plan);
    // fftw_free(p_top_in);
    // fftw_free(p_in);
    // fftw_free(f_in);
    // fftw_free(p_top_out);
    // fftw_free(f_out);
    // fftw_free(p_out);
}