#include "../include/poisson.h"
#include <stdio.h>

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
// void thomas_algorithm(double *a, double *b, double *c, double complex *d, double complex *x, double complex *l, double complex *u, double complex *y, int N) {
// void thomas_algorithm(double *a, double *b, double *c, fftw_complex *d, fftw_complex *x, fftw_complex *l, fftw_complex *u, fftw_complex *y, int N) {
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

void thomas_algorithm_debug(fftw_complex *a, fftw_complex *b, fftw_complex *c, fftw_complex *d, fftw_complex *x, fftw_complex *l, fftw_complex *u, 
    // fftw_complex *y, fftw_complex *U, fftw_complex *L, fftw_complex *Y, fftw_complex *D, int r, int s, int N) {
    fftw_complex *y, fftw_complex *X, fftw_complex *V, fftw_complex *L, fftw_complex *Y, fftw_complex *D, int r, int s, int Nx, int Ny, int N) {
    // Create L and U
    u[0] = b[0];
    V[IDX(r, s, 0, Nx, Ny, N)] = u[0];
    for (int i = 1; i < N; i++) {
        l[i-1] = a[i-1] / u[i - 1];
        u[i] = b[i] - l[i-1] * c[i - 1];
        V[IDX(r, s, i, Nx, Ny, N)] = u[i];
        L[IDX(r, s, i-1, Nx, Ny, N - 1)] = l[i-1];
    }
    // Solve Ly = d
    y[0] = d[0];
    Y[IDX(r, s, 0, Nx, Ny, N)] = y[0];
    D[IDX(r, s, 0, Nx, Ny, N)] = d[0];
    for (int i = 1; i < N; i++) {
        y[i] = d[i] - l[i-1] * y[i - 1];
        Y[IDX(r, s, i, Nx, Ny, N)] = y[i];
        D[IDX(r, s, i, Nx, Ny, N)] = d[i];
    }
    // Solve Ux = y
    x[N - 1] = y[N - 1] / u[N - 1];
    X[IDX(r, s, N-1, Nx, Ny, N)] = x[N - 1];
    for (int i = N - 2; i >= 0; i--) {
        x[i] = (y[i] - c[i] * x[i+1]) / u[i];
        X[IDX(r, s, i, Nx, Ny, N)] = x[i];
    }
    // for (int i = N - 1; i >= 0; i--) {
    //     x[i-1] = (y[i-1] - c[i] * x[i]) / u[i-1];
    //     X[IDX(r, s, i, Nx, Ny, N)] = x[i];
    // }
}

// void solve_pressure(double *U, double *p, double *a, double *b, double *c, double complex *d, double complex *l, double complex *u, double complex *y, double complex *pk, fftw_plan p_plan, fftw_plan f_plan, fftw_plan p_top_plan, fftw_complex *f_in, fftw_complex *f_out, fftw_complex *p_top_in, fftw_complex *p_top_out, fftw_complex *p_in, fftw_complex *p_out, Parameters *parameters) {
// void solve_pressure(double *U, double *p, double *a, double *b, double *c, fftw_complex *d, fftw_complex *l, fftw_complex *u, fftw_complex *y, fftw_complex *pk, fftw_plan p_plan, fftw_plan f_plan, fftw_plan p_top_plan, fftw_complex *f_in, fftw_complex *f_out, fftw_complex *p_top_in, fftw_complex *p_top_out, fftw_complex *p_in, fftw_complex *p_out, Parameters *parameters) {
// void solve_pressure(double *U, double *p, double complex *a, double complex *b, double complex *c, double complex *d, double complex *l, double complex *u, double complex *y, double complex *pk, Parameters *parameters) {
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
    // double u_ip1jk = 0.0, u_im1jk = 0.0, u_ijk = 0.0, u_iphjk = 0.0, u_imhjk = 0.0;
    // double v_ijp1k = 0.0, v_ijm1k = 0.0, v_ijk = 0.0, v_ijphk = 0.0, v_ijmhk = 0.0;
    // double w_ijk = 0.0, w_ijkp1 = 0.0, w_ijkm1 = 0.0, w_ijkp2 = 0.0, w_ijkm2 = 0.0, w_ijkph = 0.0, w_ijkmh = 0.0;
    // double ux = 0.0, vy = 0.0, wz = 0.0;
    double u_ijk, u_ip1jk, u_im1jk, u_iphjk, u_imhjk;
    double v_ijk, v_ijp1k, v_ijm1k, v_ijphk, v_ijmhk;
    double w_ijk, w_ijkp1, w_ijkm1, w_ijkp2, w_ijkm2, w_ijkph, w_ijkmh;
    double ux, vy, wz;
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
    // FILE *debug_pressure = fopen("data/output/debug_pressure.csv", "w");
    // fprintf(debug_pressure, "i, j, k, f_in_real, f_in_imag, f_out_real, f_out_imag\n");
    // double *A = (double *) malloc((Nx - 1) * (Ny - 1) * (Nz - 2) * sizeof(double));
    // double *C = (double *) malloc((Nx - 1) * (Ny - 1) * (Nz - 2) * sizeof(double));
    // double *B = (double *) malloc((Nx - 1) * (Ny - 1) * (Nz - 1) * sizeof(double));
    // double *FF = (double *) malloc((Nx - 1) * (Ny - 1) * (Nz - 1) * sizeof(double));
    // double *UX = (double *) malloc((Nx - 1) * (Ny - 1) * Nz * sizeof(double));
    // double *VY = (double *) malloc((Nx - 1) * (Ny - 1) * Nz * sizeof(double));
    // double *WZ = (double *) malloc((Nx - 1) * (Ny - 1) * Nz * sizeof(double));
    // double *UU = (double *) malloc((Nx - 1) * (Ny - 1) * (Nz - 1) * sizeof(double));
    // double *L = (double *) malloc((Nx - 1) * (Ny - 1) * (Nz - 2) * sizeof(double));
    // double *Y = (double *) malloc((Nx - 1) * (Ny - 1) * (Nz - 1) * sizeof(double));
    // fftw_complex *A = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 2));
    // fftw_complex *C = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 2));
    // // fftw_complex *B = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    // fftw_complex *V = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    // fftw_complex *L = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 2));
    // fftw_complex *Y = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    // fftw_complex *X = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    // fftw_complex *D = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    // FILE *debug_F = fopen("data/output/f.bin", "wb");
    // FILE *debug_F_K = fopen("data/output/F_k.bin", "wb");
    // // FILE *debug_P_top = fopen("data/output/p_top.bin", "wb");
    // // FILE *debug_p_in = fopen("data/output/p_in.bin", "wb");
    // FILE *debug_A = fopen("data/output/A.bin", "wb");
    // FILE *debug_C = fopen("data/output/C.bin", "wb");
    // // FILE *debug_B = fopen("data/output/B.bin", "wb");
    // // FILE *debug_ux = fopen("data/output/ux.bin", "wb");
    // // FILE *debug_vy = fopen("data/output/vy.bin", "wb");
    // // FILE *debug_wz = fopen("data/output/wz.bin", "wb");
    // FILE *debug_U = fopen("data/output/U.bin", "wb");
    // FILE *debug_L = fopen("data/output/L.bin", "wb");
    // FILE *debug_Y = fopen("data/output/Y.bin", "wb");
    // FILE *debug_X = fopen("data/output/X.bin", "wb");
    // FILE *debug_D = fopen("data/output/D.bin", "wb");
    // // FILE *debug_F = fopen("data/output/f.bin", "wb");
    // // FILE *debug_ux = fopen("data/output/ux.bin", "wb");
    // // FILE *debug_vy = fopen("data/output/vy.bin", "wb");
    // // FILE *debug_wz = fopen("data/output/wz.bin", "wb");

    double f;
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
                // Debugging div(U)
                // if ((i < Nx - 1) && (j < Ny - 1)) {
                //     UX[IDX(i, j, k, Nx - 1, Ny - 1, Nz)] = ux;
                //     VY[IDX(i, j, k, Nx - 1, Ny - 1, Nz)] = vy;
                //     WZ[IDX(i, j, k, Nx - 1, Ny - 1, Nz)] = wz;
                // }
                if (i < Nx - 1 && j < Ny - 1 && k < Nz - 1) {
                    // Compute rho / dt * div(U) and store it for many DFT (contiguous z slices)
                    f = rho * (ux + vy + wz) / dt;
                    // FF[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)] = f;
                    // f_in[j + (Ny - 1) * i + (Nx - 1) * (Ny - 1) * k] = f; //f[IDX(i, j, k, Nx, Ny, Nz)] + I * 0.0;
                    f_in[FFTWIDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)] = f;
                }
                // Fill a and c
                if (i == 0 && j == 0 && k < Nz - 2) {
                    a[k] = 1.0 / (dz * dz) + 0.0 * I;
                    c[k] = 1.0 / (dz * dz) + 0.0 * I;
                }
                // Fill p_top
                if (k == Nz - 1 && j < Ny - 1 && i < Nx - 1) {
                    // p_top_in[j + Ny * i] = p[IDX(i, j, k, Nx, Ny, Nz)];
                    // p_top_in[j + (Ny - 1) * i] = p[IDX(i, j, k, Nx, Ny, Nz)];
                    p_top_in[FFTWIDX(i, j, 0, Nx - 1, Ny - 1, Nz - 1)] = p[IDX(i, j, k, Nx, Ny, Nz)];
                }
                // // Debug
                // if (k < Nz - 2) {
                //     A[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 2)] = a[k];
                //     C[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 2)] = c[k];
                // }
            }
        }
    }
    // fwrite(f_in, sizeof(fftw_complex), (Nx - 1) * (Ny - 1) * (Nz - 1), debug_F);
    // fwrite(UX, sizeof(double), (Nx - 1) * (Ny - 1) * (Nz), debug_ux);
    // fwrite(VY, sizeof(double), (Nx - 1) * (Ny - 1) * (Nz), debug_vy);
    // fwrite(WZ, sizeof(double), (Nx - 1) * (Ny - 1) * (Nz), debug_wz);
        // fwrite(L, sizeof(fftw_complex), (Nx - 1) * (Ny - 1) * (Nz - 2), debug_L);
    // fwrite(Y, sizeof(fftw_complex), (Nx - 1) * (Ny - 1) * (Nz - 1), debug_Y);
    // fwrite(D, sizeof(fftw_complex), (Nx - 1) * (Ny - 1) * (Nz - 1), debug_D);
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
    // fwrite(p_top_out, sizeof(fftw_complex), (Nx - 1) * (Ny - 1), debug_P_top);

    
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
                // d[k] = f_out[s + (Ny - 1) * r + (Nx - 1) * (Ny - 1) * k];
                d[k] = f_out[FFTWIDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)];
                // D[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = d[k];
                // B[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = b[k];
                // u[k] = 0.0 + 0.0 * I;
                // y[k] = 0.0 + 0.0 * I;
                // pk[k] = 0.0 + 0.0 * I; // To store the solution
                // if (k < Nz - 2) 
                //     l[k] = 0.0 + 0.0 * I;
            }
            // Fix first coefficient of b and c
            b[0] =  -1.0 / dz + 0.0 * I;
            c[0] = (2.0 + 0.5 * gamma_rs) / dz + 0.0 * I;
            // B[IDX(r, s, 0, Nx - 1, Ny - 1, Nz - 1)] = b[0];
            // C[IDX(r, s, 0, Nx - 1, Ny - 1, Nz - 2)] = c[0];

            // Solve tridiagonal system
            thomas_algorithm(a, b, c, d, pk, l, u, y, Nz - 1);
            // thomas_algorithm_debug(a, b, c, d, pk, l, u, y, V, L, Y, D, r, s, Nz - 1);
            // thomas_algorithm_debug(a, b, c, d, pk, l, u, y, X, V, L, Y, D, r, s, Nx - 1, Ny - 1, Nz - 1);
            
            // Fill p_in with solution 
            for (int k = 0; k < Nz - 1; k++) {
                // p_in[s + (Ny - 1) * r + (Nx - 1) * (Ny - 1) * k] = pk[k];
                p_in[FFTWIDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = pk[k];                
            }
        }
    }
    // fwrite(f_out, sizeof(fftw_complex), (Nx - 1) * (Ny - 1) * (Nz - 1), debug_F_K);
    // // fwrite(p_in, sizeof(fftw_complex), (Nx - 1) * (Ny - 1) * (Nz - 1), debug_p_in);
    // fwrite(A, sizeof(fftw_complex), (Nx - 1) * (Ny - 1) * (Nz - 2), debug_A);
    // fwrite(C, sizeof(fftw_complex), (Nx - 1) * (Ny - 1) * (Nz - 2), debug_C);
    // // fwrite(B, sizeof(fftw_complex), (Nx - 1) * (Ny - 1) * (Nz - 1), debug_B);
    // fwrite(V, sizeof(fftw_complex), (Nx - 1) * (Ny - 1) * (Nz - 1), debug_U);
    // fwrite(L, sizeof(fftw_complex), (Nx - 1) * (Ny - 1) * (Nz - 2), debug_L);
    // fwrite(Y, sizeof(fftw_complex), (Nx - 1) * (Ny - 1) * (Nz - 1), debug_Y);
    // fwrite(X, sizeof(fftw_complex), (Nx - 1) * (Ny - 1) * (Nz - 1), debug_X);
    // fwrite(D, sizeof(fftw_complex), (Nx - 1) * (Ny - 1) * (Nz - 1), debug_D);
    // fclose(debug_pressure);

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

    // free(FF);
    // free(A);
    // free(C);
    // free(B);
    // free(UX);
    // free(VY);
    // free(WZ);
    // fftw_free(A);
    // fftw_free(C);
    // fftw_free(B);
    // fftw_free(UU);
    // fftw_free(L);
    // fftw_free(YY);
    // fftw_free(D);


    // fclose(debug_F);
    // // fclose(debug_ux);
    // // fclose(debug_vy);
    // // fclose(debug_wz);
    // fclose(debug_F_K);
    // // fclose(debug_P_top);
    // // fclose(debug_p_in);
    // fclose(debug_A);
    // fclose(debug_C);
    // // fclose(debug_B);
    // fclose(debug_U);
    // fclose(debug_L);
    // fclose(debug_Y);
    // fclose(debug_X);
    // fclose(debug_D);

    fftw_destroy_plan(p_top_plan);
    fftw_destroy_plan(f_plan);
    fftw_destroy_plan(p_plan);
    // fftw_free(p_top_in);
    // fftw_free(p_in);
    // fftw_free(f_in);
    // fftw_free(p_top_out);
    // fftw_free(f_out);
    // fftw_free(p_out);
}