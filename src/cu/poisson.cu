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

#include "../../include/cu/poisson.cuh"

__global__
// void thomas_algorithm(cufftDoubleComplex *a, cufftDoubleComplex *b, cufftDoubleComplex *c, cufftDoubleComplex *d, cufftDoubleComplex *x, cufftDoubleComplex *l, cufftDoubleComplex *u, cufftDoubleComplex *y, Parameters parameters) {
void thomas_algorithm(double *a, double *b, double *c, cufftDoubleComplex *d, cufftDoubleComplex *x, cufftDoubleComplex *l, cufftDoubleComplex *u, cufftDoubleComplex *y, Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int size = (Nx - 1) * (Ny - 1);
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;
    for (int ij = idx; ij < size; ij += stride) {
        int i = ij / (Ny - 1);
        int j = ij % (Ny - 1);
        u[IDX(i, j, 0, Nx - 1, Ny - 1, Nz - 1)] = make_cuDoubleComplex(b[IDX(i, j, 0, Nx - 1, Ny - 1, Nz - 1)], 0.0);
        for (int k = 1; k < Nz - 1; k++) {
            l[IDX(i, j, k - 1, Nx - 1, Ny - 1, Nz - 2)] = cuCdiv(make_cuDoubleComplex(a[IDX(i, j, k - 1, Nx - 1, Ny - 1, Nz - 2)], 0.0), u[IDX(i, j, k - 1, Nx - 1, Ny - 1, Nz - 1)]); 
            u[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)] = cuCsub(make_cuDoubleComplex(b[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)], 0.0), cuCmul(l[IDX(i, j, k - 1, Nx - 1, Ny - 1, Nz - 2)], make_cuDoubleComplex(c[IDX(i, j, k - 1, Nx - 1, Ny - 1, Nz - 2)], 0.0)));
        }
        y[IDX(i, j, 0, Nx - 1, Ny - 1, Nz - 1)] = d[IDX(i, j, 0, Nx - 1, Ny - 1, Nz - 1)];
        for (int k = 1; k < Nz - 1; k++) {
            y[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)] = cuCsub(d[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)], cuCmul(l[IDX(i, j, k - 1, Nx - 1, Ny - 1, Nz - 2)], y[IDX(i, j, k - 1, Nx - 1, Ny - 1, Nz - 1)]));
        }
        x[FFTWIDX(i, j, Nz - 2, Nx - 1, Ny - 1, Nz - 1)] = cuCdiv(y[IDX(i, j, Nz - 2, Nx - 1, Ny - 1, Nz - 1)], u[IDX(i, j, Nz - 2, Nx - 1, Ny - 1, Nz - 1)]);
        for (int k = Nz - 3; k >= 0; k--) {
            x[FFTWIDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)] = cuCdiv(cuCsub(y[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)], cuCmul(make_cuDoubleComplex(c[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 2)], 0.0), x[FFTWIDX(i, j, k + 1, Nx - 1, Ny - 1, Nz - 1)])), u[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)]);
        }
        // u[IDX(i, j, 0, Nx - 1, Ny - 1, Nz - 1)].x = b[IDX(i, j, 0, Nx - 1, Ny - 1, Nz - 1)];
        // for (int k = 1; k < Nz - 1; k++) {
        //     l[IDX(i, j, k - 1, Nx - 1, Ny - 1, Nz - 2)] = a[IDX(i, j, k - 1, Nx - 1, Ny - 1, Nz - 2)] / u[IDX(i, j, k - 1, Nx - 1, Ny - 1, Nz - 1)]; 
        //     u[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)] = b[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)] - l[IDX(i, j, k - 1, Nx - 1, Ny - 1, Nz - 2)] * c[IDX(i, j, k - 1, Nx - 1, Ny - 1, Nz - 2)];
        // }
        // y[IDX(i, j, 0, Nx - 1, Ny - 1, Nz - 1)] = d[IDX(i, j, 0, Nx - 1, Ny - 1, Nz - 1)];
        // for (int k = 1; k < Nz - 1; k++) {
        //     y[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)] = cuCsub(d[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)], cuCmul(l[IDX(i, j, k - 1, Nx - 1, Ny - 1, Nz - 2)], y[IDX(i, j, k - 1, Nx - 1, Ny - 1, Nz - 1)]));
        // }
        // x[FFTWIDX(i, j, Nz - 2, Nx - 1, Ny - 1, Nz - 1)] = cuCdiv(y[IDX(i, j, Nz - 2, Nx - 1, Ny - 1, Nz - 1)], u[IDX(i, j, Nz - 2, Nx - 1, Ny - 1, Nz - 1)]);
        // for (int k = Nz - 3; k >= 0; k--) {
        //     x[FFTWIDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)] = cuCdiv(cuCsub(y[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)], cuCmul(c[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 2)], x[FFTWIDX(i, j, k + 1, Nx - 1, Ny - 1, Nz - 1)])), u[IDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)]);
        // }
    }
}

__global__
void compute_f(double *U, cufftDoubleComplex *f_in, cufftDoubleComplex *p_top_in, double *p, Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int u_index = parameters.field_indexes.u;
    int v_index = parameters.field_indexes.v;
    int w_index = parameters.field_indexes.w;
    int im1, ip1, jm1, jp1;
    double dx = parameters.dx;
    double dy = parameters.dy;
    double dz = parameters.dz;
    double dt = parameters.dt;
    double rho = parameters.rho;
    double u_ijk, u_ip1jk, u_im1jk, u_iphjk, u_imhjk;
    double v_ijk, v_ijp1k, v_ijm1k, v_ijphk, v_ijmhk;
    double w_ijk, w_ijkp1, w_ijkm1, w_ijkp2, w_ijkm2, w_ijkph, w_ijkmh;
    double ux, vy, wz, f;
    // Loop over nodes to compute f
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;
    int size = Nx * Ny * Nz;
    for (int ijk = idx; ijk < size; ijk += stride) {
        int i = ijk / (Ny * Nz);
        int j = (ijk % (Ny * Nz)) / Nz;
        int k = ijk % Nz;               
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
            f_in[FFTWIDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)] = make_cuDoubleComplex(f, 0.0);
        }
        // Fill p_top
        if (k == Nz - 1 && j < Ny - 1 && i < Nx - 1) {
            p_top_in[FFTWIDX(i, j, 0, Nx - 1, Ny - 1, Nz - 1)] = make_cuDoubleComplex(p[IDX(i, j, k, Nx, Ny, Nz)], 0.0);
        }
    }
}

__global__
void post_fft(double *p, cufftDoubleComplex *p_out, Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    // Loop over nodes to compute f
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;
    int size = Nx * Ny * Nz;
    for (int ijk = idx; ijk < size; ijk += stride) {
        int i = ijk / (Ny * Nz);
        int j = (ijk % (Ny * Nz)) / Nz;
        int k = ijk % Nz;  
        if (i < Nx - 1 && j < Ny - 1 && k < Nz - 1) {
            p[IDX(i, j, k, Nx, Ny, Nz)] = cuCreal(cuCdiv(p_out[FFTWIDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)], make_cuDoubleComplex((Nx - 1) * (Ny - 1), 0.0)));
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

__global__
// void gammas_and_coefficients(double *kx, double *ky, double *gamma, cufftDoubleComplex *a, cufftDoubleComplex *b, cufftDoubleComplex *c, Parameters parameters) {
void gammas_and_coefficients(double *kx, double *ky, double *gamma, double *a, double *b, double *c, Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    double dz = parameters.dz;
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;
    int size = (Nx - 1) * (Ny - 1) * (Nz - 2);
    for (int rsk = idx; rsk < size; rsk += stride) {
        int r = rsk / ((Ny - 1) * (Nz - 2));
        int s = (rsk % ((Ny - 1) * (Nz - 2))) / (Nz - 2);
        int k = rsk % (Nz - 2);
        if (k == 0) { // Use k = 0 to fill gamma matrix and first coefficients of a, b and c
            gamma[FFTWIDX(r, s, k, Nx - 1, Ny - 1, 0)] = -2 - kx[r] * kx[r] - ky[s] * ky[s];
            // a[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 2)] = make_cuDoubleComplex(1.0 / (dz * dz), 0.0); 
            // b[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = make_cuDoubleComplex(-1.0 / dz, 0.0);
            // c[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 2)] = make_cuDoubleComplex((2.0 + 0.5 * gamma[FFTWIDX(r, s, k, Nx - 1, Ny - 1, 0)]) / dz, 0.0);
            a[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 2)] = 1.0 / (dz * dz), 0.0; 
            b[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = -1.0 / dz, 0.0;
            c[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 2)] = (2.0 + 0.5 * gamma[FFTWIDX(r, s, k, Nx - 1, Ny - 1, 0)]) / dz;
        } else { // The rest of the coefficients a and c
            // a[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 2)] = make_cuDoubleComplex(1.0 / (dz * dz), 0.0); 
            // c[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 2)] = make_cuDoubleComplex(1.0 / (dz * dz), 0.0);
            a[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 2)] = 1.0 / (dz * dz); 
            c[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 2)] = 1.0 / (dz * dz);
        }
    }
}

__global__ 
// void update_coefficients(double *gamma, cufftDoubleComplex *b,  cufftDoubleComplex *d, cufftDoubleComplex *f_out, cufftDoubleComplex *p_top_out, Parameters parameters) {
void update_coefficients(double *gamma, double *b,  cufftDoubleComplex *d, cufftDoubleComplex *f_out, cufftDoubleComplex *p_top_out, Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    double dz = parameters.dz;
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;
    int size = (Nx - 1) * (Ny - 1) * (Nz - 1);
    for (int rsk = idx; rsk < size; rsk += stride) {
        int r = rsk / ((Ny - 1) * (Nz - 1));
        int s = (rsk % ((Ny - 1) * (Nz - 1))) / (Nz - 1);
        int k = rsk % (Nz - 1);
        if (k == 0) {
            // First equation k=0
            // f_out[FFTWIDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = cuCmul(make_cuDoubleComplex(0.5 * dz, 0.0), f_out[FFTWIDX(r, s, 1, Nx - 1, Ny - 1, Nz - 1)]);
            // d[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = f_out[FFTWIDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)];
            d[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = cuCmul(make_cuDoubleComplex(0.5 * dz, 0.0), f_out[FFTWIDX(r, s, 1, Nx - 1, Ny - 1, Nz - 1)]);
            // Last equation k=Nz-2
            // f_out[FFTWIDX(r, s, Nz - 2, Nx - 1, Ny - 1, Nz - 1)] = cuCsub(f_out[FFTWIDX(r, s, Nz - 2, Nx - 1, Ny - 1, Nz - 1)], cuCdiv(p_top_out[FFTWIDX(r, s, 0, Nx - 1, Ny - 1, Nz - 1)], make_cuDoubleComplex(dz * dz, 0.0)));
            // d[IDX(r, s, Nz - 2, Nx - 1, Ny - 1, Nz - 1)] = f_out[FFTWIDX(r, s, Nz - 2, Nx - 1, Ny - 1, Nz - 1)];
            d[IDX(r, s, Nz - 2, Nx - 1, Ny - 1, Nz - 1)] = cuCsub(f_out[FFTWIDX(r, s, Nz - 2, Nx - 1, Ny - 1, Nz - 1)], cuCdiv(p_top_out[FFTWIDX(r, s, 0, Nx - 1, Ny - 1, Nz - 1)], make_cuDoubleComplex(dz * dz, 0.0)));
        } else {
            if (k < Nz - 1) {
                // b[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = make_cuDoubleComplex(gamma[FFTWIDX(r, s, 0, Nx - 1, Ny - 1, 0)] / (dz * dz), 0.0);
                b[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = gamma[FFTWIDX(r, s, 0, Nx - 1, Ny - 1, 0)] / (dz * dz);
                if (k < Nz - 2)
                    d[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = f_out[FFTWIDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)];
            }
        }
    }
}

// void solve_pressure(double *U, double *p, double *gamma, cufftDoubleComplex *a, cufftDoubleComplex *b, cufftDoubleComplex *c, cufftDoubleComplex *d, cufftDoubleComplex *l, cufftDoubleComplex *u, cufftDoubleComplex *y, cufftDoubleComplex *pk, cufftHandle p_plan, cufftHandle f_plan, cufftHandle p_top_plan, cufftDoubleComplex *f_in, cufftDoubleComplex *f_out, cufftDoubleComplex *p_top_in, cufftDoubleComplex *p_top_out, cufftDoubleComplex *p_in, cufftDoubleComplex *p_out, Parameters parameters) {
// void solve_pressure(double *U, double *p, double *gamma, cufftDoubleComplex *a, cufftDoubleComplex *b, cufftDoubleComplex *c, cufftDoubleComplex *d, cufftDoubleComplex *l, cufftDoubleComplex *u, cufftDoubleComplex *y, cufftDoubleComplex *pk, cufftDoubleComplex *f_in, cufftDoubleComplex *f_out, cufftDoubleComplex *p_top_in, cufftDoubleComplex *p_top_out, cufftDoubleComplex *p_in, cufftDoubleComplex *p_out, Parameters parameters) {
// void solve_pressure(double *U, double *p, double *gamma, cufftDoubleComplex *a, cufftDoubleComplex *b, cufftDoubleComplex *c, cufftDoubleComplex *d, cufftDoubleComplex *l, cufftDoubleComplex *u, cufftDoubleComplex *y, cufftDoubleComplex *data_in, cufftDoubleComplex *data_out, cufftDoubleComplex *p_top_in, cufftDoubleComplex *p_top_out, Parameters parameters) {
void solve_pressure(double *U, double *p, double *gamma, double *a, double *b, double *c, cufftDoubleComplex *d, cufftDoubleComplex *l, cufftDoubleComplex *u, cufftDoubleComplex *y, cufftDoubleComplex *data_in, cufftDoubleComplex *data_out, cufftDoubleComplex *p_top_in, cufftDoubleComplex *p_top_out, Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int n[] = {Nx - 1, Ny - 1};
    int idist = (Nx - 1) * (Ny - 1);
    int odist = (Nx - 1) * (Ny - 1);
    int istride = 1;
    int ostride = 1;
    int howmany = Nz - 1;
    int *inembed = NULL;
    int *onembed = NULL;
    cufftHandle p_plan = 0, f_plan = 0, p_top_plan = 0;

    // Compute f = rho / dt * div(U)
    // printf("Compute f");
    // compute_f<<<BLOCKS, THREADS>>>(U, f_in, p_top_in, p, parameters);
    compute_f<<<BLOCKS, THREADS>>>(U, data_in, p_top_in, p, parameters);
    // printf("Compute f ok\n");
    CHECK(cudaDeviceSynchronize());
    
    // Plans for FFT2D
    CHECK_CUFFT(cufftPlanMany(&p_plan, 2, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, howmany)); // FFT2(f_k) for each z slice
    CHECK_CUFFT(cufftPlan2d(&p_top_plan, Nx - 1, Ny - 1, CUFFT_Z2Z)); // FFT2(p_top)

    // Compute FFT2D
    CHECK_CUFFT(cufftExecZ2Z(p_top_plan, p_top_in, p_top_out, CUFFT_FORWARD)); // FFT2(p_top)
    // CHECK_CUFFT(cufftExecZ2Z(p_plan, f_in, f_out, CUFFT_FORWARD)); // FFT2D(f_k) for each z slice
    CHECK_CUFFT(cufftExecZ2Z(p_plan, data_in, data_out, CUFFT_FORWARD)); // FFT2D(f_k) for each z slice

    CHECK(cudaDeviceSynchronize());

    // Update coefficients, including f in pseudo-Fourier space
    // update_coefficients<<<BLOCKS, THREADS>>>(gamma, b, d, f_out, p_top_out, parameters);
    update_coefficients<<<BLOCKS, THREADS>>>(gamma, b, d, data_out, p_top_out, parameters);
    // printf("Update coefficients ok\n");
    CHECK(cudaDeviceSynchronize());

    // Compute r,s systems of equations using thomas algorithm
    // thomas_algorithm<<<BLOCKS, THREADS>>>(a, b, c, d, p_in, l, u, y, parameters);
    thomas_algorithm<<<BLOCKS, THREADS>>>(a, b, c, d, data_in, l, u, y, parameters);
    // printf("Thomas algorithm ok\n");
    CHECK(cudaDeviceSynchronize());

    // Plan for IFFT2D(p) for each z slice
    CHECK_CUFFT(cufftPlanMany(&p_plan, 2, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, howmany));

    // Compute IFFT2D
    // CHECK_CUFFT(cufftExecZ2Z(p_plan, p_in, p_out, CUFFT_INVERSE));
    CHECK_CUFFT(cufftExecZ2Z(p_plan, data_in, data_out, CUFFT_INVERSE));
    CHECK(cudaDeviceSynchronize());

    // Post FFT
    // printf("post fft\n");
    // post_fft<<<BLOCKS, THREADS>>>(p, p_out, parameters);
    post_fft<<<BLOCKS, THREADS>>>(p, data_out, parameters);
    // printf("Post FFT ok\n");
    CHECK(cudaDeviceSynchronize());
    // printf("post fft ok\n");

    // Destroy plans
    cufftDestroy(p_top_plan);
    cufftDestroy(p_plan);
    cufftDestroy(f_plan);
}