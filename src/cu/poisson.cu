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
void thomas_algorithm(cufftDoubleComplex *a, cufftDoubleComplex *b, cufftDoubleComplex *c, cufftDoubleComplex *d, cufftDoubleComplex *x, cufftDoubleComplex *l, cufftDoubleComplex *u, cufftDoubleComplex *y, int N) {
    int i = threadIdx.x;
    if (i == 0) {
        // Create L and U
        u[0] = b[0];
        for (int i = 1; i < N; i++) {
            l[i-1] = cuCdiv(a[i-1], u[i-1]); //l[i-1] = a[i-1] / u[i - 1];
            u[i] = cuCsub(b[i], cuCmul(l[i-1], c[i-1])); //b[i] - l[i-1] * c[i - 1];
        }
        // Solve Ly = d
        y[0] = d[0];
        for (int i = 1; i < N; i++) {
            // y[i] = d[i] - l[i-1] * y[i - 1];
            y[i] = cuCsub(d[i], cuCmul(l[i-1], y[i - 1]));
        }
        // Solve Ux = y
        // x[N - 1] = y[N - 1] / u[N - 1];
        x[N - 1] = cuCdiv(y[N - 1], u[N - 1]);
        for (int i = N - 2; i >= 0; i--) {
            // x[i] = (y[i] - c[i] * x[i+1]) / u[i];
            x[i] = cuCdiv(cuCsub(y[i], cuCmul(c[i], x[i+1])), u[i]);
        }
    }
}

void thomas_algorithm_cpu(cufftDoubleComplex *a, cufftDoubleComplex *b, cufftDoubleComplex *c, cufftDoubleComplex *d, cufftDoubleComplex *x, cufftDoubleComplex *l, cufftDoubleComplex *u, cufftDoubleComplex *y, int N) {
    // Create L and U
    u[0] = b[0];
    for (int i = 1; i < N; i++) {
        l[i-1] = cuCdiv(a[i-1], u[i-1]); //l[i-1] = a[i-1] / u[i - 1];
        u[i] = cuCsub(b[i], cuCmul(l[i-1], c[i-1])); //b[i] - l[i-1] * c[i - 1];
    }
    // Solve Ly = d
    y[0] = d[0];
    for (int i = 1; i < N; i++) {
        // y[i] = d[i] - l[i-1] * y[i - 1];
        y[i] = cuCsub(d[i], cuCmul(l[i-1], y[i - 1]));
    }
    // Solve Ux = y
    // x[N - 1] = y[N - 1] / u[N - 1];
    x[N - 1] = cuCdiv(y[N - 1], u[N - 1]);
    for (int i = N - 2; i >= 0; i--) {
        // x[i] = (y[i] - c[i] * x[i+1]) / u[i];
        x[i] = cuCdiv(cuCsub(y[i], cuCmul(c[i], x[i+1])), u[i]);
    }
}

__global__
void pre_fft(double *U, cufftDoubleComplex *f_in, cufftDoubleComplex *p_top_in, cufftDoubleComplex *a, cufftDoubleComplex *c, double *p, Parameters parameters) {
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
            f_in[FFTWIDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)].x = f;
            f_in[FFTWIDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)].y = 0.0;
        }
        // Fill a and c
        if (i == 0 && j == 0 && k < Nz - 2) {
            // a[k] = make_cuDoubleComplex(1.0 / (dz * dz), 0.0);
            // c[k] = make_cuDoubleComplex(1.0 / (dz * dz), 0.0);
            a[k].x = 1.0 / (dz * dz);
            a[k].y = 0.0;
            c[k].x = 1.0 / (dz * dz);
            c[k].y = 0.0;
        }
        // Fill p_top
        if (k == Nz - 1 && j < Ny - 1 && i < Nx - 1) {
            // p_top_in[FFTWIDX(i, j, 0, Nx - 1, Ny - 1, Nz - 1)] = make_cuDoubleComplex(p[IDX(i, j, k, Nx, Ny, Nz)], 0.0);
            p_top_in[FFTWIDX(i, j, 0, Nx - 1, Ny - 1, Nz - 1)].x = p[IDX(i, j, k, Nx, Ny, Nz)];
            p_top_in[FFTWIDX(i, j, 0, Nx - 1, Ny - 1, Nz - 1)].y = 0.0;
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
void update_values(cufftDoubleComplex *b, cufftDoubleComplex *c, cufftDoubleComplex *d, cufftDoubleComplex *f_out, cufftDoubleComplex *p_top_out, Parameters parameters, int r, int s, double gamma_rs) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int idx = threadIdx.x;
    double dz = parameters.dz;
    if (idx == 0) {
        // First equation k=0
        f_out[FFTWIDX(r, s, 0, Nx - 1, Ny - 1, Nz - 1)] = cuCmul(make_cuDoubleComplex(0.5 * dz, 0.0), f_out[FFTWIDX(r, s, 1, Nx - 1, Ny - 1, Nz - 1)]);
        // Last equation k=Nz-2
        f_out[FFTWIDX(r, s, Nz - 2, Nx - 1, Ny - 1, Nz - 1)] = cuCsub(f_out[FFTWIDX(r, s, Nz - 2, Nx - 1, Ny - 1, Nz - 1)], cuCdiv(p_top_out[FFTWIDX(r, s, 0, Nx - 1, Ny - 1, Nz - 1)], make_cuDoubleComplex(dz * dz, 0.0)));
    }
    __syncthreads();
    if (idx < Nz - 1) {
        b[idx] = make_cuDoubleComplex(gamma_rs / (dz * dz), 0.0);
        d[idx] = f_out[FFTWIDX(r, s, idx, Nx - 1, Ny - 1, Nz - 1)];
    }
    // Fix first coefficient of b and c
    __syncthreads();
    if (idx == 0) {
        b[idx] = make_cuDoubleComplex(-1.0 / dz, 0.0);
        c[idx] = make_cuDoubleComplex((2.0 + 0.5 * gamma_rs) / dz, 0.0);
    }
}

__global__
void copy_to_p_in(cufftDoubleComplex *p_in, cufftDoubleComplex *pk, Parameters parameters, int r, int s) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int idx = threadIdx.x;
    if (idx < Nz - 1) {
        p_in[FFTWIDX(r, s, idx, Nx - 1, Ny - 1, Nz - 1)] = pk[idx];
    }
}

void solve_pressure(double *U, double *p, double *kx, double *ky, cufftDoubleComplex *a, cufftDoubleComplex *b, cufftDoubleComplex *c, cufftDoubleComplex *d, cufftDoubleComplex *l, cufftDoubleComplex *u, cufftDoubleComplex *y, cufftDoubleComplex *pk, cufftHandle p_plan, cufftHandle f_plan, cufftHandle p_top_plan, cufftDoubleComplex *f_in, cufftDoubleComplex *f_out, cufftDoubleComplex *p_top_in, cufftDoubleComplex *p_top_out, cufftDoubleComplex *p_in, cufftDoubleComplex *p_out, Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    double gamma_rs;
    int n[] = {Nx - 1, Ny - 1};
    int idist = (Nx - 1) * (Ny - 1);
    int odist = (Nx - 1) * (Ny - 1);
    int istride = 1;
    int ostride = 1;
    int howmany = Nz - 1;
    int *inembed = NULL;
    int *onembed = NULL;
    // Allocate memory for temporary arrays in host
    cufftDoubleComplex *a_h, *b_h, *c_h, *d_h, *pk_h, *l_h, *u_h, *y_h;
    a_h = (cufftDoubleComplex *)malloc((Nz - 2) * sizeof(cufftDoubleComplex));
    b_h = (cufftDoubleComplex *)malloc((Nz - 1) * sizeof(cufftDoubleComplex));
    c_h = (cufftDoubleComplex *)malloc((Nz - 2) * sizeof(cufftDoubleComplex));
    d_h = (cufftDoubleComplex *)malloc((Nz - 1) * sizeof(cufftDoubleComplex));
    pk_h = (cufftDoubleComplex *)malloc((Nz - 1) * sizeof(cufftDoubleComplex));
    l_h = (cufftDoubleComplex *)malloc((Nz - 2) * sizeof(cufftDoubleComplex));
    u_h = (cufftDoubleComplex *)malloc((Nz - 1) * sizeof(cufftDoubleComplex));
    y_h = (cufftDoubleComplex *)malloc((Nz - 1) * sizeof(cufftDoubleComplex));

    // printf("pre fft\n");
    pre_fft<<<BLOCKS, THREADS>>>(U, f_in, p_top_in, a, c, p, parameters);
    CHECK(cudaDeviceSynchronize());
    // printf("pre fft ok\n");

    // Plan for FFT2(f) for each z slice
    CHECK_CUFFT(cufftPlanMany(&p_plan, 2, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, howmany));
    // Plan for FFT2(p_top)
    CHECK_CUFFT(cufftPlan2d(&p_top_plan, Nx - 1, Ny - 1, CUFFT_Z2Z));

    // Compute FFT
    CHECK_CUFFT(cufftExecZ2Z(p_top_plan, p_top_in, p_top_out, CUFFT_FORWARD)); // FFT2(p_top)
    CHECK_CUFFT(cufftExecZ2Z(p_plan, f_in, f_out, CUFFT_FORWARD)); // FFT2(f) for each z slice

    // Compute r,s systems of equations
    // printf("Compute r,s systems of equations\n");
    for (int r = 0; r < Nx - 1; r++) {
        for (int s = 0; s < Ny - 1; s++) {
            gamma_rs = -2 - kx[r] * kx[r] - ky[s] * ky[s];
            // printf("update values\n");
            update_values<<<1, Nz - 1>>>(b, c, d, f_out, p_top_out, parameters, r, s, gamma_rs);
            CHECK(cudaDeviceSynchronize());
            // printf("update values ok\n");
            // Solve tridiagonal system
            // printf("thomas algorithm\n");
            // thomas_algorithm<<<1, 1>>>(a, b, c, d, pk, l, u, y, Nz - 1);
            // CHECK(cudaDeviceSynchronize());
            // Copy to host
            cudaMemcpy(a_h, a, (Nz - 2) * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
            cudaMemcpy(b_h, b, (Nz - 1) * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
            cudaMemcpy(c_h, c, (Nz - 2) * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
            cudaMemcpy(d_h, d, (Nz - 1) * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
            cudaMemcpy(pk_h, pk, (Nz - 1) * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
            cudaMemcpy(l_h, l, (Nz - 2) * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
            cudaMemcpy(u_h, u, (Nz - 1) * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
            cudaMemcpy(y_h, y, (Nz - 1) * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
            thomas_algorithm_cpu(a_h, b_h, c_h, d_h, pk_h, l_h, u_h, y_h, Nz - 1);
            // Copy to device
            cudaMemcpy(a, a_h, (Nz - 2) * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
            cudaMemcpy(b, b_h, (Nz - 1) * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
            cudaMemcpy(c, c_h, (Nz - 2) * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
            cudaMemcpy(d, d_h, (Nz - 1) * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
            cudaMemcpy(pk, pk_h, (Nz - 1) * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
            cudaMemcpy(l, l_h, (Nz - 2) * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
            cudaMemcpy(u, u_h, (Nz - 1) * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
            cudaMemcpy(y, y_h, (Nz - 1) * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);            
            // printf("thomas algorithm ok\n");
            // Fill p_in with solution 
            // printf("copy to p in\n");
            copy_to_p_in<<<1, Nz - 1>>>(p_in, pk, parameters, r, s);
            CHECK(cudaDeviceSynchronize());
            // printf("copy to p in ok\n");
            // for (int k = 0; k < Nz - 1; k++) {
            //     p_in[CUFFTIDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = pk[k];
            // }
        }
    }
    

    // Plan for IFFT2(p) for each z slice
    CHECK_CUFFT(cufftPlanMany(&p_plan, 2, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, howmany));

    // Compute IFFT
    CHECK_CUFFT(cufftExecZ2Z(p_plan, p_in, p_out, CUFFT_INVERSE));

    // Post FFT
    // printf("post fft\n");
    post_fft<<<BLOCKS, THREADS>>>(p, p_out, parameters);
    CHECK(cudaDeviceSynchronize());
    // printf("post fft ok\n");

    // Destroy plans
    cufftDestroy(p_top_plan);
    cufftDestroy(p_plan);

    // Free memory
    free(a_h);
    free(b_h);
    free(c_h);
    free(d_h);
    free(pk_h);
    free(l_h);
    free(u_h);
    free(y_h);
}
