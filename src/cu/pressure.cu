/**
 * @file pressure.c
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Functions for solving the Poisson equation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "../../include/cu/pressure.cuh"

__global__
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
    }
}

__global__
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
            // gamma[FFTWIDX(r, s, k, Nx - 1, Ny - 1, 0)] = -2 - kx[r] * kx[r] - ky[s] * ky[s];
            gamma[IDX(r, s, 0, Nx - 1, Ny - 1, 1)] = -2 - kx[r] * kx[r] - ky[s] * ky[s];
            a[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 2)] = 1.0 / (dz * dz), 0.0; 
            b[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = -1.0 / dz, 0.0;
            // c[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 2)] = (2.0 + 0.5 * gamma[FFTWIDX(r, s, k, Nx - 1, Ny - 1, 0)]) / dz;
            c[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 2)] = (2.0 + 0.5 * gamma[IDX(r, s, 0, Nx - 1, Ny - 1, 1)]) / dz;
        } else { // The rest of the coefficients a and c
            a[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 2)] = 1.0 / (dz * dz); 
            c[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 2)] = 1.0 / (dz * dz);
        }
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
    double rho_inf = parameters.rho_inf;
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
            f = rho_inf * (ux + vy + wz) / dt;
            f_in[FFTWIDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)] = make_cuDoubleComplex(f, 0.0);
        }
        // Fill p_top
        if (k == Nz - 1 && j < Ny - 1 && i < Nx - 1) {
            p_top_in[FFTWIDX(i, j, 0, Nx - 1, Ny - 1, Nz - 1)] = make_cuDoubleComplex(p[IDX(i, j, k, Nx, Ny, Nz)], 0.0);
        }
    }
}

__global__
void compute_f_density(double *U, cufftDoubleComplex *f_in, cufftDoubleComplex *p_top_in, double *p, Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int u_index = parameters.field_indexes.u;
    int v_index = parameters.field_indexes.v;
    int w_index = parameters.field_indexes.w;
    int T_index = parameters.field_indexes.T;
    int im1, ip1, jm1, jp1;
    double dx = parameters.dx;
    double dy = parameters.dy;
    double dz = parameters.dz;
    double dt = parameters.dt;
    double T_inf = parameters.T_inf;
    double rho_inf = parameters.rho_inf;
    double u_ijk, u_ip1jk, u_im1jk, u_iphjk, u_imhjk;
    double v_ijk, v_ijp1k, v_ijm1k, v_ijphk, v_ijmhk;
    double w_ijk, w_ijkp1, w_ijkm1, w_ijkp2, w_ijkm2, w_ijkph, w_ijkmh;
    double T_ijk, T_im1jk, T_ip1jk, T_ijm1k, T_ijp1k, T_ijkp1, T_ijkp2, T_ijkm1, T_ijkm2;
    double p_ijk, p_im1jk, p_ip1jk, p_ijm1k, p_ijp1k, p_ijkp1, p_ijkp2, p_ijkm1, p_ijkm2;
    double rho_ijk, rho_im1jk, rho_ip1jk, rho_ijm1k, rho_ijp1k, rho_ijkp1, rho_ijkp2, rho_ijkm1, rho_ijkm2;
    double ux, vy, wz, f;
    double rho, rhox, rhoy, rhoz, px, py, pz;
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
        T_ijk = U[T_index + IDX(i, j, k, Nx, Ny, Nz)];
        p_ijk = p[IDX(i, j, k, Nx, Ny, Nz)];
        rho_ijk = T_inf * rho_inf / T_ijk;
        // Periodic boundary conditions on xy
        u_im1jk = U[u_index + IDX(im1, j, k, Nx, Ny, Nz)];
        u_ip1jk = U[u_index + IDX(ip1, j, k, Nx, Ny, Nz)];
        v_ijm1k = U[v_index + IDX(i, jm1, k, Nx, Ny, Nz)];
        v_ijp1k = U[v_index + IDX(i, jp1, k, Nx, Ny, Nz)];
        T_im1jk = U[T_index + IDX(im1, j, k, Nx, Ny, Nz)];
        T_ip1jk = U[T_index + IDX(ip1, j, k, Nx, Ny, Nz)];
        T_ijm1k = U[T_index + IDX(i, jm1, k, Nx, Ny, Nz)];
        T_ijp1k = U[T_index + IDX(i, jp1, k, Nx, Ny, Nz)];
        p_im1jk = p[IDX(im1, j, k, Nx, Ny, Nz)];
        p_ip1jk = p[IDX(ip1, j, k, Nx, Ny, Nz)];
        p_ijm1k = p[IDX(i, jm1, k, Nx, Ny, Nz)];
        p_ijp1k = p[IDX(i, jp1, k, Nx, Ny, Nz)];
        rho_im1jk = T_inf * rho_inf / T_im1jk;
        rho_ip1jk = T_inf * rho_inf / T_ip1jk;
        rho_ijm1k = T_inf * rho_inf / T_ijm1k;
        rho_ijp1k = T_inf * rho_inf / T_ijp1k;
        // dw/dz 
        if (k == 0) { // Bottom boundary                    
            w_ijkp1 = U[w_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
            w_ijkp2 = U[w_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
            T_ijkp1 = U[T_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
            T_ijkp2 = U[T_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
            p_ijkp1 = p[IDX(i, j, k + 1, Nx, Ny, Nz)];
            p_ijkp2 = p[IDX(i, j, k + 2, Nx, Ny, Nz)];
            rho_ijkp1 = T_inf * rho_inf / T_ijkp1;
            rho_ijkp2 = T_inf * rho_inf / T_ijkp2;
            wz = (-3 * w_ijk + 4 * w_ijkp1 - w_ijkp2) / (2 * dz); // dw/dz at z = z_min
            rhoz = (-3 * rho_ijk + 4 * rho_ijkp1 - rho_ijkp2) / (2 * dz);
            pz = (-3 * p_ijk + 4 * p_ijkp1 - p_ijkp2) / (2 * dz);
        } else if (k == Nz - 1) { // Top boundary
            w_ijkm1 = U[w_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
            w_ijkm2 = U[w_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
            T_ijkm1 = U[T_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
            T_ijkm2 = U[T_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
            p_ijkm1 = p[IDX(i, j, k - 1, Nx, Ny, Nz)];
            p_ijkm2 = p[IDX(i, j, k - 2, Nx, Ny, Nz)];
            rho_ijkm1 = T_inf * rho_inf / T_ijkm1;
            rho_ijkm2 = T_inf * rho_inf / T_ijkm2;
            wz = (3 * w_ijk - 4 * w_ijkm1 + w_ijkm2) / (2 * dz); // dw/dz at z = z_max
            rhoz = (3 * rho_ijk - 4 * rho_ijkm1 + rho_ijkm2) / (2 * dz);
            pz = (3 * p_ijk - 4 * p_ijkm1 + p_ijkm2) / (2 * dz);
        } else { // Interior
            w_ijkp1 = U[w_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
            w_ijkm1 = U[w_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
            T_ijkp1 = U[T_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
            T_ijkm1 = U[T_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
            rho_ijkp1 = T_inf * rho_inf / T_ijkp1;
            rho_ijkm1 = T_inf * rho_inf / T_ijkm1;
            p_ijkp1 = p[IDX(i, j, k + 1, Nx, Ny, Nz)];
            p_ijkm1 = p[IDX(i, j, k - 1, Nx, Ny, Nz)];
            w_ijkph = 0.5 * (w_ijk + w_ijkp1);
            w_ijkmh = 0.5 * (w_ijk + w_ijkm1);
            wz = (w_ijkph - w_ijkmh) / dz; // dw/dz at z = z_k
            rhoz = (rho_ijkp1 - rho_ijkm1) / (2 * dz); // drho/dz at z = z_k
            pz = (p_ijkp1 - p_ijkm1) / (2 * dz); // dp/dz at z = z_k
        }
        // Half derivatives
        u_iphjk = 0.5 * (u_ip1jk + u_ijk);
        u_imhjk = 0.5 * (u_ijk + u_im1jk);
        v_ijphk = 0.5 * (v_ijp1k + v_ijk);
        v_ijmhk = 0.5 * (v_ijk + v_ijm1k);
        ux = (u_iphjk - u_imhjk) / dx; // du/dx
        vy = (v_ijphk - v_ijmhk) / dy; // dv/dy
        // Density and pressure gradients using central differences
        rhox = (rho_ip1jk - rho_im1jk) / (2 * dx);
        rhoy = (rho_ijp1k - rho_ijm1k) / (2 * dy);
        px = (p_ip1jk - p_im1jk) / (2 * dx);
        py = (p_ijp1k - p_ijm1k) / (2 * dy);
        if (i < Nx - 1 && j < Ny - 1 && k < Nz - 1) {
            // Compute rho / dt * div(U) and store it for many DFT (contiguous z slices)
            rho = T_inf * rho_inf / T_ijk;
            f = rho * (ux + vy + wz) / dt + (rhox * px + rhoy * py + rhoz * pz) / rho;
            f_in[FFTWIDX(i, j, k, Nx - 1, Ny - 1, Nz - 1)] = make_cuDoubleComplex(f, 0.0);
        }
        // Fill p_top
        if (k == Nz - 1 && j < Ny - 1 && i < Nx - 1) {
            p_top_in[FFTWIDX(i, j, 0, Nx - 1, Ny - 1, Nz - 1)] = make_cuDoubleComplex(p[IDX(i, j, k, Nx, Ny, Nz)], 0.0);
        }
    }
}

__global__ 
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
            d[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = cuCmul(make_cuDoubleComplex(0.5 * dz, 0.0), f_out[FFTWIDX(r, s, 1, Nx - 1, Ny - 1, Nz - 1)]);
            // Last equation k=Nz-2
            d[IDX(r, s, Nz - 2, Nx - 1, Ny - 1, Nz - 1)] = cuCsub(f_out[FFTWIDX(r, s, Nz - 2, Nx - 1, Ny - 1, Nz - 1)], cuCdiv(p_top_out[FFTWIDX(r, s, 0, Nx - 1, Ny - 1, Nz - 1)], make_cuDoubleComplex(dz * dz, 0.0)));
        } else {
            if (k < Nz - 1) {
                // b[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = gamma[FFTWIDX(r, s, 0, Nx - 1, Ny - 1, 0)] / (dz * dz);
                b[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = gamma[IDX(r, s, 0, Nx - 1, Ny - 1, 1)] / (dz * dz);
                if (k < Nz - 2)
                    d[IDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)] = f_out[FFTWIDX(r, s, k, Nx - 1, Ny - 1, Nz - 1)];
            }
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
    if (parameters.variable_density == 0)
        compute_f<<<BLOCKS, THREADS>>>(U, data_in, p_top_in, p, parameters);
    else
        compute_f_density<<<BLOCKS, THREADS>>>(U, data_in, p_top_in, p, parameters);
    checkCuda(cudaGetLastError());
    
    // Plans for FFT2D
    CHECK_CUFFT(cufftPlanMany(&p_plan, 2, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, howmany)); // FFT2(f_k) for each z slice
    CHECK_CUFFT(cufftPlan2d(&p_top_plan, Nx - 1, Ny - 1, CUFFT_Z2Z)); // FFT2(p_top)

    // Compute FFT2D
    CHECK_CUFFT(cufftExecZ2Z(p_top_plan, p_top_in, p_top_out, CUFFT_FORWARD)); // FFT2(p_top)
    CHECK_CUFFT(cufftExecZ2Z(p_plan, data_in, data_out, CUFFT_FORWARD)); // FFT2D(f_k) for each z slice
    CHECK(cudaDeviceSynchronize());

    // Update coefficients, including f in pseudo-Fourier space
    update_coefficients<<<BLOCKS, THREADS>>>(gamma, b, d, data_out, p_top_out, parameters);
    checkCuda(cudaGetLastError());

    // Compute r,s systems of equations using thomas algorithm
    thomas_algorithm<<<BLOCKS, THREADS>>>(a, b, c, d, data_in, l, u, y, parameters);
    checkCuda(cudaGetLastError());

    // Compute IFFT2D
    CHECK_CUFFT(cufftExecZ2Z(p_plan, data_in, data_out, CUFFT_INVERSE));
    CHECK(cudaDeviceSynchronize());

    // Post FFT
    post_fft<<<BLOCKS, THREADS>>>(p, data_out, parameters);
    checkCuda(cudaGetLastError());

    // Destroy plans
    cufftDestroy(p_top_plan);
    cufftDestroy(p_plan);
    cufftDestroy(f_plan);
}

void solve_pressure_iterative(double *U, double *p, double *gamma, double *a, double *b, double *c, cufftDoubleComplex *d, cufftDoubleComplex *l, cufftDoubleComplex *u, cufftDoubleComplex *y, cufftDoubleComplex *data_in, cufftDoubleComplex *data_out, cufftDoubleComplex *p_top_in, cufftDoubleComplex *p_top_out, Parameters parameters, double *error, int *max_iter) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int size = Nx * Ny * Nz;
    double h_tol, *d_tol; 
    char solver_log_message[128];
    // Create a temporary array on GPU to store the pressure field
    double *p_tmp;
    int m;
    checkCuda(cudaMalloc((void **)&p_tmp, size * sizeof(double)));
    // Allocate memory for the tolerance
    checkCuda(cudaMalloc((void **)&d_tol, sizeof(double)));
    for (m = 0; m < parameters.pressure_solver_iter; m++) {
        // Copy the initial pressure field to the temporary array
        checkCuda(cudaMemcpy(p_tmp, p, size * sizeof(double), cudaMemcpyDeviceToDevice));
        solve_pressure(U, p, gamma, a, b, c, d, l, u, y, data_in, data_out, p_top_in, p_top_out, parameters);
        norm<<<BLOCKS, THREADS>>>(p, p_tmp, d_tol, INFINITY, size);
        checkCuda(cudaGetLastError());
        checkCuda(cudaMemcpy(&h_tol, d_tol, sizeof(double), cudaMemcpyDeviceToHost));
        if (h_tol <= parameters.pressure_solver_tol) {
            break;
        }        
    }
    if (parameters.pressure_solver_log == 1) {
        sprintf(solver_log_message, "Pressure solver: Error = %e, iterations = %d", h_tol, m);
        log_message(parameters, solver_log_message);
    }
    *max_iter = m;
    *error = h_tol;
    // Free memory
    checkCuda(cudaFree(p_tmp));
    checkCuda(cudaFree(d_tol));
}