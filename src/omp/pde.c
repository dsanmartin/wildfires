/**
 * @file pde.c
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Functions for solving the partial differential equations of the wildfire simulation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "../../include/omp/pde.h"

void RHS_v1(double t, double *R_old, double *R_new, double *R_turbulence, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->Nz_Y;
    double dx = parameters->dx;
    double dy = parameters->dy;
    double dz = parameters->dz;
    double *z = parameters->z;
    double nu = parameters->nu;
    double alpha = parameters->alpha;
    double Y_f = parameters->Y_f;
    double H_R = parameters->H_R;
    double A = parameters->A;
    double T_a = parameters->T_a;
    double T_pc = parameters->T_pc;
    double h = parameters->h;
    double a_v = parameters->a_v;
    double T_inf = parameters->T_inf;
    double c_p = parameters->c_p;
    double rho = parameters->rho;
    double Y_D = parameters->Y_D;
    double g = parameters->g;
    // Fields indexes
    int u_index = parameters->field_indexes.u;
    int v_index = parameters->field_indexes.v;
    int w_index = parameters->field_indexes.w;
    int T_index = parameters->field_indexes.T;
    int Y_index = parameters->field_indexes.Y;
    // Fields nodes
    double u_ijk, u_ip1jk, u_im1jk, u_ijp1k, u_ijm1k, u_ijkp1, u_ijkm1;
    double u_ip2jk, u_im2jk, u_ijp2k, u_ijm2k, u_ijkp2, u_ijkm2, u_ijkp3, u_ijkm3;
    double v_ijk, v_ip1jk, v_im1jk, v_ijp1k, v_ijm1k, v_ijkp1, v_ijkm1;
    double v_ip2jk, v_im2jk, v_ijp2k, v_ijm2k, v_ijkp2, v_ijkm2, v_ijkp3, v_ijkm3;
    double w_ijk, w_ip1jk, w_im1jk, w_ijp1k, w_ijm1k, w_ijkp1, w_ijkm1;
    double w_ip2jk, w_im2jk, w_ijp2k, w_ijm2k, w_ijkp2, w_ijkm2, w_ijkp3, w_ijkm3;
    double T_ijk, T_ip1jk, T_im1jk, T_ijp1k, T_ijm1k, T_ijkp1, T_ijkm1;
    double T_ijkp2, T_ijkm2, T_ijkp3, T_ijkm3;
    double Y_ijk;
    // Upwind scheme terms
    double u_plu, u_min, v_plu, v_min, w_plu, w_min;
    double u_ip, u_im, u_jp, u_jm, u_kp, u_km;
    double v_ip, v_im, v_jp, v_jm, v_kp, v_km;
    double w_ip, w_im, w_jp, w_jm, w_kp, w_km;
    // First partial derivatives
    double ux, uy, uz, vx, vy, vz, wx, wy, wz, Tx, Ty, Tz; // For central difference
    // Second partial derivatives
    double uxx, uyy, uzz, vxx, vyy, vzz, wxx, wyy, wzz, Txx, Tyy, Tzz;
    double uux, vuy, wuz, uvx, vvy, wvz, uwx, vwy, wwz;
    double lap_u, lap_v, lap_w, lap_T;
    double u_RHS, v_RHS, w_RHS, T_RHS, Y_RHS, S;    
    double u_tau, tau_p, fw;
    double mod_U;
    double F_x, F_y, F_z;
    int im1, ip1, jm1, jp1;
    int im2, ip2, jm2, jp2;
    int n_threads = parameters->n_threads;
    // Number of elements to process by each thread
    int chunk = Nx * Ny * Nz / n_threads;
    // Number of elements to process by the last thread
    int last_chunk = chunk + Nx * Ny * Nz % n_threads;
    // Process all nodes in parallel
    #pragma omp parallel num_threads(n_threads) 
    {
        int thread_id = omp_get_thread_num();
        int start = thread_id * chunk;
        int end = (thread_id == n_threads - 1) ? start + last_chunk : start + chunk;
        for (int ijk = start; ijk < end; ijk++) {
            int i = ijk / (Ny * Nz);
            int j = (ijk % (Ny * Nz)) / Nz;
            int k = ijk % Nz;
            // Indexes for periodic boundary conditions
            im1 = (i - 1 + Nx - 1) % (Nx - 1);
            im2 = (i - 2 + Nx - 1) % (Nx - 1);
            jm1 = (j - 1 + Ny - 1) % (Ny - 1);
            jm2 = (j - 2 + Ny - 1) % (Ny - 1);
            ip1 = (i + 1) % (Nx - 1);
            ip2 = (i + 2) % (Nx - 1);
            jp1 = (j + 1) % (Ny - 1);
            jp2 = (j + 2) % (Ny - 1);
            // Current nodes \phi_{i,j,k}
            u_ijk = R_old[u_index + IDX(i, j, k, Nx, Ny, Nz)];
            v_ijk = R_old[v_index + IDX(i, j, k, Nx, Ny, Nz)];
            w_ijk = R_old[w_index + IDX(i, j, k, Nx, Ny, Nz)];
            T_ijk = R_old[T_index + IDX(i, j, k, Nx, Ny, Nz)];
            // \phi_{i-1,j,k}
            u_im1jk = R_old[u_index + IDX(im1, j, k, Nx, Ny, Nz)];
            v_im1jk = R_old[v_index + IDX(im1, j, k, Nx, Ny, Nz)];
            w_im1jk = R_old[w_index + IDX(im1, j, k, Nx, Ny, Nz)];
            T_im1jk = R_old[T_index + IDX(im1, j, k, Nx, Ny, Nz)];
            // \phi_{i+1,j,k}
            u_ip1jk = R_old[u_index + IDX(ip1, j, k, Nx, Ny, Nz)];
            v_ip1jk = R_old[v_index + IDX(ip1, j, k, Nx, Ny, Nz)];
            w_ip1jk = R_old[w_index + IDX(ip1, j, k, Nx, Ny, Nz)];
            T_ip1jk = R_old[T_index + IDX(ip1, j, k, Nx, Ny, Nz)];
            // \phi_{i-2,j,k}
            u_im2jk = R_old[u_index + IDX(im2, j, k, Nx, Ny, Nz)];
            v_im2jk = R_old[v_index + IDX(im2, j, k, Nx, Ny, Nz)];
            w_im2jk = R_old[w_index + IDX(im2, j, k, Nx, Ny, Nz)];
            // \phi_{i+2,j,k}
            u_ip2jk = R_old[u_index + IDX(ip2, j, k, Nx, Ny, Nz)];
            v_ip2jk = R_old[v_index + IDX(ip2, j, k, Nx, Ny, Nz)];
            w_ip2jk = R_old[w_index + IDX(ip2, j, k, Nx, Ny, Nz)];
            // \phi_{i,j-1,k}
            u_ijm1k = R_old[u_index + IDX(i, jm1, k, Nx, Ny, Nz)];
            v_ijm1k = R_old[v_index + IDX(i, jm1, k, Nx, Ny, Nz)];
            w_ijm1k = R_old[w_index + IDX(i, jm1, k, Nx, Ny, Nz)];
            T_ijm1k = R_old[T_index + IDX(i, jm1, k, Nx, Ny, Nz)];
            // \phi_{i,j+1,k}
            u_ijp1k = R_old[u_index + IDX(i, jp1, k, Nx, Ny, Nz)];
            v_ijp1k = R_old[v_index + IDX(i, jp1, k, Nx, Ny, Nz)];
            w_ijp1k = R_old[w_index + IDX(i, jp1, k, Nx, Ny, Nz)];
            T_ijp1k = R_old[T_index + IDX(i, jp1, k, Nx, Ny, Nz)];
            // \phi_{i,j-2,k}
            u_ijm2k = R_old[u_index + IDX(i, jm2, k, Nx, Ny, Nz)];
            v_ijm2k = R_old[v_index + IDX(i, jm2, k, Nx, Ny, Nz)];
            w_ijm2k = R_old[w_index + IDX(i, jm2, k, Nx, Ny, Nz)];
            // \phi_{i,j+2,k}
            u_ijp2k = R_old[u_index + IDX(i, jp2, k, Nx, Ny, Nz)];
            v_ijp2k = R_old[v_index + IDX(i, jp2, k, Nx, Ny, Nz)];
            w_ijp2k = R_old[w_index + IDX(i, jp2, k, Nx, Ny, Nz)];
            // \phi_{i,j,k-1}
            if (k > 0) {
                u_ijkm1 = R_old[u_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
                v_ijkm1 = R_old[v_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
                w_ijkm1 = R_old[w_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
                T_ijkm1 = R_old[T_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
            } else {
                u_ijkm1 = 0.0;
                v_ijkm1 = 0.0;
                w_ijkm1 = 0.0;
                T_ijkm1 = 0.0;
            }
            // \phi_{i,j,k+1}
            if (k < Nz - 1) {
                u_ijkp1 = R_old[u_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
                v_ijkp1 = R_old[v_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
                w_ijkp1 = R_old[w_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
                T_ijkp1 = R_old[T_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
            } else {
                u_ijkp1 = 0.0;
                v_ijkp1 = 0.0;
                w_ijkp1 = 0.0;
                T_ijkp1 = 0.0;
            }
            // \phi_{i,j,k-2}
            if (k > 1) {
                u_ijkm2 = R_old[u_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
                v_ijkm2 = R_old[v_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
                w_ijkm2 = R_old[w_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
                T_ijkm2 = R_old[T_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
            } else {
                u_ijkm2 = 0.0;
                v_ijkm2 = 0.0;
                w_ijkm2 = 0.0;
                T_ijkm2 = 0.0;
            }
            // \phi_{i,j,k+2}
            if (k < Nz - 2) {
                u_ijkp2 = R_old[u_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
                v_ijkp2 = R_old[v_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
                w_ijkp2 = R_old[w_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
                T_ijkp2 = R_old[T_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
            } else {
                u_ijkp2 = 0.0;
                v_ijkp2 = 0.0;
                w_ijkp2 = 0.0;
                T_ijkp2 = 0.0;
            }
            // \phi_{i,j,k-3}
            if (k > 2) {
                u_ijkm3 = R_old[u_index + IDX(i, j, k - 3, Nx, Ny, Nz)];
                v_ijkm3 = R_old[v_index + IDX(i, j, k - 3, Nx, Ny, Nz)];
                w_ijkm3 = R_old[w_index + IDX(i, j, k - 3, Nx, Ny, Nz)];
                T_ijkm3 = R_old[T_index + IDX(i, j, k - 3, Nx, Ny, Nz)];
            } else {
                u_ijkm3 = 0.0;
                v_ijkm3 = 0.0;
                w_ijkm3 = 0.0;
                T_ijkm3 = 0.0;
            }
            // \phi_{i,j,k+3}
            if (k < Nz - 3) {
                u_ijkp3 = R_old[u_index + IDX(i, j, k + 3, Nx, Ny, Nz)];
                v_ijkp3 = R_old[v_index + IDX(i, j, k + 3, Nx, Ny, Nz)];
                w_ijkp3 = R_old[w_index + IDX(i, j, k + 3, Nx, Ny, Nz)];
                T_ijkp3 = R_old[T_index + IDX(i, j, k + 3, Nx, Ny, Nz)];
            } else {
                u_ijkp3 = 0.0;
                v_ijkp3 = 0.0;
                w_ijkp3 = 0.0;
                T_ijkp3 = 0.0;
            }
            /* Computing upwind scheme terms */
            u_plu = MAX(u_ijk, 0.0);
            u_min = MIN(u_ijk, 0.0);
            v_plu = MAX(v_ijk, 0.0);
            v_min = MIN(v_ijk, 0.0);
            w_plu = MAX(w_ijk, 0.0);
            w_min = MIN(w_ijk, 0.0);
            u_im = (3 * u_ijk - 4 * u_im1jk + u_im2jk) / (2 * dx);
            v_im = (3 * v_ijk - 4 * v_im1jk + v_im2jk) / (2 * dx);
            w_im = (3 * w_ijk - 4 * w_im1jk + w_im2jk) / (2 * dx);
            u_ip = (-3 * u_ijk + 4 * u_ip1jk - u_ip2jk) / (2 * dx);
            v_ip = (-3 * v_ijk + 4 * v_ip1jk - v_ip2jk) / (2 * dx);
            w_ip = (-3 * w_ijk + 4 * w_ip1jk - w_ip2jk) / (2 * dx);
            u_jm = (3 * u_ijk - 4 * u_ijm1k + u_ijm2k) / (2 * dy);
            v_jm = (3 * v_ijk - 4 * v_ijm1k + v_ijm2k) / (2 * dy);
            w_jm = (3 * w_ijk - 4 * w_ijm1k + w_ijm2k) / (2 * dy);
            u_jp = (-3 * u_ijk + 4 * u_ijp1k - u_ijp2k) / (2 * dy);
            v_jp = (-3 * v_ijk + 4 * v_ijp1k - v_ijp2k) / (2 * dy);
            w_jp = (-3 * w_ijk + 4 * w_ijp1k - w_ijp2k) / (2 * dy);
            // Second order forward difference at k=0 and k=1
            if (k <= 1) {
                u_km = (-3 * u_ijk + 4 * u_ijkp1 - u_ijkp2) / (2 * dz);
                v_km = (-3 * v_ijk + 4 * v_ijkp1 - v_ijkp2) / (2 * dz);
                w_km = (-3 * w_ijk + 4 * w_ijkp1 - w_ijkp2) / (2 * dz);
            } else {
                u_km = (3 * u_ijk - 4 * u_ijkm1 + u_ijkm2) / (2 * dz);
                v_km = (3 * v_ijk - 4 * v_ijkm1 + v_ijkm2) / (2 * dz);
                w_km = (3 * w_ijk - 4 * w_ijkm1 + w_ijkm2) / (2 * dz);
            }
            // Second order backward difference at k=Nz-2 and k=Nz-1
            if (k >= Nz - 2) {
                u_kp = (3 * u_ijk - 4 * u_ijkm1 + u_ijkm2) / (2 * dz);
                v_kp = (3 * v_ijk - 4 * v_ijkm1 + v_ijkm2) / (2 * dz);
                w_kp = (3 * w_ijk - 4 * w_ijkm1 + w_ijkm2) / (2 * dz);
            } else {
                u_kp = (-3 * u_ijk + 4 * u_ijkp1 - u_ijkp2) / (2 * dz);
                v_kp = (-3 * v_ijk + 4 * v_ijkp1 - v_ijkp2) / (2 * dz);
                w_kp = (-3 * w_ijk + 4 * w_ijkp1 - w_ijkp2) / (2 * dz);
            }
            /* Compute first partial derivatives */
            // Upwind scheme for velocity
            uux= u_plu * u_im + u_min * u_ip; // u * du/dx
            vuy = v_plu * u_jm + v_min * u_jp; // v * du/dy
            wuz = w_plu * u_km + w_min * u_kp; // w * du/dz
            uvx = u_plu * v_im + u_min * v_ip; // u * dv/dx
            vvy = v_plu * v_jm + v_min * v_jp; // v * dv/dy
            wvz = w_plu * v_km + w_min * v_kp; // w * dv/dz
            uwx = u_plu * w_im + u_min * w_ip; // u * dw/dx
            vwy = v_plu * w_jm + v_min * w_jp; // v * dw/dy
            wwz = w_plu * w_km + w_min * w_kp; // w * dw/dz
            // Central difference for temperature and turbulence velocity
            ux  = (u_ip1jk - u_im1jk) / (2 * dx); // du/dx
            uy  = (u_ijp1k - u_ijm1k) / (2 * dy); // du/dy                
            vx  = (v_ip1jk - v_im1jk) / (2 * dx); // dv/dx
            vy  = (v_ijp1k - v_ijm1k) / (2 * dy); // dv/dy                
            wx  = (w_ip1jk - w_im1jk) / (2 * dx); // dw/dx
            wy  = (w_ijp1k - w_ijm1k) / (2 * dy); // dw/dy                
            Tx  = (T_ip1jk - T_im1jk) / (2.0 * dx); // dT/dx
            Ty  = (T_ijp1k - T_ijm1k) / (2.0 * dy); // dT/dy                
            /* Compute second partial derivatives */
            // Central difference for velocity and temperature
            uxx = (u_ip1jk - 2.0 * u_ijk + u_im1jk) / (dx * dx); // d^2u/dx^2
            uyy = (u_ijp1k - 2.0 * u_ijk + u_ijm1k) / (dy * dy); // d^2u/dy^2
            vxx = (v_ip1jk - 2.0 * v_ijk + v_im1jk) / (dx * dx); // d^2v/dx^2
            vyy = (v_ijp1k - 2.0 * v_ijk + v_ijm1k) / (dy * dy); // d^2v/dy^2
            wxx = (w_ip1jk - 2.0 * w_ijk + w_im1jk) / (dx * dx); // d^2w/dx^2
            wyy = (w_ijp1k - 2.0 * w_ijk + w_ijm1k) / (dy * dy); // d^2w/dy^2
            Txx = (T_ip1jk - 2.0 * T_ijk + T_im1jk) / (dx * dx); // d^2T/dx^2 
            Tyy = (T_ijp1k - 2.0 * T_ijk + T_ijm1k) / (dy * dy); // d^2T/dy^2
            if (k == 0) { // Forward difference
                uz = (-3 * u_ijk + 4 * u_ijkp1 - u_ijkp2) / (2 * dz); // du/dz
                vz = (-3 * v_ijk + 4 * v_ijkp1 - v_ijkp2) / (2 * dz); // dv/dz
                wz = (-3 * w_ijk + 4 * w_ijkp1 - w_ijkp2) / (2 * dz); // dw/dz
                Tz = (-3.0 * T_ijk + 4.0 * T_ijkp1 - T_ijkp2) / (2.0 * dz); // dT/dz
                uzz = (2 * u_ijk - 5 * u_ijkp1 + 4 * u_ijkp2 - u_ijkp3) / (dz * dz); // d^2u/dz^2
                vzz = (2 * v_ijk - 5 * v_ijkp1 + 4 * v_ijkp2 - v_ijkp3) / (dz * dz); // d^2v/dz^2
                wzz = (2 * w_ijk - 5 * w_ijkp1 + 4 * w_ijkp2 - w_ijkp3) / (dz * dz); // d^2w/dz^2
                Tzz = (2 * T_ijk - 5 * T_ijkp1 + 4 * T_ijkp2 - T_ijkp3) / (dz * dz); // d^2T/dz^2
            } else if (k == Nz - 1) {
                uz = (3 * u_ijk - 4 * u_ijkm1 + u_ijkm2) / (2 * dz); // du/dz
                vz = (3 * v_ijk - 4 * v_ijkm1 + v_ijkm2) / (2 * dz); // dv/dz
                wz = (3 * w_ijk - 4 * w_ijkm1 + w_ijkm2) / (2 * dz); // dw/dz
                Tz = (3.0 * T_ijk - 4.0 * T_ijkm1 + T_ijkm2) / (2.0 * dz); // dT/dz
                uzz = (2 * u_ijk - 5 * u_ijkm1 + 4 * u_ijkm2 - u_ijkm3) / (dz * dz); // d^2u/dz^2
                vzz = (2 * v_ijk - 5 * v_ijkm1 + 4 * v_ijkm2 - v_ijkm3) / (dz * dz); // d^2v/dz^2
                wzz = (2 * w_ijk - 5 * w_ijkm1 + 4 * w_ijkm2 - w_ijkm3) / (dz * dz); // d^2w/dz^2
                Tzz = (2 * T_ijk - 5 * T_ijkm1 + 4 * T_ijkm2 - T_ijkm3) / (dz * dz); // d^2T/dz^2
            } else {
                uz  = (u_ijkp1 - u_ijkm1) / (2 * dz); // du/dz
                vz  = (v_ijkp1 - v_ijkm1) / (2 * dz); // dv/dz
                wz  = (w_ijkp1 - w_ijkm1) / (2 * dz); // dw/dz
                Tz  = (T_ijkp1 - T_ijkm1) / (2.0 * dz); // dT/dz
                uzz = (u_ijkp1 - 2.0 * u_ijk + u_ijkm1) / (dz * dz); // d^2u/dz^2
                vzz = (v_ijkp1 - 2.0 * v_ijk + v_ijkm1) / (dz * dz); // d^2v/dz^2
                wzz = (w_ijkp1 - 2.0 * w_ijk + w_ijkm1) / (dz * dz); // d^2w/dz^2
                Tzz = (T_ijkp1 - 2.0 * T_ijk + T_ijkm1) / (dz * dz); // d^2T/dz^2
            }
            // Damping function
            if (k == 0) {
                tau_p = sqrt((nu * 0.5 * (uz + wx)) * (nu * 0.5 * (uz + wx)) + (nu * 0.5 * (vz + wy)) * (nu * 0.5 * (vz + wy)));
            } else {
                tau_p = 0.0;
            }
            u_tau = sqrt(tau_p);
            fw = f_damping(z[k], u_tau, nu);
            // Copy to turbulence array
            R_turbulence[parameters->turbulence_indexes.ux  + IDX(i, j, k, Nx, Ny, Nz)] = ux;
            R_turbulence[parameters->turbulence_indexes.uy  + IDX(i, j, k, Nx, Ny, Nz)] = uy;
            R_turbulence[parameters->turbulence_indexes.uz  + IDX(i, j, k, Nx, Ny, Nz)] = uz;
            R_turbulence[parameters->turbulence_indexes.vx  + IDX(i, j, k, Nx, Ny, Nz)] = vx;
            R_turbulence[parameters->turbulence_indexes.vy  + IDX(i, j, k, Nx, Ny, Nz)] = vy;
            R_turbulence[parameters->turbulence_indexes.vz  + IDX(i, j, k, Nx, Ny, Nz)] = vz;
            R_turbulence[parameters->turbulence_indexes.wx  + IDX(i, j, k, Nx, Ny, Nz)] = wx;
            R_turbulence[parameters->turbulence_indexes.wy  + IDX(i, j, k, Nx, Ny, Nz)] = wy;
            R_turbulence[parameters->turbulence_indexes.wz  + IDX(i, j, k, Nx, Ny, Nz)] = wz;
            R_turbulence[parameters->turbulence_indexes.Tx  + IDX(i, j, k, Nx, Ny, Nz)] = Tx;
            R_turbulence[parameters->turbulence_indexes.Ty  + IDX(i, j, k, Nx, Ny, Nz)] = Ty;
            R_turbulence[parameters->turbulence_indexes.Tz  + IDX(i, j, k, Nx, Ny, Nz)] = Tz;
            R_turbulence[parameters->turbulence_indexes.uxx + IDX(i, j, k, Nx, Ny, Nz)] = uxx;
            R_turbulence[parameters->turbulence_indexes.uyy + IDX(i, j, k, Nx, Ny, Nz)] = uyy;
            R_turbulence[parameters->turbulence_indexes.uzz + IDX(i, j, k, Nx, Ny, Nz)] = uzz;
            R_turbulence[parameters->turbulence_indexes.vxx + IDX(i, j, k, Nx, Ny, Nz)] = vxx;
            R_turbulence[parameters->turbulence_indexes.vyy + IDX(i, j, k, Nx, Ny, Nz)] = vyy;
            R_turbulence[parameters->turbulence_indexes.vzz + IDX(i, j, k, Nx, Ny, Nz)] = vzz;
            R_turbulence[parameters->turbulence_indexes.wxx + IDX(i, j, k, Nx, Ny, Nz)] = wxx;
            R_turbulence[parameters->turbulence_indexes.wyy + IDX(i, j, k, Nx, Ny, Nz)] = wyy;
            R_turbulence[parameters->turbulence_indexes.wzz + IDX(i, j, k, Nx, Ny, Nz)] = wzz;
            R_turbulence[parameters->turbulence_indexes.Txx + IDX(i, j, k, Nx, Ny, Nz)] = Txx;
            R_turbulence[parameters->turbulence_indexes.Tyy + IDX(i, j, k, Nx, Ny, Nz)] = Tyy;
            R_turbulence[parameters->turbulence_indexes.Tzz + IDX(i, j, k, Nx, Ny, Nz)] = Tzz;
            R_turbulence[parameters->turbulence_indexes.fw  + IDX(i, j, k, Nx, Ny, Nz)] = fw;
            /* Compute fuel and source term */
            if (k < Nz_Y) {
                Y_ijk = R_old[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)];
                Y_RHS = -Y_f * K(T_ijk, A, T_a) * H(T_ijk, T_pc) * Y_ijk;
                R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] = Y_RHS;
            } else {
                Y_ijk = 0.0;
            }
            // Compute source and force terms
            S = source(T_ijk, Y_ijk, H_R, A, T_a, h, a_v, T_inf, c_p, rho, T_pc);
            mod_U = sqrt(u_ijk * u_ijk + v_ijk * v_ijk + w_ijk * w_ijk);
            // Force terms
            F_x = - Y_D * a_v * Y_ijk * mod_U * u_ijk;
            F_y = - Y_D * a_v * Y_ijk * mod_U * v_ijk;                             
            F_z = - g * (T_ijk - T_inf) / T_ijk - Y_D * a_v * Y_ijk * mod_U * w_ijk;
            // Compute Laplacian terms
            lap_u = uxx + uyy + uzz;
            lap_v = vxx + vyy + vzz;
            lap_w = wxx + wyy + wzz;
            lap_T = Txx + Tyy + Tzz;
            // Compute RHS
            u_RHS = nu * lap_u - (uux + vuy + wuz) + F_x;
            v_RHS = nu * lap_v - (uvx + vvy + wvz) + F_y;
            w_RHS = nu * lap_w - (uwx + vwy + wwz) + F_z;
            T_RHS = alpha * lap_T - (u_ijk * Tx + v_ijk * Ty + w_ijk * Tz) + S;
            // Save RHS into R_new
            R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] = u_RHS;
            R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] = v_RHS;
            R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] = w_RHS;
            R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = T_RHS;  
        }
    }
}

void RHS(double t, double *R_old, double *R_new, double *R_turbulence, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->Nz_Y;
    double dx = parameters->dx;
    double dy = parameters->dy;
    double dz = parameters->dz;
    double *z = parameters->z;
    double nu = parameters->nu;
    double alpha = parameters->alpha;
    double Y_f = parameters->Y_f;
    double H_R = parameters->H_R;
    double A = parameters->A;
    double T_a = parameters->T_a;
    double T_pc = parameters->T_pc;
    double h = parameters->h;
    double a_v = parameters->a_v;
    double T_inf = parameters->T_inf;
    double c_p = parameters->c_p;
    double rho = parameters->rho;
    double Y_D = parameters->Y_D;
    double g = parameters->g;
    // Fields indexes
    int u_index = parameters->field_indexes.u;
    int v_index = parameters->field_indexes.v;
    int w_index = parameters->field_indexes.w;
    int T_index = parameters->field_indexes.T;
    int Y_index = parameters->field_indexes.Y;
    // Fields nodes
    double u_ijk, u_ip1jk, u_im1jk, u_ijp1k, u_ijm1k, u_ijkp1, u_ijkm1;
    double u_ip2jk, u_im2jk, u_ijp2k, u_ijm2k, u_ijkp2, u_ijkm2, u_ijkp3, u_ijkm3;
    double v_ijk, v_ip1jk, v_im1jk, v_ijp1k, v_ijm1k, v_ijkp1, v_ijkm1;
    double v_ip2jk, v_im2jk, v_ijp2k, v_ijm2k, v_ijkp2, v_ijkm2, v_ijkp3, v_ijkm3;
    double w_ijk, w_ip1jk, w_im1jk, w_ijp1k, w_ijm1k, w_ijkp1, w_ijkm1;
    double w_ip2jk, w_im2jk, w_ijp2k, w_ijm2k, w_ijkp2, w_ijkm2, w_ijkp3, w_ijkm3;
    double T_ijk, T_ip1jk, T_im1jk, T_ijp1k, T_ijm1k, T_ijkp1, T_ijkm1;
    double T_ijkp2, T_ijkm2, T_ijkp3, T_ijkm3;
    double Y_ijk;
    // Upwind scheme terms
    double u_plu, u_min, v_plu, v_min, w_plu, w_min;
    double u_ip, u_im, u_jp, u_jm, u_kp, u_km;
    double v_ip, v_im, v_jp, v_jm, v_kp, v_km;
    double w_ip, w_im, w_jp, w_jm, w_kp, w_km;
    // First partial derivatives
    double ux, uy, uz, vx, vy, vz, wx, wy, wz, Tx, Ty, Tz; // For central difference
    // Second partial derivatives
    double uxx, uyy, uzz, vxx, vyy, vzz, wxx, wyy, wzz, Txx, Tyy, Tzz;
    double uux, vuy, wuz, uvx, vvy, wvz, uwx, vwy, wwz;
    double lap_u, lap_v, lap_w, lap_T;
    double u_RHS, v_RHS, w_RHS, T_RHS, Y_RHS, S;    
    double u_tau, tau_p, fw;
    double mod_U;
    double F_x, F_y, F_z;
    int im1, ip1, jm1, jp1;
    int im2, ip2, jm2, jp2;
    // int n_threads = parameters->n_threads;
    // Process all nodes in parallel
    #pragma omp parallel for collapse(3) schedule(dynamic) 
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                // Indexes for periodic boundary conditions
                im1 = (i - 1 + Nx - 1) % (Nx - 1);
                im2 = (i - 2 + Nx - 1) % (Nx - 1);
                jm1 = (j - 1 + Ny - 1) % (Ny - 1);
                jm2 = (j - 2 + Ny - 1) % (Ny - 1);
                ip1 = (i + 1) % (Nx - 1);
                ip2 = (i + 2) % (Nx - 1);
                jp1 = (j + 1) % (Ny - 1);
                jp2 = (j + 2) % (Ny - 1);
                // Current nodes \phi_{i,j,k}
                u_ijk = R_old[u_index + IDX(i, j, k, Nx, Ny, Nz)];
                v_ijk = R_old[v_index + IDX(i, j, k, Nx, Ny, Nz)];
                w_ijk = R_old[w_index + IDX(i, j, k, Nx, Ny, Nz)];
                T_ijk = R_old[T_index + IDX(i, j, k, Nx, Ny, Nz)];
                // \phi_{i-1,j,k}
                u_im1jk = R_old[u_index + IDX(im1, j, k, Nx, Ny, Nz)];
                v_im1jk = R_old[v_index + IDX(im1, j, k, Nx, Ny, Nz)];
                w_im1jk = R_old[w_index + IDX(im1, j, k, Nx, Ny, Nz)];
                T_im1jk = R_old[T_index + IDX(im1, j, k, Nx, Ny, Nz)];
                // \phi_{i+1,j,k}
                u_ip1jk = R_old[u_index + IDX(ip1, j, k, Nx, Ny, Nz)];
                v_ip1jk = R_old[v_index + IDX(ip1, j, k, Nx, Ny, Nz)];
                w_ip1jk = R_old[w_index + IDX(ip1, j, k, Nx, Ny, Nz)];
                T_ip1jk = R_old[T_index + IDX(ip1, j, k, Nx, Ny, Nz)];
                // \phi_{i-2,j,k}
                u_im2jk = R_old[u_index + IDX(im2, j, k, Nx, Ny, Nz)];
                v_im2jk = R_old[v_index + IDX(im2, j, k, Nx, Ny, Nz)];
                w_im2jk = R_old[w_index + IDX(im2, j, k, Nx, Ny, Nz)];
                // \phi_{i+2,j,k}
                u_ip2jk = R_old[u_index + IDX(ip2, j, k, Nx, Ny, Nz)];
                v_ip2jk = R_old[v_index + IDX(ip2, j, k, Nx, Ny, Nz)];
                w_ip2jk = R_old[w_index + IDX(ip2, j, k, Nx, Ny, Nz)];
                // \phi_{i,j-1,k}
                u_ijm1k = R_old[u_index + IDX(i, jm1, k, Nx, Ny, Nz)];
                v_ijm1k = R_old[v_index + IDX(i, jm1, k, Nx, Ny, Nz)];
                w_ijm1k = R_old[w_index + IDX(i, jm1, k, Nx, Ny, Nz)];
                T_ijm1k = R_old[T_index + IDX(i, jm1, k, Nx, Ny, Nz)];
                // \phi_{i,j+1,k}
                u_ijp1k = R_old[u_index + IDX(i, jp1, k, Nx, Ny, Nz)];
                v_ijp1k = R_old[v_index + IDX(i, jp1, k, Nx, Ny, Nz)];
                w_ijp1k = R_old[w_index + IDX(i, jp1, k, Nx, Ny, Nz)];
                T_ijp1k = R_old[T_index + IDX(i, jp1, k, Nx, Ny, Nz)];
                // \phi_{i,j-2,k}
                u_ijm2k = R_old[u_index + IDX(i, jm2, k, Nx, Ny, Nz)];
                v_ijm2k = R_old[v_index + IDX(i, jm2, k, Nx, Ny, Nz)];
                w_ijm2k = R_old[w_index + IDX(i, jm2, k, Nx, Ny, Nz)];
                // \phi_{i,j+2,k}
                u_ijp2k = R_old[u_index + IDX(i, jp2, k, Nx, Ny, Nz)];
                v_ijp2k = R_old[v_index + IDX(i, jp2, k, Nx, Ny, Nz)];
                w_ijp2k = R_old[w_index + IDX(i, jp2, k, Nx, Ny, Nz)];
                // \phi_{i,j,k-1}
                if (k > 0) {
                    u_ijkm1 = R_old[u_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
                    v_ijkm1 = R_old[v_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
                    w_ijkm1 = R_old[w_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
                    T_ijkm1 = R_old[T_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
                } else {
                    u_ijkm1 = 0.0;
                    v_ijkm1 = 0.0;
                    w_ijkm1 = 0.0;
                    T_ijkm1 = 0.0;
                }
                // \phi_{i,j,k+1}
                if (k < Nz - 1) {
                    u_ijkp1 = R_old[u_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
                    v_ijkp1 = R_old[v_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
                    w_ijkp1 = R_old[w_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
                    T_ijkp1 = R_old[T_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
                } else {
                    u_ijkp1 = 0.0;
                    v_ijkp1 = 0.0;
                    w_ijkp1 = 0.0;
                    T_ijkp1 = 0.0;
                }
                // \phi_{i,j,k-2}
                if (k > 1) {
                    u_ijkm2 = R_old[u_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
                    v_ijkm2 = R_old[v_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
                    w_ijkm2 = R_old[w_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
                    T_ijkm2 = R_old[T_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
                } else {
                    u_ijkm2 = 0.0;
                    v_ijkm2 = 0.0;
                    w_ijkm2 = 0.0;
                    T_ijkm2 = 0.0;
                }
                // \phi_{i,j,k+2}
                if (k < Nz - 2) {
                    u_ijkp2 = R_old[u_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
                    v_ijkp2 = R_old[v_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
                    w_ijkp2 = R_old[w_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
                    T_ijkp2 = R_old[T_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
                } else {
                    u_ijkp2 = 0.0;
                    v_ijkp2 = 0.0;
                    w_ijkp2 = 0.0;
                    T_ijkp2 = 0.0;
                }
                // \phi_{i,j,k-3}
                if (k > 2) {
                    u_ijkm3 = R_old[u_index + IDX(i, j, k - 3, Nx, Ny, Nz)];
                    v_ijkm3 = R_old[v_index + IDX(i, j, k - 3, Nx, Ny, Nz)];
                    w_ijkm3 = R_old[w_index + IDX(i, j, k - 3, Nx, Ny, Nz)];
                    T_ijkm3 = R_old[T_index + IDX(i, j, k - 3, Nx, Ny, Nz)];
                } else {
                    u_ijkm3 = 0.0;
                    v_ijkm3 = 0.0;
                    w_ijkm3 = 0.0;
                    T_ijkm3 = 0.0;
                }
                // \phi_{i,j,k+3}
                if (k < Nz - 3) {
                    u_ijkp3 = R_old[u_index + IDX(i, j, k + 3, Nx, Ny, Nz)];
                    v_ijkp3 = R_old[v_index + IDX(i, j, k + 3, Nx, Ny, Nz)];
                    w_ijkp3 = R_old[w_index + IDX(i, j, k + 3, Nx, Ny, Nz)];
                    T_ijkp3 = R_old[T_index + IDX(i, j, k + 3, Nx, Ny, Nz)];
                } else {
                    u_ijkp3 = 0.0;
                    v_ijkp3 = 0.0;
                    w_ijkp3 = 0.0;
                    T_ijkp3 = 0.0;
                }
                /* Computing upwind scheme terms */
                u_plu = MAX(u_ijk, 0.0);
                u_min = MIN(u_ijk, 0.0);
                v_plu = MAX(v_ijk, 0.0);
                v_min = MIN(v_ijk, 0.0);
                w_plu = MAX(w_ijk, 0.0);
                w_min = MIN(w_ijk, 0.0);
                u_im = (3 * u_ijk - 4 * u_im1jk + u_im2jk) / (2 * dx);
                v_im = (3 * v_ijk - 4 * v_im1jk + v_im2jk) / (2 * dx);
                w_im = (3 * w_ijk - 4 * w_im1jk + w_im2jk) / (2 * dx);
                u_ip = (-3 * u_ijk + 4 * u_ip1jk - u_ip2jk) / (2 * dx);
                v_ip = (-3 * v_ijk + 4 * v_ip1jk - v_ip2jk) / (2 * dx);
                w_ip = (-3 * w_ijk + 4 * w_ip1jk - w_ip2jk) / (2 * dx);
                u_jm = (3 * u_ijk - 4 * u_ijm1k + u_ijm2k) / (2 * dy);
                v_jm = (3 * v_ijk - 4 * v_ijm1k + v_ijm2k) / (2 * dy);
                w_jm = (3 * w_ijk - 4 * w_ijm1k + w_ijm2k) / (2 * dy);
                u_jp = (-3 * u_ijk + 4 * u_ijp1k - u_ijp2k) / (2 * dy);
                v_jp = (-3 * v_ijk + 4 * v_ijp1k - v_ijp2k) / (2 * dy);
                w_jp = (-3 * w_ijk + 4 * w_ijp1k - w_ijp2k) / (2 * dy);
                // Second order forward difference at k=0 and k=1
                if (k <= 1) {
                    u_km = (-3 * u_ijk + 4 * u_ijkp1 - u_ijkp2) / (2 * dz);
                    v_km = (-3 * v_ijk + 4 * v_ijkp1 - v_ijkp2) / (2 * dz);
                    w_km = (-3 * w_ijk + 4 * w_ijkp1 - w_ijkp2) / (2 * dz);
                } else {
                    u_km = (3 * u_ijk - 4 * u_ijkm1 + u_ijkm2) / (2 * dz);
                    v_km = (3 * v_ijk - 4 * v_ijkm1 + v_ijkm2) / (2 * dz);
                    w_km = (3 * w_ijk - 4 * w_ijkm1 + w_ijkm2) / (2 * dz);
                }
                // Second order backward difference at k=Nz-2 and k=Nz-1
                if (k >= Nz - 2) {
                    u_kp = (3 * u_ijk - 4 * u_ijkm1 + u_ijkm2) / (2 * dz);
                    v_kp = (3 * v_ijk - 4 * v_ijkm1 + v_ijkm2) / (2 * dz);
                    w_kp = (3 * w_ijk - 4 * w_ijkm1 + w_ijkm2) / (2 * dz);
                } else {
                    u_kp = (-3 * u_ijk + 4 * u_ijkp1 - u_ijkp2) / (2 * dz);
                    v_kp = (-3 * v_ijk + 4 * v_ijkp1 - v_ijkp2) / (2 * dz);
                    w_kp = (-3 * w_ijk + 4 * w_ijkp1 - w_ijkp2) / (2 * dz);
                }
                /* Compute first partial derivatives */
                // Upwind scheme for velocity
                uux = u_plu * u_im + u_min * u_ip; // u * du/dx
                vuy = v_plu * u_jm + v_min * u_jp; // v * du/dy
                wuz = w_plu * u_km + w_min * u_kp; // w * du/dz
                uvx = u_plu * v_im + u_min * v_ip; // u * dv/dx
                vvy = v_plu * v_jm + v_min * v_jp; // v * dv/dy
                wvz = w_plu * v_km + w_min * v_kp; // w * dv/dz
                uwx = u_plu * w_im + u_min * w_ip; // u * dw/dx
                vwy = v_plu * w_jm + v_min * w_jp; // v * dw/dy
                wwz = w_plu * w_km + w_min * w_kp; // w * dw/dz
                // Central difference for temperature and turbulence velocity
                ux  = (u_ip1jk - u_im1jk) / (2 * dx); // du/dx
                uy  = (u_ijp1k - u_ijm1k) / (2 * dy); // du/dy                
                vx  = (v_ip1jk - v_im1jk) / (2 * dx); // dv/dx
                vy  = (v_ijp1k - v_ijm1k) / (2 * dy); // dv/dy                
                wx  = (w_ip1jk - w_im1jk) / (2 * dx); // dw/dx
                wy  = (w_ijp1k - w_ijm1k) / (2 * dy); // dw/dy                
                Tx  = (T_ip1jk - T_im1jk) / (2.0 * dx); // dT/dx
                Ty  = (T_ijp1k - T_ijm1k) / (2.0 * dy); // dT/dy                
                /* Compute second partial derivatives */
                // Central difference for velocity and temperature
                uxx = (u_ip1jk - 2.0 * u_ijk + u_im1jk) / (dx * dx); // d^2u/dx^2
                uyy = (u_ijp1k - 2.0 * u_ijk + u_ijm1k) / (dy * dy); // d^2u/dy^2
                vxx = (v_ip1jk - 2.0 * v_ijk + v_im1jk) / (dx * dx); // d^2v/dx^2
                vyy = (v_ijp1k - 2.0 * v_ijk + v_ijm1k) / (dy * dy); // d^2v/dy^2
                wxx = (w_ip1jk - 2.0 * w_ijk + w_im1jk) / (dx * dx); // d^2w/dx^2
                wyy = (w_ijp1k - 2.0 * w_ijk + w_ijm1k) / (dy * dy); // d^2w/dy^2
                Txx = (T_ip1jk - 2.0 * T_ijk + T_im1jk) / (dx * dx); // d^2T/dx^2 
                Tyy = (T_ijp1k - 2.0 * T_ijk + T_ijm1k) / (dy * dy); // d^2T/dy^2
                if (k == 0) { // Forward difference
                    uz = (-3 * u_ijk + 4 * u_ijkp1 - u_ijkp2) / (2 * dz); // du/dz
                    vz = (-3 * v_ijk + 4 * v_ijkp1 - v_ijkp2) / (2 * dz); // dv/dz
                    wz = (-3 * w_ijk + 4 * w_ijkp1 - w_ijkp2) / (2 * dz); // dw/dz
                    Tz = (-3.0 * T_ijk + 4.0 * T_ijkp1 - T_ijkp2) / (2.0 * dz); // dT/dz
                    uzz = (2 * u_ijk - 5 * u_ijkp1 + 4 * u_ijkp2 - u_ijkp3) / (dz * dz); // d^2u/dz^2
                    vzz = (2 * v_ijk - 5 * v_ijkp1 + 4 * v_ijkp2 - v_ijkp3) / (dz * dz); // d^2v/dz^2
                    wzz = (2 * w_ijk - 5 * w_ijkp1 + 4 * w_ijkp2 - w_ijkp3) / (dz * dz); // d^2w/dz^2
                    Tzz = (2 * T_ijk - 5 * T_ijkp1 + 4 * T_ijkp2 - T_ijkp3) / (dz * dz); // d^2T/dz^2
                } else if (k == Nz - 1) {
                    uz = (3 * u_ijk - 4 * u_ijkm1 + u_ijkm2) / (2 * dz); // du/dz
                    vz = (3 * v_ijk - 4 * v_ijkm1 + v_ijkm2) / (2 * dz); // dv/dz
                    wz = (3 * w_ijk - 4 * w_ijkm1 + w_ijkm2) / (2 * dz); // dw/dz
                    Tz = (3.0 * T_ijk - 4.0 * T_ijkm1 + T_ijkm2) / (2.0 * dz); // dT/dz
                    uzz = (2 * u_ijk - 5 * u_ijkm1 + 4 * u_ijkm2 - u_ijkm3) / (dz * dz); // d^2u/dz^2
                    vzz = (2 * v_ijk - 5 * v_ijkm1 + 4 * v_ijkm2 - v_ijkm3) / (dz * dz); // d^2v/dz^2
                    wzz = (2 * w_ijk - 5 * w_ijkm1 + 4 * w_ijkm2 - w_ijkm3) / (dz * dz); // d^2w/dz^2
                    Tzz = (2 * T_ijk - 5 * T_ijkm1 + 4 * T_ijkm2 - T_ijkm3) / (dz * dz); // d^2T/dz^2
                } else {
                    uz  = (u_ijkp1 - u_ijkm1) / (2 * dz); // du/dz
                    vz  = (v_ijkp1 - v_ijkm1) / (2 * dz); // dv/dz
                    wz  = (w_ijkp1 - w_ijkm1) / (2 * dz); // dw/dz
                    Tz  = (T_ijkp1 - T_ijkm1) / (2.0 * dz); // dT/dz
                    uzz = (u_ijkp1 - 2.0 * u_ijk + u_ijkm1) / (dz * dz); // d^2u/dz^2
                    vzz = (v_ijkp1 - 2.0 * v_ijk + v_ijkm1) / (dz * dz); // d^2v/dz^2
                    wzz = (w_ijkp1 - 2.0 * w_ijk + w_ijkm1) / (dz * dz); // d^2w/dz^2
                    Tzz = (T_ijkp1 - 2.0 * T_ijk + T_ijkm1) / (dz * dz); // d^2T/dz^2
                }
                // Damping function
                if (k == 0) {
                    tau_p = sqrt((nu * 0.5 * (uz + wx)) * (nu * 0.5 * (uz + wx)) + (nu * 0.5 * (vz + wy)) * (nu * 0.5 * (vz + wy)));
                } else {
                    tau_p = 0.0;
                }
                u_tau = sqrt(tau_p);
                fw = f_damping(z[k], u_tau, nu);
                // Copy to turbulence array
                R_turbulence[parameters->turbulence_indexes.ux  + IDX(i, j, k, Nx, Ny, Nz)] = ux;
                R_turbulence[parameters->turbulence_indexes.uy  + IDX(i, j, k, Nx, Ny, Nz)] = uy;
                R_turbulence[parameters->turbulence_indexes.uz  + IDX(i, j, k, Nx, Ny, Nz)] = uz;
                R_turbulence[parameters->turbulence_indexes.vx  + IDX(i, j, k, Nx, Ny, Nz)] = vx;
                R_turbulence[parameters->turbulence_indexes.vy  + IDX(i, j, k, Nx, Ny, Nz)] = vy;
                R_turbulence[parameters->turbulence_indexes.vz  + IDX(i, j, k, Nx, Ny, Nz)] = vz;
                R_turbulence[parameters->turbulence_indexes.wx  + IDX(i, j, k, Nx, Ny, Nz)] = wx;
                R_turbulence[parameters->turbulence_indexes.wy  + IDX(i, j, k, Nx, Ny, Nz)] = wy;
                R_turbulence[parameters->turbulence_indexes.wz  + IDX(i, j, k, Nx, Ny, Nz)] = wz;
                R_turbulence[parameters->turbulence_indexes.Tx  + IDX(i, j, k, Nx, Ny, Nz)] = Tx;
                R_turbulence[parameters->turbulence_indexes.Ty  + IDX(i, j, k, Nx, Ny, Nz)] = Ty;
                R_turbulence[parameters->turbulence_indexes.Tz  + IDX(i, j, k, Nx, Ny, Nz)] = Tz;
                R_turbulence[parameters->turbulence_indexes.uxx + IDX(i, j, k, Nx, Ny, Nz)] = uxx;
                R_turbulence[parameters->turbulence_indexes.uyy + IDX(i, j, k, Nx, Ny, Nz)] = uyy;
                R_turbulence[parameters->turbulence_indexes.uzz + IDX(i, j, k, Nx, Ny, Nz)] = uzz;
                R_turbulence[parameters->turbulence_indexes.vxx + IDX(i, j, k, Nx, Ny, Nz)] = vxx;
                R_turbulence[parameters->turbulence_indexes.vyy + IDX(i, j, k, Nx, Ny, Nz)] = vyy;
                R_turbulence[parameters->turbulence_indexes.vzz + IDX(i, j, k, Nx, Ny, Nz)] = vzz;
                R_turbulence[parameters->turbulence_indexes.wxx + IDX(i, j, k, Nx, Ny, Nz)] = wxx;
                R_turbulence[parameters->turbulence_indexes.wyy + IDX(i, j, k, Nx, Ny, Nz)] = wyy;
                R_turbulence[parameters->turbulence_indexes.wzz + IDX(i, j, k, Nx, Ny, Nz)] = wzz;
                R_turbulence[parameters->turbulence_indexes.Txx + IDX(i, j, k, Nx, Ny, Nz)] = Txx;
                R_turbulence[parameters->turbulence_indexes.Tyy + IDX(i, j, k, Nx, Ny, Nz)] = Tyy;
                R_turbulence[parameters->turbulence_indexes.Tzz + IDX(i, j, k, Nx, Ny, Nz)] = Tzz;
                R_turbulence[parameters->turbulence_indexes.fw  + IDX(i, j, k, Nx, Ny, Nz)] = fw;
                /* Compute fuel and source term */
                if (k < Nz_Y) {
                    Y_ijk = R_old[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)];
                    Y_RHS = -Y_f * K(T_ijk, A, T_a) * H(T_ijk, T_pc) * Y_ijk;
                    R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] = Y_RHS;
                } else {
                    Y_ijk = 0.0;
                }
                // Compute source and force terms
                S = source(T_ijk, Y_ijk, H_R, A, T_a, h, a_v, T_inf, c_p, rho, T_pc);
                mod_U = sqrt(u_ijk * u_ijk + v_ijk * v_ijk + w_ijk * w_ijk);
                // Force terms
                F_x = - Y_D * a_v * Y_ijk * mod_U * u_ijk;
                F_y = - Y_D * a_v * Y_ijk * mod_U * v_ijk;                             
                F_z = - g * (T_ijk - T_inf) / T_ijk - Y_D * a_v * Y_ijk * mod_U * w_ijk;
                // Compute Laplacian terms
                lap_u = uxx + uyy + uzz;
                lap_v = vxx + vyy + vzz;
                lap_w = wxx + wyy + wzz;
                lap_T = Txx + Tyy + Tzz;
                // Compute RHS
                u_RHS = nu * lap_u - (uux + vuy + wuz) + F_x;
                v_RHS = nu * lap_v - (uvx + vvy + wvz) + F_y;
                w_RHS = nu * lap_w - (uwx + vwy + wwz) + F_z;
                T_RHS = alpha * lap_T - (u_ijk * Tx + v_ijk * Ty + w_ijk * Tz) + S;
                // Save RHS into R_new
                R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] = u_RHS;
                R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] = v_RHS;
                R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] = w_RHS;
                R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = T_RHS;  
            }
        }
    }
}

void Phi(double t, double *R_old, double *R_new, double *R_turbulence, Parameters *parameters) {
    RHS_v1(t, R_old, R_new, R_turbulence, parameters);
    turbulence(R_turbulence, R_new, parameters);
    boundary_conditions(R_new, parameters);
}

void boundary_conditions(double *R_new, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->Nz_Y;
    int u_index = parameters->field_indexes.u;
    int v_index = parameters->field_indexes.v;
    int w_index = parameters->field_indexes.w;
    int T_index = parameters->field_indexes.T;
    int Y_index = parameters->field_indexes.Y;
    int k;
    double u_ijkm1, u_ijkm2;
    double v_ijkm1, v_ijkm2;
    double w_ijkm1, w_ijkm2;
    double T_ijkp1, T_ijkm1, T_ijkm2, T_ijkp2;
    // double Y_ijkp1, Y_ijkm1, Y_ijkm2, Y_ijkp2;
    double Y_ijkp1, Y_ijkp2;
    // Periodic boundary conditions are included in RHS computation, we only compute in top and bottom boundaries (z=z_min and z=z_max)
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            // Bottom boundary z_k = z_min, k=0
            k = 0; 
            T_ijkp1 = R_new[T_index + IDX(i, j, k + 1, Nx, Ny, Nz)]; // T_{i,j,k+1}
            T_ijkp2 = R_new[T_index + IDX(i, j, k + 2, Nx, Ny, Nz)]; // T_{i,j,k+2}
            if (k + 1 < Nz_Y) // Check if Y_ijkp1 is part of the fuel
                Y_ijkp1 = R_new[Y_index + IDX(i, j, k + 1, Nx, Ny, Nz_Y)]; // Y_{i,j,k+1}
            else // If not, set Y_ijkp1 = 0 (because we don't store Y when it's 0)
                Y_ijkp1 = 0.0;
            if (k + 2 < Nz_Y) // Same as above
                Y_ijkp2 = R_new[Y_index + IDX(i, j, k + 2, Nx, Ny, Nz_Y)]; // Y_{i,j,k+2}
            else
                Y_ijkp2 = 0.0;
            R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] = 0; // u = 0
            R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] = 0; // v = 0
            R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] = 0; // w = 0
            R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = (4 * T_ijkp1 - T_ijkp2) / 3; // dT/dz = 0
            R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] = (4 * Y_ijkp1 - Y_ijkp2) / 3; // dY/dz = 0
            // Top boundary z_k = z_max, k=Nz-1
            k = Nz - 1;
            u_ijkm1 = R_new[u_index + IDX(i, j, k - 1, Nx, Ny, Nz)]; // u_{i,j,k-1}
            u_ijkm2 = R_new[u_index + IDX(i, j, k - 2, Nx, Ny, Nz)]; // u_{i,j,k-2}
            v_ijkm1 = R_new[v_index + IDX(i, j, k - 1, Nx, Ny, Nz)]; // v_{i,j,k-1}
            v_ijkm2 = R_new[v_index + IDX(i, j, k - 2, Nx, Ny, Nz)]; // v_{i,j,k-2}
            w_ijkm1 = R_new[w_index + IDX(i, j, k - 1, Nx, Ny, Nz)]; // w_{i,j,k-1}
            w_ijkm2 = R_new[w_index + IDX(i, j, k - 2, Nx, Ny, Nz)]; // w_{i,j,k-2}
            T_ijkm1 = R_new[T_index + IDX(i, j, k - 1, Nx, Ny, Nz)]; // T_{i,j,k-1}
            T_ijkm2 = R_new[T_index + IDX(i, j, k - 2, Nx, Ny, Nz)]; // T_{i,j,k-2}
            R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] = (4 * u_ijkm1 - u_ijkm2) / 3; // du/dz = 0
            R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] = (4 * v_ijkm1 - v_ijkm2) / 3; // dv/dz = 0
            R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] = (4 * w_ijkm1 - w_ijkm2) / 3; // dw/dz = 0
            R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = (4 * T_ijkm1 - T_ijkm2) / 3; // dT/dz = 0   
            // Actually we don't need to set Y at the top boundary, because it's not used in the computation of the RHS 
        }
    }
}

void bounds(double *R_new, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->Nz_Y;
    int T_index = parameters->field_indexes.T;
    int Y_index = parameters->field_indexes.Y;
    double T_min = parameters->T_min;
    double T_max = parameters->T_max;
    double Y_min = parameters->Y_min;
    double Y_max = parameters->Y_max;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                // Check if Y_ijk is not valid
                if (k < Nz_Y) {
                    if (R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] < Y_min) {
                        R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] = Y_min;
                    }
                    if (R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] > Y_max) {
                        R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] = Y_max;
                    }
                }
                // Check if T_ijk is less than T_inf and higher than T_max
                if (R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] < T_min) {
                    R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = T_min;
                }
                if (R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] > T_max) {
                    R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = T_max;
                }
            }
        }
    }
}

void velocity_correction(double *R_new, double *p, double dt, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int u_index = parameters->field_indexes.u;
    int v_index = parameters->field_indexes.v;
    int w_index = parameters->field_indexes.w;
    double rho = parameters->rho;
    double dx = parameters->dx;
    double dy = parameters->dy;
    double dz = parameters->dz;
    double p_ijk, p_im1jk, p_ip1jk, p_ijm1k, p_ijp1k, p_ijkp1, p_ijkm1, p_ijkp2, p_ijkm2;
    double px, py, pz;
    int im1, ip1, jm1, jp1;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                // Indexes for periodic boundary conditions
                im1 = (i - 1 + Nx - 1) % (Nx - 1); // i-1
                jm1 = (j - 1 + Ny - 1) % (Ny - 1); // j-1
                ip1 = (i + 1) % (Nx - 1); // i+1
                jp1 = (j + 1) % (Ny - 1); // j+1
                p_im1jk = p[IDX(im1, j, k, Nx, Ny, Nz)]; // p_{i-1,j,k}
                p_ip1jk = p[IDX(ip1, j, k, Nx, Ny, Nz)]; // p_{i+1,j,k}
                p_ijm1k = p[IDX(i, jm1, k, Nx, Ny, Nz)]; // p_{i,j-1,k}
                p_ijp1k = p[IDX(i, jp1, k, Nx, Ny, Nz)]; // p_{i,j+1,k}
                px = (p_ip1jk - p_im1jk) / (2 * dx);
                py = (p_ijp1k - p_ijm1k) / (2 * dy);
                if (k == 0) { // Second order forward difference
                    p_ijk = p[IDX(i, j, k, Nx, Ny, Nz)]; // p_{i,j,k}
                    p_ijkp1 = p[IDX(i, j, k + 1, Nx, Ny, Nz)]; // p_{i,j,k+1}
                    p_ijkp2 = p[IDX(i, j, k + 2, Nx, Ny, Nz)]; // p_{i,j,k+2}
                    pz = (-3 * p_ijk + 4 * p_ijkp1 - p_ijkp2) / (2 * dz);
                } else if (k == Nz - 1) { // Second order backward difference
                    p_ijk = p[IDX(i, j, k, Nx, Ny, Nz)]; // p_{i,j,k}
                    p_ijkm1 = p[IDX(i, j, k - 1, Nx, Ny, Nz)]; // p_{i,j,k-1}
                    p_ijkm2 = p[IDX(i, j, k - 2, Nx, Ny, Nz)]; // p_{i,j,k-2}
                    pz = (3 * p_ijk - 4 * p_ijkm1 + p_ijkm2) / (2 * dz);
                }
                else { // Central difference
                    p_ijkp1 = p[IDX(i, j, k + 1, Nx, Ny, Nz)]; // p_{i,j,k+1}
                    p_ijkm1 = p[IDX(i, j, k - 1, Nx, Ny, Nz)]; // p_{i,j,k-1}
                    pz = (p_ijkp1 - p_ijkm1) / (2 * dz);
                }
                // Correct velocity
                R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * px;
                R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * py;
                R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * pz;
            }
        }
    }
}

void velocity_correction_fw(double *R_new, double *p, double dt, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int u_index = parameters->field_indexes.u;
    int v_index = parameters->field_indexes.v;
    int w_index = parameters->field_indexes.w;
    double rho = parameters->rho;
    double dx = parameters->dx;
    double dy = parameters->dy;
    double dz = parameters->dz;
    double p_ijk, p_im1jk, p_ip1jk, p_ijm1k, p_ijp1k, p_ijkp1, p_ijkm1, p_ijkp2, p_ijkm2;
    double px, py, pz;
    int im1, ip1, jm1, jp1;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                // Indexes for periodic boundary conditions
                im1 = (i - 1 + Nx - 1) % (Nx - 1); // i-1
                jm1 = (j - 1 + Ny - 1) % (Ny - 1); // j-1
                ip1 = (i + 1) % (Nx - 1); // i+1
                jp1 = (j + 1) % (Ny - 1); // j+1
                p_im1jk = p[IDX(im1, j, k, Nx, Ny, Nz)]; // p_{i-1,j,k}
                p_ip1jk = p[IDX(ip1, j, k, Nx, Ny, Nz)]; // p_{i+1,j,k}
                p_ijm1k = p[IDX(i, jm1, k, Nx, Ny, Nz)]; // p_{i,j-1,k}
                p_ijp1k = p[IDX(i, jp1, k, Nx, Ny, Nz)]; // p_{i,j+1,k}
                px = (p_ip1jk - p_im1jk) / (2 * dx);
                py = (p_ijp1k - p_ijm1k) / (2 * dy);
                if (k < Nz - 2) { // Second order forward difference
                    p_ijk = p[IDX(i, j, k, Nx, Ny, Nz)]; // p_{i,j,k}
                    p_ijkp1 = p[IDX(i, j, k + 1, Nx, Ny, Nz)]; // p_{i,j,k+1}
                    p_ijkp2 = p[IDX(i, j, k + 2, Nx, Ny, Nz)]; // p_{i,j,k+2}
                    pz = (-3 * p_ijk + 4 * p_ijkp1 - p_ijkp2) / (2 * dz);
                } else { // Second order backward difference
                    p_ijk = p[IDX(i, j, k, Nx, Ny, Nz)]; // p_{i,j,k}
                    p_ijkm1 = p[IDX(i, j, k - 1, Nx, Ny, Nz)]; // p_{i,j,k-1}
                    p_ijkm2 = p[IDX(i, j, k - 2, Nx, Ny, Nz)]; // p_{i,j,k-2}
                    pz = (3 * p_ijk - 4 * p_ijkm1 + p_ijkm2) / (2 * dz);
                }
                // Correct velocity
                R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * px;
                R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * py;
                R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * pz;
            }
        }
    }
}

void velocity_correction_bw(double *R_new, double *p, double dt, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int u_index = parameters->field_indexes.u;
    int v_index = parameters->field_indexes.v;
    int w_index = parameters->field_indexes.w;
    double rho = parameters->rho;
    double dx = parameters->dx;
    double dy = parameters->dy;
    double dz = parameters->dz;
    double p_ijk, p_im1jk, p_ip1jk, p_ijm1k, p_ijp1k, p_ijkp1, p_ijkm1, p_ijkp2, p_ijkm2;
    double px, py, pz;
    int im1, ip1, jm1, jp1;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                // Indexes for periodic boundary conditions
                im1 = (i - 1 + Nx - 1) % (Nx - 1); // i-1
                jm1 = (j - 1 + Ny - 1) % (Ny - 1); // j-1
                ip1 = (i + 1) % (Nx - 1); // i+1
                jp1 = (j + 1) % (Ny - 1); // j+1
                p_im1jk = p[IDX(im1, j, k, Nx, Ny, Nz)]; // p_{i-1,j,k}
                p_ip1jk = p[IDX(ip1, j, k, Nx, Ny, Nz)]; // p_{i+1,j,k}
                p_ijm1k = p[IDX(i, jm1, k, Nx, Ny, Nz)]; // p_{i,j-1,k}
                p_ijp1k = p[IDX(i, jp1, k, Nx, Ny, Nz)]; // p_{i,j+1,k}
                px = (p_ip1jk - p_im1jk) / (2 * dx);
                py = (p_ijp1k - p_ijm1k) / (2 * dy);
                if (k < 2) { // Second order forward difference
                    p_ijk = p[IDX(i, j, k, Nx, Ny, Nz)]; // p_{i,j,k}
                    p_ijkp1 = p[IDX(i, j, k + 1, Nx, Ny, Nz)]; // p_{i,j,k+1}
                    p_ijkp2 = p[IDX(i, j, k + 2, Nx, Ny, Nz)]; // p_{i,j,k+2}
                    pz = (-3 * p_ijk + 4 * p_ijkp1 - p_ijkp2) / (2 * dz);
                } else { // Second order backward difference
                    p_ijk = p[IDX(i, j, k, Nx, Ny, Nz)]; // p_{i,j,k}
                    p_ijkm1 = p[IDX(i, j, k - 1, Nx, Ny, Nz)]; // p_{i,j,k-1}
                    p_ijkm2 = p[IDX(i, j, k - 2, Nx, Ny, Nz)]; // p_{i,j,k-2}
                    pz = (3 * p_ijk - 4 * p_ijkm1 + p_ijkm2) / (2 * dz);
                }
                // Correct velocity
                R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * px;
                R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * py;
                R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * pz;
            }
        }
    }
}