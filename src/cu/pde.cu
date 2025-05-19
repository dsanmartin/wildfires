/**
 * @file pde.cu
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Functions for solving the partial differential equations of the wildfire simulation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "../../include/cu/pde.cuh"

__global__
void boundary_conditions(double *R_new, double *z, int *Nz_Y, int *cut_nodes, Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int Nz_Y_max = parameters.Nz_Y_max;
    int u_index = parameters.field_indexes.u;
    int v_index = parameters.field_indexes.v;
    int w_index = parameters.field_indexes.w;
    int T_index = parameters.field_indexes.T;
    int Y_index = parameters.field_indexes.Y;
    int k;
    double u_dead_nodes = parameters.u_dead_nodes;
    double v_dead_nodes = parameters.v_dead_nodes;
    double w_dead_nodes = parameters.w_dead_nodes;
    double T_dead_nodes = parameters.T_dead_nodes;
    double Y_dead_nodes = parameters.Y_dead_nodes;
    double u_ijkm1, u_ijkm2;
    double v_ijkm1, v_ijkm2;
    double w_ijkm1, w_ijkm2;
    double T_ijkp1, T_ijkm1, T_ijkm2, T_ijkp2;
    double Y_ijkp1, Y_ijkp2;
    double dz_k, dz_kp1, dz_Nzm2, dz_Nzm1;
    // Periodic boundary conditions are included in RHS computation, we only compute in top and bottom boundaries (z=z_min and z=z_max)
    int size = Nx * Ny;
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;
    for (int ij = idx; ij < size; ij += stride) {
        int i = ij / Ny;
        int j = ij % Ny;        
        // Bottom boundary z_k = z_min, k=0
        k = cut_nodes[IDX(i, j, 0, Nx, Ny, 1)]; 
        dz_k = z[k + 1] - z[k];
        dz_kp1 = z[k + 2] - z[k + 1];
        T_ijkp1 = R_new[T_index + IDX(i, j, k + 1, Nx, Ny, Nz)]; // T_{i,j,k+1}
        T_ijkp2 = R_new[T_index + IDX(i, j, k + 2, Nx, Ny, Nz)]; // T_{i,j,k+2}
        if (k + 1 < Nz_Y[IDX(i, j, 0, Nx, Ny, 1)]) // Check if Y_ijkp1 is part of the fuel
            Y_ijkp1 = R_new[Y_index + IDX(i, j, k + 1, Nx, Ny, Nz_Y_max)]; // Y_{i,j,k+1}
        else // If not, set Y_ijkp1 = 0 (because we don't store Y when it's 0)
            Y_ijkp1 = 0.0;
        if (k + 2 < Nz_Y[IDX(i, j, 0, Nx, Ny, 1)]) // Same as above
            Y_ijkp2 = R_new[Y_index + IDX(i, j, k + 2, Nx, Ny, Nz_Y_max)]; // Y_{i,j,k+2}
        else
            Y_ijkp2 = 0.0;
        R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] = 0; // u = 0
        R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] = 0; // v = 0
        R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] = 0; // w = 0
        // R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = (4 * T_ijkp1 - T_ijkp2) / 3; // dT/dz = 0 (equispaced nodes)
        R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = ((dz_k + dz_kp1) * (dz_k + dz_kp1) * T_ijkp1 - dz_k * dz_k * T_ijkp2) / (dz_kp1 * (2 * dz_k + dz_kp1)); // dT/dz = 0 (non-equispaced nodes)
        if (k < Nz_Y[IDX(i, j, 0, Nx, Ny, 1)])
            // R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y_max)] = (4 * Y_ijkp1 - Y_ijkp2) / 3; // dY/dz = 0 (equispaced nodes)
            R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y_max)] = ((dz_k + dz_kp1) * (dz_k + dz_kp1) * Y_ijkp1 - dz_k * dz_k * Y_ijkp2) / (dz_kp1 * (2 * dz_k + dz_kp1)); // dY/dz = 0 (non-equispaced nodes)
        // Top boundary z_k = z_max, k=Nz-1
        k = Nz - 1;
        dz_Nzm1 = z[k] - z[k - 1];
        dz_Nzm2 = z[k - 1] - z[k - 2];
        u_ijkm1 = R_new[u_index + IDX(i, j, k - 1, Nx, Ny, Nz)]; // u_{i,j,k-1}
        u_ijkm2 = R_new[u_index + IDX(i, j, k - 2, Nx, Ny, Nz)]; // u_{i,j,k-2}
        v_ijkm1 = R_new[v_index + IDX(i, j, k - 1, Nx, Ny, Nz)]; // v_{i,j,k-1}
        v_ijkm2 = R_new[v_index + IDX(i, j, k - 2, Nx, Ny, Nz)]; // v_{i,j,k-2}
        w_ijkm1 = R_new[w_index + IDX(i, j, k - 1, Nx, Ny, Nz)]; // w_{i,j,k-1}
        w_ijkm2 = R_new[w_index + IDX(i, j, k - 2, Nx, Ny, Nz)]; // w_{i,j,k-2}
        T_ijkm1 = R_new[T_index + IDX(i, j, k - 1, Nx, Ny, Nz)]; // T_{i,j,k-1}
        T_ijkm2 = R_new[T_index + IDX(i, j, k - 2, Nx, Ny, Nz)]; // T_{i,j,k-2}
        // Equispaced nodes
        // R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] = (4 * u_ijkm1 - u_ijkm2) / 3; // du/dz = 0
        // R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] = (4 * v_ijkm1 - v_ijkm2) / 3; // dv/dz = 0
        // R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] = (4 * w_ijkm1 - w_ijkm2) / 3; // dw/dz = 0
        // R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = (4 * T_ijkm1 - T_ijkm2) / 3; // dT/dz = 0   
        // Non-equispaced nodes
        R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] = ((dz_Nzm2 + dz_Nzm1) * (dz_Nzm2 + dz_Nzm1) * u_ijkm1 - dz_Nzm1 * dz_Nzm1 * u_ijkm2) / (dz_Nzm2 * (dz_Nzm2 + 2 * dz_Nzm1)); // du/dz = 0
        R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] = ((dz_Nzm2 + dz_Nzm1) * (dz_Nzm2 + dz_Nzm1) * v_ijkm1 - dz_Nzm1 * dz_Nzm1 * v_ijkm2) / (dz_Nzm2 * (dz_Nzm2 + 2 * dz_Nzm1)); // dv/dz = 0
        R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] = ((dz_Nzm2 + dz_Nzm1) * (dz_Nzm2 + dz_Nzm1) * w_ijkm1 - dz_Nzm1 * dz_Nzm1 * w_ijkm2) / (dz_Nzm2 * (dz_Nzm2 + 2 * dz_Nzm1)); // dw/dz = 0
        R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = ((dz_Nzm2 + dz_Nzm1) * (dz_Nzm2 + dz_Nzm1) * T_ijkm1 - dz_Nzm1 * dz_Nzm1 * T_ijkm2) / (dz_Nzm2 * (dz_Nzm2 + 2 * dz_Nzm1)); // dT/dz = 0   
        // Actually we don't need to set Y at the top boundary, because it's not used in the computation of the RHS 
        // Set dead nodes values
        for (k = 0; k < cut_nodes[IDX(i, j, 0, Nx, Ny, 1)]; k++) {
            R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] = u_dead_nodes;
            R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] = v_dead_nodes;
            R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] = w_dead_nodes;
            R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = T_dead_nodes;
            if (k < Nz_Y[IDX(i, j, 0, Nx, Ny, 1)])
                R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y_max)] = Y_dead_nodes;
        }
    }
}

__global__
void bounds(double *R_new, Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int Nz_Y_max = parameters.Nz_Y_max;
    int T_index = parameters.field_indexes.T;
    int Y_index = parameters.field_indexes.Y;
    double T_min = parameters.T_min;
    double T_max = parameters.T_max;
    double Y_min = parameters.Y_min;
    double Y_max = parameters.Y_max;
    int size = Nx * Ny * Nz;
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;
    for (int ijk = idx; ijk < size; ijk += stride) {
        int i = ijk / (Ny * Nz);
        int j = (ijk % (Ny * Nz)) / Nz;
        int k = ijk % Nz;
        // Check if Y_ijk is not valid
        if (k < Nz_Y_max) {
            if (R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y_max)] < Y_min) {
                R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y_max)] = Y_min;
            }
            if (R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y_max)] > Y_max) {
                R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y_max)] = Y_max;
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

__global__
void RHS(double t, double *R_old, double *R_new, double *R_turbulence, double *z, double *z_ibm, int *Nz_Y, Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int Nz_Y_max = parameters.Nz_Y_max;
    double dx = parameters.dx;
    double dy = parameters.dy;
    // double dz = parameters.dz;
    // Model parameters
    double mu = parameters.mu;
    // double alpha = parameters.alpha;
    double kappa = parameters.kappa;
    double Y_f = parameters.Y_f;
    double H_C = parameters.H_C;
    double A = parameters.A;
    double T_act = parameters.T_act;
    double T_pc = parameters.T_pc;
    double h_c = parameters.h_c;
    // double a_v = parameters.a_v;
    double alpha_s = parameters.alpha_s;
    double sigma_s = parameters.sigma_s;
    double T_inf = parameters.T_inf;
    double c_p = parameters.c_p;
    double rho_inf = parameters.rho_inf;
    // double Y_D = parameters.Y_D;
    double C_d = parameters.C_d;
    double g = parameters.g;
    double delta = parameters.delta;
    double C_s = parameters.C_s;
    // Fields indexes
    int u_index = parameters.field_indexes.u;
    int v_index = parameters.field_indexes.v;
    int w_index = parameters.field_indexes.w;
    int T_index = parameters.field_indexes.T;
    int Y_index = parameters.field_indexes.Y;
    // Fields nodes
    double u_ijk, u_ip1jk, u_im1jk, u_ijp1k, u_ijm1k, u_ijkp1, u_ijkm1;
    double u_ip2jk, u_im2jk, u_ijp2k, u_ijm2k, u_ijkp2, u_ijkm2, u_ijkp3, u_ijkm3;
    double v_ijk, v_ip1jk, v_im1jk, v_ijp1k, v_ijm1k, v_ijkp1, v_ijkm1;
    double v_ip2jk, v_im2jk, v_ijp2k, v_ijm2k, v_ijkp2, v_ijkm2, v_ijkp3, v_ijkm3;
    double w_ijk, w_ip1jk, w_im1jk, w_ijp1k, w_ijm1k, w_ijkp1, w_ijkm1;
    double w_ip2jk, w_im2jk, w_ijp2k, w_ijm2k, w_ijkp2, w_ijkm2, w_ijkp3, w_ijkm3;
    double T_ijk, T_ip1jk, T_im1jk, T_ijp1k, T_ijm1k, T_ijkp1, T_ijkm1;
    double T_ijkp2, T_ijkm2, T_ijkp3, T_ijkm3;
    double rho_ijk;
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
    double u_RHS, v_RHS, w_RHS, T_RHS, Y_RHS, q, T_RHS_tmp;
    double mod_U, nu;
    double F_x, F_y, F_z;
    double H_step, K_T;
    double dz_km3, dz_km2, dz_km1, dz_k, dz_kp1, dz_kp2;
    double u_tau, tau_w, fw;
    double S_11, S_12, S_13, S_21, S_22, S_23, S_31, S_32, S_33;
    double mod_S, mu_sgs;
    double T_gas, T_solid;
    int i, j, k;
    int im1, ip1, jm1, jp1;
    int im2, ip2, jm2, jp2;
    int size = Nx * Ny * Nz;
    // Loop over interior nodes. Periodic boundary conditions in x and y.
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;
    for (int ijk = idx; ijk < size; ijk += stride) {
        i = ijk / (Ny * Nz);
        j = (ijk % (Ny * Nz)) / Nz;
        k = ijk % Nz; 
        /* Get fields nodes */
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
        rho_ijk = T_inf * rho_inf / T_ijk; // Density
        Y_ijk = 0.0; // Solid fuel fraction
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
        // Get dz values
        // \phi_{i,j,k-1}
        u_ijkm1 = 0, v_ijkm1 = 0, w_ijkm1 = 0, T_ijkm1 = 0;
        if (k > 0) {
            u_ijkm1 = R_old[u_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
            v_ijkm1 = R_old[v_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
            w_ijkm1 = R_old[w_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
            T_ijkm1 = R_old[T_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
        } 
        // \phi_{i,j,k+1}
        u_ijkp1 = 0, v_ijkp1 = 0, w_ijkp1 = 0, T_ijkp1 = 0;
        if (k < Nz - 1) {
            u_ijkp1 = R_old[u_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
            v_ijkp1 = R_old[v_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
            w_ijkp1 = R_old[w_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
            T_ijkp1 = R_old[T_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
        }
        // \phi_{i,j,k-2}
        u_ijkm2 = 0, v_ijkm2 = 0, w_ijkm2 = 0, T_ijkm2 = 0;
        if (k > 1) {
            u_ijkm2 = R_old[u_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
            v_ijkm2 = R_old[v_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
            w_ijkm2 = R_old[w_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
            T_ijkm2 = R_old[T_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
        }
        // \phi_{i,j,k+2}
        u_ijkp2 = 0, v_ijkp2 = 0, w_ijkp2 = 0, T_ijkp2 = 0;
        if (k < Nz - 2) {
            u_ijkp2 = R_old[u_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
            v_ijkp2 = R_old[v_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
            w_ijkp2 = R_old[w_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
            T_ijkp2 = R_old[T_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
        }
        // \phi_{i,j,k-3}
        u_ijkm3 = 0, v_ijkm3 = 0, w_ijkm3 = 0, T_ijkm3 = 0;
        if (k > 2) {
            u_ijkm3 = R_old[u_index + IDX(i, j, k - 3, Nx, Ny, Nz)];
            v_ijkm3 = R_old[v_index + IDX(i, j, k - 3, Nx, Ny, Nz)];
            w_ijkm3 = R_old[w_index + IDX(i, j, k - 3, Nx, Ny, Nz)];
            T_ijkm3 = R_old[T_index + IDX(i, j, k - 3, Nx, Ny, Nz)];
        }
        // \phi_{i,j,k+3}
        u_ijkp3 = 0, v_ijkp3 = 0, w_ijkp3 = 0, T_ijkp3 = 0;
        if (k < Nz - 3) {
            u_ijkp3 = R_old[u_index + IDX(i, j, k + 3, Nx, Ny, Nz)];
            v_ijkp3 = R_old[v_index + IDX(i, j, k + 3, Nx, Ny, Nz)];
            w_ijkp3 = R_old[w_index + IDX(i, j, k + 3, Nx, Ny, Nz)];
            T_ijkp3 = R_old[T_index + IDX(i, j, k + 3, Nx, Ny, Nz)];
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
            dz_k = z[k + 1] - z[k];
            dz_kp1 = z[k + 2] - z[k + 1]; 
            // Equispaced
            // u_km = (-3 * u_ijk + 4 * u_ijkp1 - u_ijkp2) / (2 * dz);
            // v_km = (-3 * v_ijk + 4 * v_ijkp1 - v_ijkp2) / (2 * dz);
            // w_km = (-3 * w_ijk + 4 * w_ijkp1 - w_ijkp2) / (2 * dz);
            // Non-equispaced grid
            u_km = - (2 * dz_k + dz_kp1) * u_ijk / (dz_k * (dz_k + dz_kp1)) 
                + (dz_k + dz_kp1) * u_ijkp1 / (dz_k * dz_kp1) 
                - dz_k * u_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            v_km = - (2 * dz_k + dz_kp1) * v_ijk / (dz_k * (dz_k + dz_kp1)) 
                + (dz_k + dz_kp1) * v_ijkp1 / (dz_k * dz_kp1) 
                - dz_k * v_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            w_km = - (2 * dz_k + dz_kp1) * w_ijk / (dz_k * (dz_k + dz_kp1)) 
                + (dz_k + dz_kp1) * w_ijkp1 / (dz_k * dz_kp1) 
                - dz_k * w_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
        } else {
            dz_km1 = z[k] - z[k - 1];
            dz_km2 = z[k - 1] - z[k - 2];
            // Equispaced
            // u_km = (3 * u_ijk - 4 * u_ijkm1 + u_ijkm2) / (2 * dz);
            // v_km = (3 * v_ijk - 4 * v_ijkm1 + v_ijkm2) / (2 * dz);
            // w_km = (3 * w_ijk - 4 * w_ijkm1 + w_ijkm2) / (2 * dz);
            // Non-equispaced grid
            u_km = dz_km1 * u_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) 
                - (dz_km2 + dz_km1) * u_ijkm1 / (dz_km2 * dz_km1) 
                + (dz_km2 + 2 * dz_km1) * u_ijk / (dz_km1 * (dz_km2 + dz_km1));
            v_km = dz_km1 * v_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) 
                - (dz_km2 + dz_km1) * v_ijkm1 / (dz_km2 * dz_km1) 
                + (dz_km2 + 2 * dz_km1) * v_ijk / (dz_km1 * (dz_km2 + dz_km1));
            w_km = dz_km1 * w_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) 
                - (dz_km2 + dz_km1) * w_ijkm1 / (dz_km2 * dz_km1) 
                + (dz_km2 + 2 * dz_km1) * w_ijk / (dz_km1 * (dz_km2 + dz_km1));
        }
        // Second order backward difference at k=Nz-2 and k=Nz-1
        if (k >= Nz - 2) {
            dz_km1 = z[k] - z[k - 1];
            dz_km2 = z[k - 1] - z[k - 2];
            // Equispaced
            // u_kp = (3 * u_ijk - 4 * u_ijkm1 + u_ijkm2) / (2 * dz);
            // v_kp = (3 * v_ijk - 4 * v_ijkm1 + v_ijkm2) / (2 * dz);
            // w_kp = (3 * w_ijk - 4 * w_ijkm1 + w_ijkm2) / (2 * dz);
            // Non-equispaced grid
            u_kp = dz_km1 * u_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) 
                - (dz_km2 + dz_km1) * u_ijkm1 / (dz_km2 * dz_km1) 
                + (dz_km2 + 2 * dz_km1) * u_ijk / (dz_km1 * (dz_km2 + dz_km1));
            v_kp = dz_km1 * v_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) 
                - (dz_km2 + dz_km1) * v_ijkm1 / (dz_km2 * dz_km1) 
                + (dz_km2 + 2 * dz_km1) * v_ijk / (dz_km1 * (dz_km2 + dz_km1));
            w_kp = dz_km1 * w_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) 
                - (dz_km2 + dz_km1) * w_ijkm1 / (dz_km2 * dz_km1) 
                + (dz_km2 + 2 * dz_km1) * w_ijk / (dz_km1 * (dz_km2 + dz_km1));
        } else {
            dz_k = z[k + 1] - z[k];
            dz_kp1 = z[k + 2] - z[k + 1]; 
            // Equispaced
            // u_kp = (-3 * u_ijk + 4 * u_ijkp1 - u_ijkp2) / (2 * dz);
            // v_kp = (-3 * v_ijk + 4 * v_ijkp1 - v_ijkp2) / (2 * dz);
            // w_kp = (-3 * w_ijk + 4 * w_ijkp1 - w_ijkp2) / (2 * dz);
            // Non-equispaced
            u_kp = - (2 * dz_k + dz_kp1) * u_ijk / (dz_k * (dz_k + dz_kp1)) 
                + (dz_k + dz_kp1) * u_ijkp1 / (dz_k * dz_kp1) 
                - dz_k * u_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            v_kp = - (2 * dz_k + dz_kp1) * v_ijk / (dz_k * (dz_k + dz_kp1)) 
                + (dz_k + dz_kp1) * v_ijkp1 / (dz_k * dz_kp1) 
                - dz_k * v_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            w_kp = - (2 * dz_k + dz_kp1) * w_ijk / (dz_k * (dz_k + dz_kp1)) 
                + (dz_k + dz_kp1) * w_ijkp1 / (dz_k * dz_kp1) 
                - dz_k * w_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
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
            dz_k = z[k + 1] - z[k]; 
            dz_kp1 = z[k + 2] - z[k + 1]; 
            dz_kp2 = z[k + 3] - z[k + 2];
            // Equispaced grid
            // uz = (-3 * u_ijk + 4 * u_ijkp1 - u_ijkp2) / (2 * dz); // du/dz
            // vz = (-3 * v_ijk + 4 * v_ijkp1 - v_ijkp2) / (2 * dz); // dv/dz
            // wz = (-3 * w_ijk + 4 * w_ijkp1 - w_ijkp2) / (2 * dz); // dw/dz
            // Tz = (-3.0 * T_ijk + 4.0 * T_ijkp1 - T_ijkp2) / (2.0 * dz); // dT/dz
            // uzz = (2 * u_ijk - 5 * u_ijkp1 + 4 * u_ijkp2 - u_ijkp3) / (dz * dz); // d^2u/dz^2
            // vzz = (2 * v_ijk - 5 * v_ijkp1 + 4 * v_ijkp2 - v_ijkp3) / (dz * dz); // d^2v/dz^2
            // wzz = (2 * w_ijk - 5 * w_ijkp1 + 4 * w_ijkp2 - w_ijkp3) / (dz * dz); // d^2w/dz^2
            // Tzz = (2 * T_ijk - 5 * T_ijkp1 + 4 * T_ijkp2 - T_ijkp3) / (dz * dz); // d^2T/dz^2
            // Non-equispaced grid
            uz = - (2 * dz_k + dz_kp1) * u_ijk / (dz_k * (dz_k + dz_kp1)) 
                + (dz_k + dz_kp1) * u_ijkp1 / (dz_k * dz_kp1) 
                - dz_k * u_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            vz = - (2 * dz_k + dz_kp1) * v_ijk / (dz_k * (dz_k + dz_kp1)) 
                + (dz_k + dz_kp1) * v_ijkp1 / (dz_k * dz_kp1) 
                - dz_k * v_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            wz = - (2 * dz_k + dz_kp1) * w_ijk / (dz_k * (dz_k + dz_kp1)) 
                + (dz_k + dz_kp1) * w_ijkp1 / (dz_k * dz_kp1) 
                - dz_k * w_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            Tz = - (2 * dz_k + dz_kp1) * T_ijk / (dz_k * (dz_k + dz_kp1)) 
                + (dz_k + dz_kp1) * T_ijkp1 / (dz_k * dz_kp1) 
                - dz_k * T_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            uzz = (6 * dz_k + 4 * dz_kp1 + 2 * dz_kp2) * u_ijk / (dz_k * (dz_k + dz_kp1) * (dz_k + dz_kp1 + dz_kp2))
                - (4 * dz_k + 4 * dz_kp1 + 2 * dz_kp2) * u_ijkp1 / (dz_k * dz_kp1 * (dz_kp1 + dz_kp2))
                + (4 * dz_k + 2 * dz_kp1 + 2 * dz_kp2) * u_ijkp2 / (dz_kp1 * dz_kp2 * (dz_k + dz_kp1))
                - (4 * dz_k + 2 * dz_kp1) * u_ijkp3 / (dz_kp2 * (dz_kp1 + dz_kp2) * (dz_k + dz_kp1 + dz_kp2));
            vzz = (6 * dz_k + 4 * dz_kp1 + 2 * dz_kp2) * v_ijk / (dz_k * (dz_k + dz_kp1) * (dz_k + dz_kp1 + dz_kp2))
                - (4 * dz_k + 4 * dz_kp1 + 2 * dz_kp2) * v_ijkp1 / (dz_k * dz_kp1 * (dz_kp1 + dz_kp2))
                + (4 * dz_k + 2 * dz_kp1 + 2 * dz_kp2) * v_ijkp2 / (dz_kp1 * dz_kp2 * (dz_k + dz_kp1))
                - (4 * dz_k + 2 * dz_kp1) * v_ijkp3 / (dz_kp2 * (dz_kp1 + dz_kp2) * (dz_k + dz_kp1 + dz_kp2));
            wzz = (6 * dz_k + 4 * dz_kp1 + 2 * dz_kp2) * w_ijk / (dz_k * (dz_k + dz_kp1) * (dz_k + dz_kp1 + dz_kp2))
                - (4 * dz_k + 4 * dz_kp1 + 2 * dz_kp2) * w_ijkp1 / (dz_k * dz_kp1 * (dz_kp1 + dz_kp2))
                + (4 * dz_k + 2 * dz_kp1 + 2 * dz_kp2) * w_ijkp2 / (dz_kp1 * dz_kp2 * (dz_k + dz_kp1))
                - (4 * dz_k + 2 * dz_kp1) * w_ijkp3 / (dz_kp2 * (dz_kp1 + dz_kp2) * (dz_k + dz_kp1 + dz_kp2));
            Tzz = (6 * dz_k + 4 * dz_kp1 + 2 * dz_kp2) * T_ijk / (dz_k * (dz_k + dz_kp1) * (dz_k + dz_kp1 + dz_kp2))
                - (4 * dz_k + 4 * dz_kp1 + 2 * dz_kp2) * T_ijkp1 / (dz_k * dz_kp1 * (dz_kp1 + dz_kp2))
                + (4 * dz_k + 2 * dz_kp1 + 2 * dz_kp2) * T_ijkp2 / (dz_kp1 * dz_kp2 * (dz_k + dz_kp1))
                - (4 * dz_k + 2 * dz_kp1) * T_ijkp3 / (dz_kp2 * (dz_kp1 + dz_kp2) * (dz_k + dz_kp1 + dz_kp2));
        } else if (k == Nz - 1) {
            dz_km1 = z[k] - z[k - 1];
            dz_km2 = z[k - 1] - z[k - 2];
            dz_km3 = z[k - 2] - z[k - 3];
            // Equispaced grid
            // uz = (3 * u_ijk - 4 * u_ijkm1 + u_ijkm2) / (2 * dz); // du/dz
            // vz = (3 * v_ijk - 4 * v_ijkm1 + v_ijkm2) / (2 * dz); // dv/dz
            // wz = (3 * w_ijk - 4 * w_ijkm1 + w_ijkm2) / (2 * dz); // dw/dz
            // Tz = (3.0 * T_ijk - 4.0 * T_ijkm1 + T_ijkm2) / (2.0 * dz); // dT/dz
            // uzz = (2 * u_ijk - 5 * u_ijkm1 + 4 * u_ijkm2 - u_ijkm3) / (dz * dz); // d^2u/dz^2
            // vzz = (2 * v_ijk - 5 * v_ijkm1 + 4 * v_ijkm2 - v_ijkm3) / (dz * dz); // d^2v/dz^2
            // wzz = (2 * w_ijk - 5 * w_ijkm1 + 4 * w_ijkm2 - w_ijkm3) / (dz * dz); // d^2w/dz^2
            // Tzz = (2 * T_ijk - 5 * T_ijkm1 + 4 * T_ijkm2 - T_ijkm3) / (dz * dz); // d^2T/dz^2
            // Non-equispaced grid
            uz = dz_km1 * u_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) 
                - (dz_km2 + dz_km1) * u_ijkm1 / (dz_km2 * dz_km1) 
                + (dz_km2 + 2 * dz_km1) * u_ijk / (dz_km1 * (dz_km2 + dz_km1));
            vz = dz_km1 * v_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) 
                - (dz_km2 + dz_km1) * v_ijkm1 / (dz_km2 * dz_km1) 
                + (dz_km2 + 2 * dz_km1) * v_ijk / (dz_km1 * (dz_km2 + dz_km1));
            wz = dz_km1 * w_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) 
                - (dz_km2 + dz_km1) * w_ijkm1 / (dz_km2 * dz_km1) 
                + (dz_km2 + 2 * dz_km1) * w_ijk / (dz_km1 * (dz_km2 + dz_km1));
            Tz = dz_km1 * T_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) 
                - (dz_km2 + dz_km1) * T_ijkm1 / (dz_km2 * dz_km1) 
                + (dz_km2 + 2 * dz_km1) * T_ijk / (dz_km1 * (dz_km2 + dz_km1));
            uzz = - (2 * dz_km2 + 4 * dz_km1) * u_ijkm3 / (dz_km3 * (dz_km3 + dz_km2) * (dz_km3 + dz_km2 + dz_km1))
                + (2 * dz_km3 + 2 * dz_km2 + 4 * dz_km1) * u_ijkm2 / (dz_km3 * dz_km2 * (dz_km2 + dz_km1))
                - (2 * dz_km3 + 4 * dz_km2 + 4 * dz_km1) * u_ijkm1 / (dz_km2 * dz_km1 * (dz_km3 + dz_km2))
                + (2 * dz_km3 + 4 * dz_km2 + 6 * dz_km1) * u_ijk / (dz_km1 * (dz_km2 + dz_km1) * (dz_km3 + dz_km2 + dz_km1)); 
            vzz = - (2 * dz_km2 + 4 * dz_km1) * v_ijkm3 / (dz_km3 * (dz_km3 + dz_km2) * (dz_km3 + dz_km2 + dz_km1))
                + (2 * dz_km3 + 2 * dz_km2 + 4 * dz_km1) * v_ijkm2 / (dz_km3 * dz_km2 * (dz_km2 + dz_km1))
                - (2 * dz_km3 + 4 * dz_km2 + 4 * dz_km1) * v_ijkm1 / (dz_km2 * dz_km1 * (dz_km3 + dz_km2))
                + (2 * dz_km3 + 4 * dz_km2 + 6 * dz_km1) * v_ijk / (dz_km1 * (dz_km2 + dz_km1) * (dz_km3 + dz_km2 + dz_km1)); 
            wzz = - (2 * dz_km2 + 4 * dz_km1) * w_ijkm3 / (dz_km3 * (dz_km3 + dz_km2) * (dz_km3 + dz_km2 + dz_km1))
                + (2 * dz_km3 + 2 * dz_km2 + 4 * dz_km1) * w_ijkm2 / (dz_km3 * dz_km2 * (dz_km2 + dz_km1))
                - (2 * dz_km3 + 4 * dz_km2 + 4 * dz_km1) * w_ijkm1 / (dz_km2 * dz_km1 * (dz_km3 + dz_km2))
                + (2 * dz_km3 + 4 * dz_km2 + 6 * dz_km1) * w_ijk / (dz_km1 * (dz_km2 + dz_km1) * (dz_km3 + dz_km2 + dz_km1)); 
            Tzz = - (2 * dz_km2 + 4 * dz_km1) * T_ijkm3 / (dz_km3 * (dz_km3 + dz_km2) * (dz_km3 + dz_km2 + dz_km1))
                + (2 * dz_km3 + 2 * dz_km2 + 4 * dz_km1) * T_ijkm2 / (dz_km3 * dz_km2 * (dz_km2 + dz_km1))
                - (2 * dz_km3 + 4 * dz_km2 + 4 * dz_km1) * T_ijkm1 / (dz_km2 * dz_km1 * (dz_km3 + dz_km2))
                + (2 * dz_km3 + 4 * dz_km2 + 6 * dz_km1) * T_ijk / (dz_km1 * (dz_km2 + dz_km1) * (dz_km3 + dz_km2 + dz_km1)); 
        } else {
            dz_km1 = z[k] - z[k-1];
            dz_k = z[k + 1] - z[k];
            // Equispaced grid
            // uz  = (u_ijkp1 - u_ijkm1) / (2 * dz); // du/dz
            // vz  = (v_ijkp1 - v_ijkm1) / (2 * dz); // dv/dz
            // wz  = (w_ijkp1 - w_ijkm1) / (2 * dz); // dw/dz
            // Tz  = (T_ijkp1 - T_ijkm1) / (2.0 * dz); // dT/dz
            // uzz = (u_ijkp1 - 2.0 * u_ijk + u_ijkm1) / (dz * dz); // d^2u/dz^2
            // vzz = (v_ijkp1 - 2.0 * v_ijk + v_ijkm1) / (dz * dz); // d^2v/dz^2
            // wzz = (w_ijkp1 - 2.0 * w_ijk + w_ijkm1) / (dz * dz); // d^2w/dz^2
            // Tzz = (T_ijkp1 - 2.0 * T_ijk + T_ijkm1) / (dz * dz); // d^2T/dz^2
            // Non-equispaced grid
            uz = (u_ijkp1 - u_ijkm1) / (dz_km1 + dz_k); // du/dz
            vz = (v_ijkp1 - v_ijkm1) / (dz_km1 + dz_k); // dv/dz
            wz = (w_ijkp1 - w_ijkm1) / (dz_km1 + dz_k); // dw/dz
            Tz = (T_ijkp1 - T_ijkm1) / (dz_km1 + dz_k); // dT/dz
            uzz = 2 * u_ijkm1 / (dz_km1 * (dz_km1 + dz_k)) - 2 * u_ijk / (dz_km1 * dz_k) + 2 * u_ijkp1 / (dz_k * (dz_km1 + dz_k)); // d^2u/dz^2
            vzz = 2 * v_ijkm1 / (dz_km1 * (dz_km1 + dz_k)) - 2 * v_ijk / (dz_km1 * dz_k) + 2 * v_ijkp1 / (dz_k * (dz_km1 + dz_k)); // d^2v/dz^2
            wzz = 2 * w_ijkm1 / (dz_km1 * (dz_km1 + dz_k)) - 2 * w_ijk / (dz_km1 * dz_k) + 2 * w_ijkp1 / (dz_k * (dz_km1 + dz_k)); // d^2w/dz^2
            Tzz = 2 * T_ijkm1 / (dz_km1 * (dz_km1 + dz_k)) - 2 * T_ijk / (dz_km1 * dz_k) + 2 * T_ijkp1 / (dz_k * (dz_km1 + dz_k)); // d^2T/dz^2
        }
        /* Turbulence components */
        nu = mu / rho_ijk;
        // Damping function
        // tau_w = 0.0;
        // if (k == 0) {
        //     tau_w = sqrt((mu * 0.5 * (uz + wx)) * (mu * 0.5 * (uz + wx)) + (mu * 0.5 * (vz + wy)) * (mu * 0.5 * (vz + wy)));
        // }
        /*
        tau_w = sqrt((mu * 0.5 * (uz + wx)) * (mu * 0.5 * (uz + wx)) + (mu * 0.5 * (vz + wy)) * (mu * 0.5 * (vz + wy)));
        u_tau = sqrt(pow(tau_w / rho_ijk, 2.0));
        if (z_ibm[IDX(i, j, k, Nx, Ny, Nz)] * u_tau / nu > 5.0) {
            u_tau = 0.0;
        } */
        tau_w = 0.0;
        if (k == 0)
            tau_w = sqrt((mu * 0.5 * (uz + wx)) * (mu * 0.5 * (uz + wx)) + (mu * 0.5 * (vz + wy)) * (mu * 0.5 * (vz + wy)));
        u_tau = sqrt(pow(tau_w / rho_ijk, 2.0)); 
        fw = 1 - exp(-z_ibm[IDX(i, j, k, Nx, Ny, Nz)] * u_tau / 26 / nu);
        // Modified strain rate tensor S'
        S_11 = 2 * ux / 3 - (vy + wz) / 3;
        S_22 = 2 * vy / 3 - (uz + wx) / 3;
        S_33 = 2 * wz / 3 - (ux + vy) / 3;
        S_12 = 0.5 * (uy + vx);
        S_13 = 0.5 * (uz + wx);
        S_23 = 0.5 * (wy + vz);
        S_21 = S_12;
        S_31 = S_13;
        S_32 = S_23;
        // |S'| = |2 S_ij S_ij - 1/3 S_kk S_ij|
        mod_S = sqrt(
            2.0 * (ux * ux + vy * vy + wz * wz) 
            + (uz + wx) * (uz + wx) + (vx + uy) * (vx + uy) + (wy + vz) * (wy + vz) 
            - 2/3 * (ux + vy + wz) * (ux + vy + wz)
        );
        // Delta and l
        if (k < Nz - 1)
            dz_k = z[k + 1] - z[k];
        else
            dz_k = z[k] - z[k - 1];
        // Turbulent viscosity
        mu_sgs = rho_ijk * pow(C_s * pow(dx * dy * dz_k, 1.0 / 3.0) * fw, 2.0) * mod_S; 
        // if (isnan(fw)) {
        //     printf("NAN detected in fw at i=%d, j=%d, k=%d\n", i, j, k);
        //     printf("nu: %f\n", nu);
        //     printf("fw: %f\n", fw);
        //     printf("u_tau: %f\n", u_tau);
        //     printf("tau_w: %f\n", tau_w);
        //     break;
        // }
        /* Compute fuel and source term */
        H_step = (T_ijk > T_pc) ? 1.0 : 0.0;
        K_T = A * exp(-T_act / T_ijk);
        // Convection out of the solid zone and the next gas zone cancels the reaction because T_gas = T_solid = T_ijk and because there is no solid fuel
        T_gas = T_ijk;
        T_solid = T_ijk;
        // Solid fuel zone and gas zone
        if (k <= Nz_Y[IDX(i, j, 0, Nx, Ny, 1)]) {
            // Solid fuel zone
            if (k < Nz_Y[IDX(i, j, 0, Nx, Ny, 1)]) {
                Y_ijk = R_old[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y_max)];          
                Y_RHS = -Y_f * K_T * H_step * Y_ijk;
                R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y_max)] = Y_RHS;
            } 
            // Gas zone next to the solid zone (k+1) or inside the solid zone but with no solid fuel
            if (Y_ijk > 0.0) {
                T_gas = T_pc + (T_ijk - T_pc);
            } else {
                T_gas = T_ijk;
            }
        }
        // if (k < Nz_Y[IDX(i, j, 0, Nx, Ny, 1)]) {
        //     Y_ijk = R_old[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y_max)];        
        //     Y_RHS = -Y_f * K_T * H_step * Y_ijk;
        //     R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y_max)] = Y_RHS;
        // }
        // Compute source and force terms
        h_c = h_c * dz_k * pow((T_gas - T_solid) / dz_k, 1.0 / 4.0); // Old FDS heat transfer coefficient
        // q = H_R * Y_ijk * K_T * H_step / c_p - h_c * alpha_s * sigma_s * (T_ijk - T_inf) / (c_p * rho_inf);
        q = H_C * Y_ijk * K_T * H_step / c_p - h_c * alpha_s * sigma_s * (T_gas - T_solid) / (c_p * T_inf * rho_inf / T_gas);
        mod_U = sqrt(u_ijk * u_ijk + v_ijk * v_ijk + w_ijk * w_ijk);
        // Force terms
        // F_x = - Y_D * a_v * Y_ijk * mod_U * u_ijk;
        // F_y = - Y_D * a_v * Y_ijk * mod_U * v_ijk;
        // F_z = - g * (T_ijk - T_inf) / (T_ijk + 1e-16) - Y_D * a_v * Y_ijk * mod_U * w_ijk;
        F_x = - 0.5 * C_d * alpha_s * sigma_s * Y_ijk * mod_U * u_ijk;
        F_y = - 0.5 * C_d * alpha_s * sigma_s * Y_ijk * mod_U * v_ijk;
        F_z = g * (rho_ijk - rho_inf) / rho_ijk - 0.5 * C_d * alpha_s * sigma_s * Y_ijk * mod_U * w_ijk;
        // Compute Laplacian terms
        lap_u = uxx + uyy + uzz;
        lap_v = vxx + vyy + vzz;
        lap_w = wxx + wyy + wzz;
        lap_T = Txx + Tyy + Tzz;
        // Compute RHS        
        u_RHS = nu * lap_u - (uux + vuy + wuz) + F_x;
        v_RHS = nu * lap_v - (uvx + vvy + wvz) + F_y;
        w_RHS = nu * lap_w - (uwx + vwy + wwz) + F_z;
        // T_RHS = alpha * lap_T - (u_ijk * Tx + v_ijk * Ty + w_ijk * Tz) + S;
        // Temperature RHS with radiation and conduction
        T_RHS_tmp = (12 * SIGMA * delta * pow(T_ijk, 2) * (Tx * Tx + Ty * Ty + Tz * Tz) + (kappa + 4 * SIGMA * delta * pow(T_ijk, 3)) * lap_T) / (rho_ijk * c_p) + q;
        T_RHS = T_RHS_tmp - (u_ijk * Tx + v_ijk * Ty + w_ijk * Tz);
        // Copy to turbulence array. These calculations will be used in the next step to include the turbulence model
        R_turbulence[parameters.turbulence_indexes.rho + IDX(i, j, k, Nx, Ny, Nz)] = rho_ijk;
        R_turbulence[parameters.turbulence_indexes.Tx  + IDX(i, j, k, Nx, Ny, Nz)] = Tx;
        R_turbulence[parameters.turbulence_indexes.Ty  + IDX(i, j, k, Nx, Ny, Nz)] = Ty;
        R_turbulence[parameters.turbulence_indexes.Tz  + IDX(i, j, k, Nx, Ny, Nz)] = Tz;
        R_turbulence[parameters.turbulence_indexes.Txx + IDX(i, j, k, Nx, Ny, Nz)] = Txx;
        R_turbulence[parameters.turbulence_indexes.Tyy + IDX(i, j, k, Nx, Ny, Nz)] = Tyy;
        R_turbulence[parameters.turbulence_indexes.Tzz + IDX(i, j, k, Nx, Ny, Nz)] = Tzz;
        R_turbulence[parameters.turbulence_indexes.mu_sgs + IDX(i, j, k, Nx, Ny, Nz)] = mu_sgs;
        R_turbulence[parameters.turbulence_indexes.S_11 + IDX(i, j, k, Nx, Ny, Nz)] = S_11;
        R_turbulence[parameters.turbulence_indexes.S_12 + IDX(i, j, k, Nx, Ny, Nz)] = S_12;
        R_turbulence[parameters.turbulence_indexes.S_13 + IDX(i, j, k, Nx, Ny, Nz)] = S_13;
        R_turbulence[parameters.turbulence_indexes.S_21 + IDX(i, j, k, Nx, Ny, Nz)] = S_21;
        R_turbulence[parameters.turbulence_indexes.S_22 + IDX(i, j, k, Nx, Ny, Nz)] = S_22;
        R_turbulence[parameters.turbulence_indexes.S_23 + IDX(i, j, k, Nx, Ny, Nz)] = S_23;
        R_turbulence[parameters.turbulence_indexes.S_31 + IDX(i, j, k, Nx, Ny, Nz)] = S_31;
        R_turbulence[parameters.turbulence_indexes.S_32 + IDX(i, j, k, Nx, Ny, Nz)] = S_32;
        R_turbulence[parameters.turbulence_indexes.S_33 + IDX(i, j, k, Nx, Ny, Nz)] = S_33;
        R_turbulence[parameters.turbulence_indexes.div_U + IDX(i, j, k, Nx, Ny, Nz)] = nu * T_RHS_tmp / T_ijk;
        // Save RHS into R_new
        R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] = u_RHS;
        R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] = v_RHS;
        R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] = w_RHS;
        R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = T_RHS;  
    }
}

__global__
void velocity_correction(double *R_new, double *R_old, double *p, double *z, int fd_z, Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int u_index = parameters.field_indexes.u;
    int v_index = parameters.field_indexes.v;
    int w_index = parameters.field_indexes.w;
    int T_index = parameters.field_indexes.T;
    double rho_inf = parameters.rho_inf;
    double T_inf = parameters.T_inf;
    double dx = parameters.dx;
    double dy = parameters.dy;
    // double dz = parameters.dz;
    double dt = parameters.dt;
    double p_ijk, p_im1jk, p_ip1jk, p_ijm1k, p_ijp1k, p_ijkp1, p_ijkm1, p_ijkp2, p_ijkm2, rho_ijk;
    double px, py, pz;
    double dz_km2, dz_km1, dz_k, dz_kp1;
    int im1, ip1, jm1, jp1;
    int size = Nx * Ny * Nz;
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;
    for (int ijk = idx; ijk < size; ijk += stride) {
        int i = ijk / (Ny * Nz);
        int j = (ijk % (Ny * Nz)) / Nz;
        int k = ijk % Nz;
        // Indexes for periodic boundary conditions
        im1 = (i - 1 + Nx - 1) % (Nx - 1); // i-1
        jm1 = (j - 1 + Ny - 1) % (Ny - 1); // j-1
        ip1 = (i + 1) % (Nx - 1); // i+1
        jp1 = (j + 1) % (Ny - 1); // j+1
        p_ijk = p[IDX(i, j, k, Nx, Ny, Nz)]; // p_{i,j,k}
        rho_ijk = T_inf * rho_inf / R_old[T_index + IDX(i, j, k, Nx, Ny, Nz)]; // \rho_{i,j,k}
        p_im1jk = p[IDX(im1, j, k, Nx, Ny, Nz)]; // p_{i-1,j,k}
        p_ip1jk = p[IDX(ip1, j, k, Nx, Ny, Nz)]; // p_{i+1,j,k}
        p_ijm1k = p[IDX(i, jm1, k, Nx, Ny, Nz)]; // p_{i,j-1,k}
        p_ijp1k = p[IDX(i, jp1, k, Nx, Ny, Nz)]; // p_{i,j+1,k}
        // Central difference in x and y directions
        px = (p_ip1jk - p_im1jk) / (2 * dx);
        py = (p_ijp1k - p_ijm1k) / (2 * dy);
        // Depending on the finite difference scheme used for the z direction
        if (fd_z == 0) {  // Central difference 
            if (k == 0) { // Second order forward difference at the bottom boundary
                dz_k = z[k + 1] - z[k];
                dz_kp1 = z[k + 2] - z[k + 1];
                p_ijkp1 = p[IDX(i, j, k + 1, Nx, Ny, Nz)]; // p_{i,j,k+1}
                p_ijkp2 = p[IDX(i, j, k + 2, Nx, Ny, Nz)]; // p_{i,j,k+2}
                // pz = (-3 * p_ijk + 4 * p_ijkp1 - p_ijkp2) / (2 * dz); // Equispaced grid
                pz = - (2 * dz_k + dz_kp1) * p_ijk / (dz_k * (dz_k + dz_kp1)) 
                    + (dz_k + dz_kp1) * p_ijkp1 / (dz_k * dz_kp1) 
                    - dz_k * p_ijkp2 / (dz_kp1 * (dz_k + dz_kp1)); // Non equispaced grid
            } else if (k == Nz - 1) { // Second order backward difference at the top boundary
                dz_km1 = z[k] - z[k - 1];
                dz_km2 = z[k - 1] - z[k - 2];
                p_ijkm1 = p[IDX(i, j, k - 1, Nx, Ny, Nz)]; // p_{i,j,k-1}
                p_ijkm2 = p[IDX(i, j, k - 2, Nx, Ny, Nz)]; // p_{i,j,k-2}
                //pz = (3 * p_ijk - 4 * p_ijkm1 + p_ijkm2) / (2 * dz); // Equispaced grid
                pz = dz_km1 * p_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) 
                    - (dz_km2 + dz_km1) * p_ijkm1 / (dz_km2 * dz_km1) 
                    + (dz_km2 + 2 * dz_km1) * p_ijk / (dz_km1 * (dz_km2 + dz_km1)); // Non-equispaced grid
            } else { // Inside the domain
                dz_km1 = z[k] - z[k - 1];
                dz_k = z[k + 1] - z[k];
                p_ijkp1 = p[IDX(i, j, k + 1, Nx, Ny, Nz)]; // p_{i,j,k+1}
                p_ijkm1 = p[IDX(i, j, k - 1, Nx, Ny, Nz)]; // p_{i,j,k-1}
                // pz = (p_ijkp1 - p_ijkm1) / (2 * dz); // Equispaced grid
                pz = (p_ijkp1 - p_ijkm1) / (dz_km1 + dz_k); // Non equispaced grid
            }
        } else if (fd_z == 1) { // Forward difference
            if (k < Nz - 2) { // Second order forward difference for all nodes except the 2 last nodes
                dz_k = z[k + 1] - z[k];
                dz_kp1 = z[k + 2] - z[k + 1];
                p_ijkp1 = p[IDX(i, j, k + 1, Nx, Ny, Nz)]; // p_{i,j,k+1}
                p_ijkp2 = p[IDX(i, j, k + 2, Nx, Ny, Nz)]; // p_{i,j,k+2}
                // pz = (-3 * p_ijk + 4 * p_ijkp1 - p_ijkp2) / (2 * dz); // Equispaced grid
                pz = - (2 * dz_k + dz_kp1) * p_ijk / (dz_k * (dz_k + dz_kp1)) 
                    + (dz_k + dz_kp1) * p_ijkp1 / (dz_k * dz_kp1) 
                    - dz_k * p_ijkp2 / (dz_kp1 * (dz_k + dz_kp1)); // Non equispaced grid
            } else { // Second order backward difference for the 2 last nodes
                dz_km1 = z[k] - z[k - 1];
                dz_km2 = z[k - 1] - z[k - 2];
                p_ijkm1 = p[IDX(i, j, k - 1, Nx, Ny, Nz)]; // p_{i,j,k-1}
                p_ijkm2 = p[IDX(i, j, k - 2, Nx, Ny, Nz)]; // p_{i,j,k-2}
                // pz = (3 * p_ijk - 4 * p_ijkm1 + p_ijkm2) / (2 * dz); // Equispaced grid
                pz = dz_km1 * p_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) 
                    - (dz_km2 + dz_km1) * p_ijkm1 / (dz_km2 * dz_km1) 
                    + (dz_km2 + 2 * dz_km1) * p_ijk / (dz_km1 * (dz_km2 + dz_km1)); // Non equispaced grid
            }
        } else if (fd_z == -1) { // Backward difference
            if (k < 2) { // Second order forward difference for the 2 first nodes
                dz_k = z[k + 1] - z[k];
                dz_kp1 = z[k + 2] - z[k + 1];
                p_ijkp1 = p[IDX(i, j, k + 1, Nx, Ny, Nz)]; // p_{i,j,k+1}
                p_ijkp2 = p[IDX(i, j, k + 2, Nx, Ny, Nz)]; // p_{i,j,k+2}
                // pz = (-3 * p_ijk + 4 * p_ijkp1 - p_ijkp2) / (2 * dz); // Equispaced grid
                pz = - (2 * dz_k + dz_kp1) * p_ijk / (dz_k * (dz_k + dz_kp1)) 
                    + (dz_k + dz_kp1) * p_ijkp1 / (dz_k * dz_kp1) 
                    - dz_k * p_ijkp2 / (dz_kp1 * (dz_k + dz_kp1)); // Non equispaced grid
            } else { // Second order backward difference for all nodes except the 2 first nodes
                dz_km1 = z[k] - z[k - 1];
                dz_km2 = z[k - 1] - z[k - 2];
                p_ijkm1 = p[IDX(i, j, k - 1, Nx, Ny, Nz)]; // p_{i,j,k-1}
                p_ijkm2 = p[IDX(i, j, k - 2, Nx, Ny, Nz)]; // p_{i,j,k-2}
                // pz = (3 * p_ijk - 4 * p_ijkm1 + p_ijkm2) / (2 * dz); // Equispaced grid
                pz = dz_km1 * p_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) 
                    - (dz_km2 + dz_km1) * p_ijkm1 / (dz_km2 * dz_km1) 
                    + (dz_km2 + 2 * dz_km1) * p_ijk / (dz_km1 * (dz_km2 + dz_km1)); // Non equispaced grid
            }
        }
        // Get density
        if (parameters.variable_density == 0) {
            rho_ijk = rho_inf;
        } 
        // Velocity correction
        R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt * px / rho_ijk;
        R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt * py / rho_ijk;
        R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt * pz / rho_ijk;
    }
}

void Phi(double t, double *R_old, double *R_new, double *R_turbulence, double *z, double *z_ibm, int *Nz_Y, int *cut_nodes, Parameters parameters) {
    RHS<<<BLOCKS, THREADS>>>(t, R_old, R_new, R_turbulence, z, z_ibm, Nz_Y, parameters);
    checkCuda(cudaGetLastError());
    turbulence<<<BLOCKS, THREADS>>>(R_turbulence, R_new, z, parameters);
    checkCuda(cudaGetLastError());
    boundary_conditions<<<BLOCKS, THREADS>>>(R_new, z, Nz_Y, cut_nodes, parameters);
    checkCuda(cudaGetLastError());    
}