/**
 * @file turbulence.cu
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Functions for adding turbulence to the wildfire simulation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "../../include/cu/turbulence.cuh"

__global__ 
void turbulence(double *R_turbulence, double *R_new, double *z, double *z_ibm, Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int size = Nx * Ny * Nz;
    double dx = parameters.dx;
    double dy = parameters.dy;
    // double dz = parameters.dz;
    // double Delta = pow(dx * dy * dz, 1.0 / 3.0);
    // double Delta;
    // double C_s = parameters.C_s;
    double mu = parameters.mu;
    double nu;
    double Pr = parameters.Pr;
    double Tx, Ty, Tz, Txx, Tyy, Tzz;
    double rho;
    double mu_sgs_ijk, mu_sgs_ip1jk, mu_sgs_im1jk, mu_sgs_ijp1k, mu_sgs_ijm1k, mu_sgs_ijkp1, mu_sgs_ijkm1;
    double mu_sgs_ijkp2, mu_sgs_ijkm2;
    double S_11_ijk, S_12_ijk, S_13_ijk, S_21_ijk, S_22_ijk, S_23_ijk, S_31_ijk, S_32_ijk, S_33_ijk;
    double S_11_ip1jk, S_11_im1jk, S_12_ip1jk, S_12_im1jk, S_13_ip1jk, S_13_im1jk;
    double S_21_ijp1k, S_21_ijm1k, S_22_ijp1k, S_22_ijm1k, S_23_ijp1k, S_23_ijm1k;
    double S_31_ijkp1, S_31_ijkm1, S_31_ijkp2, S_31_ijkm2, S_32_ijkp1, S_32_ijkm1, S_32_ijkp2, S_32_ijkm2, S_33_ijkp1, S_33_ijkm1, S_33_ijkp2, S_33_ijkm2;
    double div_U_ijk, div_U_ip1jk, div_U_im1jk, div_U_ijp1k, div_U_ijm1k, div_U_ijkp1, div_U_ijkm1, div_U_ijkp2, div_U_ijkm2;
    double mu_sgs_x, mu_sgs_y, mu_sgs_z;
    double S_11_x, S_12_x, S_13_x, S_21_y, S_22_y, S_23_y, S_31_z, S_32_z, S_33_z;
    double div_S_x, div_S_y, div_S_z;
    double sgs_x, sgs_y, sgs_z, sgs_q;
    double dz_km2, dz_km1, dz_k, dz_kp1;
    double div_U_x, div_U_y, div_U_z;
    double u_tau_wall_ij0, u_tau_wall_ip1j0, u_tau_wall_im1j0, u_tau_wall_ijp10, u_tau_wall_ijm10;
    double fw_ijk, fw_im1jk, fw_ip1jk, fw_ijm1k, fw_ijp1k, fw_ijkm1, fw_ijkp1, fw_ijkm2, fw_ijkp2;
    int im1, ip1, jm1, jp1;
    // Loop over nodes
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;
    for (int ijk = idx; ijk < size; ijk += stride) {
        int i = ijk / (Ny * Nz);
        int j = (ijk % (Ny * Nz)) / Nz;
        int k = ijk % Nz;
        // Indexes for periodic boundary conditions
        im1 = (i - 1 + Nx - 1) % (Nx - 1);
        jm1 = (j - 1 + Ny - 1) % (Ny - 1);
        ip1 = (i + 1) % (Nx - 1);
        jp1 = (j + 1) % (Ny - 1);
        // Velocity gradients from R
        rho = R_turbulence[parameters.turbulence_indexes.rho + IDX(i, j, k, Nx, Ny, Nz)];
        Tx = R_turbulence[parameters.turbulence_indexes.Tx + IDX(i, j, k, Nx, Ny, Nz)];
        Ty = R_turbulence[parameters.turbulence_indexes.Ty + IDX(i, j, k, Nx, Ny, Nz)];
        Tz = R_turbulence[parameters.turbulence_indexes.Tz + IDX(i, j, k, Nx, Ny, Nz)];
        Txx = R_turbulence[parameters.turbulence_indexes.Txx + IDX(i, j, k, Nx, Ny, Nz)];
        Tyy = R_turbulence[parameters.turbulence_indexes.Tyy + IDX(i, j, k, Nx, Ny, Nz)];
        Tzz = R_turbulence[parameters.turbulence_indexes.Tzz + IDX(i, j, k, Nx, Ny, Nz)];
        mu_sgs_ijk = R_turbulence[parameters.turbulence_indexes.mu_sgs + IDX(i, j, k, Nx, Ny, Nz)];
        S_11_ijk = R_turbulence[parameters.turbulence_indexes.S_11 + IDX(i, j, k, Nx, Ny, Nz)];
        S_12_ijk = R_turbulence[parameters.turbulence_indexes.S_12 + IDX(i, j, k, Nx, Ny, Nz)];
        S_13_ijk = R_turbulence[parameters.turbulence_indexes.S_13 + IDX(i, j, k, Nx, Ny, Nz)];
        S_21_ijk = R_turbulence[parameters.turbulence_indexes.S_21 + IDX(i, j, k, Nx, Ny, Nz)];
        S_22_ijk = R_turbulence[parameters.turbulence_indexes.S_22 + IDX(i, j, k, Nx, Ny, Nz)];
        S_23_ijk = R_turbulence[parameters.turbulence_indexes.S_23 + IDX(i, j, k, Nx, Ny, Nz)];
        S_31_ijk = R_turbulence[parameters.turbulence_indexes.S_31 + IDX(i, j, k, Nx, Ny, Nz)];
        S_32_ijk = R_turbulence[parameters.turbulence_indexes.S_32 + IDX(i, j, k, Nx, Ny, Nz)];
        S_33_ijk = R_turbulence[parameters.turbulence_indexes.S_33 + IDX(i, j, k, Nx, Ny, Nz)];
        mu_sgs_ip1jk = R_turbulence[parameters.turbulence_indexes.mu_sgs + IDX(ip1, j, k, Nx, Ny, Nz)];
        mu_sgs_im1jk = R_turbulence[parameters.turbulence_indexes.mu_sgs + IDX(im1, j, k, Nx, Ny, Nz)];
        mu_sgs_ijp1k = R_turbulence[parameters.turbulence_indexes.mu_sgs + IDX(i, jp1, k, Nx, Ny, Nz)];
        mu_sgs_ijm1k = R_turbulence[parameters.turbulence_indexes.mu_sgs + IDX(i, jm1, k, Nx, Ny, Nz)];
        S_11_ip1jk = R_turbulence[parameters.turbulence_indexes.S_11 + IDX(ip1, j, k, Nx, Ny, Nz)];
        S_11_im1jk = R_turbulence[parameters.turbulence_indexes.S_11 + IDX(im1, j, k, Nx, Ny, Nz)];
        S_12_ip1jk = R_turbulence[parameters.turbulence_indexes.S_12 + IDX(ip1, j, k, Nx, Ny, Nz)];
        S_12_im1jk = R_turbulence[parameters.turbulence_indexes.S_12 + IDX(im1, j, k, Nx, Ny, Nz)];
        S_13_ip1jk = R_turbulence[parameters.turbulence_indexes.S_13 + IDX(ip1, j, k, Nx, Ny, Nz)];
        S_13_im1jk = R_turbulence[parameters.turbulence_indexes.S_13 + IDX(im1, j, k, Nx, Ny, Nz)];
        S_21_ijp1k = R_turbulence[parameters.turbulence_indexes.S_21 + IDX(i, jp1, k, Nx, Ny, Nz)];
        S_21_ijm1k = R_turbulence[parameters.turbulence_indexes.S_21 + IDX(i, jm1, k, Nx, Ny, Nz)];
        S_22_ijp1k = R_turbulence[parameters.turbulence_indexes.S_22 + IDX(i, jp1, k, Nx, Ny, Nz)];
        S_22_ijm1k = R_turbulence[parameters.turbulence_indexes.S_22 + IDX(i, jm1, k, Nx, Ny, Nz)];
        S_23_ijp1k = R_turbulence[parameters.turbulence_indexes.S_23 + IDX(i, jp1, k, Nx, Ny, Nz)];
        S_23_ijm1k = R_turbulence[parameters.turbulence_indexes.S_23 + IDX(i, jm1, k, Nx, Ny, Nz)];
        div_U_ip1jk = R_turbulence[parameters.turbulence_indexes.div_U + IDX(ip1, j, k, Nx, Ny, Nz)];
        div_U_im1jk = R_turbulence[parameters.turbulence_indexes.div_U + IDX(im1, j, k, Nx, Ny, Nz)];
        div_U_ijp1k = R_turbulence[parameters.turbulence_indexes.div_U + IDX(i, jp1, k, Nx, Ny, Nz)];
        div_U_ijm1k = R_turbulence[parameters.turbulence_indexes.div_U + IDX(i, jm1, k, Nx, Ny, Nz)];
        u_tau_wall_ij0 = R_turbulence[parameters.turbulence_indexes.u_tau_wall + IDX(i, j, 0, Nx, Ny, 1)];
        u_tau_wall_ip1j0 = R_turbulence[parameters.turbulence_indexes.u_tau_wall + IDX(ip1, j, 0, Nx, Ny, 1)];
        u_tau_wall_im1j0 = R_turbulence[parameters.turbulence_indexes.u_tau_wall + IDX(im1, j, 0, Nx, Ny, 1)];
        u_tau_wall_ijp10 = R_turbulence[parameters.turbulence_indexes.u_tau_wall + IDX(i, jp1, 0, Nx, Ny, 1)];
        u_tau_wall_ijm10 = R_turbulence[parameters.turbulence_indexes.u_tau_wall + IDX(i, jm1, 0, Nx, Ny, 1)];
        mu_sgs_ijkm1 = 0.0, S_31_ijkm1 = 0.0, S_32_ijkm1 = 0.0, S_33_ijkm1 = 0.0, div_U_ijkm1 = 0.0;
        mu_sgs_ijkp1 = 0.0, S_31_ijkp1 = 0.0, S_32_ijkp1 = 0.0, S_33_ijkp1 = 0.0, div_U_ijkp1 = 0.0;
        nu = mu / rho;
        fw_ijk = 1 - exp(-z_ibm[IDX(i, j, k, Nx, Ny, Nz)] * u_tau_wall_ij0 / A_VD / nu);
        fw_ip1jk = 1 - exp(-z_ibm[IDX(ip1, j, k, Nx, Ny, Nz)] * u_tau_wall_ip1j0 / A_VD / nu);
        fw_im1jk = 1 - exp(-z_ibm[IDX(im1, j, k, Nx, Ny, Nz)] * u_tau_wall_im1j0 / A_VD / nu);
        fw_ijp1k = 1 - exp(-z_ibm[IDX(i, jp1, k, Nx, Ny, Nz)] * u_tau_wall_ijp10 / A_VD / nu);
        fw_ijm1k = 1 - exp(-z_ibm[IDX(i, jm1, k, Nx, Ny, Nz)] * u_tau_wall_ijm10 / A_VD / nu);
        // I need to scale mu_sgs by fw^2
        mu_sgs_ijk *= fw_ijk * fw_ijk;
        mu_sgs_im1jk *= fw_im1jk * fw_im1jk;
        mu_sgs_ip1jk *= fw_ip1jk * fw_ip1jk;
        mu_sgs_ijm1k *= fw_ijm1k * fw_ijm1k;
        mu_sgs_ijp1k *= fw_ijp1k * fw_ijp1k;
        if (k > 0) {
            fw_ijkm1 = 1 - exp(-z_ibm[IDX(i, j, k - 1, Nx, Ny, Nz)] * u_tau_wall_ij0 / A_VD / nu);
            mu_sgs_ijkm1 = R_turbulence[parameters.turbulence_indexes.mu_sgs + IDX(i, j, k - 1, Nx, Ny, Nz)];
            S_31_ijkm1 = R_turbulence[parameters.turbulence_indexes.S_31 + IDX(i, j, k - 1, Nx, Ny, Nz)];
            S_32_ijkm1 = R_turbulence[parameters.turbulence_indexes.S_32 + IDX(i, j, k - 1, Nx, Ny, Nz)];
            S_33_ijkm1 = R_turbulence[parameters.turbulence_indexes.S_33 + IDX(i, j, k - 1, Nx, Ny, Nz)];
            div_U_ijkm1 = R_turbulence[parameters.turbulence_indexes.div_U + IDX(i, j, k - 1, Nx, Ny, Nz)];
            mu_sgs_ijkm1 *= fw_ijkm1 * fw_ijkm1;
        } 
        if (k < Nz - 1) {
            fw_ijkp1 = 1 - exp(-z_ibm[IDX(i, j, k + 1, Nx, Ny, Nz)] * u_tau_wall_ij0 / A_VD / nu);
            mu_sgs_ijkp1 = R_turbulence[parameters.turbulence_indexes.mu_sgs + IDX(i, j, k + 1, Nx, Ny, Nz)];
            S_31_ijkp1 = R_turbulence[parameters.turbulence_indexes.S_31 + IDX(i, j, k + 1, Nx, Ny, Nz)];
            S_32_ijkp1 = R_turbulence[parameters.turbulence_indexes.S_32 + IDX(i, j, k + 1, Nx, Ny, Nz)];
            S_33_ijkp1 = R_turbulence[parameters.turbulence_indexes.S_33 + IDX(i, j, k + 1, Nx, Ny, Nz)];
            div_U_ijkp1 = R_turbulence[parameters.turbulence_indexes.div_U + IDX(i, j, k + 1, Nx, Ny, Nz)];
            mu_sgs_ijkp1 *= fw_ijkp1 * fw_ijkp1;
        } 
        // Compute grad(mu_sgs)
        mu_sgs_x = (mu_sgs_ip1jk - mu_sgs_im1jk) / (2 * dx);
        mu_sgs_y = (mu_sgs_ijp1k - mu_sgs_ijm1k) / (2 * dy);
        // div(S_ij)
        S_11_x = (S_11_ip1jk - S_11_im1jk) / (2.0 * dx);
        S_12_x = (S_12_ip1jk - S_12_im1jk) / (2.0 * dx);
        S_13_x = (S_13_ip1jk - S_13_im1jk) / (2.0 * dx);
        S_21_y = (S_21_ijp1k - S_21_ijm1k) / (2.0 * dy);
        S_22_y = (S_22_ijp1k - S_22_ijm1k) / (2.0 * dy);
        S_23_y = (S_23_ijp1k - S_23_ijm1k) / (2.0 * dy);
        // grad(div(U))
        div_U_x = (div_U_ip1jk - div_U_im1jk) / (2.0 * dx);
        div_U_y = (div_U_ijp1k - div_U_ijm1k) / (2.0 * dy);
        if (k == 0) { // Second-order forward difference    
            dz_k = z[k + 1] - z[k];
            dz_kp1 = z[k + 2] - z[k + 1];
            fw_ijkp2 = 1 - exp(-z_ibm[IDX(i, j, k + 2, Nx, Ny, Nz)] * u_tau_wall_ij0 / A_VD / nu);
            mu_sgs_ijkp2 = R_turbulence[parameters.turbulence_indexes.mu_sgs + IDX(i, j, k + 2, Nx, Ny, Nz)];
            S_31_ijkp2 = R_turbulence[parameters.turbulence_indexes.S_31 + IDX(i, j, k + 2, Nx, Ny, Nz)];
            S_32_ijkp2 = R_turbulence[parameters.turbulence_indexes.S_32 + IDX(i, j, k + 2, Nx, Ny, Nz)];
            S_33_ijkp2 = R_turbulence[parameters.turbulence_indexes.S_33 + IDX(i, j, k + 2, Nx, Ny, Nz)];
            mu_sgs_ijkp2 *= fw_ijkp2 * fw_ijkp2;
            // Non-equispaced grid
            mu_sgs_z = - (2 * dz_k + dz_kp1) * mu_sgs_ijk / (dz_k * (dz_k + dz_kp1)) + (dz_k + dz_kp1) * mu_sgs_ijkp1 / (dz_k * dz_kp1) - dz_k * mu_sgs_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            S_31_z = - (2 * dz_k + dz_kp1) * S_31_ijk / (dz_k * (dz_k + dz_kp1)) + (dz_k + dz_kp1) * S_31_ijkp1 / (dz_k * dz_kp1) - dz_k * S_31_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            S_32_z = - (2 * dz_k + dz_kp1) * S_32_ijk / (dz_k * (dz_k + dz_kp1)) + (dz_k + dz_kp1) * S_32_ijkp1 / (dz_k * dz_kp1) - dz_k * S_32_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            S_33_z = - (2 * dz_k + dz_kp1) * S_33_ijk / (dz_k * (dz_k + dz_kp1)) + (dz_k + dz_kp1) * S_33_ijkp1 / (dz_k * dz_kp1) - dz_k * S_33_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            div_U_z = - (2 * dz_k + dz_kp1) * div_U_ijk / (dz_k * (dz_k + dz_kp1)) + (dz_k + dz_kp1) * div_U_ijkp1 / (dz_k * dz_kp1) - dz_k * div_U_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
        } else if (k == Nz - 1) { // Second-order backward difference
            dz_km1 = z[k] - z[k - 1];
            dz_km2 = z[k - 1] - z[k - 2];
            fw_ijkm2 = 1 - exp(-z_ibm[IDX(i, j, k - 2, Nx, Ny, Nz)] * u_tau_wall_ij0 / A_VD / nu);
            mu_sgs_ijkm2 = R_turbulence[parameters.turbulence_indexes.mu_sgs + IDX(i, j, k - 2, Nx, Ny, Nz)];
            S_31_ijkm2 = R_turbulence[parameters.turbulence_indexes.S_31 + IDX(i, j, k - 2, Nx, Ny, Nz)];
            S_32_ijkm2 = R_turbulence[parameters.turbulence_indexes.S_32 + IDX(i, j, k - 2, Nx, Ny, Nz)];
            S_33_ijkm2 = R_turbulence[parameters.turbulence_indexes.S_33 + IDX(i, j, k - 2, Nx, Ny, Nz)];
            mu_sgs_ijkm2 *= fw_ijkm2 * fw_ijkm2;
            // Non-equispaced grid
            mu_sgs_z = dz_km1 * mu_sgs_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) - (dz_km2 + dz_km1) * mu_sgs_ijkm1 / (dz_km2 * dz_km1) + (dz_km2 + 2 * dz_km1) * mu_sgs_ijk / (dz_km1 * (dz_km2 + dz_km1));
            S_31_z = dz_km1 * S_31_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) - (dz_km2 + dz_km1) * S_31_ijkm1 / (dz_km2 * dz_km1) + (dz_km2 + 2 * dz_km1) * S_31_ijk / (dz_km1 * (dz_km2 + dz_km1));
            S_32_z = dz_km1 * S_32_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) - (dz_km2 + dz_km1) * S_32_ijkm1 / (dz_km2 * dz_km1) + (dz_km2 + 2 * dz_km1) * S_32_ijk / (dz_km1 * (dz_km2 + dz_km1));
            S_33_z = dz_km1 * S_33_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) - (dz_km2 + dz_km1) * S_33_ijkm1 / (dz_km2 * dz_km1) + (dz_km2 + 2 * dz_km1) * S_33_ijk / (dz_km1 * (dz_km2 + dz_km1));
            div_U_z = dz_km1 * div_U_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) - (dz_km2 + dz_km1) * div_U_ijkm1 / (dz_km2 * dz_km1) + (dz_km2 + 2 * dz_km1) * div_U_ijk / (dz_km1 * (dz_km2 + dz_km1));
        } else {
            dz_km1 = z[k] - z[k - 1];
            dz_k = z[k + 1] - z[k];
            // Non-equispaced grid
            mu_sgs_z = (mu_sgs_ijkp1 - mu_sgs_ijkm1) / (dz_k + dz_km1);
            S_31_z = (S_31_ijkp1 - S_31_ijkm1) / (dz_k + dz_km1);
            S_32_z = (S_32_ijkp1 - S_32_ijkm1) / (dz_k + dz_km1);
            S_33_z = (S_33_ijkp1 - S_33_ijkm1) / (dz_k + dz_km1);
            div_U_z = (div_U_ijkp1 - div_U_ijkm1) / (dz_k + dz_km1);
        }
        // Divergence of S
        div_S_x = S_11_x + S_21_y + S_31_z;
        div_S_y = S_12_x + S_22_y + S_32_z;
        div_S_z = S_13_x + S_23_y + S_33_z;
        // Compute SGS model
        sgs_x = -2 * (
            mu_sgs_x * S_11_ijk + mu_sgs_y * S_21_ijk + mu_sgs_z * S_31_ijk + 
            mu_sgs_ijk * div_S_x 
        ) / rho;
        sgs_y = -2 * (
            mu_sgs_x * S_12_ijk + mu_sgs_y * S_22_ijk + mu_sgs_z * S_32_ijk +
            mu_sgs_ijk * div_S_y
        ) / rho;
        sgs_z = -2 * (
            mu_sgs_x * S_13_ijk + mu_sgs_y * S_23_ijk + mu_sgs_z * S_33_ijk +
            mu_sgs_ijk *  div_S_z
        ) / rho;
        // q^{sgs} = -c_p/Pr * (div(mu_sgs)\dot\grad(T) + mu_sgs * \nabla^2(T))
        // c_p is ommited because it divides the whole term
        sgs_q = -(mu_sgs_x * Tx + mu_sgs_y * Ty + mu_sgs_z * Tz + mu_sgs_ijk * (Txx + Tyy + Tzz)) / (rho * Pr);
        // Add SGS model to R
        R_new[parameters.field_indexes.u + IDX(i, j, k, Nx, Ny, Nz)] += -sgs_x + div_U_x / 3.0;
        R_new[parameters.field_indexes.v + IDX(i, j, k, Nx, Ny, Nz)] += -sgs_y + div_U_y / 3.0;
        R_new[parameters.field_indexes.w + IDX(i, j, k, Nx, Ny, Nz)] += -sgs_z + div_U_z / 3.0;
        R_new[parameters.field_indexes.T + IDX(i, j, k, Nx, Ny, Nz)] -= sgs_q;
    }
}

/*
__global__ 
void turbulence_old(double *R_turbulence, double *R_new, double *z, Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int size = Nx * Ny * Nz;
    double dx = parameters.dx;
    double dy = parameters.dy;
    // double dz = parameters.dz;
    // double Delta = pow(dx * dy * dz, 1.0 / 3.0);
    double Delta;
    double C_s = parameters.C_s;
    double Pr = parameters.Pr;
    // double l = C_s * Delta;
    double l;
    double ux, uy, uz, vx, vy, vz, wx, wy, wz, Tx, Ty, Tz;
    double uxx, uyy, uzz, vxx, vyy, vzz, wxx, wyy, wzz, Txx, Tyy, Tzz;
    double uyx, uzx, uxy, uzy, uxz, uyz;
    double vyx, vzx, vxy, vzy, vxz, vyz;
    double wyx, wzx, wxy, wzy, wxz, wyz;
    double ux_ijp1k, ux_ijm1k, ux_ijkp1, ux_ijkm1, uy_ip1jk, uy_im1jk, uz_ip1jk, uz_im1jk, uz_ijp1k, uz_ijm1k, uy_ijkp1, uy_ijkm1;
    double vx_ijp1k, vx_ijm1k, vx_ijkp1, vx_ijkm1, vy_ip1jk, vy_im1jk, vy_ijkp1, vy_ijkm1, vz_ip1jk, vz_im1jk, vz_ijp1k, vz_ijm1k;
    double wx_ijp1k, wx_ijm1k, wx_ijkp1, wx_ijkm1, wy_ip1jk, wy_im1jk, wy_ijkp1, wy_ijkm1, wz_ip1jk, wz_im1jk, wz_ijp1k, wz_ijm1k;
    double ux_ijkp2, ux_ijkm2, uy_ijkp2, uy_ijkm2;
    double vx_ijkp2, vx_ijkm2, vy_ijkp2, vy_ijkm2;
    double wx_ijkp2, wx_ijkm2, wy_ijkp2, wy_ijkm2;
    double fw_ip1jk, fw_im1jk, fw_ijp1k, fw_ijm1k, fw_ijkp1, fw_ijkm1, fw_ijkp2, fw_ijkm2;
    double sgs_x, sgs_y, sgs_z, sgs_q;
    double mod_S, psi_x, psi_y, psi_z;
    double sgs_x_no_damp, sgs_y_no_damp, sgs_z_no_damp, sgs_q_no_damp;
    double sgs_x_damp, sgs_y_damp, sgs_z_damp, sgs_q_damp;
    double fw, fwx, fwy, fwz;
    double dz_km2, dz_km1, dz_k, dz_kp1;
    int im1, ip1, jm1, jp1;
    // Loop over nodes
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;
    for (int ijk = idx; ijk < size; ijk += stride) {
        int i = ijk / (Ny * Nz);
        int j = (ijk % (Ny * Nz)) / Nz;
        int k = ijk % Nz;
        // Indexes for periodic boundary conditions
        im1 = (i - 1 + Nx - 1) % (Nx - 1);
        jm1 = (j - 1 + Ny - 1) % (Ny - 1);
        ip1 = (i + 1) % (Nx - 1);
        jp1 = (j + 1) % (Ny - 1);
        // Velocity gradients from R
        ux = R_turbulence[parameters.turbulence_indexes.ux + IDX(i, j, k, Nx, Ny, Nz)];
        uy = R_turbulence[parameters.turbulence_indexes.uy + IDX(i, j, k, Nx, Ny, Nz)];
        uz = R_turbulence[parameters.turbulence_indexes.uz + IDX(i, j, k, Nx, Ny, Nz)];
        vx = R_turbulence[parameters.turbulence_indexes.vx + IDX(i, j, k, Nx, Ny, Nz)];
        vy = R_turbulence[parameters.turbulence_indexes.vy + IDX(i, j, k, Nx, Ny, Nz)];
        vz = R_turbulence[parameters.turbulence_indexes.vz + IDX(i, j, k, Nx, Ny, Nz)];
        wx = R_turbulence[parameters.turbulence_indexes.wx + IDX(i, j, k, Nx, Ny, Nz)];
        wy = R_turbulence[parameters.turbulence_indexes.wy + IDX(i, j, k, Nx, Ny, Nz)];
        wz = R_turbulence[parameters.turbulence_indexes.wz + IDX(i, j, k, Nx, Ny, Nz)];
        Tx = R_turbulence[parameters.turbulence_indexes.Tx + IDX(i, j, k, Nx, Ny, Nz)];
        Ty = R_turbulence[parameters.turbulence_indexes.Ty + IDX(i, j, k, Nx, Ny, Nz)];
        Tz = R_turbulence[parameters.turbulence_indexes.Tz + IDX(i, j, k, Nx, Ny, Nz)];
        uxx = R_turbulence[parameters.turbulence_indexes.uxx + IDX(i, j, k, Nx, Ny, Nz)];
        uyy = R_turbulence[parameters.turbulence_indexes.uyy + IDX(i, j, k, Nx, Ny, Nz)];
        uzz = R_turbulence[parameters.turbulence_indexes.uzz + IDX(i, j, k, Nx, Ny, Nz)];
        vxx = R_turbulence[parameters.turbulence_indexes.vxx + IDX(i, j, k, Nx, Ny, Nz)];
        vyy = R_turbulence[parameters.turbulence_indexes.vyy + IDX(i, j, k, Nx, Ny, Nz)];
        vzz = R_turbulence[parameters.turbulence_indexes.vzz + IDX(i, j, k, Nx, Ny, Nz)];
        wxx = R_turbulence[parameters.turbulence_indexes.wxx + IDX(i, j, k, Nx, Ny, Nz)];
        wyy = R_turbulence[parameters.turbulence_indexes.wyy + IDX(i, j, k, Nx, Ny, Nz)];
        wzz = R_turbulence[parameters.turbulence_indexes.wzz + IDX(i, j, k, Nx, Ny, Nz)];
        Txx = R_turbulence[parameters.turbulence_indexes.Txx + IDX(i, j, k, Nx, Ny, Nz)];
        Tyy = R_turbulence[parameters.turbulence_indexes.Tyy + IDX(i, j, k, Nx, Ny, Nz)];
        Tzz = R_turbulence[parameters.turbulence_indexes.Tzz + IDX(i, j, k, Nx, Ny, Nz)];
        fw  = R_turbulence[parameters.turbulence_indexes.fw + IDX(i, j, k, Nx, Ny, Nz)];
        // Get the values of the derivatives
        ux_ijp1k = R_turbulence[parameters.turbulence_indexes.ux + IDX(i, jp1, k, Nx, Ny, Nz)];
        ux_ijm1k = R_turbulence[parameters.turbulence_indexes.ux + IDX(i, jm1, k, Nx, Ny, Nz)];
        uy_ip1jk = R_turbulence[parameters.turbulence_indexes.uy + IDX(ip1, j, k, Nx, Ny, Nz)];
        uy_im1jk = R_turbulence[parameters.turbulence_indexes.uy + IDX(im1, j, k, Nx, Ny, Nz)];
        // uy_ijkp1 = R_turbulence[parameters.turbulence_indexes.uy + IDX(i, j, k + 1, Nx, Ny, Nz)];
        // uy_ijkm1 = R_turbulence[parameters.turbulence_indexes.uy + IDX(i, j, k - 1, Nx, Ny, Nz)];
        uz_ip1jk = R_turbulence[parameters.turbulence_indexes.uz + IDX(ip1, j, k, Nx, Ny, Nz)];
        uz_im1jk = R_turbulence[parameters.turbulence_indexes.uz + IDX(im1, j, k, Nx, Ny, Nz)];
        uz_ijp1k = R_turbulence[parameters.turbulence_indexes.uz + IDX(i, jp1, k, Nx, Ny, Nz)];
        uz_ijm1k = R_turbulence[parameters.turbulence_indexes.uz + IDX(i, jm1, k, Nx, Ny, Nz)];
        vx_ijp1k = R_turbulence[parameters.turbulence_indexes.vx + IDX(i, jp1, k, Nx, Ny, Nz)];
        vx_ijm1k = R_turbulence[parameters.turbulence_indexes.vx + IDX(i, jm1, k, Nx, Ny, Nz)];
        // vx_ijkp1 = R_turbulence[parameters.turbulence_indexes.vx + IDX(i, j, k + 1, Nx, Ny, Nz)];
        // vx_ijkm1 = R_turbulence[parameters.turbulence_indexes.vx + IDX(i, j, k - 1, Nx, Ny, Nz)];
        vy_ip1jk = R_turbulence[parameters.turbulence_indexes.vy + IDX(ip1, j, k, Nx, Ny, Nz)];
        vy_im1jk = R_turbulence[parameters.turbulence_indexes.vy + IDX(im1, j, k, Nx, Ny, Nz)];
        // vy_ijkp1 = R_turbulence[parameters.turbulence_indexes.vy + IDX(i, j, k + 1, Nx, Ny, Nz)];
        // vy_ijkm1 = R_turbulence[parameters.turbulence_indexes.vy + IDX(i, j, k - 1, Nx, Ny, Nz)];
        vz_ip1jk = R_turbulence[parameters.turbulence_indexes.vz + IDX(ip1, j, k, Nx, Ny, Nz)];
        vz_im1jk = R_turbulence[parameters.turbulence_indexes.vz + IDX(im1, j, k, Nx, Ny, Nz)];
        vz_ijp1k = R_turbulence[parameters.turbulence_indexes.vz + IDX(i, jp1, k, Nx, Ny, Nz)];
        vz_ijm1k = R_turbulence[parameters.turbulence_indexes.vz + IDX(i, jm1, k, Nx, Ny, Nz)];
        wx_ijp1k = R_turbulence[parameters.turbulence_indexes.wx + IDX(i, jp1, k, Nx, Ny, Nz)];
        wx_ijm1k = R_turbulence[parameters.turbulence_indexes.wx + IDX(i, jm1, k, Nx, Ny, Nz)];
        // wx_ijkp1 = R_turbulence[parameters.turbulence_indexes.wx + IDX(i, j, k + 1, Nx, Ny, Nz)];
        // wx_ijkm1 = R_turbulence[parameters.turbulence_indexes.wx + IDX(i, j, k - 1, Nx, Ny, Nz)];
        wy_ip1jk = R_turbulence[parameters.turbulence_indexes.wy + IDX(ip1, j, k, Nx, Ny, Nz)];
        wy_im1jk = R_turbulence[parameters.turbulence_indexes.wy + IDX(im1, j, k, Nx, Ny, Nz)];
        // wy_ijkp1 = R_turbulence[parameters.turbulence_indexes.wy + IDX(i, j, k + 1, Nx, Ny, Nz)];
        // wy_ijkm1 = R_turbulence[parameters.turbulence_indexes.wy + IDX(i, j, k - 1, Nx, Ny, Nz)];
        wz_ip1jk = R_turbulence[parameters.turbulence_indexes.wz + IDX(ip1, j, k, Nx, Ny, Nz)];
        wz_im1jk = R_turbulence[parameters.turbulence_indexes.wz + IDX(im1, j, k, Nx, Ny, Nz)];
        wz_ijp1k = R_turbulence[parameters.turbulence_indexes.wz + IDX(i, jp1, k, Nx, Ny, Nz)];
        wz_ijm1k = R_turbulence[parameters.turbulence_indexes.wz + IDX(i, jm1, k, Nx, Ny, Nz)];
        fw_ip1jk = R_turbulence[parameters.turbulence_indexes.fw + IDX(ip1, j, k, Nx, Ny, Nz)];
        fw_im1jk = R_turbulence[parameters.turbulence_indexes.fw + IDX(im1, j, k, Nx, Ny, Nz)];
        fw_ijp1k = R_turbulence[parameters.turbulence_indexes.fw + IDX(i, jp1, k, Nx, Ny, Nz)];
        fw_ijm1k = R_turbulence[parameters.turbulence_indexes.fw + IDX(i, jm1, k, Nx, Ny, Nz)];
        // fw_ijkp1 = R_turbulence[parameters.turbulence_indexes.fw + IDX(i, j, k + 1, Nx, Ny, Nz)];
        // fw_ijkm1 = R_turbulence[parameters.turbulence_indexes.fw + IDX(i, j, k - 1, Nx, Ny, Nz)];
        if (k > 0) {
            ux_ijkm1 = R_turbulence[parameters.turbulence_indexes.ux + IDX(i, j, k - 1, Nx, Ny, Nz)];
            uy_ijkm1 = R_turbulence[parameters.turbulence_indexes.uy + IDX(i, j, k - 1, Nx, Ny, Nz)];
            vx_ijkm1 = R_turbulence[parameters.turbulence_indexes.vx + IDX(i, j, k - 1, Nx, Ny, Nz)];
            vy_ijkm1 = R_turbulence[parameters.turbulence_indexes.vy + IDX(i, j, k - 1, Nx, Ny, Nz)];
            wx_ijkm1 = R_turbulence[parameters.turbulence_indexes.wx + IDX(i, j, k - 1, Nx, Ny, Nz)];
            wy_ijkm1 = R_turbulence[parameters.turbulence_indexes.wy + IDX(i, j, k - 1, Nx, Ny, Nz)];
            fw_ijkm1 = R_turbulence[parameters.turbulence_indexes.fw + IDX(i, j, k - 1, Nx, Ny, Nz)];

        } else {
            ux_ijkm1 = 0.0;
            uy_ijkm1 = 0.0;
            vx_ijkm1 = 0.0;
            vy_ijkm1 = 0.0;
            wx_ijkm1 = 0.0;
            wy_ijkm1 = 0.0;
            fw_ijkm1 = 0.0;
        }
        if (k < Nz - 1) {
            ux_ijkp1 = R_turbulence[parameters.turbulence_indexes.ux + IDX(i, j, k + 1, Nx, Ny, Nz)];
            uy_ijkp1 = R_turbulence[parameters.turbulence_indexes.uy + IDX(i, j, k + 1, Nx, Ny, Nz)];
            vx_ijkp1 = R_turbulence[parameters.turbulence_indexes.vx + IDX(i, j, k + 1, Nx, Ny, Nz)];
            vy_ijkp1 = R_turbulence[parameters.turbulence_indexes.vy + IDX(i, j, k + 1, Nx, Ny, Nz)];
            wx_ijkp1 = R_turbulence[parameters.turbulence_indexes.wx + IDX(i, j, k + 1, Nx, Ny, Nz)];
            wy_ijkp1 = R_turbulence[parameters.turbulence_indexes.wy + IDX(i, j, k + 1, Nx, Ny, Nz)];
            fw_ijkp1 = R_turbulence[parameters.turbulence_indexes.fw + IDX(i, j, k + 1, Nx, Ny, Nz)];
        } else {
            ux_ijkp1 = 0.0;
            uy_ijkp1 = 0.0;
            vx_ijkp1 = 0.0;
            vy_ijkp1 = 0.0;
            wx_ijkp1 = 0.0;
            wy_ijkp1 = 0.0;
            fw_ijkp1 = 0.0;
        }
        // Mixed partial derivatives d/dx and d/dy
        uyx = (uy_ip1jk - uy_im1jk) / (2.0 * dx);
        uzx = (uz_ip1jk - uz_im1jk) / (2.0 * dx);
        vyx = (vy_ip1jk - vy_im1jk) / (2.0 * dx);
        vzx = (vz_ip1jk - vz_im1jk) / (2.0 * dx);
        wyx = (wy_ip1jk - wy_im1jk) / (2.0 * dx);
        wzx = (wz_ip1jk - wz_im1jk) / (2.0 * dx);
        uzy = (uz_ijp1k - uz_ijm1k) / (2.0 * dy);
        uxy = (ux_ijp1k - ux_ijm1k) / (2.0 * dy);
        vzy = (vz_ijp1k - vz_ijm1k) / (2.0 * dy);
        vxy = (vx_ijp1k - vx_ijm1k) / (2.0 * dy);
        wzy = (wz_ijp1k - wz_ijm1k) / (2.0 * dy);
        wxy = (wx_ijp1k - wx_ijm1k) / (2.0 * dy);
        fwx = (fw_ip1jk - fw_im1jk) / (2 * dx);
        fwy = (fw_ijp1k - fw_ijm1k) / (2 * dy);
        if (k == 0) { // Second-order forward difference    
            dz_k = z[k + 1] - z[k];
            dz_kp1 = z[k + 2] - z[k + 1];  
            ux_ijkp2 = R_turbulence[parameters.turbulence_indexes.ux + IDX(i, j, k + 2, Nx, Ny, Nz)];
            uy_ijkp2 = R_turbulence[parameters.turbulence_indexes.uy + IDX(i, j, k + 2, Nx, Ny, Nz)];
            vx_ijkp2 = R_turbulence[parameters.turbulence_indexes.vx + IDX(i, j, k + 2, Nx, Ny, Nz)];
            vy_ijkp2 = R_turbulence[parameters.turbulence_indexes.vy + IDX(i, j, k + 2, Nx, Ny, Nz)];
            wx_ijkp2 = R_turbulence[parameters.turbulence_indexes.wx + IDX(i, j, k + 2, Nx, Ny, Nz)];
            wy_ijkp2 = R_turbulence[parameters.turbulence_indexes.wy + IDX(i, j, k + 2, Nx, Ny, Nz)];
            fw_ijkp2 = R_turbulence[parameters.turbulence_indexes.fw + IDX(i, j, k + 2, Nx, Ny, Nz)];
            // Equispaced grid
            // uxz = (-3 * ux + 4 * ux_ijkp1 - ux_ijkp2) / (2 * dz);
            // uyz = (-3 * uy + 4 * uy_ijkp1 - uy_ijkp2) / (2 * dz);
            // vxz = (-3 * vx + 4 * vx_ijkp1 - vx_ijkp2) / (2 * dz);
            // vyz = (-3 * vy + 4 * vy_ijkp1 - vy_ijkp2) / (2 * dz);
            // wxz = (-3 * wx + 4 * wx_ijkp1 - wx_ijkp2) / (2 * dz);
            // wyz = (-3 * wy + 4 * wy_ijkp1 - wy_ijkp2) / (2 * dz);
            // fwz = (-3 * fw + 4 * fw_ijkp1 - fw_ijkp2) / (2 * dz);
            // Non-equispaced grid
            uxz = - (2 * dz_k + dz_kp1) * ux / (dz_k * (dz_k + dz_kp1)) + (dz_k + dz_kp1) * ux_ijkp1 / (dz_k * dz_kp1) - dz_k * ux_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            uyz = - (2 * dz_k + dz_kp1) * uy / (dz_k * (dz_k + dz_kp1)) + (dz_k + dz_kp1) * uy_ijkp1 / (dz_k * dz_kp1) - dz_k * uy_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            vxz = - (2 * dz_k + dz_kp1) * vx / (dz_k * (dz_k + dz_kp1)) + (dz_k + dz_kp1) * vx_ijkp1 / (dz_k * dz_kp1) - dz_k * vx_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            vyz = - (2 * dz_k + dz_kp1) * vy / (dz_k * (dz_k + dz_kp1)) + (dz_k + dz_kp1) * vy_ijkp1 / (dz_k * dz_kp1) - dz_k * vy_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            wxz = - (2 * dz_k + dz_kp1) * wx / (dz_k * (dz_k + dz_kp1)) + (dz_k + dz_kp1) * wx_ijkp1 / (dz_k * dz_kp1) - dz_k * wx_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            wyz = - (2 * dz_k + dz_kp1) * wy / (dz_k * (dz_k + dz_kp1)) + (dz_k + dz_kp1) * wy_ijkp1 / (dz_k * dz_kp1) - dz_k * wy_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
            fwz = - (2 * dz_k + dz_kp1) * fw / (dz_k * (dz_k + dz_kp1)) + (dz_k + dz_kp1) * fw_ijkp1 / (dz_k * dz_kp1) - dz_k * fw_ijkp2 / (dz_kp1 * (dz_k + dz_kp1));
        } else if (k == Nz - 1) { // Second-order backward difference
            dz_km1 = z[k] - z[k - 1];
            dz_km2 = z[k - 1] - z[k - 2];
            ux_ijkm2 = R_turbulence[parameters.turbulence_indexes.ux + IDX(i, j, k - 2, Nx, Ny, Nz)];
            uy_ijkm2 = R_turbulence[parameters.turbulence_indexes.uy + IDX(i, j, k - 2, Nx, Ny, Nz)];
            vx_ijkm2 = R_turbulence[parameters.turbulence_indexes.vx + IDX(i, j, k - 2, Nx, Ny, Nz)];
            vy_ijkm2 = R_turbulence[parameters.turbulence_indexes.vy + IDX(i, j, k - 2, Nx, Ny, Nz)];
            wx_ijkm2 = R_turbulence[parameters.turbulence_indexes.wx + IDX(i, j, k - 2, Nx, Ny, Nz)];
            wy_ijkm2 = R_turbulence[parameters.turbulence_indexes.wy + IDX(i, j, k - 2, Nx, Ny, Nz)];
            fw_ijkm2 = R_turbulence[parameters.turbulence_indexes.fw + IDX(i, j, k - 2, Nx, Ny, Nz)];
            // Equispaced grid
            // uxz = (3 * ux - 4 * ux_ijkm1 + ux_ijkm2) / (2 * dz);
            // uyz = (3 * uy - 4 * uy_ijkm1 + uy_ijkm2) / (2 * dz);
            // vxz = (3 * vx - 4 * vx_ijkm1 + vx_ijkm2) / (2 * dz);
            // vyz = (3 * vy - 4 * vy_ijkm1 + vy_ijkm2) / (2 * dz);
            // wxz = (3 * wx - 4 * wx_ijkm1 + wx_ijkm2) / (2 * dz);
            // wyz = (3 * wy - 4 * wy_ijkm1 + wy_ijkm2) / (2 * dz);
            // fwz = (3 * fw - 4 * fw_ijkm1 + fw_ijkm2) / (2 * dz);
            // Non-equispaced grid
            uxz = dz_km1 * ux_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) - (dz_km2 + dz_km1) * ux_ijkm1 / (dz_km2 * dz_km1) + (dz_km2 + 2 * dz_km1) * ux / (dz_km1 * (dz_km2 + dz_km1));
            uyz = dz_km1 * uy_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) - (dz_km2 + dz_km1) * uy_ijkm1 / (dz_km2 * dz_km1) + (dz_km2 + 2 * dz_km1) * uy / (dz_km1 * (dz_km2 + dz_km1));
            vxz = dz_km1 * vx_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) - (dz_km2 + dz_km1) * vx_ijkm1 / (dz_km2 * dz_km1) + (dz_km2 + 2 * dz_km1) * vx / (dz_km1 * (dz_km2 + dz_km1));
            vyz = dz_km1 * vy_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) - (dz_km2 + dz_km1) * vy_ijkm1 / (dz_km2 * dz_km1) + (dz_km2 + 2 * dz_km1) * vy / (dz_km1 * (dz_km2 + dz_km1));
            wxz = dz_km1 * wx_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) - (dz_km2 + dz_km1) * wx_ijkm1 / (dz_km2 * dz_km1) + (dz_km2 + 2 * dz_km1) * wx / (dz_km1 * (dz_km2 + dz_km1));
            wyz = dz_km1 * wy_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) - (dz_km2 + dz_km1) * wy_ijkm1 / (dz_km2 * dz_km1) + (dz_km2 + 2 * dz_km1) * wy / (dz_km1 * (dz_km2 + dz_km1));
            fwz = dz_km1 * fw_ijkm2 / (dz_km2 * (dz_km2 + dz_km1)) - (dz_km2 + dz_km1) * fw_ijkm1 / (dz_km2 * dz_km1) + (dz_km2 + 2 * dz_km1) * fw / (dz_km1 * (dz_km2 + dz_km1));
        } else {
            dz_km1 = z[k] - z[k - 1];
            dz_k = z[k + 1] - z[k];
            // Equispaced grid
            // uxz = (ux_ijkp1 - ux_ijkm1) / (2.0 * dz);
            // uyz = (uy_ijkp1 - uy_ijkm1) / (2.0 * dz);
            // vxz = (vx_ijkp1 - vx_ijkm1) / (2.0 * dz);
            // vyz = (vy_ijkp1 - vy_ijkm1) / (2.0 * dz);
            // wxz = (wx_ijkp1 - wx_ijkm1) / (2.0 * dz);
            // wyz = (wy_ijkp1 - wy_ijkm1) / (2.0 * dz);
            // fwz = (fw_ijkp1 - fw_ijkm1) / (2.0 * dz);
            // Non-equispaced grid
            uxz = (ux_ijkp1 - ux_ijkm1) / (dz_k + dz_km1);
            uyz = (uy_ijkp1 - uy_ijkm1) / (dz_k + dz_km1);
            vxz = (vx_ijkp1 - vx_ijkm1) / (dz_k + dz_km1);
            vyz = (vy_ijkp1 - vy_ijkm1) / (dz_k + dz_km1);
            wxz = (wx_ijkp1 - wx_ijkm1) / (dz_k + dz_km1);
            wyz = (wy_ijkp1 - wy_ijkm1) / (dz_k + dz_km1);
            fwz = (fw_ijkp1 - fw_ijkm1) / (dz_k + dz_km1);
        }
        // Delta and l
        if (k < Nz - 1)
            dz_k = z[k + 1] - z[k];
        else
            dz_k = z[k] - z[k - 1];
        Delta = pow(dx * dy * dz_k, 1.0 / 3.0);
        l = C_s * Delta;
        // |S| = sqrt(2 * S_ij * S_ij)
        mod_S = sqrt(2.0 * (ux * ux + vy * vy + wz * wz) + (uz + wx) * (uz + wx) + (vx + uy) * (vx + uy) + (wy + vz) * (wy + vz)) + 1e-16;
        psi_x = 4 * (ux * uxx + vy * vyx + wz * wzx) + 2 * (uz + wx) * (uzx + wxx) + 2 * (vx + uy) * (vxx + uyx) + 2 * (wy + vz) * (wyx + vzx);
        psi_y = 4 * (ux * uxy + vy * vyy + wz * wzy) + 2 * (uz + wx) * (uzy + wxy) + 2 * (vx + uy) * (vxy + uyy) + 2 * (wy + vz) * (wyy + vzy);
        psi_z = 4 * (ux * uxz + vy * vyz + wz * wzz) + 2 * (uz + wx) * (uzz + wxz) + 2 * (vx + uy) * (vxz + uyz) + 2 * (wy + vz) * (wyz + vzz);
        // SGS model
        sgs_x_damp = 2 * mod_S * fw * (fwx * ux + 0.5 * fwy * (vx + uy) + 0.5 * fwz * (wx + uz));
        sgs_y_damp = 2 * mod_S * fw * (fwy * vy + 0.5 * fwx * (uy + vx) + 0.5 * fwz * (wy + vz));
        sgs_z_damp = 2 * mod_S * fw * (fwz * wz + 0.5 * fwx * (wx + uz) + 0.5 * fwy * (wy + vz));
        sgs_q_damp = 2 * mod_S * fw * (fwx * Tx + fwy * Ty + fwz * Tz);
        sgs_x_no_damp = 1 / (2 * mod_S) * (psi_x * ux + 0.5 * psi_y * (uy + vx) + 0.5 * psi_z * (wx + uz)) + mod_S * (uxx + 0.5 * (vxy + uyy) + 0.5 * (wxz + uzz));
        sgs_y_no_damp = 1 / (2 * mod_S) * (psi_y * vy + 0.5 * psi_x * (vx + uy) + 0.5 * psi_z * (wy + vz)) + mod_S * (vyy + 0.5 * (uyx + vxx) + 0.5 * (wyz + vzz));
        sgs_z_no_damp = 1 / (2 * mod_S) * (psi_z * wz + 0.5 * psi_x * (wx + uz) + 0.5 * psi_y * (wy + vz)) + mod_S * (wzz + 0.5 * (uzx + wxx) + 0.5 * (vzx + wyx));
        sgs_q_no_damp = 1 / (2 * mod_S) * (psi_x * Tx  + psi_y * Ty + psi_z * Tz) + mod_S * (Txx + Tyy + Tzz);
        sgs_x = -2 * l * l * (sgs_x_no_damp * fw * fw + sgs_x_damp);
        sgs_y = -2 * l * l * (sgs_y_no_damp * fw * fw + sgs_y_damp);
        sgs_z = -2 * l * l * (sgs_z_no_damp * fw * fw + sgs_z_damp);
        sgs_q = -l * l / Pr * (sgs_q_no_damp * fw * fw + sgs_q_damp);
        // Add SGS model to R
        R_new[parameters.field_indexes.u + IDX(i, j, k, Nx, Ny, Nz)] -= sgs_x;
        R_new[parameters.field_indexes.v + IDX(i, j, k, Nx, Ny, Nz)] -= sgs_y;
        R_new[parameters.field_indexes.w + IDX(i, j, k, Nx, Ny, Nz)] -= sgs_z;
        R_new[parameters.field_indexes.T + IDX(i, j, k, Nx, Ny, Nz)] -= sgs_q;
    }
}
*/