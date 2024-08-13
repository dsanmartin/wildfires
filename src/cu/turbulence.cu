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
void turbulence(double *R_turbulence, double *R_new, Parameters parameters) {
    // printf("Turbulence\n");
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int size = Nx * Ny * Nz;
    double dx = parameters.dx;
    double dy = parameters.dy;
    double dz = parameters.dz;
    double Delta = pow(dx * dy * dz, 1.0 / 3.0);
    double C_s = parameters.C_s;
    double Pr = parameters.Pr;
    double l = C_s * Delta;
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
            ux_ijkp2 = R_turbulence[parameters.turbulence_indexes.ux + IDX(i, j, k + 2, Nx, Ny, Nz)];
            uy_ijkp2 = R_turbulence[parameters.turbulence_indexes.uy + IDX(i, j, k + 2, Nx, Ny, Nz)];
            vx_ijkp2 = R_turbulence[parameters.turbulence_indexes.vx + IDX(i, j, k + 2, Nx, Ny, Nz)];
            vy_ijkp2 = R_turbulence[parameters.turbulence_indexes.vy + IDX(i, j, k + 2, Nx, Ny, Nz)];
            wx_ijkp2 = R_turbulence[parameters.turbulence_indexes.wx + IDX(i, j, k + 2, Nx, Ny, Nz)];
            wy_ijkp2 = R_turbulence[parameters.turbulence_indexes.wy + IDX(i, j, k + 2, Nx, Ny, Nz)];
            fw_ijkp2 = R_turbulence[parameters.turbulence_indexes.fw + IDX(i, j, k + 2, Nx, Ny, Nz)];
            uxz = (-3 * ux + 4 * ux_ijkp1 - ux_ijkp2) / (2 * dz);
            uyz = (-3 * uy + 4 * uy_ijkp1 - uy_ijkp2) / (2 * dz);
            vxz = (-3 * vx + 4 * vx_ijkp1 - vx_ijkp2) / (2 * dz);
            vyz = (-3 * vy + 4 * vy_ijkp1 - vy_ijkp2) / (2 * dz);
            wxz = (-3 * wx + 4 * wx_ijkp1 - wx_ijkp2) / (2 * dz);
            wyz = (-3 * wy + 4 * wy_ijkp1 - wy_ijkp2) / (2 * dz);
            fwz = (-3 * fw + 4 * fw_ijkp1 - fw_ijkp2) / (2 * dz);
        } else if (k == Nz - 1) { // Second-order backward difference
            ux_ijkm2 = R_turbulence[parameters.turbulence_indexes.ux + IDX(i, j, k - 2, Nx, Ny, Nz)];
            uy_ijkm2 = R_turbulence[parameters.turbulence_indexes.uy + IDX(i, j, k - 2, Nx, Ny, Nz)];
            vx_ijkm2 = R_turbulence[parameters.turbulence_indexes.vx + IDX(i, j, k - 2, Nx, Ny, Nz)];
            vy_ijkm2 = R_turbulence[parameters.turbulence_indexes.vy + IDX(i, j, k - 2, Nx, Ny, Nz)];
            wx_ijkm2 = R_turbulence[parameters.turbulence_indexes.wx + IDX(i, j, k - 2, Nx, Ny, Nz)];
            wy_ijkm2 = R_turbulence[parameters.turbulence_indexes.wy + IDX(i, j, k - 2, Nx, Ny, Nz)];
            fw_ijkm2 = R_turbulence[parameters.turbulence_indexes.fw + IDX(i, j, k - 2, Nx, Ny, Nz)];
            uxz = (3 * ux - 4 * ux_ijkm1 + ux_ijkm2) / (2 * dz);
            uyz = (3 * uy - 4 * uy_ijkm1 + uy_ijkm2) / (2 * dz);
            vxz = (3 * vx - 4 * vx_ijkm1 + vx_ijkm2) / (2 * dz);
            vyz = (3 * vy - 4 * vy_ijkm1 + vy_ijkm2) / (2 * dz);
            wxz = (3 * wx - 4 * wx_ijkm1 + wx_ijkm2) / (2 * dz);
            wyz = (3 * wy - 4 * wy_ijkm1 + wy_ijkm2) / (2 * dz);
            fwz = (3 * fw - 4 * fw_ijkm1 + fw_ijkm2) / (2 * dz);
        } else {
            uxz = (ux_ijkp1 - ux_ijkm1) / (2.0 * dz);
            uyz = (uy_ijkp1 - uy_ijkm1) / (2.0 * dz);
            vxz = (vx_ijkp1 - vx_ijkm1) / (2.0 * dz);
            vyz = (vy_ijkp1 - vy_ijkm1) / (2.0 * dz);
            wxz = (wx_ijkp1 - wx_ijkm1) / (2.0 * dz);
            wyz = (wy_ijkp1 - wy_ijkm1) / (2.0 * dz);
            fwz = (fw_ijkp1 - fw_ijkm1) / (2.0 * dz);
        }
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