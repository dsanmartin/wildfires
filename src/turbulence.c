#include "../include/turbulence.h"

void turbulence(double *R_tubulence, double *R, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    double dx = parameters->dx;
    double dy = parameters->dy;
    double dz = parameters->dz;
    double Delta = pow(dx * dy * dz, 1.0 / 3.0);
    double C_s = parameters->C_s;
    double Pr = parameters->Pr;
    // double nu = parameters->nu;
    double l = C_s * Delta;
    double ux, uy, uz, vx, vy, vz, wx, wy, wz, Tx, Ty, Tz;
    double uxx, uyy, uzz, vxx, vyy, vzz, wxx, wyy, wzz, Txx, Tyy, Tzz;
    double uyx, uzx, uxy, uzy, uxz, uyz;
    double vyx, vzx, vxy, vzy, vxz, vyz;
    double wyx, wzx, wxy, wzy, wxz, wyz;
    double uy_ip1jk, uy_im1jk, uz_ip1jk, uz_im1jk, uz_ijp1k, uz_ijm1k;
    double vy_ip1jk, vy_im1jk, vz_ip1jk, vz_im1jk, vz_ijp1k, vz_ijm1k;
    double wy_ip1jk, wy_im1jk, wz_ip1jk, wz_im1jk, wz_ijp1k, wz_ijm1k;
    double fw_ip1jk, fw_im1jk, fw_ijp1k, fw_ijm1k, fw_ijkp1, fw_ijkm1;
    double sgs_x, sgs_y, sgs_z, sgs_q;
    double mod_S, psi_x, psi_y, psi_z;
    double sgs_x_no_damp, sgs_y_no_damp, sgs_z_no_damp, sgs_q_no_damp;
    double sgs_x_damp, sgs_y_damp, sgs_z_damp, sgs_q_damp;
    double fw, fwx, fwy, fwz;
    // double u_tau, tau_p;
    // double *z = parameters->z;

    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int k = 1; k < Nz - 1; k++) {
                uy_ip1jk = R_tubulence[parameters->turbulence_indexes.uy + IDX(i + 1, j, k, Nx, Ny, Nz)];
                uy_im1jk = R_tubulence[parameters->turbulence_indexes.uy + IDX(i - 1, j, k, Nx, Ny, Nz)];
                uz_ip1jk = R_tubulence[parameters->turbulence_indexes.uz + IDX(i + 1, j, k, Nx, Ny, Nz)];
                uz_im1jk = R_tubulence[parameters->turbulence_indexes.uz + IDX(i - 1, j, k, Nx, Ny, Nz)];
                uz_ijp1k = R_tubulence[parameters->turbulence_indexes.uz + IDX(i, j + 1, k, Nx, Ny, Nz)];
                uz_ijm1k = R_tubulence[parameters->turbulence_indexes.uz + IDX(i, j - 1, k, Nx, Ny, Nz)];
                vy_ip1jk = R_tubulence[parameters->turbulence_indexes.vy + IDX(i + 1, j, k, Nx, Ny, Nz)];
                vy_im1jk = R_tubulence[parameters->turbulence_indexes.vy + IDX(i - 1, j, k, Nx, Ny, Nz)];
                vz_ip1jk = R_tubulence[parameters->turbulence_indexes.vz + IDX(i + 1, j, k, Nx, Ny, Nz)];
                vz_im1jk = R_tubulence[parameters->turbulence_indexes.vz + IDX(i - 1, j, k, Nx, Ny, Nz)];
                vz_ijp1k = R_tubulence[parameters->turbulence_indexes.vz + IDX(i, j + 1, k, Nx, Ny, Nz)];
                vz_ijm1k = R_tubulence[parameters->turbulence_indexes.vz + IDX(i, j - 1, k, Nx, Ny, Nz)];
                wy_ip1jk = R_tubulence[parameters->turbulence_indexes.wy + IDX(i + 1, j, k, Nx, Ny, Nz)];
                wy_im1jk = R_tubulence[parameters->turbulence_indexes.wy + IDX(i - 1, j, k, Nx, Ny, Nz)];
                wz_ip1jk = R_tubulence[parameters->turbulence_indexes.wz + IDX(i + 1, j, k, Nx, Ny, Nz)];
                wz_im1jk = R_tubulence[parameters->turbulence_indexes.wz + IDX(i - 1, j, k, Nx, Ny, Nz)];
                wz_ijp1k = R_tubulence[parameters->turbulence_indexes.wz + IDX(i, j + 1, k, Nx, Ny, Nz)];
                wz_ijm1k = R_tubulence[parameters->turbulence_indexes.wz + IDX(i, j - 1, k, Nx, Ny, Nz)];
                fw_ip1jk = R_tubulence[parameters->turbulence_indexes.fw + IDX(i + 1, j, k, Nx, Ny, Nz)];
                fw_im1jk = R_tubulence[parameters->turbulence_indexes.fw + IDX(i - 1, j, k, Nx, Ny, Nz)];
                fw_ijp1k = R_tubulence[parameters->turbulence_indexes.fw + IDX(i, j + 1, k, Nx, Ny, Nz)];
                fw_ijm1k = R_tubulence[parameters->turbulence_indexes.fw + IDX(i, j - 1, k, Nx, Ny, Nz)];
                fw_ijkp1 = R_tubulence[parameters->turbulence_indexes.fw + IDX(i, j, k + 1, Nx, Ny, Nz)];
                fw_ijkm1 = R_tubulence[parameters->turbulence_indexes.fw + IDX(i, j, k - 1, Nx, Ny, Nz)];
                fw       = R_tubulence[parameters->turbulence_indexes.fw + IDX(i, j, k, Nx, Ny, Nz)];
                // Mixed partial derivatives
                uyx = (uy_ip1jk - uy_im1jk) / (2.0 * dx);
                uzx = (uz_ip1jk - uz_im1jk) / (2.0 * dx);
                uzy = (uz_ijp1k - uz_ijm1k) / (2.0 * dy);
                uxy = uyx;
                uxz = uzx;
                uyz = uzy;
                vyx = (vy_ip1jk - vy_im1jk) / (2.0 * dx);
                vzx = (vz_ip1jk - vz_im1jk) / (2.0 * dx);
                vzy = (vz_ijp1k - vz_ijm1k) / (2.0 * dy);
                vxy = vyx;
                vxz = vzx;
                vyz = vzy;
                wyx = (wy_ip1jk - wy_im1jk) / (2.0 * dx);
                wzx = (wz_ip1jk - wz_im1jk) / (2.0 * dx);
                wzy = (wz_ijp1k - wz_ijm1k) / (2.0 * dy);
                wxy = wyx;
                wxz = wzx;
                wyz = wzy;
                // Velocity gradients from R
                ux = R_tubulence[parameters->turbulence_indexes.ux + IDX(i, j, k, Nx, Ny, Nz)];
                uy = R_tubulence[parameters->turbulence_indexes.uy + IDX(i, j, k, Nx, Ny, Nz)];
                uz = R_tubulence[parameters->turbulence_indexes.uz + IDX(i, j, k, Nx, Ny, Nz)];
                vx = R_tubulence[parameters->turbulence_indexes.vx + IDX(i, j, k, Nx, Ny, Nz)];
                vy = R_tubulence[parameters->turbulence_indexes.vy + IDX(i, j, k, Nx, Ny, Nz)];
                vz = R_tubulence[parameters->turbulence_indexes.vz + IDX(i, j, k, Nx, Ny, Nz)];
                wx = R_tubulence[parameters->turbulence_indexes.wx + IDX(i, j, k, Nx, Ny, Nz)];
                wy = R_tubulence[parameters->turbulence_indexes.wy + IDX(i, j, k, Nx, Ny, Nz)];
                wz = R_tubulence[parameters->turbulence_indexes.wz + IDX(i, j, k, Nx, Ny, Nz)];
                Tx = R_tubulence[parameters->turbulence_indexes.Tx + IDX(i, j, k, Nx, Ny, Nz)];
                Ty = R_tubulence[parameters->turbulence_indexes.Ty + IDX(i, j, k, Nx, Ny, Nz)];
                Tz = R_tubulence[parameters->turbulence_indexes.Tz + IDX(i, j, k, Nx, Ny, Nz)];
                uxx = R_tubulence[parameters->turbulence_indexes.uxx + IDX(i, j, k, Nx, Ny, Nz)];
                uyy = R_tubulence[parameters->turbulence_indexes.uyy + IDX(i, j, k, Nx, Ny, Nz)];
                uzz = R_tubulence[parameters->turbulence_indexes.uzz + IDX(i, j, k, Nx, Ny, Nz)];
                vxx = R_tubulence[parameters->turbulence_indexes.vxx + IDX(i, j, k, Nx, Ny, Nz)];
                vyy = R_tubulence[parameters->turbulence_indexes.vyy + IDX(i, j, k, Nx, Ny, Nz)];
                vzz = R_tubulence[parameters->turbulence_indexes.vzz + IDX(i, j, k, Nx, Ny, Nz)];
                wxx = R_tubulence[parameters->turbulence_indexes.wxx + IDX(i, j, k, Nx, Ny, Nz)];
                wyy = R_tubulence[parameters->turbulence_indexes.wyy + IDX(i, j, k, Nx, Ny, Nz)];
                wzz = R_tubulence[parameters->turbulence_indexes.wzz + IDX(i, j, k, Nx, Ny, Nz)];
                Txx = R_tubulence[parameters->turbulence_indexes.Txx + IDX(i, j, k, Nx, Ny, Nz)];
                Tyy = R_tubulence[parameters->turbulence_indexes.Tyy + IDX(i, j, k, Nx, Ny, Nz)];
                Tzz = R_tubulence[parameters->turbulence_indexes.Tzz + IDX(i, j, k, Nx, Ny, Nz)];

                // tau_p = sqrt((nu * 0.5 * (uz + wx)) * (nu * 0.5 * (uz + wx)) + (nu * 0.5 * (vz + wy)) * (nu * 0.5 * (vz + wy)));
                // u_tau = sqrt(tau_p);
                // fw = f_damping(z[k], u_tau, nu);
                fwx = (fw_ip1jk - fw_im1jk) / (2 * dx);
                fwy = (fw_ijp1k - fw_ijm1k) / (2 * dy);
                fwz = (fw_ijkp1 - fw_ijkm1) / (2 * dz);
    
                // |S| = sqrt(2 * S_ij * S_ij)
                mod_S = sqrt(2.0 * (ux * ux + vy * vy + wz * wz + (uz + wx) * (uz + wx) + (vx + uy) * (vx + uy) + (wy + vz) * (wy + vz))) + 1e-16;
                psi_x = 4 * (ux * uxx + vy * vyx + wz * wzx) + 2 * (uz + wx) * (uzx + wxx) + 2 * (vx + uy) * (vxx + uyx) + 2 * (wy + vz) * (wyx + vzx);
                psi_y = 4 * (ux * uxy + vy * vyy + wz * wzy) + 2 * (uz + wx) * (uzy + wxy) + 2 * (vx + uy) * (vxy + uyy) + 2 * (wy + vz) * (wyy + vzy);
                psi_z = 4 * (ux * uxz + vy * vyz + wz * wzz) + 2 * (uz + wx) * (uzz + wxz) + 2 * (vx + uy) * (vxz + uyz) + 2 * (wy + vz) * (wyz + vzz);

                // SGS model
                sgs_x_damp = 2 * mod_S * fw * (fwx * ux + 0.5 * fwy * (vx + uy) + 0.5 * fwz * (wx + uz));
                sgs_y_damp = 2 * mod_S * fw * (fwy * vy + 0.5 * fwx * (uy + vx) + 0.5 * fwz * (wy + vz));
                sgs_z_damp = 2 * mod_S * fw * (fwz * wz + 0.5 * fwx * (wx + uz) + 0.5 * fwy * (wy + vz));
                sgs_q_damp = 2 * fw * mod_S * (fwx * Tx + fwy * Ty + fwz * Tz);
                sgs_x_no_damp = 1 / (2 * mod_S) * (psi_x * ux + 0.5 * psi_y * (uy + vx) + 0.5 * psi_z * (wx + uz)) + mod_S * (uxx + 0.5 * (vxy + uyy) + 0.5 * (wxz + uzz));
                sgs_y_no_damp = 1 / (2 * mod_S) * (psi_y * vy + 0.5 * psi_x * (vx + uy) + 0.5 * psi_z * (wy + vz)) + mod_S * (vyy + 0.5 * (uyx + vxx) + 0.5 * (wyz + vzz));
                sgs_z_no_damp = 1 / (2 * mod_S) * (psi_z * wz + 0.5 * psi_x * (wx + uz) + 0.5 * psi_y * (wy + vz)) + mod_S * (wzz + 0.5 * (uzx + wxx) + 0.5 * (vzx + wyx));
                sgs_q_no_damp = 1 / (2 * mod_S) * (psi_x * Tx  + psi_y * Ty + psi_z * Tz) + mod_S * (Txx + Tyy + Tzz);
                sgs_x = -2 * l * l * (sgs_x_no_damp * fw * fw + sgs_x_damp);
                sgs_y = -2 * l * l * (sgs_y_no_damp * fw * fw + sgs_y_damp);
                sgs_z = -2 * l * l * (sgs_z_no_damp * fw * fw + sgs_z_damp);
                sgs_q = -l * l / Pr * (sgs_q_no_damp * fw * fw + sgs_q_damp);

                // Add SGS model to R
                R[parameters->field_indexes.u + IDX(i, j, k, Nx, Ny, Nz)] -= sgs_x;
                R[parameters->field_indexes.v + IDX(i, j, k, Nx, Ny, Nz)] -= sgs_y;
                R[parameters->field_indexes.w + IDX(i, j, k, Nx, Ny, Nz)] -= sgs_z;
                R[parameters->field_indexes.T + IDX(i, j, k, Nx, Ny, Nz)] -= sgs_q;
            }
        }
    }
}