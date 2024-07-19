#include "../include/pde.h"

void Phi(double t, double *R_old, double *R_new, double *R_turbulence, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->Nz_Y;
    double dx = parameters->dx;
    double dy = parameters->dy;
    double dz = parameters->dz;
    // double dt = parameters->dt;
    // double *x = parameters->x;
    // double *y = parameters->y;
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
    double u_ip2jk, u_im2jk, u_ijp2k, u_ijm2k, u_ijkp2, u_ijkm2;
    double v_ijk, v_ip1jk, v_im1jk, v_ijp1k, v_ijm1k, v_ijkp1, v_ijkm1;
    double v_ip2jk, v_im2jk, v_ijp2k, v_ijm2k, v_ijkp2, v_ijkm2;
    double w_ijk, w_ip1jk, w_im1jk, w_ijp1k, w_ijm1k, w_ijkp1, w_ijkm1;
    double w_ip2jk, w_im2jk, w_ijp2k, w_ijm2k, w_ijkp2, w_ijkm2;
    double T_ijk, T_ip1jk, T_im1jk, T_ijp1k, T_ijm1k, T_ijkp1, T_ijkm1;
    double Y_ijk;
    // Upwind scheme terms
    double u_plu, u_min, v_plu, v_min, w_plu, w_min;
    double u_ip, u_im, u_jp, u_jm, u_kp, u_km;
    double v_ip, v_im, v_jp, v_jm, v_kp, v_km;
    double w_ip, w_im, w_jp, w_jm, w_kp, w_km;
    // First partial derivatives
    double ux_uw, uy_uw, uz_uw, vx_uw, vy_uw, vz_uw, wx_uw, wy_uw, wz_uw; // For upwind scheme
    double ux, uy, uz, vx, vy, vz, wx, wy, wz, Tx, Ty, Tz; // For central difference
    // Second partial derivatives
    double uxx, uyy, uzz, vxx, vyy, vzz, wxx, wyy, wzz, Txx, Tyy, Tzz;
    double u_RHS, v_RHS, w_RHS, T_RHS, Y_RHS, S;    
    double u_tau, tau_p, fw;
    double mod_U;
    double F_x, F_y, F_z;
    int im1, ip1, jm1, jp1;
    int im2, ip2, jm2, jp2;
    // Loop over interior nodes. Periodic boundary conditions in x and y.
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 1; k < Nz - 1; k++) {   
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
                u_ijkm1 = R_old[u_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
                v_ijkm1 = R_old[v_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
                w_ijkm1 = R_old[w_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
                T_ijkm1 = R_old[T_index + IDX(i, j, k - 1, Nx, Ny, Nz)];
                // \phi_{i,j,k+1}
                u_ijkp1 = R_old[u_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
                v_ijkp1 = R_old[v_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
                w_ijkp1 = R_old[w_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
                T_ijkp1 = R_old[T_index + IDX(i, j, k + 1, Nx, Ny, Nz)];
                // \phi_{i,j,k-2}
                if (k > 1) {
                    u_ijkm2 = R_old[u_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
                    v_ijkm2 = R_old[v_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
                    w_ijkm2 = R_old[w_index + IDX(i, j, k - 2, Nx, Ny, Nz)];
                } else {
                    u_ijkm2 = 0.0;
                    v_ijkm2 = 0.0;
                    w_ijkm2 = 0.0;
                }
                // \phi_{i,j,k+2}
                if (k < Nz - 2) {
                    u_ijkp2 = R_old[u_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
                    v_ijkp2 = R_old[v_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
                    w_ijkp2 = R_old[w_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
                } else {
                    u_ijkp2 = 0.0;
                    v_ijkp2 = 0.0;
                    w_ijkp2 = 0.0;
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
                // Second order forward difference at k=1
                if (k == 1) {
                    u_km = (-3 * u_ijk + 4 * u_ijkp1 - u_ijkp2) / (2 * dz);
                    v_km = (-3 * v_ijk + 4 * v_ijkp1 - v_ijkp2) / (2 * dz);
                    w_km = (-3 * w_ijk + 4 * w_ijkp1 - w_ijkp2) / (2 * dz);
                } else {
                    u_km = (3 * u_ijk - 4 * u_ijkm1 + u_ijkm2) / (2 * dz);
                    v_km = (3 * v_ijk - 4 * v_ijkm1 + v_ijkm2) / (2 * dz);
                    w_km = (3 * w_ijk - 4 * w_ijkm1 + w_ijkm2) / (2 * dz);
                }
                // Second order backward difference at k=Nz-2
                if (k == Nz - 2) {
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
                ux_uw = u_plu * u_im + u_min * u_ip; // du/dx
                uy_uw = u_plu * u_jm + u_min * u_jp; // du/dy
                uz_uw = u_plu * u_km + u_min * u_kp; // du/dz
                vx_uw = v_plu * v_im + v_min * v_ip; // dv/dx
                vy_uw = v_plu * v_jm + v_min * v_jp; // dv/dy
                vz_uw = v_plu * v_km + v_min * v_kp; // dv/dz
                wx_uw = w_plu * w_im + w_min * w_ip; // dw/dx
                wy_uw = w_plu * w_jm + w_min * w_jp; // dw/dy
                wz_uw = w_plu * w_km + w_min * w_kp; // dw/dz
                // Central difference for temperature and turbulence velocity
                ux  = (u_ip1jk - u_im1jk) / (2 * dx); // du/dx
                uy  = (u_ijp1k - u_ijm1k) / (2 * dy); // du/dy
                uz  = (u_ijkp1 - u_ijkm1) / (2 * dz); // du/dz
                vx  = (v_ip1jk - v_im1jk) / (2 * dx); // dv/dx
                vy  = (v_ijp1k - v_ijm1k) / (2 * dy); // dv/dy
                vz  = (v_ijkp1 - v_ijkm1) / (2 * dz); // dv/dz
                wx  = (w_ip1jk - w_im1jk) / (2 * dx); // dw/dx
                wy  = (w_ijp1k - w_ijm1k) / (2 * dy); // dw/dy
                wz  = (w_ijkp1 - w_ijkm1) / (2 * dz); // dw/dz
                Tx  = (T_ip1jk - T_im1jk) / (2.0 * dx); // dT/dx
                Ty  = (T_ijp1k - T_ijm1k) / (2.0 * dy); // dT/dy
                Tz  = (T_ijkp1 - T_ijkm1) / (2.0 * dz); // dT/dz
                /* Compute second partial derivatives */
                // Central difference for velocity and temperature
                uxx = (u_ip1jk - 2 * u_ijk + u_im1jk) / (dx * dx); // d^2u/dx^2
                uyy = (u_ijp1k - 2 * u_ijk + u_ijm1k) / (dy * dy); // d^2u/dy^2
                uzz = (u_ijkp1 - 2 * u_ijk + u_ijkm1) / (dz * dz); // d^2u/dz^2
                vxx = (v_ip1jk - 2 * v_ijk + v_im1jk) / (dx * dx); // d^2v/dx^2
                vyy = (v_ijp1k - 2 * v_ijk + v_ijm1k) / (dy * dy); // d^2v/dy^2
                vzz = (v_ijkp1 - 2 * v_ijk + v_ijkm1) / (dz * dz); // d^2v/dz^2
                wxx = (w_ip1jk - 2 * w_ijk + w_im1jk) / (dx * dx); // d^2w/dx^2
                wyy = (w_ijp1k - 2 * w_ijk + w_ijm1k) / (dy * dy); // d^2w/dy^2
                wzz = (w_ijkp1 - 2 * w_ijk + w_ijkm1) / (dz * dz); // d^2w/dz^2
                Txx = (T_ip1jk - 2.0 * T_ijk + T_im1jk) / (dx * dx); // d^2T/dx^2 
                Tyy = (T_ijp1k - 2.0 * T_ijk + T_ijm1k) / (dy * dy); // d^2T/dy^2
                Tzz = (T_ijkp1 - 2.0 * T_ijk + T_ijkm1) / (dz * dz); // d^2T/dz^2
                // Damping function
                tau_p = sqrt((nu * 0.5 * (uz + wx)) * (nu * 0.5 * (uz + wx)) + (nu * 0.5 * (vz + wy)) * (nu * 0.5 * (vz + wy)));
                u_tau = sqrt(tau_p);
                fw = f_damping(z[k], u_tau, nu);
                // fw = 0;
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
                // mod_U = sqrt(pow(u_ijk, 2) + pow(v_ijk, 2) + pow(w_ijk, 2));
                F_x = - Y_D * a_v * Y_ijk * mod_U * u_ijk;
                F_y = - Y_D * a_v * Y_ijk * mod_U * v_ijk;                             
                F_z = - g * (T_ijk - T_inf) / T_ijk - Y_D * a_v * Y_ijk * mod_U * w_ijk;
                // Compute RHS for velocity and temperature
                u_RHS = nu * (uxx + uyy + uzz) - (u_ijk * ux_uw + v_ijk * uy_uw + w_ijk * uz_uw);// + F_x;
                v_RHS = nu * (vxx + vyy + vzz) - (u_ijk * vx_uw + v_ijk * vy_uw + w_ijk * vz_uw);// + F_y*0;
                w_RHS = nu * (wxx + wyy + wzz) - (u_ijk * wx_uw + v_ijk * wy_uw + w_ijk * wz_uw);// + F_z*0;
                T_RHS = alpha * (Txx + Tyy + Tzz) - (u_ijk * Tx + v_ijk * Ty + w_ijk * Tz);// + S*0;
                // Save RHS into R_new
                R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] = u_RHS;
                R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] = v_RHS;
                R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] = w_RHS;
                R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = T_RHS;  
            }
        }
    }
}

void Phi_2(double t, double *R_old, double *R_new, double *R_turbulence, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->Nz_Y;
    double dx = parameters->dx;
    double dy = parameters->dy;
    double dz = parameters->dz;
    // double dt = parameters->dt;
    // double *x = parameters->x;
    // double *y = parameters->y;
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
    double ux_uw, uy_uw, uz_uw, vx_uw, vy_uw, vz_uw, wx_uw, wy_uw, wz_uw; // For upwind scheme
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
    // Loop over interior nodes. Periodic boundary conditions in x and y.
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {   
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
                // ux_uw = u_plu * u_im + u_min * u_ip; // u * du/dx
                // uy_uw = u_plu * u_jm + u_min * u_jp; // u * du/dy
                // uz_uw = u_plu * u_km + u_min * u_kp; // u * du/dz
                // vx_uw = v_plu * v_im + v_min * v_ip; // v * dv/dx
                // vy_uw = v_plu * v_jm + v_min * v_jp; // v * dv/dy
                // vz_uw = v_plu * v_km + v_min * v_kp; // v * dv/dz
                // wx_uw = w_plu * w_im + w_min * w_ip; // w * dw/dx
                // wy_uw = w_plu * w_jm + w_min * w_jp; // w * dw/dy
                // wz_uw = w_plu * w_km + w_min * w_kp; // w * dw/dz
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
                // fw = 0.0;
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
                // Y_D = 1.0;
                // a_v = 1.0;
                F_x = - Y_D * a_v * Y_ijk * mod_U * u_ijk;
                F_y = - Y_D * a_v * Y_ijk * mod_U * v_ijk;                             
                F_z = - g * (T_ijk - T_inf) / T_ijk - Y_D * a_v * Y_ijk * mod_U * w_ijk;
                // double T_tmp = (T_ijk - T_inf) / T_ijk;
                // F_z = (T_inf / T_ijk - 1.0) * g - Y_D * a_v * Y_ijk * mod_U * w_ijk;
                // Compute RHS for velocity and temperature
                // uux = ux_uw;
                // vuy = v_ijk * uy_uw;
                // wuz = w_ijk * uz_uw;
                // uvx = u_ijk * vx_uw;
                // vvy = v_ijk * vy_uw;
                // wvz = w_ijk * vz_uw;
                // uwx = u_ijk * wx_uw;
                // vwy = v_ijk * wy_uw;
                // wwz = w_ijk * wz_uw;
                lap_u = uxx + uyy + uzz;
                lap_v = vxx + vyy + vzz;
                lap_w = wxx + wyy + wzz;
                lap_T = Txx + Tyy + Tzz;
                u_RHS = nu * lap_u - (uux + vuy + wuz) + F_x;
                v_RHS = nu * lap_v - (uvx + vvy + wvz) + F_y;
                w_RHS = nu * lap_w - (uwx + vwy + wwz) + F_z;
                // T_RHS = alpha * lap_T - (u_ijk * Tx + v_ijk * Ty + w_ijk * Tz) + S;
                // u_RHS = nu * lap_u - (uux + vuy + wuz);
                // v_RHS = nu * lap_v - (uvx + vvy + wvz);
                // w_RHS = nu * lap_w - (uwx + vwy + wwz);
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
    double T_inf = parameters->T_inf;
    double u_ijkm1, u_ijkm2;
    double v_ijkm1, v_ijkm2;
    double w_ijkm1, w_ijkm2;
    double T_ijkp1, T_ijkm1, T_ijkm2, T_ijkp2;
    double Y_ijkp1, Y_ijkm1, Y_ijkm2, Y_ijkp2;
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

void bound(double *R_new, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->Nz_Y;
    int T_index = parameters->field_indexes.T;
    int Y_index = parameters->field_indexes.Y;
    double T_inf = parameters->T_inf;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                // Check if Y_ijk is not valid
                if (k < Nz_Y) {
                    if (R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] < 0) {
                        R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] = 0;
                    }
                    if (R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] > 1) {
                        R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] = 1;
                    }
                }
                // Check if T_ijk is less than T_inf and higher than 1500
                if (R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] < T_inf) {
                    R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = T_inf;
                }
                if (R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] > 1500) {
                    R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = 1500;
                }
            }
        }
    }
}

void velocity_correction(double *R_new, double *p, double dt, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->Nz_Y;
    int u_index = parameters->field_indexes.u;
    int v_index = parameters->field_indexes.v;
    int w_index = parameters->field_indexes.w;
    int T_index = parameters->field_indexes.T;
    int Y_index = parameters->field_indexes.Y;
    double rho = parameters->rho;
    double T_inf = parameters->T_inf;
    double dx = parameters->dx;
    double dy = parameters->dy;
    double dz = parameters->dz;
    // double u_ijkm1, u_ijkm2, v_ijkm1, v_ijkm2, w_ijkm1, w_ijkm2, T_ijkm1, T_ijkm2, T_ijkp1, T_ijkp2, Y_ijkm1, Y_ijkm2, Y_ijkp1, Y_ijkp2;
    double pim1jk, pip1jk, pijm1k, pijp1k, pijkp1, pijkm1;
    double px, py, pz;
    int im1, ip1, jm1, jp1;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 1; k < Nz-1; k++) {
                // Indexes for periodic boundary conditions
                im1 = (i - 1 + Nx - 1) % (Nx - 1); // i-1
                jm1 = (j - 1 + Ny - 1) % (Ny - 1); // j-1
                ip1 = (i + 1) % (Nx - 1); // i+1
                jp1 = (j + 1) % (Ny - 1); // j+1
                pim1jk = p[IDX(im1, j, k, Nx, Ny, Nz)]; // p_{i-1,j,k}
                pip1jk = p[IDX(ip1, j, k, Nx, Ny, Nz)]; // p_{i+1,j,k}
                pijm1k = p[IDX(i, jm1, k, Nx, Ny, Nz)]; // p_{i,j-1,k}
                pijp1k = p[IDX(i, jp1, k, Nx, Ny, Nz)]; // p_{i,j+1,k}
                // if (k > 0 && k < Nz - 1) { // Interior points
                    pijkp1 = p[IDX(i, j, k + 1, Nx, Ny, Nz)]; // p_{i,j,k+1}
                    pijkm1 = p[IDX(i, j, k - 1, Nx, Ny, Nz)]; // p_{i,j,k-1}
                    px = (pip1jk - pim1jk) / (2 * dx);
                    py = (pijp1k - pijm1k) / (2 * dy);
                    pz = (pijkp1 - pijkm1) / (2 * dz);
                    // px = 0;
                    // py = 0;
                    // pz = 0;
                    R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * px;
                    R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * py;
                    R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * pz;
            }
        }
    }
}

void velocity_correction_boundaries_bounding(double *R_new, double *p, double dt, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->Nz_Y;
    int u_index = parameters->field_indexes.u;
    int v_index = parameters->field_indexes.v;
    int w_index = parameters->field_indexes.w;
    int T_index = parameters->field_indexes.T;
    int Y_index = parameters->field_indexes.Y;
    double rho = parameters->rho;
    double T_inf = parameters->T_inf;
    double dx = parameters->dx;
    double dy = parameters->dy;
    double dz = parameters->dz;
    double u_ijkm1, u_ijkm2, v_ijkm1, v_ijkm2, w_ijkm1, w_ijkm2, T_ijkm1, T_ijkm2, T_ijkp1, T_ijkp2, Y_ijkm1, Y_ijkm2, Y_ijkp1, Y_ijkp2;
    double pim1jk, pip1jk, pijm1k, pijp1k, pijkp1, pijkm1;
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
                pim1jk = p[IDX(im1, j, k, Nx, Ny, Nz)]; // p_{i-1,j,k}
                pip1jk = p[IDX(ip1, j, k, Nx, Ny, Nz)]; // p_{i+1,j,k}
                pijm1k = p[IDX(i, jm1, k, Nx, Ny, Nz)]; // p_{i,j-1,k}
                pijp1k = p[IDX(i, jp1, k, Nx, Ny, Nz)]; // p_{i,j+1,k}
                if (k > 0 && k < Nz - 1) { // Interior points
                    pijkp1 = p[IDX(i, j, k + 1, Nx, Ny, Nz)]; // p_{i,j,k+1}
                    pijkm1 = p[IDX(i, j, k - 1, Nx, Ny, Nz)]; // p_{i,j,k-1}
                    px = (pip1jk - pim1jk) / (2 * dx);
                    py = (pijp1k - pijm1k) / (2 * dy);
                    pz = (pijkp1 - pijkm1) / (2 * dz);
                    R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * px;
                    R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * py;
                    R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * pz;
                }
                // Bottom boundary z_k = z_min, k=0
                if (k == 0) {
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
                }
                // Top boundary z_k = z_max, k=Nz-1
                if (k == Nz - 1) {
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
                // Check if Y_ijk is not valid
                if (k < Nz_Y) {
                    if (R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] < 0) {
                        R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] = 0;
                    }
                    if (R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] > 1) {
                        R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] = 1;
                    }
                }
                // Check if T_ijk is less than T_inf and higher than 1500
                if (R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] < T_inf) {
                    R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = T_inf;
                }
                if (R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] > 1500) {
                    R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = 1500;
                }
            }
        }
    }
}

void euler(double t_n, double *y_n, double *y_np1, double *F, double *U_turbulence, double dt, int size, Parameters *parameters) {
    Phi_2(t_n, y_n, F, U_turbulence, parameters);
    turbulence(U_turbulence, F, parameters);
    boundary_conditions(F, parameters);
    for (int i = 0; i < size; i++) {        
        y_np1[i] = y_n[i] + dt * F[i];
    }
    // boundary_conditions(y_np1, parameters);
    // bound(y_np1, parameters);
}

void rk2(double t_n, double *y_n, double *y_np1, double *k, double *F, double *U_turbulence, double dt, int size, Parameters *parameters) {
    int k1_index = 0;
    int k2_index = size;
    Phi_2(t_n, y_n, k + k1_index, U_turbulence, parameters);
    turbulence(U_turbulence, k + k1_index, parameters);
    boundary_conditions(k + k1_index, parameters);
    caxpy(F, k + k1_index, y_n, dt, size);
    Phi_2(t_n + dt, F, k + k2_index, U_turbulence, parameters);
    turbulence(U_turbulence, k + k2_index, parameters);
    boundary_conditions(k + k2_index, parameters);
    for (int i = 0; i < size; i++) {
        y_np1[i] = y_n[i] + 0.5 * dt * (k[k1_index + i] + k[k2_index + i]);
    }
}

void rk4(double t_n, double *y_n, double *y_np1, double *k, double *F, double *U_turbulence,double dt, int size, Parameters *parameters) {
    int k1_index = 0;
    int k2_index = size;
    int k3_index = 2 * size;
    int k4_index = 3 * size;
    Phi_2(t_n, y_n, k + k1_index, U_turbulence, parameters);
    caxpy(F, k + k1_index, y_n, dt * 0.5, size);
    Phi_2(t_n + 0.5 * dt, F, k + k2_index, U_turbulence, parameters);
    caxpy(F, k + k2_index, y_n, dt * 0.5, size);
    Phi_2(t_n + 0.5 * dt, F, k + k3_index, U_turbulence, parameters);
    caxpy(F, k + k3_index, y_n, dt, size);
    Phi_2(t_n + dt, F, k + k4_index, U_turbulence, parameters);
    for (int i = 0; i < size; i++) {
        y_np1[i] = y_n[i] + (dt / 6) * (k[k1_index + i] + 2 * k[k2_index + i] + 2 * k[k3_index + i] + k[k4_index + i]);
    }
}

void create_y_0(double *u, double *v, double *w, double *T, double *Y, double *y_0, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->Nz_Y;
    int size = Nx * Ny * Nz;
    int size_Y = Nx * Ny * Nz_Y;
    for (int i = 0; i < size; i++) {
        y_0[parameters->field_indexes.u + i] = u[i];
        y_0[parameters->field_indexes.v + i] = v[i];
        y_0[parameters->field_indexes.w + i] = w[i];
        y_0[parameters->field_indexes.T + i] = T[i];
        if (i < size_Y) {
            y_0[parameters->field_indexes.Y + i] = Y[i];
        }
    }
}

void solve_PDE(double *y_n, double *p, Parameters *parameters) {
    setbuf(stdout, NULL);
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nt = parameters->Nt;
    int NT = parameters->NT;
    int Nz_Y = parameters->Nz_Y;
    int size = 4 * Nx * Ny * Nz + Nx * Ny * Nz_Y;
    int n_save;
    double step_time;
    double *t = parameters->t;
    double dt = parameters->dt;
    double *y_np1 = (double *) malloc(size * sizeof(double));
    double *F = (double *) malloc(size * sizeof(double));
    double *k = (double *) malloc(2 * size * sizeof(double));
    double *R_turbulence = (double *) malloc(27 * Nx * Ny * Nz * sizeof(double));
    // Arrays for pressure Poisson Problem
    // double *p = (double *) malloc(Nx * Ny * Nz * sizeof(double));
    // double *p_top = (double *) malloc(Nx * Ny * Nz * sizeof(double));
    // double *f = (double *) malloc(Nx * Ny * Nz * sizeof(double)); 
    // double *a = (double *) malloc((Nz - 2) * sizeof(double));
    // double *b = (double *) malloc((Nz - 1) * sizeof(double));
    // double *c = (double *) malloc((Nz - 2) * sizeof(double));
    // double complex *d = (double complex *) malloc((Nz - 1) * sizeof(double complex));
    // double complex *l = (double complex *) malloc((Nz - 2) * sizeof(double complex));
    // double complex *u = (double complex *) malloc((Nz - 1) * sizeof(double complex));
    // double complex *y = (double complex *) malloc((Nz - 1) * sizeof(double complex));
    fftw_complex *a = fftw_alloc_complex((Nz - 2));
    fftw_complex *b = fftw_alloc_complex((Nz - 1));
    fftw_complex *c = fftw_alloc_complex((Nz - 2));
    fftw_complex *d = fftw_alloc_complex((Nz - 1));
    fftw_complex *l = fftw_alloc_complex((Nz - 2));
    fftw_complex *u = fftw_alloc_complex((Nz - 1));
    fftw_complex *y = fftw_alloc_complex((Nz - 1));
    // double complex *pk = (double complex *) malloc((Nz - 1) * sizeof(double complex));
    fftw_complex *pk = fftw_alloc_complex((Nz - 1));
    // fftw_complex *pk = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    fftw_complex *f_in = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    fftw_complex *f_out = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    fftw_complex *p_top_in = fftw_alloc_complex((Nx - 1) * (Ny - 1));
    fftw_complex *p_top_out = fftw_alloc_complex((Nx - 1) * (Ny - 1));
    fftw_complex *p_in = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    fftw_complex *p_out = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    fftw_plan p_plan, f_plan, p_top_plan;
    char *save_path = "data/output/";
    clock_t start, end, step_start, step_end;
    // copy_slice(p_top, p_0, Nz - 1, Nx, Ny, Nz);
    // FILE *logs = fopen("data/output/logs.txt", "w");
    int u_index = parameters->field_indexes.u;
    int v_index = parameters->field_indexes.v;
    int w_index = parameters->field_indexes.w;

    // Solver time
    start = clock();
    // Time integration
    for (int n = 1; n < Nt; n++) { 
        step_start = clock();
        // printf("Euler step\n");
        // Euler step to compute U^*, T^{n+1}, Y^{n+1}
        // euler(t[n], y_n, y_np1, F, R_turbulence, dt, size, parameters);

        // RK2 step to compute U^*, T^{n+1}, Y^{n+1}
        rk2(t[n], y_n, y_np1, k, F, R_turbulence, dt, size, parameters);

        // RK4 step to compute U^*, T^{n+1}, Y^{n+1}
        // rk4(t[n], y_n, y_np1, k, F, R_turbulence, dt, size, parameters);

        // Boundary conditions
        // boundary_conditions(y_np1, parameters);

        // Solve Poisson problem for pressure (it only uses U^*)
        // (THE PROBLEM COULD BE HERE)
        // printf("Solving pressure Poisson problem\n");
        solve_pressure(y_np1, p, a, b, c, d, l, u, y, pk, p_plan, f_plan, p_top_plan, f_in, f_out, p_top_in, p_top_out, p_in, p_out, parameters);

        // // Velocity correction 
        // velocity_correction(y_np1, p, dt, parameters);

        // // Boundary conditions
        // boundary_conditions(y_np1, parameters); 

        // // Bound temperature and fuel
        // bound(y_np1, parameters);

        // printf("Velocity correction\n");
        velocity_correction_boundaries_bounding(y_np1, p, dt, parameters);
        
        step_end = clock();

        step_time = (double) (step_end - step_start) / CLOCKS_PER_SEC;
        
        // Save data each NT steps and at the last step
        if (n % NT == 0 || n == Nt - 1) {  
            n_save = n / NT;
            printf("n = %d, t_n = %lf\n", n, t[n]);      
            printf("Time per iteration: %lf s\n", step_time); 
            printf("Saving data...\n");
            // fprintf(logs, "n = %d, t_n = %lf\n", n, t[n]);
            // fprintf(logs, "Time per iteration: %lf s\n", step_time); 
            // fprintf(logs, "Saving data...\n");
            save_data(save_path, y_np1, p, n_save, parameters);
        }

        // Update y_n
        copy(y_n, y_np1, size);
    }
    end = clock();
    printf("Solver time: %lf s\n", (double) (end - start) / CLOCKS_PER_SEC);

    // fclose(logs);
    // Free memory
    // free(t);
    // free(y_n);
    // free(y_np1);
    // free(F);
    // free(R_turbulence);
    // free(k);
    // free(p);
    // free(p_top);
    // free(f);
    // Free memory
    // fftw_destroy_plan(p_top_plan);
    // fftw_destroy_plan(f_plan);
    // fftw_destroy_plan(p_plan);
    // fftw_free(p_top_in);
    // fftw_free(p_in);
    // fftw_free(f_in);
    // fftw_free(p_top_out);
    // fftw_free(f_out);
    // fftw_free(p_out);
    // fftw_free(pk);
    // fftw_free(d);
    // fftw_free(l);
    // fftw_free(uu);
    // fftw_free(y);
    // fftw_free(a);
    // fftw_free(b);
    // fftw_free(c);
    // free(a);
    // free(b);
    // free(c);
    // free(d);
    // free(l);
    // free(u);
    // free(y);
}