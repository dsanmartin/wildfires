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
    // double Y_h = parameters->Y_h;
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
    double u_ip2jk = 0.0, u_im2jk = 0.0, u_ijp2k = 0.0, u_ijm2k = 0.0, u_ijkp2 = 0.0, u_ijkm2 = 0.0;
    double v_ijk, v_ip1jk, v_im1jk, v_ijp1k, v_ijm1k, v_ijkp1, v_ijkm1;
    double v_ip2jk = 0.0, v_im2jk = 0.0, v_ijp2k = 0.0, v_ijm2k = 0.0, v_ijkp2 = 0.0, v_ijkm2 = 0.0;
    double w_ijk, w_ip1jk, w_im1jk, w_ijp1k, w_ijm1k, w_ijkp1, w_ijkm1;
    double w_ip2jk = 0.0, w_im2jk = 0.0, w_ijp2k = 0.0, w_ijm2k = 0.0, w_ijkp2 = 0.0, w_ijkm2 = 0.0;
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
    double u_RHS, v_RHS, w_RHS, T_RHS, Y_RHS, S = 0.0;    
    double u_tau, tau_p, fw;
    double mod_U;
    double F_x = 0.0, F_y = 0.0, F_z = 0.0;
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
                // Actual nodes \phi_{i,j,k}
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
                }
                // \phi_{i,j,k+2}
                if (k < Nz - 2) {
                    u_ijkp2 = R_old[u_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
                    v_ijkp2 = R_old[v_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
                    w_ijkp2 = R_old[w_index + IDX(i, j, k + 2, Nx, Ny, Nz)];
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
                Tx  = (T_ip1jk - T_im1jk) / (2 * dx); // dT/dx
                Ty  = (T_ijp1k - T_ijm1k) / (2 * dy); // dT/dy
                Tz  = (T_ijkp1 - T_ijkm1) / (2 * dz); // dT/dz
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
                Txx = (T_ip1jk - 2 * T_ijk + T_im1jk) / (dx * dx); // d^2T/dx^2 
                Tyy = (T_ijp1k - 2 * T_ijk + T_ijm1k) / (dy * dy); // d^2T/dy^2
                Tzz = (T_ijkp1 - 2 * T_ijk + T_ijkm1) / (dz * dz); // d^2T/dz^2
                // Damping function
                tau_p = sqrt((nu * 0.5 * (uz + wx)) * (nu * 0.5 * (uz + wx)) + (nu * 0.5 * (vz + wy)) * (nu * 0.5 * (vz + wy)));
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
                // R_turbulence[parameters->turbulence_indexes.ux  + IDX(i, j, k, Nx, Ny, Nz)] = ux_uw;
                // R_turbulence[parameters->turbulence_indexes.uy  + IDX(i, j, k, Nx, Ny, Nz)] = uy_uw;
                // R_turbulence[parameters->turbulence_indexes.uz  + IDX(i, j, k, Nx, Ny, Nz)] = uz_uw;
                // R_turbulence[parameters->turbulence_indexes.vx  + IDX(i, j, k, Nx, Ny, Nz)] = vx_uw;
                // R_turbulence[parameters->turbulence_indexes.vy  + IDX(i, j, k, Nx, Ny, Nz)] = vy_uw;
                // R_turbulence[parameters->turbulence_indexes.vz  + IDX(i, j, k, Nx, Ny, Nz)] = vz_uw;
                // R_turbulence[parameters->turbulence_indexes.wx  + IDX(i, j, k, Nx, Ny, Nz)] = wx_uw;
                // R_turbulence[parameters->turbulence_indexes.wy  + IDX(i, j, k, Nx, Ny, Nz)] = wy_uw;
                // R_turbulence[parameters->turbulence_indexes.wz  + IDX(i, j, k, Nx, Ny, Nz)] = wz_uw;
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
                // if (k < Nz_Y) {
                //     Y_ijk = R_old[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)];
                //     S = source(T_ijk, Y_ijk, H_R, A, T_a, h, a_v, T_inf, c_p, rho);
                //     Y_RHS = -Y_f * Y_ijk * K(T_ijk, A, T_a) * H(T_ijk - T_pc);
                //     R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] = Y_RHS;
                //     // Compute force terms
                //     mod_U = sqrt(u_ijk * u_ijk + v_ijk * v_ijk + w_ijk * w_ijk);
                //     F_x = - Y_D * a_v * Y_ijk * mod_U * u_ijk;
                //     F_y = - Y_D * a_v * Y_ijk * mod_U * v_ijk;                             
                //     F_z = - g * (T_ijk - T_inf) / T_ijk - Y_D * a_v * Y_ijk * mod_U * w_ijk;
                // }
                if (k < Nz_Y) {
                    Y_ijk = R_old[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)];
                    Y_RHS = -Y_f * Y_ijk * K(T_ijk, A, T_a) * H(T_ijk - T_pc);
                    R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] = Y_RHS;
                } else {
                    Y_ijk = 0.0;
                }
                // Compute source and force terms
                S = source(T_ijk, Y_ijk, H_R, A, T_a, h, a_v, T_inf, c_p, rho);
                mod_U = sqrt(u_ijk * u_ijk + v_ijk * v_ijk + w_ijk * w_ijk);
                F_x = - Y_D * a_v * Y_ijk * mod_U * u_ijk;
                F_y = - Y_D * a_v * Y_ijk * mod_U * v_ijk;                             
                F_z = - g * (T_ijk - T_inf) / T_ijk - Y_D * a_v * Y_ijk * mod_U * w_ijk;
                // Compute RHS for velocity and temperature
                u_RHS = nu * (uxx + uyy + uzz) - (u_ijk * ux_uw + v_ijk * uy_uw + w_ijk * uz_uw) + F_x;
                v_RHS = nu * (vxx + vyy + vzz) - (u_ijk * vx_uw + v_ijk * vy_uw + w_ijk * vz_uw) + F_y;
                w_RHS = nu * (wxx + wyy + wzz) - (u_ijk * wx_uw + v_ijk * wy_uw + w_ijk * wz_uw) + F_z;
                T_RHS = alpha * (Txx + Tyy + Tzz) - (u_ijk * Tx + v_ijk * Ty + wz_uw * Tz) + S;
                // Save RHS into R_new
                R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] = u_RHS;
                R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] = v_RHS;
                R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] = w_RHS;
                R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = T_RHS;  
                // printf("u[%d,%d,%d] = %e\n", i, j, k, u_RHS);              
            }
        }
    }
    // Add turbulence terms
    turbulence(R_turbulence, R_new, parameters);
}

void boundary_conditions(double *R_old, double *R_new, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->Nz_Y;
    int u_index = parameters->field_indexes.u;
    int v_index = parameters->field_indexes.v;
    int w_index = parameters->field_indexes.w;
    int T_index = parameters->field_indexes.T;
    int Y_index = parameters->field_indexes.Y;
    double T_inf = parameters->T_inf;
    double u_ijkm1, u_ijkm2;
    double v_ijkm1, v_ijkm2;
    double w_ijkm1, w_ijkm2;
    double T_ijkp1, T_ijkm1, T_ijkm2, T_ijkp2;
    double Y_ijkp1, Y_ijkm1, Y_ijkm2, Y_ijkp2;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                // Periodic boundary conditions in x and y
                // Left boundary
                // if (i == 0) {
                //     R_new[u_index + IDX(0, j, k, Nx, Ny, Nz)] = R_old[u_index + IDX(Nx - 2, j, k, Nx, Ny, Nz)];
                //     R_new[v_index + IDX(0, j, k, Nx, Ny, Nz)] = R_old[v_index + IDX(Nx - 2, j, k, Nx, Ny, Nz)];
                //     R_new[w_index + IDX(0, j, k, Nx, Ny, Nz)] = R_old[w_index + IDX(Nx - 2, j, k, Nx, Ny, Nz)];
                //     R_new[T_index + IDX(0, j, k, Nx, Ny, Nz)] = R_old[T_index + IDX(Nx - 2, j, k, Nx, Ny, Nz)];
                //     if (k < Nz_Y) {
                //         R_new[Y_index + IDX(0, j, k, Nx, Ny, Nz_Y)] = R_old[Y_index + IDX(Nx - 2, j, k, Nx, Ny, Nz_Y)];
                //     }
                // }
                // // Right boundary
                // if (i == Nx - 1) {
                //     R_new[u_index + IDX(Nx - 1, j, k, Nx, Ny, Nz)] = R_old[u_index + IDX(1, j, k, Nx, Ny, Nz)];
                //     R_new[v_index + IDX(Nx - 1, j, k, Nx, Ny, Nz)] = R_old[v_index + IDX(1, j, k, Nx, Ny, Nz)];
                //     R_new[w_index + IDX(Nx - 1, j, k, Nx, Ny, Nz)] = R_old[w_index + IDX(1, j, k, Nx, Ny, Nz)];
                //     R_new[T_index + IDX(Nx - 1, j, k, Nx, Ny, Nz)] = R_old[T_index + IDX(1, j, k, Nx, Ny, Nz)];
                //     if (k < Nz_Y) {
                //         R_new[Y_index + IDX(Nx - 1, j, k, Nx, Ny, Nz_Y)] = R_old[Y_index + IDX(1, j, k, Nx, Ny, Nz_Y)];
                //     }
                // }
                // // Back boundary
                // if (j == 0) {
                //     R_new[u_index + IDX(i, 0, k, Nx, Ny, Nz)] = R_old[u_index + IDX(i, Ny - 2, k, Nx, Ny, Nz)];
                //     R_new[v_index + IDX(i, 0, k, Nx, Ny, Nz)] = R_old[v_index + IDX(i, Ny - 2, k, Nx, Ny, Nz)];
                //     R_new[w_index + IDX(i, 0, k, Nx, Ny, Nz)] = R_old[w_index + IDX(i, Ny - 2, k, Nx, Ny, Nz)];
                //     R_new[T_index + IDX(i, 0, k, Nx, Ny, Nz)] = R_old[T_index + IDX(i, Ny - 2, k, Nx, Ny, Nz)];
                //     if (k < Nz_Y) {
                //         R_new[Y_index + IDX(i, 0, k, Nx, Ny, Nz_Y)] = R_old[Y_index + IDX(i, Ny - 2, k, Nx, Ny, Nz_Y)];
                //     }
                // }
                // // Front boundary
                // if (j == Ny - 1) {
                //     R_new[u_index + IDX(i, Ny - 1, k, Nx, Ny, Nz)] = R_old[u_index + IDX(i, 1, k, Nx, Ny, Nz)];
                //     R_new[v_index + IDX(i, Ny - 1, k, Nx, Ny, Nz)] = R_old[v_index + IDX(i, 1, k, Nx, Ny, Nz)];
                //     R_new[w_index + IDX(i, Ny - 1, k, Nx, Ny, Nz)] = R_old[w_index + IDX(i, 1, k, Nx, Ny, Nz)];
                //     R_new[T_index + IDX(i, Ny - 1, k, Nx, Ny, Nz)] = R_old[T_index + IDX(i, 1, k, Nx, Ny, Nz)];
                //     if (k < Nz_Y) {
                //         R_new[Y_index + IDX(i, Ny - 1, k, Nx, Ny, Nz_Y)] = R_old[Y_index + IDX(i, 1, k, Nx, Ny, Nz_Y)];
                //     }
                // }
                // u = v = w = dT/dz = 0 at z = z_min and z = z_max using a second order approximation
                // du/dz = dv/dz = dw/dz = dY/dz = 0 at z = z_min and z = Y_height using a second order approximation
                // Bottom boundary
                if (k == 0) {
                    T_ijkp1 = R_old[T_index + IDX(i, j, 1, Nx, Ny, Nz)]; // T_{i,j,k+1}
                    T_ijkp2 = R_old[T_index + IDX(i, j, 2, Nx, Ny, Nz)]; // T_{i,j,k+2}
                    Y_ijkp1 = R_old[Y_index + IDX(i, j, 1, Nx, Ny, Nz_Y)]; // Y_{i,j,k+1}
                    Y_ijkp2 = R_old[Y_index + IDX(i, j, 2, Nx, Ny, Nz_Y)]; // Y_{i,j,k+2}
                    R_new[u_index + IDX(i, j, 0, Nx, Ny, Nz)] = 0; // u = 0
                    R_new[v_index + IDX(i, j, 0, Nx, Ny, Nz)] = 0; // v = 0
                    R_new[w_index + IDX(i, j, 0, Nx, Ny, Nz)] = 0; // w = 0
                    R_new[T_index + IDX(i, j, 0, Nx, Ny, Nz)] = (4 * T_ijkp1 - T_ijkp2) / 3; // dT/dz = 0
                    R_new[Y_index + IDX(i, j, 0, Nx, Ny, Nz_Y)] = (4 * Y_ijkp1 - Y_ijkp2) / 3; // dY/dz = 0
                }
                // Top boundary
                if (k == Nz - 1) {
                    u_ijkm1 = R_old[u_index + IDX(i, j, Nz - 2, Nx, Ny, Nz)]; // u_{i,j,k-1}
                    u_ijkm2 = R_old[u_index + IDX(i, j, Nz - 3, Nx, Ny, Nz)]; // u_{i,j,k-2}
                    v_ijkm1 = R_old[v_index + IDX(i, j, Nz - 2, Nx, Ny, Nz)]; // v_{i,j,k-1}
                    v_ijkm2 = R_old[v_index + IDX(i, j, Nz - 3, Nx, Ny, Nz)]; // v_{i,j,k-2}
                    w_ijkm1 = R_old[w_index + IDX(i, j, Nz - 2, Nx, Ny, Nz)]; // w_{i,j,k-1}
                    w_ijkm2 = R_old[w_index + IDX(i, j, Nz - 3, Nx, Ny, Nz)]; // w_{i,j,k-2}
                    T_ijkm1 = R_old[T_index + IDX(i, j, Nz - 2, Nx, Ny, Nz)]; // T_{i,j,k-1}
                    T_ijkm2 = R_old[T_index + IDX(i, j, Nz - 3, Nx, Ny, Nz)]; // T_{i,j,k-2}
                    R_new[u_index + IDX(i, j, Nz - 1, Nx, Ny, Nz)] = (4 * u_ijkm1 - u_ijkm2) / 3; // du/dz = 0
                    R_new[v_index + IDX(i, j, Nz - 1, Nx, Ny, Nz)] = (4 * v_ijkm1 - v_ijkm2) / 3; // dv/dz = 0
                    R_new[w_index + IDX(i, j, Nz - 1, Nx, Ny, Nz)] = (4 * w_ijkm1 - w_ijkm2) / 3; // dw/dz = 0
                    R_new[T_index + IDX(i, j, Nz - 1, Nx, Ny, Nz)] = (4 * T_ijkm1 - T_ijkm2) / 3; // dT/dz = 0
                    if (k < Nz_Y) {
                        Y_ijkm1 = R_old[Y_index + IDX(i, j, Nz - 2, Nx, Ny, Nz_Y)]; // Y_{i,j,k-1}
                        Y_ijkm2 = R_old[Y_index + IDX(i, j, Nz - 3, Nx, Ny, Nz_Y)]; // Y_{i,j,k-2}
                        R_new[Y_index + IDX(i, j, Nz - 1, Nx, Ny, Nz_Y)] = (4 * Y_ijkm1 - Y_ijkm2) / 3; // dY/dz = 0
                    }
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

void velocity_correction(double *R_new, double *R_old, double *p, double dt, Parameters *parameters) {
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
                    pijkp1 = p[IDX(i, j, k + 1, Nx, Ny, Nz)];
                    pijkm1 = p[IDX(i, j, k - 1, Nx, Ny, Nz)];
                    px = (pim1jk - pip1jk) / (2 * dx);
                    py = (pijm1k - pijp1k) / (2 * dy);
                    pz = (pijkm1 - pijkp1) / (2 * dz);
                    R_new[parameters->field_indexes.u + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * px;
                    R_new[parameters->field_indexes.v + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * py;
                    R_new[parameters->field_indexes.w + IDX(i, j, k, Nx, Ny, Nz)] -= dt / rho * pz;
                    // Temperature bounds
                    if (R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] < T_inf) {
                        R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = T_inf;
                    }
                    if (R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] > 1500) {
                        R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = 1500;
                    }
                    // Fuel bounds
                    if (k < Nz_Y) {
                        if (R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] < 0) {
                            R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] = 0;
                        }
                        if (R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] > 1) {
                            R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] = 1;
                        }
                    }
                } else { // Boundary points
                    if (k == 0) {
                        T_ijkp1 = R_old[T_index + IDX(i, j, 1, Nx, Ny, Nz)]; // T_{i,j,k+1}
                        T_ijkp2 = R_old[T_index + IDX(i, j, 2, Nx, Ny, Nz)]; // T_{i,j,k+2}
                        Y_ijkp1 = R_old[Y_index + IDX(i, j, 1, Nx, Ny, Nz_Y)]; // Y_{i,j,k+1}
                        Y_ijkp2 = R_old[Y_index + IDX(i, j, 2, Nx, Ny, Nz_Y)]; // Y_{i,j,k+2}
                        R_new[u_index + IDX(i, j, 0, Nx, Ny, Nz)] = 0; // u = 0
                        R_new[v_index + IDX(i, j, 0, Nx, Ny, Nz)] = 0; // v = 0
                        R_new[w_index + IDX(i, j, 0, Nx, Ny, Nz)] = 0; // w = 0
                        R_new[T_index + IDX(i, j, 0, Nx, Ny, Nz)] = (4 * T_ijkp1 - T_ijkp2) / 3; // dT/dz = 0
                        R_new[Y_index + IDX(i, j, 0, Nx, Ny, Nz_Y)] = (4 * Y_ijkp1 - Y_ijkp2) / 3; // dY/dz = 0
                    }
                    // Top boundary
                    if (k == Nz - 1) {
                        u_ijkm1 = R_old[u_index + IDX(i, j, Nz - 2, Nx, Ny, Nz)]; // u_{i,j,k-1}
                        u_ijkm2 = R_old[u_index + IDX(i, j, Nz - 3, Nx, Ny, Nz)]; // u_{i,j,k-2}
                        v_ijkm1 = R_old[v_index + IDX(i, j, Nz - 2, Nx, Ny, Nz)]; // v_{i,j,k-1}
                        v_ijkm2 = R_old[v_index + IDX(i, j, Nz - 3, Nx, Ny, Nz)]; // v_{i,j,k-2}
                        w_ijkm1 = R_old[w_index + IDX(i, j, Nz - 2, Nx, Ny, Nz)]; // w_{i,j,k-1}
                        w_ijkm2 = R_old[w_index + IDX(i, j, Nz - 3, Nx, Ny, Nz)]; // w_{i,j,k-2}
                        T_ijkm1 = R_old[T_index + IDX(i, j, Nz - 2, Nx, Ny, Nz)]; // T_{i,j,k-1}
                        T_ijkm2 = R_old[T_index + IDX(i, j, Nz - 3, Nx, Ny, Nz)]; // T_{i,j,k-2}
                        R_new[u_index + IDX(i, j, Nz - 1, Nx, Ny, Nz)] = (4 * u_ijkm1 - u_ijkm2) / 3; // du/dz = 0
                        R_new[v_index + IDX(i, j, Nz - 1, Nx, Ny, Nz)] = (4 * v_ijkm1 - v_ijkm2) / 3; // dv/dz = 0
                        R_new[w_index + IDX(i, j, Nz - 1, Nx, Ny, Nz)] = (4 * w_ijkm1 - w_ijkm2) / 3; // dw/dz = 0
                        R_new[T_index + IDX(i, j, Nz - 1, Nx, Ny, Nz)] = (4 * T_ijkm1 - T_ijkm2) / 3; // dT/dz = 0
                        if (k < Nz_Y) {
                            Y_ijkm1 = R_old[Y_index + IDX(i, j, Nz - 2, Nx, Ny, Nz_Y)]; // Y_{i,j,k-1}
                            Y_ijkm2 = R_old[Y_index + IDX(i, j, Nz - 3, Nx, Ny, Nz_Y)]; // Y_{i,j,k-2}
                            R_new[Y_index + IDX(i, j, Nz - 1, Nx, Ny, Nz_Y)] = (4 * Y_ijkm1 - Y_ijkm2) / 3; // dY/dz = 0
                        }
                    }
                }
            }
        }
    }
}

void euler(double t_n, double *y_n, double *y_np1, double *F, double *U_turbulence, double dt, int size, Parameters *parameters) {
    Phi(t_n, y_n, F, U_turbulence, parameters);
    for (int i = 0; i < size; i++) {        
        y_np1[i] = y_n[i] + dt * F[i];
    }
}

void rk4(double t_n, double *y_n, double *y_np1, double *k, double *F, double *U_turbulence,double dt, int size, Parameters *parameters) {
    int k1_index = 0;
    int k2_index = size;
    int k3_index = 2 * size;
    int k4_index = 3 * size;
    Phi(t_n, y_n, k + k1_index, U_turbulence, parameters);
    caxpy(F, k + k1_index, y_n, dt * 0.5, size);
    Phi(t_n + 0.5 * dt, F, k + k2_index, U_turbulence, parameters);
    caxpy(F, k + k2_index, y_n, dt * 0.5, size);
    Phi(t_n + 0.5 * dt, F, k + k3_index, U_turbulence, parameters);
    caxpy(F, k + k3_index, y_n, dt, size);
    Phi(t_n + dt, F, k + k4_index, U_turbulence, parameters);
    for (int i = 0; i < size; i++) {
        y_np1[i] = y_n[i] + dt / 6 * (k[k1_index + i] + 2 * k[k2_index + i] + 2 * k[k3_index + i] + k[k4_index + i]);
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
    double *k = (double *) malloc(4 * size * sizeof(double));
    double *R_turbulence = (double *) malloc(28 * Nx * Ny * Nz * sizeof(double));
    // Arrays for pressure Poisson Problem
    // double *p = (double *) malloc(Nx * Ny * Nz * sizeof(double));
    // double *p_top = (double *) malloc(Nx * Ny * Nz * sizeof(double));
    // double *f = (double *) malloc(Nx * Ny * Nz * sizeof(double)); 
    double complex *a = malloc((Nz - 2) * sizeof(double complex));
    double complex *b = malloc((Nz - 1) * sizeof(double complex));
    double complex *c = malloc((Nz - 2) * sizeof(double complex));
    double complex *d = malloc((Nz - 1) * sizeof(double complex));
    double complex *l = malloc((Nz - 2) * sizeof(double complex));
    double complex *u = malloc((Nz - 1) * sizeof(double complex));
    double complex *y = malloc((Nz - 1) * sizeof(double complex));
    double complex *pk = malloc((Nz - 1) * sizeof(double complex));
    fftw_complex *f_in = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    fftw_complex *f_out = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    fftw_complex *p_top_in = fftw_alloc_complex((Nx - 1) * (Ny - 1));
    fftw_complex *p_top_out = fftw_alloc_complex((Nx - 1) * (Ny - 1));
    fftw_complex *p_in = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    fftw_complex *p_out = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    fftw_plan p_plan, f_plan, p_top_plan;
    char *save_path = "data/output/";
    clock_t start, end;
    // copy_slice(p_top, p_0, Nz - 1, Nx, Ny, Nz);
    // FILE *logs = fopen("data/output/logs.txt", "w");

    // Time integration
    for (int n = 0; n < Nt; n++) { 
        start = clock();
        // Euler step to compute U^*, T^{n+1}, Y^{n+1}
        euler(t[n], y_n, y_np1, F, R_turbulence, dt, size, parameters);

        // RK4 step to compute U^{n+1}, T^{n+1}, Y^{n+1}
        // rk4(t[n], y_n, y_np1, k, F, R_turbulence, dt, size, parameters);

        // Boundary conditions
        // boundary_conditions(y_n, y_np1, parameters);

        // Solve Poisson problem for pressure (it only uses U^*)
        // solve_pressure_v1(y_np1, p, parameters);
        // solve_pressure(y_np1, p, a, b, c, d, l, u, y, pk, f_in, f_out, p_top_in, p_top_out, p_in, p_out, parameters);
        // solve_pressure(y_np1, p, a, b, c, d, l, u, y, pk, parameters);
        solve_pressure(y_np1, p, a, b, c, d, l, u, y, pk, p_plan, f_plan, p_top_plan, f_in, f_out, p_top_in, p_top_out, p_in, p_out, parameters);

        // Velocity correction
        velocity_correction(y_np1, y_n, p, dt, parameters);

        // Boundary conditions
        // boundary_conditions(y_n, y_np1, parameters); 
        
        end = clock();

        step_time = (double) (end - start) / CLOCKS_PER_SEC;
        
        // Save data each NT steps and at the last step
        if (n > 0 && (n % NT == 0 || n == Nt - 1)) {  
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

    // fclose(logs);
    // Free memory
    free(t);
    free(y_n);
    free(y_np1);
    free(F);
    free(R_turbulence);
    // free(p);
    // free(p_top);
    // free(f);
    // Free memory
    fftw_destroy_plan(p_top_plan);
    fftw_destroy_plan(f_plan);
    fftw_destroy_plan(p_plan);
    fftw_free(p_top_in);
    fftw_free(p_in);
    fftw_free(f_in);
    fftw_free(p_top_out);
    fftw_free(f_out);
    fftw_free(p_out);
    free(a);
    free(b);
    free(c);
    free(d);
    free(l);
    free(u);
    free(y);
}