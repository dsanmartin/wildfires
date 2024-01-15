#include "../include/pde.h"
#include "../include/poisson.h"

void Phi(double t, double *R_old, double *R_new, double *U_turbulence, double *P_RHS, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->k_Y_h + 1;
    double dx = parameters->dx;
    double dy = parameters->dy;
    double dz = parameters->dz;
    double dt = parameters->dt;
    // double *x = parameters->x;
    // double *y = parameters->y;
    // double *z = parameters->z;
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
    // Fields indexes
    int u_index = 0;
    int v_index = Nx * Ny * Nz;
    int w_index = 2 * Nx * Ny * Nz;
    int T_index = 3 * Nx * Ny * Nz;
    int Y_index = 4 * Nx * Ny * Nz;
    // Fields nodes
    double u_ijk, u_ip1jk, u_im1jk, u_ip2jk, u_im2jk, u_ijp1k, u_ijm1k, u_ijp2k, u_ijm2k, u_ijkp1, u_ijkm1, u_ijkp2, u_ijkm2;
    double v_ijk, v_ip1jk, v_im1jk, v_ip2jk, v_im2jk, v_ijp1k, v_ijm1k, v_ijp2k, v_ijm2k, v_ijkp1, v_ijkm1, v_ijkp2, v_ijkm2;
    double w_ijk, w_ip1jk, w_im1jk, w_ip2jk, w_im2jk, w_ijp1k, w_ijm1k, w_ijp2k, w_ijm2k, w_ijkp1, w_ijkm1, w_ijkp2, w_ijkm2;
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
    double fx = 0.0, fy = 0.0, fz = 0.0;
    // u_ijk = 1;
    // v_ijk = 1;
    // w_ijk = 0;
    // Loop over interior nodes
    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int k = 1; k < Nz - 1; k++) {
                /* Get fields nodes */
                u_ijk   = R_old[u_index + IDX(i, j, k, Nx, Ny, Nz)]; // u_{i,j,k}
                u_ip1jk = R_old[u_index + IDX(i + 1, j, k, Nx, Ny, Nz)]; // u_{i+1,j,k}
                u_im1jk = R_old[u_index + IDX(i - 1, j, k, Nx, Ny, Nz)]; // u_{i-1,j,k}
                u_ip2jk = R_old[u_index + IDX(i + 2, j, k, Nx, Ny, Nz)]; // u_{i+2,j,k}
                u_im2jk = R_old[u_index + IDX(i - 2, j, k, Nx, Ny, Nz)]; // u_{i-2,j,k}
                u_ijp1k = R_old[u_index + IDX(i, j + 1, k, Nx, Ny, Nz)]; // u_{i,j+1,k}
                u_ijm1k = R_old[u_index + IDX(i, j - 1, k, Nx, Ny, Nz)]; // u_{i,j-1,k}
                u_ijp2k = R_old[u_index + IDX(i, j + 2, k, Nx, Ny, Nz)]; // u_{i,j+2,k}
                u_ijm2k = R_old[u_index + IDX(i, j - 2, k, Nx, Ny, Nz)]; // u_{i,j-2,k}
                u_ijkp1 = R_old[u_index + IDX(i, j, k + 1, Nx, Ny, Nz)]; // u_{i,j,k+1}
                u_ijkm1 = R_old[u_index + IDX(i, j, k - 1, Nx, Ny, Nz)]; // u_{i,j,k-1}
                u_ijkp2 = R_old[u_index + IDX(i, j, k + 2, Nx, Ny, Nz)]; // u_{i,j,k+2}
                u_ijkm2 = R_old[u_index + IDX(i, j, k - 2, Nx, Ny, Nz)]; // u_{i,j,k-2}
                v_ijk   = R_old[v_index + IDX(i, j, k, Nx, Ny, Nz)]; // v_{i,j,k}
                v_ip1jk = R_old[v_index + IDX(i + 1, j, k, Nx, Ny, Nz)]; // v_{i+1,j,k}
                v_im1jk = R_old[v_index + IDX(i - 1, j, k, Nx, Ny, Nz)]; // v_{i-1,j,k}
                v_ip2jk = R_old[v_index + IDX(i + 2, j, k, Nx, Ny, Nz)]; // v_{i+2,j,k}
                v_im2jk = R_old[v_index + IDX(i - 2, j, k, Nx, Ny, Nz)]; // v_{i-2,j,k}
                v_ijp1k = R_old[v_index + IDX(i, j + 1, k, Nx, Ny, Nz)]; // v_{i,j+1,k}
                v_ijm1k = R_old[v_index + IDX(i, j - 1, k, Nx, Ny, Nz)]; // v_{i,j-1,k}
                v_ijp2k = R_old[v_index + IDX(i, j + 2, k, Nx, Ny, Nz)]; // v_{i,j+2,k}
                v_ijm2k = R_old[v_index + IDX(i, j - 2, k, Nx, Ny, Nz)]; // v_{i,j-2,k}
                v_ijkp1 = R_old[v_index + IDX(i, j, k + 1, Nx, Ny, Nz)]; // v_{i,j,k+1}
                v_ijkm1 = R_old[v_index + IDX(i, j, k - 1, Nx, Ny, Nz)]; // v_{i,j,k-1}
                v_ijkp2 = R_old[v_index + IDX(i, j, k + 2, Nx, Ny, Nz)]; // v_{i,j,k+2}
                v_ijkm2 = R_old[v_index + IDX(i, j, k - 2, Nx, Ny, Nz)]; // v_{i,j,k-2}
                w_ijk   = R_old[w_index + IDX(i, j, k, Nx, Ny, Nz)]; // w_{i,j,k}
                w_ip1jk = R_old[w_index + IDX(i + 1, j, k, Nx, Ny, Nz)]; // w_{i+1,j,k}
                w_im1jk = R_old[w_index + IDX(i - 1, j, k, Nx, Ny, Nz)]; // w_{i-1,j,k}
                w_ip2jk = R_old[w_index + IDX(i + 2, j, k, Nx, Ny, Nz)]; // w_{i+2,j,k}
                w_im2jk = R_old[w_index + IDX(i - 2, j, k, Nx, Ny, Nz)]; // w_{i-2,j,k}
                w_ijp1k = R_old[w_index + IDX(i, j + 1, k, Nx, Ny, Nz)]; // w_{i,j+1,k}
                w_ijm1k = R_old[w_index + IDX(i, j - 1, k, Nx, Ny, Nz)]; // w_{i,j-1,k}
                w_ijp2k = R_old[w_index + IDX(i, j + 2, k, Nx, Ny, Nz)]; // w_{i,j+2,k}
                w_ijm2k = R_old[w_index + IDX(i, j - 2, k, Nx, Ny, Nz)]; // w_{i,j-2,k}
                w_ijkp1 = R_old[w_index + IDX(i, j, k + 1, Nx, Ny, Nz)]; // w_{i,j,k+1}
                w_ijkm1 = R_old[w_index + IDX(i, j, k - 1, Nx, Ny, Nz)]; // w_{i,j,k-1}
                w_ijkp2 = R_old[w_index + IDX(i, j, k + 2, Nx, Ny, Nz)]; // w_{i,j,k+2}
                w_ijkm2 = R_old[w_index + IDX(i, j, k - 2, Nx, Ny, Nz)]; // w_{i,j,k-2} 
                T_ijk   = R_old[T_index + IDX(i, j, k, Nx, Ny, Nz)]; // T_{i,j,k}
                T_ip1jk = R_old[T_index + IDX(i + 1, j, k, Nx, Ny, Nz)]; // T_{i+1,j,k}
                T_im1jk = R_old[T_index + IDX(i - 1, j, k, Nx, Ny, Nz)]; // T_{i-1,j,k}
                T_ijp1k = R_old[T_index + IDX(i, j + 1, k, Nx, Ny, Nz)]; // T_{i,j+1,k}
                T_ijm1k = R_old[T_index + IDX(i, j - 1, k, Nx, Ny, Nz)]; // T_{i,j-1,k}
                T_ijkp1 = R_old[T_index + IDX(i, j, k + 1, Nx, Ny, Nz)]; // T_{i,j,k+1}
                T_ijkm1 = R_old[T_index + IDX(i, j, k - 1, Nx, Ny, Nz)]; // T_{i,j,k-1}
                /* Computing upwind scheme terms */
                u_plu = MAX(u_ijk, 0.0);
                u_min = MIN(u_ijk, 0.0);
                v_plu = MAX(v_ijk, 0.0);
                v_min = MIN(v_ijk, 0.0);
                w_plu = MAX(w_ijk, 0.0);
                w_min = MIN(w_ijk, 0.0);
                // Second order forward difference at i=1
                if (i == 1) {
                    u_im = (-3 * u_ijk + 4 * u_ip1jk - u_ip2jk) / (2 * dx);
                    v_im = (-3 * v_ijk + 4 * v_ip1jk - v_ip2jk) / (2 * dx);
                    w_im = (-3 * w_ijk + 4 * w_ip1jk - w_ip2jk) / (2 * dx);
                } else {
                    u_im = (-3 * u_ijk + 4 * u_im1jk - u_im2jk) / (2 * dx);
                    v_im = (-3 * v_ijk + 4 * v_im1jk - v_im2jk) / (2 * dx);
                    w_im = (-3 * w_ijk + 4 * w_im1jk - w_im2jk) / (2 * dx);
                }
                // Second order backward difference at i=Nx-2
                if (i == Nx - 2) {
                    u_ip = (3 * u_ijk - 4 * u_im1jk + u_im2jk) / (2 * dx);
                    v_ip = (3 * v_ijk - 4 * v_im1jk + v_im2jk) / (2 * dx);
                    w_ip = (3 * w_ijk - 4 * w_im1jk + w_im2jk) / (2 * dx);
                } else {
                    u_ip = (3 * u_ijk - 4 * u_ip1jk + u_ip2jk) / (2 * dx);
                    v_ip = (3 * v_ijk - 4 * v_ip1jk + v_ip2jk) / (2 * dx);
                    w_ip = (3 * w_ijk - 4 * w_ip1jk + w_ip2jk) / (2 * dx);
                }
                // Second order forward difference at j=1
                if (j == 1) {
                    u_jm = (-3 * u_ijk + 4 * u_ijp1k - u_ijp2k) / (2 * dy);
                    v_jm = (-3 * v_ijk + 4 * v_ijp1k - v_ijp2k) / (2 * dy);
                    w_jm = (-3 * w_ijk + 4 * w_ijp1k - w_ijp2k) / (2 * dy);
                } else {
                    u_jm = (-3 * u_ijk + 4 * u_ijm1k - u_ijm2k) / (2 * dy);
                    v_jm = (-3 * v_ijk + 4 * v_ijm1k - v_ijm2k) / (2 * dy);
                    w_jm = (-3 * w_ijk + 4 * w_ijm1k - w_ijm2k) / (2 * dy);
                }
                // Second order backward difference at j=Ny-2
                if (j == Ny - 2) {
                    u_jp = (3 * u_ijk - 4 * u_ijm1k + u_ijm2k) / (2 * dy);
                    v_jp = (3 * v_ijk - 4 * v_ijm1k + v_ijm2k) / (2 * dy);
                    w_jp = (3 * w_ijk - 4 * w_ijm1k + w_ijm2k) / (2 * dy);
                } else {
                    u_jp = (3 * u_ijk - 4 * u_ijp1k + u_ijp2k) / (2 * dy);
                    v_jp = (3 * v_ijk - 4 * v_ijp1k + v_ijp2k) / (2 * dy);
                    w_jp = (3 * w_ijk - 4 * w_ijp1k + w_ijp2k) / (2 * dy);
                }
                // Second order forward difference at k=1
                if (k == 1) {
                    u_km = (-3 * u_ijk + 4 * u_ijkp1 - u_ijkp2) / (2 * dz);
                    v_km = (-3 * v_ijk + 4 * v_ijkp1 - v_ijkp2) / (2 * dz);
                    w_km = (-3 * w_ijk + 4 * w_ijkp1 - w_ijkp2) / (2 * dz);
                } else {
                    u_km = (-3 * u_ijk + 4 * u_ijkm1 - u_ijkm2) / (2 * dz);
                    v_km = (-3 * v_ijk + 4 * v_ijkm1 - v_ijkm2) / (2 * dz);
                    w_km = (-3 * w_ijk + 4 * w_ijkm1 - w_ijkm2) / (2 * dz);
                }
                // Second order backward difference at k=Nz-2
                if (k == Nz - 2) {
                    u_kp = (3 * u_ijk - 4 * u_ijkm1 + u_ijkm2) / (2 * dz);
                    v_kp = (3 * v_ijk - 4 * v_ijkm1 + v_ijkm2) / (2 * dz);
                    w_kp = (3 * w_ijk - 4 * w_ijkm1 + w_ijkm2) / (2 * dz);
                } else {
                    u_kp = (3 * u_ijk - 4 * u_ijkp1 + u_ijkp2) / (2 * dz);
                    v_kp = (3 * v_ijk - 4 * v_ijkp1 + v_ijkp2) / (2 * dz);
                    w_kp = (3 * w_ijk - 4 * w_ijkp1 + w_ijkp2) / (2 * dz);
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
                /* Compute fuel and source term */
                if (k < Nz_Y) {
                    Y_ijk = R_old[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)];
                    S = source(T_ijk, Y_ijk, H_R, A, T_a, h, a_v, T_inf, c_p, rho);
                    Y_RHS = -Y_f * Y_ijk * K(T_ijk, A, T_a) * H(T_ijk - T_pc);
                    R_new[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)] = Y_RHS;
                }
                // Compute RHS for velocity and temperature
                u_RHS = nu * (uxx + uyy + uzz) + fx; // - (u_ijk * ux + v_ijk * uy + w_ijk * uz) + fx;
                v_RHS = nu * (vxx + vyy + vzz) + fy; // - (u_ijk * vx + v_ijk * vy + w_ijk * vz) + fy;
                w_RHS = nu * (wxx + wyy + wzz) + fz; // - (u_ijk * wx + v_ijk * wy + w_ijk * wz) + fz;
                T_RHS = alpha * (Txx + Tyy + Tzz) - (u_ijk * Tx + v_ijk * Ty + w_ijk * Tz) + S;
                // Save RHS into R_new
                R_new[u_index + IDX(i, j, k, Nx, Ny, Nz)] = u_RHS;
                R_new[v_index + IDX(i, j, k, Nx, Ny, Nz)] = v_RHS;
                R_new[w_index + IDX(i, j, k, Nx, Ny, Nz)] = w_RHS;
                R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = T_RHS;
            }
        }
    }
}

void euler(double t_n, double *y_n, double *y_np1, double *F, double *U_turbulence, double *P_RHS, double dt, int size, Parameters *parameters) {
    Phi(t_n, y_n, F, U_turbulence, P_RHS, parameters);
    for (int i = 0; i < size; i++) {        
        y_np1[i] = y_n[i] + dt * F[i];
    }
}

void boundary_conditions(double *R_old, double *R_new, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->k_Y_h + 1;
    int u_index = 0;
    int v_index = Nx * Ny * Nz;
    int w_index = 2 * Nx * Ny * Nz;
    int T_index = 3 * Nx * Ny * Nz;
    int Y_index = 4 * Nx * Ny * Nz;
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
                if (i == 0) {
                    R_new[u_index + IDX(0, j, k, Nx, Ny, Nz)] = R_old[u_index + IDX(Nx - 2, j, k, Nx, Ny, Nz)];
                    R_new[v_index + IDX(0, j, k, Nx, Ny, Nz)] = R_old[v_index + IDX(Nx - 2, j, k, Nx, Ny, Nz)];
                    R_new[w_index + IDX(0, j, k, Nx, Ny, Nz)] = R_old[w_index + IDX(Nx - 2, j, k, Nx, Ny, Nz)];
                    R_new[T_index + IDX(0, j, k, Nx, Ny, Nz)] = R_old[T_index + IDX(Nx - 2, j, k, Nx, Ny, Nz)];
                    if (k < Nz_Y) {
                        R_new[Y_index + IDX(0, j, k, Nx, Ny, Nz_Y)] = R_old[Y_index + IDX(Nx - 2, j, k, Nx, Ny, Nz_Y)];
                    }
                }
                // Right boundary
                if (i == Nx - 1) {
                    R_new[u_index + IDX(Nx - 1, j, k, Nx, Ny, Nz)] = R_old[u_index + IDX(1, j, k, Nx, Ny, Nz)];
                    R_new[v_index + IDX(Nx - 1, j, k, Nx, Ny, Nz)] = R_old[v_index + IDX(1, j, k, Nx, Ny, Nz)];
                    R_new[w_index + IDX(Nx - 1, j, k, Nx, Ny, Nz)] = R_old[w_index + IDX(1, j, k, Nx, Ny, Nz)];
                    R_new[T_index + IDX(Nx - 1, j, k, Nx, Ny, Nz)] = R_old[T_index + IDX(1, j, k, Nx, Ny, Nz)];
                    if (k < Nz_Y) {
                        R_new[Y_index + IDX(Nx - 1, j, k, Nx, Ny, Nz_Y)] = R_old[Y_index + IDX(1, j, k, Nx, Ny, Nz_Y)];
                    }
                }
                // Back boundary
                if (j == 0) {
                    R_new[u_index + IDX(i, 0, k, Nx, Ny, Nz)] = R_old[u_index + IDX(i, Ny - 2, k, Nx, Ny, Nz)];
                    R_new[v_index + IDX(i, 0, k, Nx, Ny, Nz)] = R_old[v_index + IDX(i, Ny - 2, k, Nx, Ny, Nz)];
                    R_new[w_index + IDX(i, 0, k, Nx, Ny, Nz)] = R_old[w_index + IDX(i, Ny - 2, k, Nx, Ny, Nz)];
                    R_new[T_index + IDX(i, 0, k, Nx, Ny, Nz)] = R_old[T_index + IDX(i, Ny - 2, k, Nx, Ny, Nz)];
                    if (k < Nz_Y) {
                        R_new[Y_index + IDX(i, 0, k, Nx, Ny, Nz_Y)] = R_old[Y_index + IDX(i, Ny - 2, k, Nx, Ny, Nz_Y)];
                    }
                }
                // Front boundary
                if (j == Ny - 1) {
                    R_new[u_index + IDX(i, Ny - 1, k, Nx, Ny, Nz)] = R_old[u_index + IDX(i, 1, k, Nx, Ny, Nz)];
                    R_new[v_index + IDX(i, Ny - 1, k, Nx, Ny, Nz)] = R_old[v_index + IDX(i, 1, k, Nx, Ny, Nz)];
                    R_new[w_index + IDX(i, Ny - 1, k, Nx, Ny, Nz)] = R_old[w_index + IDX(i, 1, k, Nx, Ny, Nz)];
                    R_new[T_index + IDX(i, Ny - 1, k, Nx, Ny, Nz)] = R_old[T_index + IDX(i, 1, k, Nx, Ny, Nz)];
                    if (k < Nz_Y) {
                        R_new[Y_index + IDX(i, Ny - 1, k, Nx, Ny, Nz_Y)] = R_old[Y_index + IDX(i, 1, k, Nx, Ny, Nz_Y)];
                    }
                }
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
                // Check if T_ijk is less than T_inf
                if (R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] < T_inf) {
                    R_new[T_index + IDX(i, j, k, Nx, Ny, Nz)] = T_inf;
                }
            }        
        }
    }
}

// void rk4(float t_n, float *y_n, float *y_np1, float *k1, float *k2, float *k3, float *k4, float dt, int size, Parameters *parameters) {
//     Phi(t_n, y_n, k1, parameters);
//     axpy(k2, )
//     Phi(t_n + dt / 2, y_n + dt / 2 * k1, k2, parameters);
//     Phi(t_n + dt / 2, y_n + dt / 2 * k2, k3, parameters);
//     Phi(t_n + dt, y_n + dt * k3, k4, parameters);
//     for (int i = 0; i < size; i++) {
//         y_np1[i] = y_n[i] + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
//     }
// }

void create_y_0(double *u, double *v, double *w, double *T, double *Y, double *y_0, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->k_Y_h + 1;
    int size = Nx * Ny * Nz;
    int size_Y = Nx * Ny * Nz_Y;
    for (int i = 0; i < size; i++) {
        y_0[i] = u[i];
        y_0[size + i] = v[i];
        y_0[2 * size + i] = w[i];
        y_0[3 * size + i] = T[i];
        if (i < size_Y) {
            y_0[4 * size + i] = Y[i];
        }
    }
}

void solve_PDE(double *y_n, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nt = parameters->Nt;
    int NT = parameters->NT;
    int Nz_Y = parameters->k_Y_h + 1;
    int size = 4 * Nx * Ny * Nz + Nx * Ny * Nz_Y;
    int n_save;
    double *t = parameters->t;
    double dt = parameters->dt;
    double *y_np1 = (double *) malloc(size * sizeof(double));
    double *F = (double *) malloc(size * sizeof(double));
    double *U_turbulence = (double *) malloc(45 * Nx * Ny * Nz * sizeof(double));
    // Pressure RHS for Poisson Problem
    double *P_RHS = (double *) malloc(Nx * Ny * Nz * sizeof(double)); 
    char *save_path = "data/output/";

    // Time integration
    for (int n = 0; n < Nt; n++) { 
        // Euler step to compute U^*, T^{n+1}, Y^{n+1}
        euler(t[n], y_n, y_np1, F, U_turbulence, P_RHS, dt, size, parameters);

        // Solve Poisson problem for pressure (it only uses U^*)
        solve_pressure(y_np1, P_RHS, parameters);

        // Boundary conditions
        boundary_conditions(y_n, y_np1, parameters);
        
        // Save data each NT steps and at the last step
        if (n > 0 && (n % NT == 0 || n == Nt - 1)) {  
            n_save = n / NT;
            printf("n = %d, t_n = %lf\n", n, t[n]);          
            printf("Saving files...\n");
            // save_data(filename, parameters->x, parameters->y, parameters->z, y_np1, Nx, Ny, Nz);   
            save_data(save_path, y_np1, n_save, parameters);
        }

        // Update y_n
        copy(y_n, y_np1, size);

    }

    // Free memory
    free(t);
    free(y_np1);
    free(F);
    free(U_turbulence);
    free(P_RHS);
}