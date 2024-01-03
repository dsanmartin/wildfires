#include "../include/pde.h"

void Phi(double t, double *R_old, double *R_new, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->k_Y_h + 1;
    double dx = parameters->dx;
    double dy = parameters->dy;
    double dz = parameters->dz;
    double *x = parameters->x;
    double *y = parameters->y;
    double *z = parameters->z;
    double alpha = parameters->alpha;
    double Y_h = parameters->Y_h;
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
    int u_index = 0;
    int v_index = Nx * Ny * Nz;
    int w_index = 2 * Nx * Ny * Nz;
    int T_index = 0;//3 * Nx * Ny * Nz * 0;
    int Y_index = Nx * Ny * Nz;
    double u_ijk, u_ip1jk, u_im1jk, u_ijp1k, u_ijm1k, u_ijkp1, u_ijkm1;
    double v_ijk, v_ip1jk, v_im1jk, v_ijp1k, v_ijm1k, v_ijkp1, v_ijkm1;
    double w_ijk, w_ip1jk, w_im1jk, w_ijp1k, w_ijm1k, w_ijkp1, w_ijkm1;
    double T_ijk, T_ip1jk, T_im1jk, T_ijp1k, T_ijm1k, T_ijkp1, T_ijkm1, T_ijkm2, T_ijkp2;
    double Y_ijk;
    double Tx, Ty, Tz, Txx, Tyy, Tzz;
    double T_RHS, Y_RHS, S = 0.0;    
    u_ijk = 1;
    v_ijk = 1;
    w_ijk = 0;
    // Loop over interior nodes
    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int k = 1; k < Nz - 1; k++) {
                // Get nodes
                T_ijk = R_old[T_index + idx(i, j, k, Nx, Ny, Nz)]; // T_{i,j,k}
                T_ip1jk = R_old[T_index + idx(i + 1, j, k, Nx, Ny, Nz)]; // T_{i+1,j,k}
                T_im1jk = R_old[T_index + idx(i - 1, j, k, Nx, Ny, Nz)]; // T_{i-1,j,k}
                T_ijp1k = R_old[T_index + idx(i, j + 1, k, Nx, Ny, Nz)]; // T_{i,j+1,k}
                T_ijm1k = R_old[T_index + idx(i, j - 1, k, Nx, Ny, Nz)]; // T_{i,j-1,k}
                T_ijkp1 = R_old[T_index + idx(i, j, k + 1, Nx, Ny, Nz)]; // T_{i,j,k+1}
                T_ijkm1 = R_old[T_index + idx(i, j, k - 1, Nx, Ny, Nz)]; // T_{i,j,k-1}
                // Get derivatives
                Tx = (T_ip1jk - T_im1jk) / (2 * dx); // dT/dx
                Ty = (T_ijp1k - T_ijm1k) / (2 * dy); // dT/dy
                Tz = (T_ijkp1 - T_ijkm1) / (2 * dz); // dT/dz
                Txx = (T_ip1jk - 2 * T_ijk + T_im1jk) / (dx * dx); // d^2T/dx^2 
                Tyy = (T_ijp1k - 2 * T_ijk + T_ijm1k) / (dy * dy); // d^2T/dy^2
                Tzz = (T_ijkp1 - 2 * T_ijk + T_ijkm1) / (dz * dz); // d^2T/dz^2
                if (k < Nz_Y) {
                    Y_ijk = R_old[Y_index + idx(i, j, k, Nx, Ny, Nz_Y)];
                    S = source(T_ijk, Y_ijk, H_R, A, T_a, h, a_v, T_inf, c_p, rho);
                    Y_RHS = -Y_f * Y_ijk * K(T_ijk, A, T_a) * H(T_ijk - T_pc);
                    R_new[Y_index + idx(i, j, k, Nx, Ny, Nz_Y)] = Y_RHS;
                }
                T_RHS = alpha * (Txx + Tyy + Tzz) - (u_ijk * Tx + v_ijk * Ty + w_ijk * Tz) + S;
                // Save RHS into R_new
                R_new[T_index + idx(i, j, k, Nx, Ny, Nz)] = T_RHS;
            }
        }
    }
}

void euler(double t_n, double *y_n, double *y_np1, double *F, double dt, int size, Parameters *parameters) {
    Phi(t_n, y_n, F, parameters);
    for (int i = 0; i < size; i++) {        
        y_np1[i] = y_n[i] + dt * F[i];
    }
}

void boundary_conditions(double *R_old, double *R_new, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->k_Y_h + 1;
    int T_index = 0;
    int Y_index = Nx * Ny * Nz;
    double T_inf = parameters->T_inf;
    double T_ijkp1, T_ijkm1, T_ijkm2, T_ijkp2;
    double Y_ijkp1, Y_ijkm1, Y_ijkm2, Y_ijkp2;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                // Periodic boundary conditions in x and y
                // Left boundary
                if (i == 0) {
                    R_new[T_index + idx(0, j, k, Nx, Ny, Nz)] = R_old[T_index + idx(Nx - 2, j, k, Nx, Ny, Nz)];
                    if (k < Nz_Y) {
                        R_new[Y_index + idx(0, j, k, Nx, Ny, Nz_Y)] = R_old[Y_index + idx(Nx - 2, j, k, Nx, Ny, Nz_Y)];
                    }
                }
                // Right boundary
                if (i == Nx - 1) {
                    R_new[T_index + idx(Nx - 1, j, k, Nx, Ny, Nz)] = R_old[T_index + idx(1, j, k, Nx, Ny, Nz)];
                    if (k < Nz_Y) {
                        R_new[Y_index + idx(Nx - 1, j, k, Nx, Ny, Nz_Y)] = R_old[Y_index + idx(1, j, k, Nx, Ny, Nz_Y)];
                    }
                }
                // Back boundary
                if (j == 0) {
                    R_new[T_index + idx(i, 0, k, Nx, Ny, Nz)] = R_old[T_index + idx(i, Ny - 2, k, Nx, Ny, Nz)];
                    if (k < Nz_Y) {
                        R_new[Y_index + idx(i, 0, k, Nx, Ny, Nz_Y)] = R_old[Y_index + idx(i, Ny - 2, k, Nx, Ny, Nz_Y)];
                    }
                }
                // Front boundary
                if (j == Ny - 1) {
                    R_new[T_index + idx(i, Ny - 1, k, Nx, Ny, Nz)] = R_old[T_index + idx(i, 1, k, Nx, Ny, Nz)];
                    if (k < Nz_Y) {
                        R_new[Y_index + idx(i, Ny - 1, k, Nx, Ny, Nz_Y)] = R_old[Y_index + idx(i, 1, k, Nx, Ny, Nz_Y)];
                    }
                }
                // dT/dz = 0 at z = z_min and z = z_max using a second order approximation
                // dY/dz = 0 at z = z_min and z = Y_height using a second order approximation
                // Bottom boundary
                if (k == 0) {
                    T_ijkp1 = R_old[T_index + idx(i, j, 1, Nx, Ny, Nz)]; // T_{i,j,k+1}
                    T_ijkp2 = R_old[T_index + idx(i, j, 2, Nx, Ny, Nz)]; // T_{i,j,k+2}
                    Y_ijkp1 = R_old[Y_index + idx(i, j, 1, Nx, Ny, Nz_Y)]; // Y_{i,j,k+1}
                    Y_ijkp2 = R_old[Y_index + idx(i, j, 2, Nx, Ny, Nz_Y)]; // Y_{i,j,k+2}
                    R_new[T_index + idx(i, j, 0, Nx, Ny, Nz)] = (4 * T_ijkp1 - T_ijkp2) / 3; // dT/dz = 0
                    R_new[Y_index + idx(i, j, 0, Nx, Ny, Nz_Y)] = (4 * Y_ijkp1 - Y_ijkp2) / 3; // dY/dz = 0
                }
                // Top boundary
                if (k == Nz - 1) {
                    T_ijkm1 = R_old[T_index + idx(i, j, Nz - 2, Nx, Ny, Nz)]; // T_{i,j,k-1}
                    T_ijkm2 = R_old[T_index + idx(i, j, Nz - 3, Nx, Ny, Nz)]; // T_{i,j,k-2}
                    R_new[T_index + idx(i, j, Nz - 1, Nx, Ny, Nz)] = (4 * T_ijkm1 - T_ijkm2) / 3; // dT/dz = 0
                    if (k < Nz_Y) {
                        Y_ijkm1 = R_old[Y_index + idx(i, j, Nz - 2, Nx, Ny, Nz_Y)]; // Y_{i,j,k-1}
                        Y_ijkm2 = R_old[Y_index + idx(i, j, Nz - 3, Nx, Ny, Nz_Y)]; // Y_{i,j,k-2}
                        R_new[Y_index + idx(i, j, Nz - 1, Nx, Ny, Nz_Y)] = (4 * Y_ijkm1 - Y_ijkm2) / 3; // dY/dz = 0
                    }
                }
                // Check if Y_ijk is negative
                if (k < Nz_Y) {
                    if (R_new[Y_index + idx(i, j, k, Nx, Ny, Nz_Y)] < 0) {
                        R_new[Y_index + idx(i, j, k, Nx, Ny, Nz_Y)] = 0;
                    }
                }
                // Check if T_ijk is less than T_inf
                if (R_new[T_index + idx(i, j, k, Nx, Ny, Nz)] < T_inf) {
                    R_new[T_index + idx(i, j, k, Nx, Ny, Nz)] = T_inf;
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

void create_y_0(double *T, double *Y, double *y_0, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->k_Y_h + 1;
    int size = Nx * Ny * Nz;
    int size_Y = Nx * Ny * Nz_Y;
    for (int i = 0; i < size; i++) {
        y_0[i] = T[i];
    }
    for (int i = 0; i < size_Y; i++) {
        y_0[size + i] = Y[i];
    }

}

void solve_PDE(double *y_n, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nt = parameters->Nt;
    int NT = parameters->NT;
    int Nz_Y = parameters->k_Y_h + 1;
    int size = Nx * Ny * Nz + Nx * Ny * Nz_Y;
    int n_save;
    double *t = parameters->t;
    double dt = parameters->dt;
    double *y_np1 = (double *) malloc(size * sizeof(double));
    double *F = (double *) malloc(size * sizeof(double));
    char *filename = (char *) malloc(100 * sizeof(char));
    char *save_path = "data/output/";


    // Time integration
    for (int n = 0; n < Nt; n++) { 
        // Euler step
        euler(t[n], y_n, y_np1, F, dt, size, parameters);

        // Boundary conditions
        boundary_conditions(y_n, y_np1, parameters);
        
        // Save data each NT steps and at the last step
        if (n > 0 && (n % NT == 0 || n == Nt - 1)) {  
            n_save = n / NT;
            printf("n = %d, t_n = %lf\n", n, t[n]);          
            sprintf(filename, "data/output/T.csv.%d", n_save);
            printf("Saving files...\n");
            // save_data(filename, parameters->x, parameters->y, parameters->z, y_np1, Nx, Ny, Nz);   
            save_data(save_path, y_np1, n_save, parameters);
        }

        // Update
        copy(y_n, y_np1, size);

    }

    // Free memory
    free(t);
    free(y_np1);
    free(F);
    free(filename);
}