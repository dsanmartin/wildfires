/**
 * @file solver.c
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Functions for solving the partial differential equations of the wildfire simulation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "../../include/c/solver.h"

void euler_step(double dt, double *y_n, double *y_np1, double *F, int size) {
    for (int i = 0; i < size; i++) {
        y_np1[i] = y_n[i] + dt * F[i];
    }
}

void euler(double t_n, double *y_n, double *y_np1, double *F, double *U_turbulence, double dt, int size, Parameters *parameters) {
    Phi(t_n, y_n, F, U_turbulence, parameters);
    euler_step(dt, y_n, y_np1, F, size);
}

void RK2_step(double dt, double *y_n, double *y_np1, double *k1, double *k2, int size) {
    for (int i = 0; i < size; i++) {
        y_np1[i] = y_n[i] + 0.5 * dt * (k1[i] + k2[i]);
    }
}

void RK2(double t_n, double *y_n, double *y_np1, double *k, double *F, double *U_turbulence, double dt, int size, Parameters *parameters) {
    int k1_index = 0;
    int k2_index = size;
    Phi(t_n, y_n, k + k1_index, U_turbulence, parameters);
    caxpy(F, k + k1_index, y_n, dt, size);
    Phi(t_n + dt, F, k + k2_index, U_turbulence, parameters);
    RK2_step(dt, y_n, y_np1, k + k1_index, k + k2_index, size);
}

void RK4_step(double dt, double *y_n, double *y_np1, double *k1, double *k2, double *k3, double *k4, int size) {
    for (int i = 0; i < size; i++) {
        y_np1[i] = y_n[i] +  (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
}

void RK4(double t_n, double *y_n, double *y_np1, double *k, double *F, double *U_turbulence, double dt, int size, Parameters *parameters) {
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
    RK4_step(dt, y_n, y_np1, k + k1_index, k + k2_index, k + k3_index, k + k4_index, size);
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
    // char log_path[100];
    int n_save;
    int k_size = 2;
    if (strncmp(parameters->method, "RK4", 3) == 0) 
        k_size = 4;
    double step_time, solver_time;
    double CFL = 0.0, T_min = 1e9, T_max = -1e9, Y_min = 1e9, Y_max = -1e9;
    // double *CFL, *T_min, *T_max, *Y_min, *Y_max;
    double *t = parameters->t;
    double dt = parameters->dt;
    double *y_np1 = (double *) malloc(size * sizeof(double));
    double *F = (double *) malloc(size * sizeof(double));
    double *k = (double *) malloc(k_size * size * sizeof(double));
    double *R_turbulence = (double *) malloc(25 * Nx * Ny * Nz * sizeof(double));
    // clock_t start, end, step_start, step_end; // Timers
    struct timeval start_solver, end_solver, start_ts, end_ts;
    char solver_time_message[256];
    // char min_max_values_message[256];
    char formatted_time[64];
    // Arrays for pressure Poisson Problem
    fftw_complex *a = fftw_alloc_complex((Nz - 2));
    fftw_complex *b = fftw_alloc_complex((Nz - 1));
    fftw_complex *c = fftw_alloc_complex((Nz - 2));
    fftw_complex *d = fftw_alloc_complex((Nz - 1));
    fftw_complex *l = fftw_alloc_complex((Nz - 2));
    fftw_complex *u = fftw_alloc_complex((Nz - 1));
    fftw_complex *y = fftw_alloc_complex((Nz - 1));
    fftw_complex *pk = fftw_alloc_complex((Nz - 1));
    fftw_complex *f_in = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    fftw_complex *f_out = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    fftw_complex *p_top_in = fftw_alloc_complex((Nx - 1) * (Ny - 1));
    fftw_complex *p_top_out = fftw_alloc_complex((Nx - 1) * (Ny - 1));
    fftw_complex *p_in = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    fftw_complex *p_out = fftw_alloc_complex((Nx - 1) * (Ny - 1) * (Nz - 1));
    fftw_plan p_plan, f_plan, p_top_plan;
    p_plan = NULL;
    f_plan = NULL;
    p_top_plan = NULL;
    // Solver time
    // start = clock();
    gettimeofday(&start_solver, NULL);
    // Time integration
    for (int n = 1; n <= Nt; n++) { 
        // step_start = clock(); // Start step timer
        gettimeofday(&start_ts, NULL);
        // Compute U^*, T^{n+1}, Y^{n+1}
        // Check time integration method
        if (strncmp(parameters->method, "Euler", 5) == 0) {
            euler(t[n], y_n, y_np1, F, R_turbulence, dt, size, parameters);
        } else if (strncmp(parameters->method, "RK2", 3) == 0) {
            RK2(t[n], y_n, y_np1, k, F, R_turbulence, dt, size, parameters);
        } else if (strncmp(parameters->method, "RK4", 3) == 0) {
            RK4(t[n], y_n, y_np1, k, F, R_turbulence, dt, size, parameters);
        } else {
            log_message(parameters, "Time integration method not found.");
            exit(1);
        }
        // Solve Poisson problem for pressure (it only uses U^*)
        solve_pressure(y_np1, p, a, b, c, d, l, u, y, pk, p_plan, f_plan, p_top_plan, f_in, f_out, p_top_in, p_top_out, p_in, p_out, parameters);
        // Chorin's projection method
        velocity_correction_fw(y_np1, p, dt, parameters);
        // Boundary conditions
        boundary_conditions(y_np1, parameters);
        // Bounds
        bounds(y_np1, parameters);
        // End step timer
        // step_end = clock(); 
        gettimeofday(&end_ts, NULL);
        // Compute step time
        // step_time = (double) (step_end - step_start) / CLOCKS_PER_SEC;
        step_time = ((end_ts.tv_sec  - start_ts.tv_sec) * 1000000u + end_ts.tv_usec - start_ts.tv_usec) / 1.e6;
        // Save data each NT steps and at the last step
        if (n % NT == 0 || n == Nt) {  
            n_save = n / NT;
            timestep_reports(y_np1, &CFL, &Y_min, &Y_max, &T_min, &T_max, parameters);
            log_timestep(parameters, n, t[n], step_time, CFL, T_min, T_max, Y_min, Y_max);
            save_data(y_np1, p, n_save, parameters);
        }
        // Update y_n
        copy(y_n, y_np1, size);
    }
    // end = clock();
    gettimeofday(&end_solver, NULL);
    // printf("Solver time: %lf s\n", (double) (end - start) / CLOCKS_PER_SEC);
    solver_time = ((end_solver.tv_sec  - start_solver.tv_sec) * 1000000u + end_solver.tv_usec - start_solver.tv_usec) / 1.e6;
    format_seconds(solver_time, formatted_time);
    // Create the solver time using format hh:mm:ss s
    sprintf(solver_time_message, "\nSolver time: %s", formatted_time);
    log_message(parameters, solver_time_message);
    // Free memory
    free(y_np1);
    free(F);
    free(R_turbulence);
    free(k);
    // Free memory for Poisson problem
    fftw_free(a);
    fftw_free(b);
    fftw_free(c);
    fftw_free(d);
    fftw_free(l);
    fftw_free(u);
    fftw_free(y);
    fftw_free(pk);
    fftw_free(f_in);
    fftw_free(f_out);
    fftw_free(p_top_in);
    fftw_free(p_top_out);
    fftw_free(p_in);
    fftw_free(p_out);
}