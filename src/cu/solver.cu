/**
 * @file solver.cu
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Functions for solving the partial differential equations of the wildfire simulation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "../../include/cu/solver.cuh"

__global__ 
void euler_step(double dt, double *y_n, double *y_np1, double *F, int size) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;    
    for (int i = idx; i < size; i += stride) {
        y_np1[i] = y_n[i] + dt * F[i];
    }
}

__global__
void RK2_step(double dt, double *y_n, double *y_np1, double *k1, double *k2, int size) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;
    for (int i = idx; i < size; i += stride) {
        y_np1[i] = y_n[i] + 0.5 * dt * (k1[i] + k2[i]);
    }
}

__global__
void RK4_step(double dt, double *y_n, double *y_np1, double *k1, double *k2, double *k3, double *k4, int size) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;
    for (int i = idx; i < size; i += stride) {
        y_np1[i] = y_n[i] +  (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
}

void create_y_0(double *u, double *v, double *w, double *T, double *Y, double *y_0, Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int Nz_Y_max = parameters.Nz_Y_max;
    int size = Nx * Ny * Nz;
    int size_Y = Nx * Ny * Nz_Y_max;
    for (int i = 0; i < size; i++) {
        y_0[parameters.field_indexes.u + i] = u[i];
        y_0[parameters.field_indexes.v + i] = v[i];
        y_0[parameters.field_indexes.w + i] = w[i];
        y_0[parameters.field_indexes.T + i] = T[i];
        if (i < size_Y) {
            y_0[parameters.field_indexes.Y + i] = Y[i];
        }
    }
}

void euler(double t_n, double *y_n, double *y_np1, double *F, double *U_turbulence, double *z, double *z_ibm, int *Nz_Y, int *cut_nodes, double dt, int size, Parameters parameters) {   
    Phi(t_n, y_n, F, U_turbulence, z, z_ibm, Nz_Y, cut_nodes, parameters);
    euler_step<<<BLOCKS, THREADS>>>(dt, y_n, y_np1, F, size);
    CHECK(cudaDeviceSynchronize());
}

void RK2(double t_n, double *y_n, double *y_np1, double *k, double *F, double *U_turbulence, double *z, double *z_ibm, int *Nz_Y, int *cut_nodes, double dt, int size, Parameters parameters) {
    int k1_index = 0;
    int k2_index = size;
    Phi(t_n, y_n, k + k1_index, U_turbulence, z, z_ibm, Nz_Y, cut_nodes, parameters);
    caxpy<<<BLOCKS, THREADS>>>(F, k + k1_index, y_n, dt, size);
    CHECK(cudaDeviceSynchronize());
    Phi(t_n + dt, F, k + k2_index, U_turbulence, z, z_ibm, Nz_Y, cut_nodes, parameters);
    RK2_step<<<BLOCKS, THREADS>>>(dt, y_n, y_np1, k + k1_index, k + k2_index, size);
    CHECK(cudaDeviceSynchronize());
}

void RK4(double t_n, double *y_n, double *y_np1, double *k, double *F, double *U_turbulence, double *z, double *z_ibm, int *Nz_Y, int *cut_nodes, double dt, int size, Parameters parameters) {
    int k1_index = 0;
    int k2_index = size;
    int k3_index = 2 * size;
    int k4_index = 3 * size;
    Phi(t_n, y_n, k + k1_index, U_turbulence, z, z_ibm, Nz_Y, cut_nodes, parameters);
    caxpy<<<BLOCKS, THREADS>>>(F, k + k1_index, y_n, dt * 0.5, size);
    CHECK(cudaDeviceSynchronize());
    Phi(t_n + 0.5 * dt, F, k + k2_index, U_turbulence, z, z_ibm, Nz_Y, cut_nodes, parameters);
    caxpy<<<BLOCKS, THREADS>>>(F, k + k2_index, y_n, dt * 0.5, size);
    CHECK(cudaDeviceSynchronize());
    Phi(t_n + 0.5 * dt, F, k + k3_index, U_turbulence, z, z_ibm, Nz_Y, cut_nodes, parameters);
    caxpy<<<BLOCKS, THREADS>>>(F, k + k3_index, y_n, dt, size);
    CHECK(cudaDeviceSynchronize());
    Phi(t_n + dt, F, k + k4_index, U_turbulence, z, z_ibm, Nz_Y, cut_nodes, parameters);
    RK4_step<<<BLOCKS, THREADS>>>(dt, y_n, y_np1, k + k1_index, k + k2_index, k + k3_index, k + k4_index, size);
    CHECK(cudaDeviceSynchronize());
}

void solve_PDE(double *y_n, double *p, Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int Nt = parameters.Nt;
    int NT = parameters.NT;
    int Nz_Y_max = parameters.Nz_Y_max;
    int size = 4 * Nx * Ny * Nz + Nx * Ny * Nz_Y_max;
    int n_save;
    int time_step_log = parameters.time_step_log;
    int k_size = (strncmp(parameters.method, "RK4", 3) == 0) ? 4 : 2;
    double step_time, solver_time, cum_step_time = 0.0;
    double CFL = 0.0, T_min = 1e9, T_max = -1e9, Y_min = 1e9, Y_max = -1e9;
    double error;
    int max_iter;
    double *t = parameters.t;
    double dt = parameters.dt;
    double t_source_start = parameters.t_source_start;
    double t_source_end = parameters.t_source_end;
    // Host data
    double *y_np1_host, *p_host;
    // Device data
    double *d_x, *d_y, *d_z, *y_np1, *F, *k, *R_turbulence, *z_ibm, *kx, *ky, *gamma, *d_T_source;
    int *Nz_Y, *cut_nodes;
    // Arrays for pressure problem
    double *a, *b, *c;
    cufftDoubleComplex *d, *l, *u, *y;
    cufftDoubleComplex *data_in, *data_out, *p_top_in, *p_top_out;
    // clock_t start, end, step_start, step_end; // Timers
    struct timeval start_solver, end_solver, start_ts, end_ts; // Timers
    // Messages
    char solver_time_message[128];
    char pressure_log_message[128];
    char formatted_time[64];
    // Host memory allocation
    y_np1_host = (double *) malloc(size * sizeof(double));
    p_host = (double *) malloc(Nx * Ny * Nz * sizeof(double));
    // Memory allocation for device data
    CHECK(cudaMalloc((void **)&d_x, Nx * sizeof(double)));
    CHECK(cudaMalloc((void **)&d_y, Ny * sizeof(double)));
    CHECK(cudaMalloc((void **)&d_z, Nz * sizeof(double)));
    CHECK(cudaMalloc((void **)&y_np1, size * sizeof(double)));
    CHECK(cudaMalloc((void **)&F, size * sizeof(double)));
    CHECK(cudaMalloc((void **)&k, k_size * size * sizeof(double)));
    // CHECK(cudaMalloc((void **)&R_turbulence, 25 * Nx * Ny * Nz * sizeof(double)));
    CHECK(cudaMalloc((void **)&R_turbulence, (18 * Nx * Ny * Nz + Nx * Ny) * sizeof(double)));
    CHECK(cudaMalloc((void **)&z_ibm, Nx * Ny * Nz * sizeof(double)));
    CHECK(cudaMalloc((void **)&Nz_Y, Nx * Ny * sizeof(int)));
    CHECK(cudaMalloc((void **)&cut_nodes, Nx * Ny * sizeof(int)));
    CHECK(cudaMalloc((void **)&kx, (Nx - 1) * sizeof(double)));
    CHECK(cudaMalloc((void **)&ky, (Ny - 1) * sizeof(double)));
    CHECK(cudaMalloc((void **)&gamma, (Nx - 1) * (Ny - 1) * (Nz - 1) * sizeof(double)));
    // Allocate memory for pressure problem
    CHECK(cudaMalloc((void **)&a, (Nx - 1) * (Ny - 1) * (Nz - 2) * sizeof(double)));
    CHECK(cudaMalloc((void **)&b, (Nx - 1) * (Ny - 1) * (Nz - 1) * sizeof(double)));
    CHECK(cudaMalloc((void **)&c, (Nx - 1) * (Ny - 1) * (Nz - 2) * sizeof(double)));
    CHECK(cudaMalloc((void **)&d, (Nx - 1) * (Ny - 1) * (Nz - 1) * sizeof(cufftDoubleComplex)));
    CHECK(cudaMalloc((void **)&l, (Nx - 1) * (Ny - 1) * (Nz - 2) * sizeof(cufftDoubleComplex)));
    CHECK(cudaMalloc((void **)&u, (Nx - 1) * (Ny - 1) * (Nz - 1) * sizeof(cufftDoubleComplex)));
    CHECK(cudaMalloc((void **)&y, (Nx - 1) * (Ny - 1) * (Nz - 1) * sizeof(cufftDoubleComplex)));
    CHECK(cudaMalloc((void **)&p_top_in, (Nx - 1) * (Ny - 1) * sizeof(cufftDoubleComplex)));
    CHECK(cudaMalloc((void **)&p_top_out, (Nx - 1) * (Ny - 1) * sizeof(cufftDoubleComplex)));
    CHECK(cudaMalloc((void **)&data_in, (Nx - 1) * (Ny - 1) * (Nz - 1) * sizeof(cufftDoubleComplex)));
    CHECK(cudaMalloc((void **)&data_out, (Nx - 1) * (Ny - 1) * (Nz - 1) * sizeof(cufftDoubleComplex)));
    CHECK(cudaMalloc((void **)&d_T_source, Nx * Ny * Nz * sizeof(double)));
    // Copy host data to device
    CHECK(cudaMemcpy(d_x, parameters.x, Nx * sizeof(double), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_y, parameters.y, Ny * sizeof(double), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_z, parameters.z, Nz * sizeof(double), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(z_ibm, parameters.z_ibm, Nx * Ny * Nz * sizeof(double), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(kx, parameters.kx, (Nx - 1) * sizeof(double), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(ky, parameters.ky, (Ny - 1) * sizeof(double), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(Nz_Y, parameters.Nz_Y, Nx * Ny * sizeof(int), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(cut_nodes, parameters.cut_nodes, Nx * Ny * sizeof(int), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_T_source, y_n + parameters.field_indexes.T, Nx * Ny * Nz * sizeof(double), cudaMemcpyDeviceToDevice));
    // Fill gamma and coefficients
    gammas_and_coefficients<<<BLOCKS, THREADS>>>(kx, ky, gamma, a, b, c, d_z, parameters);
    checkCuda(cudaGetLastError());
    CHECK(cudaDeviceSynchronize());
    // Solver time
    // start = clock();
    gettimeofday(&start_solver, NULL);
    // Time integration
    for (int n = 1; n <= Nt; n++) { 
        // step_start = clock(); // Start step timer
        gettimeofday(&start_ts, NULL);
        // Compute U^*, T^{n+1}, Y^{n+1}
        // Check time integration method
        if (strncmp(parameters.method, "Euler", 5) == 0) {
            euler(t[n - 1], y_n, y_np1, F, R_turbulence, d_z, z_ibm, Nz_Y, cut_nodes, dt, size, parameters);
        } else if (strncmp(parameters.method, "RK2", 3) == 0) {
            RK2(t[n - 1], y_n, y_np1, k, F, R_turbulence, d_z, z_ibm, Nz_Y, cut_nodes, dt, size, parameters);
        } else if (strncmp(parameters.method, "RK4", 3) == 0) {
            RK4(t[n - 1], y_n, y_np1, k, F, R_turbulence, d_z, z_ibm, Nz_Y, cut_nodes, dt, size, parameters);
        } else {
            log_message(parameters, "Time integration method not found.");
            exit(1);
        }  
        CHECK(cudaDeviceSynchronize());
        // Pressure solver
        if (parameters.variable_density == 0) { // Constant density, direct solver
            // solve_pressure(y_np1, p, d_z, gamma, a, b, c, d, l, u, y, data_in, data_out, p_top_in, p_top_out, parameters);  
            solve_pressure(y_np1, y_n, p, d_z, gamma, a, b, c, d, l, u, y, data_in, data_out, p_top_in, p_top_out, Nz_Y, parameters);            
        } else { // Variable density, iterative solver
            // solve_pressure_iterative(y_np1, p, d_z, gamma, a, b, c, d, l, u, y, data_in, data_out, p_top_in, p_top_out, parameters, &error, &max_iter);
            solve_pressure_iterative(y_np1, y_n, p, d_z, gamma, a, b, c, d, l, u, y, data_in, data_out, p_top_in, p_top_out, Nz_Y, parameters, &error, &max_iter);
        }
        checkCuda(cudaGetLastError());
        CHECK(cudaDeviceSynchronize());
        // Chorin's projection method        
        velocity_correction<<<BLOCKS, THREADS>>>(y_np1, y_n, p, d_z, parameters);
        checkCuda(cudaGetLastError());
        CHECK(cudaDeviceSynchronize());
        // Boundary conditions
        boundary_conditions<<<BLOCKS, THREADS>>>(y_np1, d_z, Nz_Y, cut_nodes, parameters);
        checkCuda(cudaGetLastError());
        CHECK(cudaDeviceSynchronize());
        // Bounds
        bounds<<<BLOCKS, THREADS>>>(y_np1, parameters);
        checkCuda(cudaGetLastError());
        CHECK(cudaDeviceSynchronize());
        // Temperature source
        if (t[n] <= t_source_end) {
            if (t_source_start <= 0) {
                temperature_source<<<BLOCKS, THREADS>>>(d_x, d_y, d_z, y_np1, d_T_source, parameters);
            } else {
                temperature_source_delay<<<BLOCKS, THREADS>>>(d_x, d_y, d_z, y_np1, t[n], parameters);
            }
        }
        // if (t[n] <= parameters.t_source_end) {
        //     // temperature_source<<<BLOCKS, THREADS>>>(d_x, d_y, d_z, y_np1, parameters);
        //     temperature_source<<<BLOCKS, THREADS>>>(d_x, d_y, d_z, y_np1, d_T_source, parameters);
        //     checkCuda(cudaGetLastError());
        //     CHECK(cudaDeviceSynchronize());
        // }
        // temperature_source_delay<<<BLOCKS, THREADS>>>(d_x, d_y, d_z, y_np1, t[n], parameters);
        checkCuda(cudaGetLastError());
        CHECK(cudaDeviceSynchronize());
        // End step timer
        // step_end = clock(); 
        gettimeofday(&end_ts, NULL);
        // Compute step time
        // step_time = (double) (step_end - step_start) / CLOCKS_PER_SEC;
        step_time = ((end_ts.tv_sec  - start_ts.tv_sec) * 1000000u + end_ts.tv_usec - start_ts.tv_usec) / 1.e6;        
        // Show time step and average time each 100 steps
        cum_step_time += step_time;
        if (n % time_step_log == 0) {
            printf("Time step: %d, Average time: %lf s\n", n, cum_step_time / n);
        }
        // Save data each NT steps and at the last step
        if (n % NT == 0 || n == Nt) {  
            n_save = n / NT;
            // Copy y_np1 and p to host
            cudaMemcpy(y_np1_host, y_np1, size * sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(p_host, p, Nx * Ny * Nz * sizeof(double), cudaMemcpyDeviceToHost);
            timestep_reports(y_np1_host, &CFL, &Y_min, &Y_max, &T_min, &T_max, n, parameters);
            log_timestep(parameters, n, t[n], step_time, CFL, T_min, T_max, Y_min, Y_max);
            if (parameters.variable_density == 1) {
                sprintf(pressure_log_message, "Pressure solver: Error = %e, iterations = %d", error, max_iter);
                log_message(parameters, pressure_log_message);
            }            
            save_data(y_np1_host, p_host, n_save, parameters);
            printf("\n");
            // Wait for the device to finish
            CHECK(cudaDeviceSynchronize());
        }
        // Update y_n
        copy<<<BLOCKS, THREADS>>>(y_n, y_np1, size);
        checkCuda(cudaGetLastError());
        CHECK(cudaDeviceSynchronize());
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
    // Host memory
    free(y_np1_host);
    free(p_host);
    // Device memory
    cudaFree(d_x);
    cudaFree(d_y);
    cudaFree(d_z);
    cudaFree(y_np1);
    cudaFree(F);
    cudaFree(R_turbulence);
    cudaFree(k);
    // Free memory for pressure problem
    cudaFree(a);
    cudaFree(b);
    cudaFree(c);
    cudaFree(d);
    cudaFree(l);
    cudaFree(u);
    cudaFree(y);
    cudaFree(p_top_in);
    cudaFree(p_top_out);
    cudaFree(data_in);
    cudaFree(data_out);
}
