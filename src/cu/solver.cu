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
     // printf("Euler step\n");
     int idx = threadIdx.x + blockIdx.x * blockDim.x;
     int stride = gridDim.x * blockDim.x;    
     for (int i = idx; i < size; i += stride) {
         y_np1[i] = y_n[i] + dt * F[i];
     }
 }
 
 void euler(double t_n, double *y_n, double *y_np1, double *F, double *U_turbulence, double *z, double *z_ibm, int *Nz_Y, int *cut_nodes, double dt, int size, Parameters parameters) {   
     Phi(t_n, y_n, F, U_turbulence, z, z_ibm, Nz_Y, cut_nodes, parameters);
     euler_step<<<BLOCKS, THREADS>>>(dt, y_n, y_np1, F, size);
     cudaDeviceSynchronize();
 }
 
 __global__
 void RK2_step(double dt, double *y_n, double *y_np1, double *k1, double *k2, int size) {
     int idx = threadIdx.x + blockIdx.x * blockDim.x;
     int stride = gridDim.x * blockDim.x;
     for (int i = idx; i < size; i += stride) {
         y_np1[i] = y_n[i] + 0.5 * dt * (k1[i] + k2[i]);
     }
 }
 
 void RK2(double t_n, double *y_n, double *y_np1, double *k, double *F, double *U_turbulence, double *z, double *z_ibm, int *Nz_Y, int *cut_nodes, double dt, int size, Parameters parameters) {
     int k1_index = 0;
     int k2_index = size;
     Phi(t_n, y_n, k + k1_index, U_turbulence, z, z_ibm, Nz_Y, cut_nodes, parameters);
     caxpy<<<BLOCKS, THREADS>>>(F, k + k1_index, y_n, dt, size);
     cudaDeviceSynchronize();
     Phi(t_n + dt, F, k + k2_index, U_turbulence, z, z_ibm, Nz_Y, cut_nodes, parameters);
     RK2_step<<<BLOCKS, THREADS>>>(dt, y_n, y_np1, k + k1_index, k + k2_index, size);
     cudaDeviceSynchronize();
 }
 
 __global__
 void RK4_step(double dt, double *y_n, double *y_np1, double *k1, double *k2, double *k3, double *k4, int size) {
     int idx = threadIdx.x + blockIdx.x * blockDim.x;
     int stride = gridDim.x * blockDim.x;
     for (int i = idx; i < size; i += stride) {
         y_np1[i] = y_n[i] +  (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
     }
 }
 
 void RK4(double t_n, double *y_n, double *y_np1, double *k, double *F, double *U_turbulence, double *z, double *z_ibm, int *Nz_Y, int *cut_nodes, double dt, int size, Parameters parameters) {
     int k1_index = 0;
     int k2_index = size;
     int k3_index = 2 * size;
     int k4_index = 3 * size;
     Phi(t_n, y_n, k + k1_index, U_turbulence, z, z_ibm, Nz_Y, cut_nodes, parameters);
     caxpy<<<BLOCKS, THREADS>>>(F, k + k1_index, y_n, dt * 0.5, size);
     cudaDeviceSynchronize();
     Phi(t_n + 0.5 * dt, F, k + k2_index, U_turbulence, z, z_ibm, Nz_Y, cut_nodes, parameters);
     caxpy<<<BLOCKS, THREADS>>>(F, k + k2_index, y_n, dt * 0.5, size);
     cudaDeviceSynchronize();
     Phi(t_n + 0.5 * dt, F, k + k3_index, U_turbulence, z, z_ibm, Nz_Y, cut_nodes, parameters);
     caxpy<<<BLOCKS, THREADS>>>(F, k + k3_index, y_n, dt, size);
     cudaDeviceSynchronize();
     Phi(t_n + dt, F, k + k4_index, U_turbulence, z, z_ibm, Nz_Y, cut_nodes, parameters);
     RK4_step<<<BLOCKS, THREADS>>>(dt, y_n, y_np1, k + k1_index, k + k2_index, k + k3_index, k + k4_index, size);
     cudaDeviceSynchronize();
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
 
 void solve_PDE(double *y_n, double *p, Parameters parameters) {
     int Nx = parameters.Nx;
     int Ny = parameters.Ny;
     int Nz = parameters.Nz;
     int Nt = parameters.Nt;
     int NT = parameters.NT;
     int Nz_Y_max = parameters.Nz_Y_max;
     int size = 4 * Nx * Ny * Nz + Nx * Ny * Nz_Y_max;
     int n_save;
     int k_size = (strncmp(parameters.method, "RK4", 3) == 0) ? 4 : 2;
     double step_time, solver_time;
     double CFL = 0.0, T_min = 1e9, T_max = -1e9, Y_min = 1e9, Y_max = -1e9;
     double *t = parameters.t;
     double dt = parameters.dt;
     // Host data
     double *y_np1_host, *p_host;
     // Device data
     double *d_x, *d_y, *d_z, *y_np1, *F, *k, *R_turbulence, *z_ibm, *kx, *ky, *gamma;
     int *Nz_Y, *cut_nodes;
     // Arrays for pressure Poisson Problema
     double *a, *b, *c;
     cufftDoubleComplex *d, *l, *u, *y;
     cufftDoubleComplex *data_in, *data_out, *p_top_in, *p_top_out;
     // clock_t start, end, step_start, step_end; // Timers
     struct timeval start_solver, end_solver, start_ts, end_ts; // Timers
     // Messages
     char solver_time_message[128];
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
     CHECK(cudaMalloc((void **)&R_turbulence, 25 * Nx * Ny * Nz * sizeof(double)));
     CHECK(cudaMalloc((void **)&z_ibm, Nx * Ny * Nz * sizeof(double)));
     CHECK(cudaMalloc((void **)&Nz_Y, Nx * Ny * sizeof(int)));
     CHECK(cudaMalloc((void **)&cut_nodes, Nx * Ny * sizeof(int)));
     CHECK(cudaMalloc((void **)&kx, (Nx - 1) * sizeof(double)));
     CHECK(cudaMalloc((void **)&ky, (Ny - 1) * sizeof(double)));
     CHECK(cudaMalloc((void **)&gamma, (Nx - 1) * (Ny - 1) * (Nz - 1) * sizeof(double)));
     // Allocate memory for Poisson problem
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
     // Copy host data to device
     CHECK(cudaMemcpy(d_x, parameters.x, Nx * sizeof(double), cudaMemcpyHostToDevice));
     CHECK(cudaMemcpy(d_y, parameters.y, Ny * sizeof(double), cudaMemcpyHostToDevice));
     CHECK(cudaMemcpy(d_z, parameters.z, Nz * sizeof(double), cudaMemcpyHostToDevice));
     CHECK(cudaMemcpy(z_ibm, parameters.z_ibm, Nx * Ny * Nz * sizeof(double), cudaMemcpyHostToDevice));
     CHECK(cudaMemcpy(kx, parameters.kx, (Nx - 1) * sizeof(double), cudaMemcpyHostToDevice));
     CHECK(cudaMemcpy(ky, parameters.ky, (Ny - 1) * sizeof(double), cudaMemcpyHostToDevice));
     CHECK(cudaMemcpy(Nz_Y, parameters.Nz_Y, Nx * Ny * sizeof(int), cudaMemcpyHostToDevice));
     CHECK(cudaMemcpy(cut_nodes, parameters.cut_nodes, Nx * Ny * sizeof(int), cudaMemcpyHostToDevice));
     // Fill gamma and coefficients
     gammas_and_coefficients<<<BLOCKS, THREADS>>>(kx, ky, gamma, a, b, c, parameters);
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
             euler(t[n], y_n, y_np1, F, R_turbulence, d_z, z_ibm, Nz_Y, cut_nodes, dt, size, parameters);
         } else if (strncmp(parameters.method, "RK2", 3) == 0) {
             RK2(t[n], y_n, y_np1, k, F, R_turbulence, d_z, z_ibm, Nz_Y, cut_nodes, dt, size, parameters);
         } else if (strncmp(parameters.method, "RK4", 3) == 0) {
             RK4(t[n], y_n, y_np1, k, F, R_turbulence, d_z, z_ibm, Nz_Y, cut_nodes, dt, size, parameters);
         } else {
             log_message(parameters, "Time integration method not found.");
             exit(1);
         }
         // Solve Poisson problem for pressure (it only uses U^*)
         solve_pressure(y_np1, p, gamma, a, b, c, d, l, u, y, data_in, data_out, p_top_in, p_top_out, parameters);
         checkCuda(cudaGetLastError());
         // Chorin's projection method
         velocity_correction_fw<<<BLOCKS, THREADS>>>(y_np1, p, dt, parameters);
         checkCuda(cudaGetLastError());
         // Boundary conditions
         boundary_conditions<<<BLOCKS, THREADS>>>(y_np1, Nz_Y, cut_nodes, parameters);
         checkCuda(cudaGetLastError());
         // Bounds
         bounds<<<BLOCKS, THREADS>>>(y_np1, parameters);
         checkCuda(cudaGetLastError());
         // Add source when t_n <= t_source
         if (t[n] <= parameters.t_source) {
             temperature_source<<<BLOCKS, THREADS>>>(d_x, d_y, d_z, y_np1, parameters);
             checkCuda(cudaGetLastError());
         }
         // End step timer
         // step_end = clock(); 
         gettimeofday(&end_ts, NULL);
         // Compute step time
         // step_time = (double) (step_end - step_start) / CLOCKS_PER_SEC;
         step_time = ((end_ts.tv_sec  - start_ts.tv_sec) * 1000000u + end_ts.tv_usec - start_ts.tv_usec) / 1.e6;
         // Save data each NT steps and at the last step
         if (n % NT == 0 || n == Nt) {  
             n_save = n / NT;
             // Copy y_np1 and p to host
             cudaMemcpy(y_np1_host, y_np1, size * sizeof(double), cudaMemcpyDeviceToHost);
             cudaMemcpy(p_host, p, Nx * Ny * Nz * sizeof(double), cudaMemcpyDeviceToHost);
             timestep_reports(y_np1_host, &CFL, &Y_min, &Y_max, &T_min, &T_max, parameters);
             log_timestep(parameters, n, t[n], step_time, CFL, T_min, T_max, Y_min, Y_max);
             save_data(y_np1_host, p_host, n_save, parameters);
         }
         // Update y_n
         copy<<<BLOCKS, THREADS>>>(y_n, y_np1, size);
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
     // Free memory for Poisson problem
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
 
 
 /*
 void solve_PDE_v1(double *y_n, double *p, Parameters parameters) {
     int Nx = parameters.Nx;
     int Ny = parameters.Ny;
     int Nz = parameters.Nz;
     int Nt = parameters.Nt;
     int NT = parameters.NT;
     int Nz_Y = parameters.Nz_Y;
     int size = 4 * Nx * Ny * Nz + Nx * Ny * Nz_Y;
     // char log_path[100];
     int n_save;
     int k_size = 2;
     if (strncmp(parameters.method, "RK4", 3) == 0) 
         k_size = 4;
     double step_time, solver_time;
     double CFL = 0.0, T_min = 1e9, T_max = -1e9, Y_min = 1e9, Y_max = -1e9;
     // double *CFL, *T_min, *T_max, *Y_min, *Y_max;
     double *t = parameters.t;
     double dt = parameters.dt;
     double *y_np1_host = (double *) malloc(size * sizeof(double));
     double *p_host = (double *) malloc(Nx * Ny * Nz * sizeof(double));
     double *kx = parameters.kx;
     double *ky = parameters.ky;
     double *y_np1, *F, *k, *R_turbulence, *z;//, *kx, *ky;
     // double *F = (double *) malloc(size * sizeof(double));
     // double *k = (double *) malloc(k_size * size * sizeof(double));
     // double *R_turbulence = (double *) malloc(27 * Nx * Ny * Nz * sizeof(double));
     // CudaMalloc
     cudaMalloc(&y_np1, size * sizeof(double));
     cudaMalloc(&F, size * sizeof(double));
     cudaMalloc(&k, k_size * size * sizeof(double));
     cudaMalloc(&R_turbulence, 25 * Nx * Ny * Nz * sizeof(double));
     cudaMalloc(&z, Nz * sizeof(double));
     // cudaMalloc(&kx, (Nx - 1) * sizeof(double));
     // cudaMalloc(&ky, (Ny - 1) * sizeof(double));
     cudaMemcpy(z, parameters.z, Nz * sizeof(double), cudaMemcpyHostToDevice);
     // cudaMemcpy(kx, parameters.kx, (Nx - 1) * sizeof(double), cudaMemcpyHostToDevice);
     // cudaMemcpy(ky, parameters.ky, (Ny - 1) * sizeof(double), cudaMemcpyHostToDevice);
     // clock_t start, end, step_start, step_end; // Timers
     struct timeval start_solver, end_solver, start_ts, end_ts;
     char solver_time_message[BLOCKS];
     // char min_max_values_message[BLOCKS];
     char formatted_time[64];
     // Arrays for pressure Poisson Problem
     // Allocate memory for cuFFT arrays
     cufftDoubleComplex *a, *b, *c, *d, *l, *u, *y, *pk;
     cufftDoubleComplex *f_in, *f_out, *p_top_in, *p_top_out, *p_in, *p_out;
     cufftHandle p_plan = 0, f_plan = 0, p_top_plan = 0;
 
     CHECK(cudaMalloc((void **)&a, (Nz - 2) * sizeof(cufftDoubleComplex)));
     CHECK(cudaMalloc((void **)&b, (Nz - 1) * sizeof(cufftDoubleComplex)));
     CHECK(cudaMalloc((void **)&c, (Nz - 2) * sizeof(cufftDoubleComplex)));
     CHECK(cudaMalloc((void **)&d, (Nz - 1) * sizeof(cufftDoubleComplex)));
     CHECK(cudaMalloc((void **)&l, (Nz - 2) * sizeof(cufftDoubleComplex)));
     CHECK(cudaMalloc((void **)&u, (Nz - 1) * sizeof(cufftDoubleComplex)));
     CHECK(cudaMalloc((void **)&y, (Nz - 1) * sizeof(cufftDoubleComplex)));
     CHECK(cudaMalloc((void **)&pk, (Nz - 1) * sizeof(cufftDoubleComplex)));
     CHECK(cudaMalloc((void **)&f_in, (Nx - 1) * (Ny - 1) * (Nz - 1) * sizeof(cufftDoubleComplex)));
     CHECK(cudaMalloc((void **)&f_out, (Nx - 1) * (Ny - 1) * (Nz - 1) * sizeof(cufftDoubleComplex)));
     CHECK(cudaMalloc((void **)&p_top_in, (Nx - 1) * (Ny - 1) * sizeof(cufftDoubleComplex)));
     CHECK(cudaMalloc((void **)&p_top_out, (Nx - 1) * (Ny - 1) * sizeof(cufftDoubleComplex)));
     CHECK(cudaMalloc((void **)&p_in, (Nx - 1) * (Ny - 1) * (Nz - 1) * sizeof(cufftDoubleComplex)));
     CHECK(cudaMalloc((void **)&p_out, (Nx - 1) * (Ny - 1) * (Nz - 1) * sizeof(cufftDoubleComplex)));
 
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
             euler(t[n], y_n, y_np1, F, R_turbulence, dt, size, parameters, z);
         } else if (strncmp(parameters.method, "RK2", 3) == 0) {
             RK2(t[n], y_n, y_np1, k, F, R_turbulence, dt, size, parameters, z);
         } else if (strncmp(parameters.method, "RK4", 3) == 0) {
             RK4(t[n], y_n, y_np1, k, F, R_turbulence, dt, size, parameters, z);
         } else {
             log_message(parameters, "Time integration method not found.");
             exit(1);
         }
         // Solve Poisson problem for pressure (it only uses U^*)
         solve_pressure(y_np1, p, kx, ky, a, b, c, d, l, u, y, pk, p_plan, f_plan, p_top_plan, f_in, f_out, p_top_in, p_top_out, p_in, p_out, parameters);
         checkCuda(cudaGetLastError());
         // Chorin's projection method
         velocity_correction_fw<<<BLOCKS, THREADS>>>(y_np1, p, dt, parameters);
         checkCuda(cudaGetLastError());
         // Boundary conditions
         boundary_conditions<<<BLOCKS, THREADS>>>(y_np1, parameters);
         checkCuda(cudaGetLastError());
         // Bounds
         bounds<<<BLOCKS, THREADS>>>(y_np1, parameters);
         checkCuda(cudaGetLastError());
         // End step timer
         // step_end = clock(); 
         gettimeofday(&end_ts, NULL);
         // Compute step time
         // step_time = (double) (step_end - step_start) / CLOCKS_PER_SEC;
         step_time = ((end_ts.tv_sec  - start_ts.tv_sec) * 1000000u + end_ts.tv_usec - start_ts.tv_usec) / 1.e6;
         // Save data each NT steps and at the last step
         if (n % NT == 0 || n == Nt) {  
             n_save = n / NT;
             // Copy y_np1 and p to host
             cudaMemcpy(y_np1_host, y_np1, size * sizeof(double), cudaMemcpyDeviceToHost);
             cudaMemcpy(p_host, p, Nx * Ny * Nz * sizeof(double), cudaMemcpyDeviceToHost);
             timestep_reports(y_np1_host, &CFL, &Y_min, &Y_max, &T_min, &T_max, parameters);
             log_timestep(parameters, n, t[n], step_time, CFL, T_min, T_max, Y_min, Y_max);
             save_data(y_np1_host, p_host, n_save, parameters);
         }
         // Update y_n
         copy<<<BLOCKS, THREADS>>>(y_n, y_np1, size);
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
     cudaFree(y_np1);
     cudaFree(F);
     cudaFree(R_turbulence);
     cudaFree(k);
     // Free memory for Poisson problem
     cudaFree(a);
     cudaFree(b);
     cudaFree(c);
     cudaFree(d);
     cudaFree(l);
     cudaFree(u);
     cudaFree(y);
     cudaFree(pk);
     cudaFree(f_in);
     cudaFree(f_out);
     cudaFree(p_top_in);
     cudaFree(p_top_out);
     cudaFree(p_in);
     cudaFree(p_out);
 }
 */
 
 /*
 void solve_PDE(double *y_n, double *p, Parameters parameters) {
     int Nx = parameters.Nx;
     int Ny = parameters.Ny;
     int Nz = parameters.Nz;
     int Nt = parameters.Nt;
     int NT = parameters.NT;
     int Nz_Y = parameters.Nz_Y;
     int size = 4 * Nx * Ny * Nz + Nx * Ny * Nz_Y;
     // char log_path[100];
     int n_save;
     int k_size = 2;
     if (strncmp(parameters.method, "RK4", 3) == 0) 
         k_size = 4;
     double step_time, solver_time;
     double CFL = 0.0, T_min = 1e9, T_max = -1e9, Y_min = 1e9, Y_max = -1e9;
     // double *CFL, *T_min, *T_max, *Y_min, *Y_max;
     double *t = parameters.t;
     double dt = parameters.dt;
     double *y_np1_host = (double *) malloc(size * sizeof(double));
     double *p_host = (double *) malloc(Nx * Ny * Nz * sizeof(double));
     double *kx = parameters.kx;
     double *ky = parameters.ky;
     double *y_np1, *F, *k, *R_turbulence, *z;//, *kx, *ky;
     // double *F = (double *) malloc(size * sizeof(double));
     // double *k = (double *) malloc(k_size * size * sizeof(double));
     // double *R_turbulence = (double *) malloc(27 * Nx * Ny * Nz * sizeof(double));
     // CudaMalloc
     cudaMalloc(&y_np1, size * sizeof(double));
     cudaMalloc(&F, size * sizeof(double));
     cudaMalloc(&k, k_size * size * sizeof(double));
     cudaMalloc(&R_turbulence, 25 * Nx * Ny * Nz * sizeof(double));
     cudaMalloc(&z, Nz * sizeof(double));
     // cudaMalloc(&kx, (Nx - 1) * sizeof(double));
     // cudaMalloc(&ky, (Ny - 1) * sizeof(double));
     cudaMemcpy(z, parameters.z, Nz * sizeof(double), cudaMemcpyHostToDevice);
     // cudaMemcpy(kx, parameters.kx, (Nx - 1) * sizeof(double), cudaMemcpyHostToDevice);
     // cudaMemcpy(ky, parameters.ky, (Ny - 1) * sizeof(double), cudaMemcpyHostToDevice);
     // clock_t start, end, step_start, step_end; // Timers
     struct timeval start_solver, end_solver, start_ts, end_ts;
     char solver_time_message[BLOCKS];
     // char min_max_values_message[BLOCKS];
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
         if (strncmp(parameters.method, "Euler", 5) == 0) {
             euler(t[n], y_n, y_np1, F, R_turbulence, dt, size, parameters, z);
         } else if (strncmp(parameters.method, "RK2", 3) == 0) {
             RK2(t[n], y_n, y_np1, k, F, R_turbulence, dt, size, parameters, z);
         } else if (strncmp(parameters.method, "RK4", 3) == 0) {
             RK4(t[n], y_n, y_np1, k, F, R_turbulence, dt, size, parameters, z);
         } else {
             log_message(parameters, "Time integration method not found.");
             exit(1);
         }
         // Solve Poisson problem for pressure (it only uses U^*)
         // solve_pressure(y_np1, p, kx, ky, a, b, c, d, l, u, y, pk, p_plan, f_plan, p_top_plan, f_in, f_out, p_top_in, p_top_out, p_in, p_out, parameters);
         // checkCuda(cudaGetLastError());
         cudaMemcpy(y_np1_host, y_np1, size * sizeof(double), cudaMemcpyDeviceToHost);
         cudaMemcpy(p_host, p, Nx * Ny * Nz * sizeof(double), cudaMemcpyDeviceToHost);
         solve_pressure(y_np1_host, p_host, a, b, c, d, l, u, y, pk, p_plan, f_plan, p_top_plan, f_in, f_out, p_top_in, p_top_out, p_in, p_out, &parameters);
         cudaMemcpy(y_np1, y_np1_host, size * sizeof(double), cudaMemcpyHostToDevice);
         cudaMemcpy(p, p_host, Nx * Ny * Nz * sizeof(double), cudaMemcpyHostToDevice);
         // Chorin's projection method
         velocity_correction_fw<<<BLOCKS, THREADS>>>(y_np1, p, dt, parameters);
         checkCuda(cudaGetLastError());
         // Boundary conditions
         boundary_conditions<<<BLOCKS, THREADS>>>(y_np1, parameters);
         checkCuda(cudaGetLastError());
         // Bounds
         bounds<<<BLOCKS, THREADS>>>(y_np1, parameters);
         checkCuda(cudaGetLastError());
         // End step timer
         // step_end = clock(); 
         gettimeofday(&end_ts, NULL);
         // Compute step time
         // step_time = (double) (step_end - step_start) / CLOCKS_PER_SEC;
         step_time = ((end_ts.tv_sec  - start_ts.tv_sec) * 1000000u + end_ts.tv_usec - start_ts.tv_usec) / 1.e6;
         // Save data each NT steps and at the last step
         if (n % NT == 0 || n == Nt) {  
             n_save = n / NT;
             // Copy y_np1 and p to host
             cudaMemcpy(y_np1_host, y_np1, size * sizeof(double), cudaMemcpyDeviceToHost);
             cudaMemcpy(p_host, p, Nx * Ny * Nz * sizeof(double), cudaMemcpyDeviceToHost);
             timestep_reports(y_np1_host, &CFL, &Y_min, &Y_max, &T_min, &T_max, parameters);
             log_timestep(parameters, n, t[n], step_time, CFL, T_min, T_max, Y_min, Y_max);
             save_data(y_np1_host, p_host, n_save, parameters);
         }
         // Update y_n
         copy<<<BLOCKS, THREADS>>>(y_n, y_np1, size);
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
     cudaFree(y_np1);
     cudaFree(F);
     cudaFree(R_turbulence);
     cudaFree(k);
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
 */