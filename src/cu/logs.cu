#include "../../include/cu/logs.cuh"

void log_parameters_report(Parameters parameters, int to_file) {
    FILE *output;
    double dz_min = INFINITY, dz_max = -INFINITY;
    setbuf(stdout, NULL);
    if (to_file) 
        output = parameters.log_files.parameters;
    else
        output = stdout;
    fprintf(output, "Simulation name: %s\n", parameters.sim_id);
    fprintf(output, "Domain:\n");
    fprintf(output, "  [%lf, %lf] x [%lf, %lf] x [%lf, %lf] x [%lf, %lf]\n",
        parameters.x_min, parameters.x_max, parameters.y_min, parameters.y_max,
        parameters.z_min, parameters.z_max, parameters.t_min, parameters.t_max);
    fprintf(output, "Grid size:\n");
    fprintf(output, "  Nx: %d, Ny: %d, Nz: %d, Nt: %d\n", parameters.Nx, parameters.Ny, parameters.Nz, parameters.Nt);
    if (parameters.z_domain == 0) {
        fprintf(output, "  dx: %lf, dy: %lf, dz: %lf, dt: %lf\n", parameters.dx, parameters.dy, parameters.dz, parameters.dt);
    } else {
        for (int k = 0; k < parameters.Nz - 1; k++) {
            dz_min = MIN(dz_min, parameters.z[k + 1] - parameters.z[k]);
            dz_max = MAX(dz_max, parameters.z[k + 1] - parameters.z[k]);
        }
        fprintf(output, "  dx: %lf, dy: %lf, dz: [%lf, %lf], dt: %lf\n", parameters.dx, parameters.dy, dz_min, dz_max, parameters.dt);
        if (parameters.z_domain == 1) {
            fprintf(output, "  Stretched grid in z-direction, k: %lf\n", parameters.k_nu_grid);
        }
    }
    fprintf(output, "  Time samples: %d\n", parameters.NT);
    fprintf(output, "Time integration: %s\n", parameters.method);
    fprintf(output, "Fluid parameters:\n");
    fprintf(output, "  mu: %e, kappa: %e, c_p: %lf\n", parameters.mu, parameters.kappa, parameters.c_p);
    fprintf(output, "  rho_inf: %lf, T_inf: %lf,  g: %lf\n", parameters.rho_inf, parameters.T_inf, parameters.g);
    fprintf(output, "  C_s: %lf, Pr: %lf\n", parameters.C_s, parameters.Pr);
    fprintf(output, "Fuel parameters:\n");
    fprintf(output, "  h_c: %lf, alpha_s: %lf, sigma_s: %lf, delta: %lf\n", parameters.h_c, parameters.alpha_s, parameters.sigma_s, parameters.delta);    
    fprintf(output, "  H_C: %e, A: %e, T_pc: %lf, T_act: %lf\n", parameters.H_C, parameters.A, parameters.T_pc, parameters.T_act);
    fprintf(output, "  C_d: %lf, Y_f: %lf\n", parameters.C_d, parameters.Y_f);
    // Check if input_path is empty
    if (parameters.input_path[0] == '\0') {
        fprintf(output, "Wind initial condition:\n");
        fprintf(output, "  Type: %s\n", parameters.U0_type);    
        if (strcmp(parameters.U0_type, "constant") == 0)
            fprintf(output, "  u_r: %lf\n", parameters.u_r);
        else if (strcmp(parameters.U0_type, "log") == 0)
            fprintf(output, "  u_z0: %lf, kappa: %lf, d: %lf, u_ast: %lf\n", parameters.u_z0, KAPPA, parameters.d, parameters.u_ast);
        else if (strcmp(parameters.U0_type, "power_law") == 0)
            fprintf(output, "  u_r: %lf, z_r: %lf, alpha_u: %lf\n", parameters.u_r, parameters.z_r, parameters.alpha_u);
        fprintf(output, "Temperature initial condition:\n");
        fprintf(output, "  Shape: %s\n", parameters.T0_shape);
        fprintf(output, "  T_source: %lf\n", parameters.T_source);
        fprintf(output, "  x: [%lf, %lf], y: [%lf, %lf], z: [%lf, %lf]\n", 
            parameters.T0_x_start, parameters.T0_x_end, 
            parameters.T0_y_start, parameters.T0_y_end, 
            parameters.T0_z_start, parameters.T0_z_end);
        fprintf(output, "  Center: [%lf, %lf, %lf], Length: %lf, Width: %lf, Height: %lf\n", 
            parameters.T0_x_center, parameters.T0_y_center, parameters.T0_z_center,
            parameters.T0_length, parameters.T0_width, parameters.T0_height);
        if (parameters.t_source > 0)
            fprintf(output, "  Source time: %lf\n", parameters.t_source);  
        // fprintf(output, "  T0_x_start: %lf, T0_x_end: %lf, T0_x_center: %lf, T0_length: %lf\n", 
        //     parameters.T0_x_start, parameters.T0_x_end, parameters.T0_x_center, parameters.T0_length);
        // fprintf(output, "  T0_y_start: %lf, T0_y_end: %lf, T0_y_center: %lf, T0_width: %lf\n", 
        //     parameters.T0_y_start, parameters.T0_y_end, parameters.T0_y_center, parameters.T0_width);
        // fprintf(output, "  T0_z_start: %lf, T0_z_end: %lf, T0_z_center: %lf, T0_height: %lf\n",
        //     parameters.T0_z_start, parameters.T0_z_end, parameters.T0_z_center, parameters.T0_height);
        fprintf(output, "Fuel initial condition:\n");
        fprintf(output, "  Fuel height: %lf\n", parameters.Y_h);
        fprintf(output, "  x: [%lf, %lf], y: [%lf, %lf]\n", 
            parameters.Y0_x_start, parameters.Y0_x_end, 
            parameters.Y0_y_start, parameters.Y0_y_end);
        fprintf(output, "  Fuel smooth x: %lf, y: %lf\n", parameters.Y0_xa, parameters.Y0_ya);
        // fprintf(output, "  Fuel relax: %d\n", parameters.fuel_relax);
    } 
    // else {
    //     fprintf(output, "  Loaded from the input path: %s\n", parameters.input_path);
    // }
    // Topography
    fprintf(output, "Topography:\n");
    fprintf(output, "  Shape: %s\n", parameters.topo_shape);
    if (strcmp(parameters.topo_shape, "hill") == 0) {
        fprintf(output, "  Hill center: [%lf, %lf], Length: %lf, Width: %lf, Height: %lf\n", 
            parameters.hill_center_x, parameters.hill_center_y, parameters.hill_length, parameters.hill_width, parameters.hill_height);
    }
    fprintf(output, "Bounds:\n");
    fprintf(output, "  Temperature: [%lf, %lf]\n", parameters.T_min, parameters.T_max);
    fprintf(output, "  Fuel: [%lf, %lf]\n", parameters.Y_min, parameters.Y_max);
    if (parameters.variable_density == 1) {
        fprintf(output, "Pressure solver for variable density:\n");
        fprintf(output, "  Tolerance: %e, Max iterations: %d\n", parameters.pressure_solver_tol, parameters.pressure_solver_iter);
    }
    if (to_file)
        fclose(output);
}

void log_parameters(Parameters parameters) {
    log_parameters_report(parameters, 0);
    log_parameters_report(parameters, 1);
}

void log_message(Parameters parameters, const char *message) {
    setbuf(stdout, NULL);
    fprintf(parameters.log_files.log, "%s\n", message);
    printf("%s\n", message);
}

void log_timestep(Parameters parameters, int n, double t_n, double step_time, double CFL, double T_min, double T_max, double Y_min, double Y_max) {
    setbuf(stdout, NULL);
    printf("\nSimulation checkpoint\n");
    printf("Time step: %6d, Simulation time: %lf s\n", n, t_n);
    printf("CFL: %lf\n", CFL);
    printf("Temperature: Min = %.2f, Max = %.2f\n", T_min, T_max);
    printf("Fuel: Min = %.2f, Max = %.2f\n", Y_min, Y_max);
    printf("Step time: %lf s\n", step_time);
    // printf("Saving data...\n");
    fprintf(parameters.log_files.log, "\nSimulation checkpoint\n");
    fprintf(parameters.log_files.log, "Time step: %6d, Simulation time: %lf s\n", n, t_n);
    fprintf(parameters.log_files.log, "CFL: %lf\n", CFL);    
    fprintf(parameters.log_files.log, "Temperature: Min = %.2f, Max = %.2f\n", T_min, T_max);
    fprintf(parameters.log_files.log, "Fuel: Min = %.2f, Max = %.2f\n", Y_min, Y_max); 
    fprintf(parameters.log_files.log, "Step time: %lf s\n", step_time); 
    // fprintf(parameters.log_files.log, "Saving data...\n");
}
