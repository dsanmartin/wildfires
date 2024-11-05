#include "../../include/cu/logs.cuh"

void log_parameters_report(Parameters parameters, int to_file) {
    FILE *output;
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
    if (parameters.z_uniform == 1) {
        fprintf(output, "  dx: %lf, dy: %lf, dz: %lf, dt: %lf\n", parameters.dx, parameters.dy, parameters.dz, parameters.dt);
    } else {
        fprintf(output, "  dx: %lf, dy: %lf, dz: [%lf, %lf], dt: %lf\n", 
            parameters.dx, parameters.dy, 
            parameters.z[1] - parameters.z[0],
            parameters.z[parameters.Nz - 1] - parameters.z[parameters.Nz - 2], 
            parameters.dt
        );
    }
    
    fprintf(output, "  Time samples: %d\n", parameters.NT);
    fprintf(output, "Time integration: %s\n", parameters.method);
    fprintf(output, "Fluid parameters:\n");
    fprintf(output, "  rho: %lf, nu: %e, alpha: %e, T_inf: %lf,  g: %lf\n", parameters.rho, parameters.nu, parameters.alpha, parameters.T_inf, parameters.g);
    fprintf(output, "  C_s: %lf, Pr: %lf\n", parameters.C_s, parameters.Pr);
    fprintf(output, "Fuel parameters:\n");
    fprintf(output, "  H_R: %e, c_p: %lf, h: %lf, a_v: %lf\n", parameters.H_R, parameters.c_p, parameters.h, parameters.a_v);
    fprintf(output, "  Y_D: %lf, Y_f: %lf\n", parameters.Y_D, parameters.Y_f);
    fprintf(output, "  A: %e, T_pc: %lf, T_a: %lf\n", parameters.A, parameters.T_pc, parameters.T_a);
    fprintf(output, "Wind initial condition:\n");
    fprintf(output, "  Type: %s\n", parameters.U0_type);    
    if (strcmp(parameters.U0_type, "constant") == 0)
        fprintf(output, "  u_r: %lf\n", parameters.u_r);
    else if (strcmp(parameters.U0_type, "log") == 0)
        fprintf(output, "  u_z0: %lf, kappa: %lf, d: %lf, u_ast: %lf\n", parameters.u_z0, parameters.kappa, parameters.d, parameters.u_ast);
    else if (strcmp(parameters.U0_type, "power_law") == 0)
        fprintf(output, "  u_r: %lf, z_r: %lf, alpha_u: %lf\n", parameters.u_r, parameters.z_r, parameters.alpha_u);
    fprintf(output, "Temperature initial condition:\n");
    fprintf(output, "  Shape: %s\n", parameters.T0_shape);
    fprintf(output, "  T_hot: %lf\n", parameters.T_hot);
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
    printf("Time step: %6d, Simulation time: %lf s\n", n, t_n);
    printf("CFL: %lf\n", CFL);
    printf("Temperature: Min = %.2f, Max = %.2f\n", T_min, T_max);
    printf("Fuel: Min = %.2f, Max = %.2f\n", Y_min, Y_max);
    printf("Step time: %lf s\n", step_time);
    // printf("Saving data...\n");
    fprintf(parameters.log_files.log, "Time step: %6d, Simulation time: %lf s\n", n, t_n);
    fprintf(parameters.log_files.log, "CFL: %lf\n", CFL);    
    fprintf(parameters.log_files.log, "Temperature: Min = %.2f, Max = %.2f\n", T_min, T_max);
    fprintf(parameters.log_files.log, "Fuel: Min = %.2f, Max = %.2f\n", Y_min, Y_max); 
    fprintf(parameters.log_files.log, "Step time: %lf s\n", step_time); 
    // fprintf(parameters.log_files.log, "Saving data...\n");
}