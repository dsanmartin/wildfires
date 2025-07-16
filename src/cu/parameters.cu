/**
 * @file parameters.c
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Functions to handle the simulation parameters from a file.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include "../../include/cu/parameters.cuh"

Parameters read_parameters_file(const char *file_path) {
    FILE *parameters_file = fopen(file_path, "r"); // Open file
    int size; // Size of the fields
    char line[MAX_LINE_LENGTH]; // Create line buffer
    char parameter_file_path[256]; // Path to the parameters file
    char log_file_path[256]; // Path to the log file
    Parameters parameters; // Create parameters struct
    memset(&parameters, 0, sizeof(Parameters)); // Initialize parameters struct
    generate_current_datetime_string(parameters.sim_id, 32); // Get YYYYMMDD for simulation ID
    strcpy(parameters.save_path, "data/output/"); // Set save path
    // Concatenate simulation ID to save path
    strcat(parameters.save_path, parameters.sim_id);
    strcat(parameters.save_path, "/");
    // Create directory if it does not exist
    mkdir(parameters.save_path, 0777);
    // Concatenate save path to parameters file path
    strcpy(parameter_file_path, parameters.save_path);
    strcat(parameter_file_path, "parameters.txt");
    // Create file for parameters.txt
    parameters.log_files.parameters = fopen(parameter_file_path, "w");
    // Concatenate save path to log file path
    strcpy(log_file_path, parameters.save_path);
    strcat(log_file_path, "log.txt");
    // Create file for log.txt
    parameters.log_files.log = fopen(log_file_path, "w");
    // Define default values for parameters
    parameters.mu = 1.802e-5; // Dynamic viscosity of air in kg/(m*s)
    parameters.nu = 1.47e-5; // Kinematic viscosity of air in m^2/s
    parameters.z_r = 2.0; // Reference height for the velocity field in m
    parameters.alpha_u = 1.0 / 7.0; // Reference height for the velocity field in m
    parameters.alpha = 2.009e-5; // Thermal diffusivity of air in m^2/s
    parameters.T_inf = 288.15; // Reference temperature in K (15 °C)
    parameters.C_d = 0.15; // Solid fuel fraction drag coefficient (C_d)
    parameters.H_C = 21e6; // Heat of combustion in J/kg (typical value for wood)
    parameters.A = 1e9; // Pre-exponential factor in 1/s (typical value for wood combustion)
    parameters.T_act = 18040.8533; // Activation temperature in K (typical value for wood combustion)
    parameters.T_pc = 523.0; // Temperature phase-change in K (typical value for wood combustion)
    parameters.h_c = 1.42; // Convective heat transfer coefficient in W/(m^2*K) (typical value for air)
    parameters.c_p = 1007.0; // Heat capacity in J/(kg*K) (typical value for air at 15 °C)
    parameters.kappa = 0.02476; // Thermal conductivity in W/(m*K) (typical value for air at 15 °C)
    parameters.delta = 0.1; // Optical path length in m (default value for FDS)
    parameters.rho_inf = 1.225; // Reference density in kg/m^3 (typical value for air at 15 °C)
    parameters.g = -9.807; // Acceleration due to gravity in m/s^2 (typical value for Earth)
    parameters.C_s = 0.2; // Smagorinsky coefficient (default value)
    parameters.Pr = 0.7329; // Prandtl number (default value)
    parameters.variable_density = 1; // Variable density (default value)
    parameters.pressure_solver_tol = 1e-10; // Pressure solver tolerance (default value)
    parameters.pressure_solver_iter = 50; // Pressure solver maximum iterations (default value)
    parameters.pressure_solver_log = 0; // Pressure solver log (default value)
    parameters.n_threads = 1; // Number of threads (default value)
    parameters.T_min = 288.15; // Minimum temperature in K (default value)
    parameters.T_max = 2500.0; // Maximum temperature in K (default value)
    parameters.Y_min = 0.0; // Minimum fuel fraction (default value)
    parameters.Y_max = 1.0; // Maximum fuel fraction (default value)
    parameters.u_dead_nodes = 0.0; // Default value for dead nodes velocity in x direction
    parameters.v_dead_nodes = 0.0; // Default value for dead nodes velocity in y direction
    parameters.w_dead_nodes = 0.0; // Default value for dead nodes velocity in z direction
    parameters.T_dead_nodes = 288.15; // Default value for dead nodes temperature in K
    parameters.Y_dead_nodes = 0.0; // Default value for dead nodes fuel fraction
    parameters.input_path[0] = '\0'; // Initialize input path to empty string
    parameters.velocity_correction_fd = 0; // Default value for velocity correction finite difference 
    // Read file line by line
    while (fgets(line, MAX_LINE_LENGTH, parameters_file) != NULL) {
        if (strncmp(line, "Nx =", 4) == 0) {
            sscanf(line + 5, "%d", &(parameters.Nx));
        } 
        if (strncmp(line, "Ny =", 4) == 0) {
            sscanf(line + 5, "%d", &(parameters.Ny));
        }
        if (strncmp(line, "Nz =", 4) == 0) {
            sscanf(line + 5, "%d", &(parameters.Nz));
        }
        if (strncmp(line, "Nt =", 4) == 0) {
            sscanf(line + 5, "%d", &(parameters.Nt));
        }
        if (strncmp(line, "NT =", 4) == 0) {
            sscanf(line + 5, "%d", &(parameters.NT));
        }
        if (strncmp(line, "x_min =", 7) == 0) {
            sscanf(line + 8, "%lf", &(parameters.x_min));
        }
        if (strncmp(line, "x_max =", 7) == 0) {
            sscanf(line + 8, "%lf", &(parameters.x_max));
        }
        if (strncmp(line, "y_min =", 7) == 0) {
            sscanf(line + 8, "%lf", &(parameters.y_min));
        }
        if (strncmp(line, "y_max =", 7) == 0) {
            sscanf(line + 8, "%lf", &(parameters.y_max));
        }
        if (strncmp(line, "z_min =", 7) == 0) {
            sscanf(line + 8, "%lf", &(parameters.z_min));
        }
        if (strncmp(line, "z_max =", 7) == 0) {
            sscanf(line + 8, "%lf", &(parameters.z_max));
        }
        if (strncmp(line, "t_min =", 7) == 0) {
            sscanf(line + 8, "%lf", &(parameters.t_min));
        }
        if (strncmp(line, "t_max =", 7) == 0) {
            sscanf(line + 8, "%lf", &(parameters.t_max));
        } 
        if (strncmp(line, "nu =", 4) == 0) {
            sscanf(line + 5, "%lf", &(parameters.nu));
        }         
        if (strncmp(line, "mu =", 4) == 0) {
            sscanf(line + 5, "%lf", &(parameters.mu));
        } 
        if (strncmp(line, "u_r =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.u_r));
        }
        if (strncmp(line, "z_r =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.z_r));
        } 
        if (strncmp(line, "alpha_u =", 9) == 0) {
            sscanf(line + 10, "%lf", &(parameters.alpha_u));
        } 
        if (strncmp(line, "alpha =", 7) == 0) {
            sscanf(line + 8, "%lf", &(parameters.alpha));
        } 
        if (strncmp(line, "T_inf =", 7) == 0) {
            sscanf(line + 8, "%lf", &(parameters.T_inf));
        } 
        if (strncmp(line, "T_source =", 10) == 0) {
            sscanf(line + 11, "%lf", &(parameters.T_source));
        }
        if (strncmp(line, "T0_shape =", 10) == 0) {
            sscanf(line + 11, "%s", parameters.T0_shape);
        }
        if (strncmp(line, "T0_x_start =", 12) == 0) {
            sscanf(line + 13, "%lf", &(parameters.T0_x_start));
        }
        if (strncmp(line, "T0_x_end =", 10) == 0) {
            sscanf(line + 11, "%lf", &(parameters.T0_x_end));
        }
        if (strncmp(line, "T0_y_start =", 12) == 0) {
            sscanf(line + 13, "%lf", &(parameters.T0_y_start));
        }
        if (strncmp(line, "T0_y_end =", 10) == 0) {
            sscanf(line + 11, "%lf", &(parameters.T0_y_end));
        }
        if (strncmp(line, "T0_z_start =", 12) == 0) {
            sscanf(line + 13, "%lf", &(parameters.T0_z_start));
        }
        if (strncmp(line, "T0_z_end =", 10) == 0) {
            sscanf(line + 11, "%lf", &(parameters.T0_z_end));
        }        
        // Fuel initial condition
        if (strncmp(line, "Y0_x_start =", 12) == 0) {
            sscanf(line + 13, "%lf", &(parameters.Y0_x_start));
        }
        if (strncmp(line, "Y0_x_end =", 10) == 0) {
            sscanf(line + 11, "%lf", &(parameters.Y0_x_end));
        }
        if (strncmp(line, "Y0_y_start =", 12) == 0) {
            sscanf(line + 13, "%lf", &(parameters.Y0_y_start));
        }
        if (strncmp(line, "Y0_y_end =", 10) == 0) {
            sscanf(line + 11, "%lf", &(parameters.Y0_y_end));
        }
        if (strncmp(line, "Y0_xa =", 7) == 0) {
            sscanf(line + 8, "%lf", &(parameters.Y0_xa));
        }
        if (strncmp(line, "Y0_ya =", 7) == 0) {
            sscanf(line + 8, "%lf", &(parameters.Y0_ya));
        }
        // if (strncmp(line, "fuel_relax =", 12) == 0) {
        //     sscanf(line + 13, "%d", &(parameters.fuel_relax));
        // }
        // Fuel height in m
        if (strncmp(line, "Y_h =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.Y_h));
        }
        if (strncmp(line, "Y_f =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.Y_f));
        }
        if (strncmp(line, "C_d =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.C_d));
        }
        // Solid volume fraction alpha_s
        if (strncmp(line, "alpha_s =", 9) == 0) {
            sscanf(line + 10, "%lf", &(parameters.alpha_s));
        }
        // Surface-to-volume ratio sigma_s
        if (strncmp(line, "sigma_s =", 9) == 0) {
            sscanf(line + 10, "%lf", &(parameters.sigma_s));
        }
        // Heat of combustion H_C in J/kg
        if (strncmp(line, "H_C =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.H_C));
        } 
        // Pre-exponential factor in 1/s
        if (strncmp(line, "A =", 3) == 0) {
            sscanf(line + 4, "%lf", &(parameters.A));
        } 
        // Activation temperature in K
        if (strncmp(line, "T_act =", 7) == 0) {
            sscanf(line + 8, "%lf", &(parameters.T_act));
        } 
        // Temperature phase-change T_pc
        if (strncmp(line, "T_pc =", 6) == 0) {
            sscanf(line + 7, "%lf", &(parameters.T_pc));
        } 
        // Convective heat transfer coefficient h
        if (strncmp(line, "h_c =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.h_c));
        } 
        // // Volumetric heat capacity a_v
        // if (strncmp(line, "a_v =", 5) == 0) {
        //     sscanf(line + 6, "%lf", &(parameters.a_v));
        // }
        // Heat capacity c_p
        if (strncmp(line, "c_p =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.c_p));
        }
        // Reference density rho_0
        if (strncmp(line, "rho_inf =", 9) == 0) {
            sscanf(line + 10, "%lf", &(parameters.rho_inf));
        } 
        // U0 Type
        if (strncmp(line, "U0_type =", 9) == 0) {
            sscanf(line + 10, "%s", (char *) &(parameters.U0_type));
        }
        // T0 Shape
        if (strncmp(line, "T0_shape =", 10) == 0) {
            sscanf(line + 11, "%s", (char *) &(parameters.T0_shape));
        }
        // Acceleration due to gravity g in m/s^2
        if (strncmp(line, "g =", 3) == 0) {
            sscanf(line + 4, "%lf", &(parameters.g));
        }
        // Smagorinsky coefficient C_s
        if (strncmp(line, "C_s =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.C_s));
        } 
        // Prandtl number Pr
        if (strncmp(line, "Pr =", 4) == 0) {
            sscanf(line + 5, "%lf", &(parameters.Pr));
        } 
        // p_top
        if (strncmp(line, "p_top =", 7) == 0) {
            sscanf(line + 8, "%lf", &(parameters.p_top));
        }
        // Time integration method
        if (strncmp(line, "method =", 8) == 0) {
            sscanf(line + 9, "%s", (char *) &(parameters.method));
        }
        // Number of threads
        if (strncmp(line, "threads =", 9) == 0) {
            sscanf(line + 10, "%d", &(parameters.n_threads));
        } 
        // Temperature bounds
        if (strncmp(line, "T_min =", 6) == 0) {
            sscanf(line + 7, "%lf", &(parameters.T_min));
        } 
        if (strncmp(line, "T_max =", 6) == 0) {
            sscanf(line + 7, "%lf", &(parameters.T_max));
        } 
        // Fuel bounds
        if (strncmp(line, "Y_min =", 6) == 0) {
            sscanf(line + 7, "%lf", &(parameters.Y_min));
        } 
        if (strncmp(line, "Y_max =", 6) == 0) {
            sscanf(line + 7, "%lf", &(parameters.Y_max));
        } 
        // Topography shape
        if (strncmp(line, "topo_shape =", 11) == 0) {
            sscanf(line + 12, "%s", (char *) &(parameters.topo_shape));
        }
        // IBM dead nodes values
        if (strncmp(line, "u_dead_nodes =", 13) == 0) {
            sscanf(line + 14, "%lf", &(parameters.u_dead_nodes));
        } 
        if (strncmp(line, "v_dead_nodes =", 13) == 0) {
            sscanf(line + 14, "%lf", &(parameters.v_dead_nodes));
        } 
        if (strncmp(line, "w_dead_nodes =", 13) == 0) {
            sscanf(line + 14, "%lf", &(parameters.w_dead_nodes));
        } 
        if (strncmp(line, "T_dead_nodes =", 13) == 0) {
            sscanf(line + 14, "%lf", &(parameters.T_dead_nodes));
        } 
        if (strncmp(line, "Y_dead_nodes =", 13) == 0) {
            sscanf(line + 14, "%lf", &(parameters.Y_dead_nodes));
        } 
        // Hill parameters
        if (strncmp(line, "hill_center_x =", 14) == 0) {
            sscanf(line + 15, "%lf", &(parameters.hill_center_x));
        }
        if (strncmp(line, "hill_center_y =", 14) == 0) {
            sscanf(line + 15, "%lf", &(parameters.hill_center_y));
        }
        if (strncmp(line, "hill_length =", 12) == 0) {
            sscanf(line + 13, "%lf", &(parameters.hill_length));
        }
        if (strncmp(line, "hill_width =", 11) == 0) {
            sscanf(line + 12, "%lf", &(parameters.hill_width));
        }
        if (strncmp(line, "hill_height =", 12) == 0) {
            sscanf(line + 13, "%lf", &(parameters.hill_height));
        }
        // Time of temperature source
        if (strncmp(line, "t_source_start =", 16) == 0) {
            sscanf(line + 17, "%lf", &(parameters.t_source_start));
        }
        if (strncmp(line, "t_source_end =", 14) == 0) {
            sscanf(line + 15, "%lf", &(parameters.t_source_end));
        }
        // Optical path length (\delta) in m
        if (strncmp(line, "delta =", 7) == 0) {
            sscanf(line + 8, "%lf", &(parameters.delta));
        } 
        // Thermal conductivity (\kappa) in W/(m*K)
        if (strncmp(line, "kappa =", 6) == 0) {
            sscanf(line + 7, "%lf", &(parameters.kappa));
        } 
        // Variable or constant density
        if (strncmp(line, "variable_density =", 17) == 0) {
            sscanf(line + 18, "%d", &(parameters.variable_density));
        } 
        // Pressure solver tolerance (pressure_solver_tol)
        if (strncmp(line, "pressure_solver_tol =", 21) == 0) {
            sscanf(line + 22, "%lf", &(parameters.pressure_solver_tol));
        } 
        // Pressure solver maximum iterations (pressure_solver_iter)
        if (strncmp(line, "pressure_solver_iter =", 22) == 0) {
            sscanf(line + 23, "%d", &(parameters.pressure_solver_iter));
        } 
        // Pressure solver log (pressure_solver_log)
        if (strncmp(line, "pressure_solver_log =", 21) == 0) {
            sscanf(line + 22, "%d", &(parameters.pressure_solver_log));
        }
        // Test
        // Check z domain
        if (strncmp(line, "z_domain =", 10) == 0) {
            sscanf(line + 11, "%d", &(parameters.z_domain));
        }
        // Find k_nu_grid parameter
        if (strncmp(line, "k_nu_grid =", 11) == 0) {
            sscanf(line + 12, "%lf", &(parameters.k_nu_grid));
        }
        // Find k_transition parameter
        if (strncmp(line, "k_transition =", 14) == 0) {
            sscanf(line + 15, "%d", &(parameters.k_transition));
        }
        // Find z_transition parameter
        if (strncmp(line, "z_transition =", 14) == 0) {
            sscanf(line + 15, "%lf", &(parameters.z_transition));
        }
        // Find input path
        if (strncmp(line, "input_path =", 12) == 0) {
            sscanf(line + 13, "%s", parameters.input_path);
        }
        // Velocity correction finite difference
        if (strncmp(line, "velocity_correction_fd =", 24) == 0) {
            sscanf(line + 25, "%d", &(parameters.velocity_correction_fd));
        }
    }
    // Initialize x, y, z, t
    parameters.x = (double *) malloc(parameters.Nx * sizeof(double));
    parameters.y = (double *) malloc(parameters.Ny * sizeof(double));
    parameters.z = (double *) malloc(parameters.Nz * sizeof(double));
    parameters.t = (double *) malloc((parameters.Nt + 1) * sizeof(double));
    // Initialize r, s
    parameters.r = (double *) malloc((parameters.Nx - 1) * sizeof(double));
    parameters.s = (double *) malloc((parameters.Ny - 1) * sizeof(double));
    parameters.kx = (double *) malloc((parameters.Nx - 1) * sizeof(double));
    parameters.ky = (double *) malloc((parameters.Ny - 1) * sizeof(double));
    // IBM malloc
    // Allocate memory
    parameters.cut_nodes = (int *)malloc(parameters.Nx * parameters.Ny * sizeof(int)); // Cut nodes
    parameters.z_ibm = (double *)malloc(parameters.Nx * parameters.Ny * parameters.Nz * sizeof(double)); // IBM z distance
    parameters.topography = (double *)malloc(parameters.Nx * parameters.Ny * sizeof(double)); // Topography
    // Nz for Y variable. To store the number of grid points in z for Y variable
    parameters.Nz_Y = (int *)malloc(parameters.Nx * parameters.Ny * sizeof(int));
    // Compute dx, dy, dz, dts
    parameters.dx = (parameters.x_max - parameters.x_min) / (parameters.Nx - 1);
    parameters.dy = (parameters.y_max - parameters.y_min) / (parameters.Ny - 1);
    parameters.dz = (parameters.z_max - parameters.z_min) / (parameters.Nz - 1);
    parameters.dt = (parameters.t_max - parameters.t_min) / parameters.Nt;
    // Frequency domain
    fft_freq(parameters.r, parameters.Nx - 1, 1.0 / (parameters.Nx - 1));
    fft_freq(parameters.s, parameters.Ny - 1, 1.0 / (parameters.Ny - 1));
    // Physical domain and frequency domain
    // t
    for (int n = 0; n <= parameters.Nt; n++) {
        parameters.t[n] = parameters.t_min + n * parameters.dt;
    }
    // x and kx
    for (int i = 0; i < parameters.Nx; i++) {
        parameters.x[i] = parameters.x_min + i * parameters.dx;
        if (i < parameters.Nx - 1) {
            // parameters.kx[i] = 2 * M_PI * parameters.r[i] * parameters.dz / (parameters.x_max - parameters.x_min);
            parameters.kx[i] = 2 * M_PI * parameters.r[i] / (parameters.x_max - parameters.x_min);
        }
    }
    // y and ky
    for (int j = 0; j < parameters.Ny; j++) {
        parameters.y[j] = parameters.y_min + j * parameters.dy;
        if (j < parameters.Ny - 1) {
            // parameters.ky[j] = 2 * M_PI * parameters.s[j] * parameters.dz / (parameters.y_max - parameters.y_min);
            parameters.ky[j] = 2 * M_PI * parameters.s[j] / (parameters.y_max - parameters.y_min);
        }
    }
    // z
    if (parameters.z_domain == 0) {
        equispaced_domain(parameters.z, parameters.Nz, parameters.z_min, parameters.dz);
    } else {
        if (parameters.z_domain == 1) {
            non_equispaced_domain(parameters.z, parameters.Nz, parameters.z_min, parameters.z_max, parameters.k_nu_grid);
        } else {
            if (parameters.z_domain == 2) {
                transition_domain(parameters.z, parameters.Nz, parameters.z_min, parameters.z_max, parameters.z_transition, parameters.k_transition);
            }
        }
    }
    // Computer T0 centers and dimensiones
    parameters.T0_x_center = (parameters.T0_x_start + parameters.T0_x_end) / 2;
    parameters.T0_y_center = (parameters.T0_y_start + parameters.T0_y_end) / 2;
    parameters.T0_z_center = 0.0;
    parameters.T0_length = parameters.T0_x_end - parameters.T0_x_start;
    parameters.T0_width = parameters.T0_y_end - parameters.T0_y_start;
    parameters.T0_height = parameters.T0_z_end - parameters.T0_z_start;
    // IBM topography
    if (strcmp(parameters.topo_shape, "hill") == 0)
        simple_hill(&parameters);
    else if (strcmp(parameters.topo_shape, "flat") == 0)
        flat_terrain(&parameters);
    ibm_parameters(&parameters);  
    // Test a cube in the domain
    // cube(&parameters);
    // Sizes
    size = parameters.Nx * parameters.Ny * parameters.Nz;
    // Fill field indexes
    parameters.field_indexes.u = 0;
    parameters.field_indexes.v = size;
    parameters.field_indexes.w = 2 * size;
    parameters.field_indexes.T = 3 * size;
    parameters.field_indexes.Y = 4 * size;
    // Fill turbulence indexes
    parameters.turbulence_indexes.rho = 0;
    parameters.turbulence_indexes.Tx = size;
    parameters.turbulence_indexes.Ty = 2 * size;
    parameters.turbulence_indexes.Tz = 3 * size;
    parameters.turbulence_indexes.Txx = 4 * size;
    parameters.turbulence_indexes.Tyy = 5 * size;
    parameters.turbulence_indexes.Tzz = 6 * size;
    parameters.turbulence_indexes.mu_sgs = 7 * size;
    parameters.turbulence_indexes.S_11 = 8 * size;
    parameters.turbulence_indexes.S_12 = 9 * size;
    parameters.turbulence_indexes.S_13 = 10 * size;
    parameters.turbulence_indexes.S_21 = 11 * size;
    parameters.turbulence_indexes.S_22 = 12 * size;
    parameters.turbulence_indexes.S_23 = 13 * size;
    parameters.turbulence_indexes.S_31 = 14 * size;
    parameters.turbulence_indexes.S_32 = 15 * size;
    parameters.turbulence_indexes.S_33 = 16 * size;
    parameters.turbulence_indexes.div_U = 17 * size;
    // parameters.turbulence_indexes.ux = 0;
    // parameters.turbulence_indexes.uy = size;
    // parameters.turbulence_indexes.uz = 2 * size;
    // parameters.turbulence_indexes.vx = 3 * size;
    // parameters.turbulence_indexes.vy = 4 * size;
    // parameters.turbulence_indexes.vz = 5 * size;
    // parameters.turbulence_indexes.wx = 6 * size;
    // parameters.turbulence_indexes.wy = 7 * size;
    // parameters.turbulence_indexes.wz = 8 * size;
    // parameters.turbulence_indexes.Tx = 9 * size;
    // parameters.turbulence_indexes.Ty = 10 * size;
    // parameters.turbulence_indexes.Tz = 11 * size;
    // parameters.turbulence_indexes.uxx = 12 * size;
    // parameters.turbulence_indexes.uyy = 13 * size;
    // parameters.turbulence_indexes.uzz = 14 * size;
    // parameters.turbulence_indexes.vxx = 15 * size;
    // parameters.turbulence_indexes.vyy = 16 * size;
    // parameters.turbulence_indexes.vzz = 17 * size;
    // parameters.turbulence_indexes.wxx = 18 * size;
    // parameters.turbulence_indexes.wyy = 19 * size;
    // parameters.turbulence_indexes.wzz = 20 * size;
    // parameters.turbulence_indexes.Txx = 21 * size;
    // parameters.turbulence_indexes.Tyy = 22 * size;
    // parameters.turbulence_indexes.Tzz = 23 * size;
    // parameters.turbulence_indexes.fw = 24 * size;
    // parameters.turbulence_indexes.mu_sgs = 25 * size;
    // parameters.turbulence_indexes.S_11 = 26 * size;
    // Close file
    fclose(parameters_file);
    return parameters;
}

void free_parameters(Parameters parameters) {
    free(parameters.x);
    free(parameters.y);
    free(parameters.z);
    free(parameters.t);
    free(parameters.r);
    free(parameters.s);
    free(parameters.kx);
    free(parameters.ky);
    free(parameters.cut_nodes);
    free(parameters.z_ibm);
    free(parameters.topography);
    free(parameters.Nz_Y);
    // fclose(parameters.log_files.parameters);
    fclose(parameters.log_files.log);
}
