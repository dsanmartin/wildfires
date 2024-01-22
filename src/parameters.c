#include "../include/parameters.h"

#define MAX_LINE_LENGTH 512

/**
 * @brief Reads the parameters from a file.
 *
 * This function reads the parameters from a file specified by the filePath parameter.
 * The file should have the following format:
 * - Nx = <value>
 * - Ny = <value>
 * - Nz = <value>
 * - Nt = <value>
 * - x_min = <value>
 * - x_max = <value>
 * - y_min = <value>
 * - y_max = <value>
 * - z_min = <value>
 * - z_max = <value>
 * - t_min = <value>
 * - t_max = <value>
 *
 * @param filePath The path to the parameters file.
 * @return The parameters read from the file.
 */
Parameters read_parameters_file(const char *filePath) {
    FILE *parameters_file = fopen(filePath, "r"); // Open file
    Parameters parameters; // Create parameters struct
    char line[MAX_LINE_LENGTH]; // Create line buffer
    int size; // Size of the fields
    memset(&parameters, 0, sizeof(Parameters)); // Initialize parameters struct
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
        if (strncmp(line, "T_hot =", 7) == 0) {
            sscanf(line + 8, "%lf", &(parameters.T_hot));
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
        // if (strncmp(line, "T0_x_center =", 13) == 0) {
        //     sscanf(line + 14, "%lf", &(parameters.T0_x_center));
        // }
        // if (strncmp(line, "T0_length =", 11) == 0) {
        //     sscanf(line + 12, "%lf", &(parameters.T0_length));
        // }
        if (strncmp(line, "T0_y_start =", 12) == 0) {
            sscanf(line + 13, "%lf", &(parameters.T0_y_start));
        }
        if (strncmp(line, "T0_y_end =", 10) == 0) {
            sscanf(line + 11, "%lf", &(parameters.T0_y_end));
        }
        // if (strncmp(line, "T0_y_center =", 13) == 0) {
        //     sscanf(line + 14, "%lf", &(parameters.T0_y_center));
        // }
        // if (strncmp(line, "T0_width =", 10) == 0) {
        //     sscanf(line + 11, "%lf", &(parameters.T0_width));
        // }
        if (strncmp(line, "T0_z_start =", 12) == 0) {
            sscanf(line + 13, "%lf", &(parameters.T0_z_start));
        }
        if (strncmp(line, "T0_z_end =", 10) == 0) {
            sscanf(line + 11, "%lf", &(parameters.T0_z_end));
        }
        // if (strncmp(line, "T0_z_center =", 13) == 0) {
        //     sscanf(line + 14, "%lf", &(parameters.T0_z_center));
        // }
        // if (strncmp(line, "T0_height =", 11) == 0) {
        //     sscanf(line + 12, "%lf", &(parameters.T0_height));
        // }
        if (strncmp(line, "Y_h =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.Y_h));
        }
        if (strncmp(line, "Y_f =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.Y_f));
        }
        // Fuel drag Y_D
        if (strncmp(line, "Y_D =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.Y_D));
        }
        // Read H_R 
        if (strncmp(line, "H_R =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.H_R));
        }
        // Read A
        if (strncmp(line, "A =", 3) == 0) {
            sscanf(line + 4, "%lf", &(parameters.A));
        }
        // Read activation temperature T_a
        if (strncmp(line, "T_a =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.T_a));
        }
        // Temperature phase-change T_pc
        if (strncmp(line, "T_pc =", 6) == 0) {
            sscanf(line + 7, "%lf", &(parameters.T_pc));
        }
        // Convective heat transfer coefficient h
        if (strncmp(line, "h =", 3) == 0) {
            sscanf(line + 4, "%lf", &(parameters.h));
        }
        // Volumetric heat capacity a_v
        if (strncmp(line, "a_v =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.a_v));
        }
        // Heat capacity c_p
        if (strncmp(line, "c_p =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.c_p));
        }
        // Density rho
        if (strncmp(line, "rho =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.rho));
        }
        // U0 Type
        if (strncmp(line, "U0_type =", 9) == 0) {
            sscanf(line + 10, "%s", &(parameters.U0_type));
        }
        // T0 Shape
        if (strncmp(line, "T0_shape =", 10) == 0) {
            sscanf(line + 11, "%s", &(parameters.T0_shape));
        }
        // Acceleration due to gravity g
        if (strncmp(line, "g =", 3) == 0) {
            sscanf(line + 4, "%lf", &(parameters.g));
        }
        // Smagorinsky constant C_s
        if (strncmp(line, "C_s =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.C_s));
        }
        // Prandtl number Pr
        if (strncmp(line, "Pr =", 4) == 0) {
            sscanf(line + 5, "%lf", &(parameters.Pr));
        }
    }
    // Compute dx, dy, dz, dts
    parameters.dx = (parameters.x_max - parameters.x_min) / (parameters.Nx - 1);
    parameters.dy = (parameters.y_max - parameters.y_min) / (parameters.Ny - 1);
    parameters.dz = (parameters.z_max - parameters.z_min) / (parameters.Nz - 1);
    parameters.dt = (parameters.t_max - parameters.t_min) / (parameters.Nt - 1);
    // Initialize x, y, z, t
    parameters.x = (double *) malloc(parameters.Nx * sizeof(double));
    parameters.y = (double *) malloc(parameters.Ny * sizeof(double));
    parameters.z = (double *) malloc(parameters.Nz * sizeof(double));
    parameters.t = (double *) malloc(parameters.Nt * sizeof(double));
    // Initialize r, s
    parameters.r = (double *) malloc((parameters.Nx - 1) * sizeof(double));
    parameters.s = (double *) malloc((parameters.Ny - 1) * sizeof(double));
    parameters.kx = (double *) malloc((parameters.Nx - 1) * sizeof(double));
    parameters.ky = (double *) malloc((parameters.Ny - 1) * sizeof(double));
    fft_freq(parameters.r, parameters.Nx - 1, 1 / (parameters.Nx - 1));
    fft_freq(parameters.s, parameters.Ny - 1, 1 / (parameters.Ny - 1));
    // Physical domain and frequency domain
    for (int i = 0; i < parameters.Nx; i++) {
        parameters.x[i] = parameters.x_min + i * parameters.dx;
        if (i < parameters.Nx - 1) {
            parameters.kx[i] = 2 * M_PI * parameters.r[i] * parameters.dz / parameters.x_max;
        }
    }
    for (int j = 0; j < parameters.Ny; j++) {
        parameters.y[j] = parameters.y_min + j * parameters.dy;
        if (j < parameters.Ny - 1) {
            parameters.ky[j] = 2 * M_PI * parameters.s[j] * parameters.dz / parameters.y_max;
        }
    }
    for (int k = 0; k < parameters.Nz; k++) {
        parameters.z[k] = parameters.z_min + k * parameters.dz;
        // Get the index of z where Y_h is located
        if (parameters.z[k] <= parameters.Y_h) {
            parameters.k_Y_h = k;
        }
    }
    for (int n = 0; n < parameters.Nt; n++) {
        parameters.t[n] = parameters.t_min + n * parameters.dt;
    }
    // Frequency domain
    
    // Computer T0 centers and dimensiones
    parameters.T0_x_center = (parameters.T0_x_start + parameters.T0_x_end) / 2;
    parameters.T0_y_center = (parameters.T0_y_start + parameters.T0_y_end) / 2;
    parameters.T0_z_center = 0.0;
    parameters.T0_length = parameters.T0_x_end - parameters.T0_x_start;
    parameters.T0_width = parameters.T0_y_end - parameters.T0_y_start;
    parameters.T0_height = parameters.T0_z_end - parameters.T0_z_start;

    // Sizes
    size = parameters.Nx * parameters.Ny * parameters.Nz;

    // Fill field indexes
    parameters.field_indexes.u = 0;
    parameters.field_indexes.v = size;
    parameters.field_indexes.w = 2 * size;
    parameters.field_indexes.T = 3 * size;
    parameters.field_indexes.Y = 4 * size;

    // Fill turbulence indexes
    parameters.turbulence_indexes.ux = 0;
    parameters.turbulence_indexes.uy = size;
    parameters.turbulence_indexes.uz = 2 * size;
    parameters.turbulence_indexes.vx = 3 * size;
    parameters.turbulence_indexes.vy = 4 * size;
    parameters.turbulence_indexes.vz = 5 * size;
    parameters.turbulence_indexes.wx = 6 * size;
    parameters.turbulence_indexes.wy = 7 * size;
    parameters.turbulence_indexes.wz = 8 * size;
    parameters.turbulence_indexes.Tx = 9 * size;
    parameters.turbulence_indexes.Ty = 10 * size;
    parameters.turbulence_indexes.Tz = 11 * size;
    parameters.turbulence_indexes.uxx = 12 * size;
    parameters.turbulence_indexes.uyy = 13 * size;
    parameters.turbulence_indexes.uzz = 14 * size;
    parameters.turbulence_indexes.vxx = 15 * size;
    parameters.turbulence_indexes.vyy = 16 * size;
    parameters.turbulence_indexes.vzz = 17 * size;
    parameters.turbulence_indexes.wxx = 18 * size;
    parameters.turbulence_indexes.wyy = 19 * size;
    parameters.turbulence_indexes.wzz = 20 * size;
    parameters.turbulence_indexes.Txx = 21 * size;
    parameters.turbulence_indexes.Tyy = 22 * size;
    parameters.turbulence_indexes.Tzz = 23 * size;
    // parameters.turbulence_indexes.uyx = 24 * size;
    // parameters.turbulence_indexes.uzx = 25 * size;
    // parameters.turbulence_indexes.uxy = 26 * size;
    // parameters.turbulence_indexes.uzy = 27 * size;
    // parameters.turbulence_indexes.uxz = 28 * size;
    // parameters.turbulence_indexes.uyz = 29 * size;
    // parameters.turbulence_indexes.vyx = 30 * size;
    // parameters.turbulence_indexes.vzx = 31 * size;
    // parameters.turbulence_indexes.vxy = 32 * size;
    // parameters.turbulence_indexes.vzy = 33 * size;
    // parameters.turbulence_indexes.vxz = 34 * size;
    // parameters.turbulence_indexes.vyz = 35 * size;
    // parameters.turbulence_indexes.wyx = 36 * size;
    // parameters.turbulence_indexes.wzx = 37 * size;
    // parameters.turbulence_indexes.wxy = 38 * size;
    // parameters.turbulence_indexes.wzy = 39 * size;
    // parameters.turbulence_indexes.wxz = 40 * size;
    // parameters.turbulence_indexes.wyz = 41 * size;
    parameters.turbulence_indexes.fw = 24 * size;
    parameters.turbulence_indexes.fwx = 25 * size;
    parameters.turbulence_indexes.fwy = 26 * size;
    parameters.turbulence_indexes.fwz = 27 * size;
    
    fclose(parameters_file);
    return parameters;
}
