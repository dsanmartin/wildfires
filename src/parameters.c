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
    memset(&parameters, 0, sizeof(Parameters)); // Initialize parameters struct
    char line[MAX_LINE_LENGTH]; // Create line buffer
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
        if (strncmp(line, "mu =", 4) == 0) {
            sscanf(line + 5, "%lf", &(parameters.mu));
        } 
        if (strncmp(line, "u_r =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.u_r));
        }
        if (strncmp(line, "z_r =", 5) == 0) {
            sscanf(line + 6, "%lf", &(parameters.z_r));
        }
        if (strncmp(line, "alpha =", 7) == 0) {
            sscanf(line + 8, "%lf", &(parameters.alpha));
        }
        if (strncmp(line, "alpha_u =", 9) == 0) {
            sscanf(line + 10, "%lf", &(parameters.alpha_u));
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
    for (int i = 0; i < parameters.Nx; i++) {
        parameters.x[i] = parameters.x_min + i * parameters.dx;
    }
    for (int j = 0; j < parameters.Ny; j++) {
        parameters.y[j] = parameters.y_min + j * parameters.dy;
    }
    for (int k = 0; k < parameters.Nz; k++) {
        parameters.z[k] = parameters.z_min + k * parameters.dz;
    }
    for (int n = 0; n < parameters.Nt; n++) {
        parameters.t[n] = parameters.t_min + n * parameters.dt;
    }
    // Computer T0 centers and dimensiones
    parameters.T0_x_center = (parameters.T0_x_start + parameters.T0_x_end) / 2;
    parameters.T0_y_center = (parameters.T0_y_start + parameters.T0_y_end) / 2;
    parameters.T0_z_center = 0.0;
    parameters.T0_length = parameters.T0_x_end - parameters.T0_x_start;
    parameters.T0_width = parameters.T0_y_end - parameters.T0_y_start;
    parameters.T0_height = parameters.T0_z_end - parameters.T0_z_start;
    fclose(parameters_file);
    return parameters;
}
