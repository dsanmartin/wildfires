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
        if (strncmp(line, "alpha_u =", 9) == 0) {
            sscanf(line + 10, "%lf", &(parameters.alpha_u));
        }
    }
    // Compute dx, dy, dz, dts
    parameters.dx = (parameters.x_max - parameters.x_min) / (parameters.Nx - 1);
    parameters.dy = (parameters.y_max - parameters.y_min) / (parameters.Ny - 1);
    parameters.dz = (parameters.z_max - parameters.z_min) / (parameters.Nz - 1);
    parameters.dt = (parameters.t_max - parameters.t_min) / (parameters.Nt - 1);
    fclose(parameters_file);
    return parameters;
}
