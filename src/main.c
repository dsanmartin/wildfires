#include <stdio.h>
#include <stdlib.h>
#include "../include/structures.h"
#include "../include/parameters.h"
#include "../include/functions.h"
#include "../include/logs.h"
#include "../include/output.h"

int main(int argc, char *argv[]) { 
    char *parameters_file_path;
    double *u, *v, *w, *T;

    // Get parameters input file path
    if (argc == 2) {
        parameters_file_path = argv[1];
    } else {
        printf("Usage: %s <parameters_file_path>\n", argv[0]);
        return 1;
    }

    Parameters parameters = read_parameters_file(parameters_file_path);
    log_parameters(&parameters, 0);
    // log_parameters(&parameters, 1);

    // Allocate memory for x, y, z, u, v, w
    u = (double *) malloc((parameters.Nx - 1) * (parameters.Ny - 1) * parameters.Nz * sizeof(double));
    v = (double *) malloc((parameters.Nx - 1) * (parameters.Ny - 1) * parameters.Nz * sizeof(double));
    w = (double *) malloc((parameters.Nx - 1) * (parameters.Ny - 1) * parameters.Nz * sizeof(double));
    T = (double *) malloc((parameters.Nx - 1) * (parameters.Ny - 1) * parameters.Nz * sizeof(double));
    // u = (double *) malloc(parameters.Nx * parameters.Ny * parameters.Nz * sizeof(double));
    // v = (double *) malloc(parameters.Nx * parameters.Ny * parameters.Nz * sizeof(double));
    // w = (double *) malloc(parameters.Nx * parameters.Ny * parameters.Nz * sizeof(double));
    // T = (double *) malloc(parameters.Nx * parameters.Ny * parameters.Nz * sizeof(double));

    // Initialize u, v, w
    power_law_initial_condition(parameters.x, parameters.y, parameters.z, u, v, w, &parameters);
    
    // Initialize T
    gaussian_temperature_initial_condition(parameters.x, parameters.y, parameters.z, T, parameters);

    // Save data
    // save_data("data/output/u.csv", x, y, z, u, parameters.Nx, parameters.Ny, parameters.Nz);
    // save_data("data/output/v.csv", x, y, z, v, parameters.Nx, parameters.Ny, parameters.Nz);
    // save_data("data/output/w.csv", x, y, z, w, parameters.Nx, parameters.Ny, parameters.Nz);
    save_data_periodic_xy("data/output/T.csv", parameters.x, parameters.y, parameters.z, T, parameters.Nx, parameters.Ny, parameters.Nz);

    // Free memory
    free(u);
    free(v);
    free(w);

    return 0;
}
