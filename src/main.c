#include <stdio.h>
#include <stdlib.h>
#include "../include/structures.h"
#include "../include/parameters.h"
#include "../include/functions.h"
#include "../include/logs.h"
#include "../include/output.h"
#include "../include/pde.h"

int main(int argc, char *argv[]) { 
    char *parameters_file_path;
    double *u, *v, *w, *T, *Y, *p, *y_0;
    int size, size_Y;

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
    // size = (parameters.Nx - 1) * (parameters.Ny - 1) * parameters.Nz;
    size = parameters.Nx * parameters.Ny * parameters.Nz;
    size_Y = parameters.Nx * parameters.Ny * (parameters.k_Y_h + 1);
    u = (double *) malloc(size * sizeof(double));
    v = (double *) malloc(size * sizeof(double));
    w = (double *) malloc(size * sizeof(double));
    T = (double *) malloc(size * sizeof(double));
    p = (double *) malloc(size * sizeof(double));
    Y = (double *) malloc(size_Y * sizeof(double));
    y_0 = (double *) malloc((size + size_Y) * sizeof(double));

    // Initialize u, v, w
    power_law_initial_condition(parameters.x, parameters.y, parameters.z, u, v, w, &parameters);
    
    // Initialize T
    gaussian_temperature_initial_condition(parameters.x, parameters.y, parameters.z, T, &parameters);

    // Initialize Y
    fuel_initial_condition(parameters.x, parameters.y, parameters.z, Y, &parameters);

    // Save initial data
    // save_data("data/output/T.csv.0", parameters.x, parameters.y, parameters.z, T, parameters.Nx, parameters.Ny, parameters.Nz);
    // save_data("data/output/Y.csv.0", parameters.x, parameters.y, parameters.z, Y, parameters.Nx, parameters.Ny, parameters.k_Y_h + 1);
    
    printf("OK0\n");

    // // Create y_0 using T and Y
    // create_y_0(T, Y, y_0, &parameters);

    // printf("OK1\n");

    // // Save initial data
    // save_data("data/output/", y_0, 0, &parameters);

    // printf("OK2\n");

    // // Solve PDE
    // solve_PDE(y_0, &parameters);

    // printf("OK3\n");

    // Free memory
    free(u);
    free(v);
    free(w);
    free(T);
    free(p);
    free(Y);
    free(y_0);

    return 0;
}
