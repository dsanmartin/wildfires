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
    log_parameters(&parameters, 1);

    // Allocate memory for x, y, z, u, v, w
    // size = (parameters.Nx - 1) * (parameters.Ny - 1) * parameters.Nz;
    size = parameters.Nx * parameters.Ny * parameters.Nz;
    size_Y = parameters.Nx * parameters.Ny * (parameters.Nz_Y);
    u = (double *) malloc(size * sizeof(double));
    v = (double *) malloc(size * sizeof(double));
    w = (double *) malloc(size * sizeof(double));
    T = (double *) malloc(size * sizeof(double));
    p = (double *) malloc(size * sizeof(double));
    Y = (double *) malloc(size_Y * sizeof(double));
    y_0 = (double *) malloc((4 * size + size_Y) * sizeof(double));

    // Initialize u, v, w, T, Y, p
    printf("Initial conditions...\n");
    initial_conditions(parameters.x, parameters.y, parameters.z, u, v, w, T, Y, p, &parameters);
    printf("Initial conditions... OK!\n");

    // Create y_0 using T and Y
    printf("Create y_0...\n");
    create_y_0(u, v, w, T, Y, y_0, &parameters);
    printf("Create y_0... OK!\n");

    // Save initial data
    printf("Saving initial data...\n");
    save_data("data/output/", y_0, p, 0, &parameters);
    save_time("data/output/", parameters.t, parameters.Nt, parameters.NT);
    printf("Initial data saved!\n");

    // Solve PDE
    printf("Solving PDE...\n");
    solve_PDE(y_0, p, &parameters);
    printf("PDE solved!\n");

    // Free memory
    printf("Freeing memory...\n");
    free(u);
    free(v);
    free(w);
    free(T);
    free(p);
    free(Y);
    // printf("Memory freed 1!\n");
    // free(parameters.x);
    // free(parameters.y);
    // free(parameters.z);
    // free(parameters.t);
    // free(parameters.r);
    // free(parameters.s);
    // free(parameters.kx);
    // free(parameters.ky);
    printf("Memory freed!\n");

    return 0;
}
