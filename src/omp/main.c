#include "../../include/omp/main.h"

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

    // Read parameters file
    Parameters parameters = read_parameters_file(parameters_file_path);
    log_parameters(&parameters);

    // Allocate memory for u, v, w, T, Y, p, y_0
    size = parameters.Nx * parameters.Ny * parameters.Nz;
    size_Y = parameters.Nx * parameters.Ny * parameters.Nz_Y;
    u = (double *) malloc(size * sizeof(double));
    v = (double *) malloc(size * sizeof(double));
    w = (double *) malloc(size * sizeof(double));
    T = (double *) malloc(size * sizeof(double));
    p = (double *) malloc(size * sizeof(double));
    Y = (double *) malloc(size_Y * sizeof(double));
    y_0 = (double *) malloc((4 * size + size_Y) * sizeof(double));

    // Initialize u, v, w, T, Y, p
    log_message(&parameters, "Initial conditions...");
    initial_conditions(u, v, w, T, Y, p, &parameters);

    // Create y_0 using T and Y
    log_message(&parameters, "Creating y_0...");
    create_y_0(u, v, w, T, Y, y_0, &parameters);

    // Save initial data
    log_message(&parameters, "Saving initial data...");
    save_data(y_0, p, 0, &parameters);
    save_domain(&parameters);

    // Solve PDE
    log_message(&parameters, "Solving PDE...");
    solve_PDE(y_0, p, &parameters);

    // Free memory
    log_message(&parameters, "Freeing memory...");
    free(u);
    free(v);
    free(w);
    free(T);
    free(p);
    free(Y);
    free_parameters(&parameters);
    printf("Memory freed!\n");

    return 0;
}
