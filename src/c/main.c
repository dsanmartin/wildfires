#include "../../include/c/main.h"

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
    size_Y = parameters.Nx * parameters.Ny * parameters.Nz_Y;
    u = (double *) malloc(size * sizeof(double));
    v = (double *) malloc(size * sizeof(double));
    w = (double *) malloc(size * sizeof(double));
    T = (double *) malloc(size * sizeof(double));
    p = (double *) malloc(size * sizeof(double));
    Y = (double *) malloc(size_Y * sizeof(double));
    y_0 = (double *) malloc((4 * size + size_Y) * sizeof(double));

    // Initialize u, v, w, T, Y, p
    // printf("Initial conditions...\n");
    log_message(&parameters, "Initial conditions...");
    initial_conditions(parameters.x, parameters.y, parameters.z, u, v, w, T, Y, p, &parameters);
    // printf("Initial conditions... OK!\n");

    // Create y_0 using T and Y
    // printf("Create y_0...\n");
    log_message(&parameters, "Creating y_0...");
    create_y_0(u, v, w, T, Y, y_0, &parameters);
    // printf("Create y_0... OK!\n");

    // Save initial data
    // printf("Saving initial data...\n");
    log_message(&parameters, "Saving initial data...");
    save_data(parameters.save_path, y_0, p, 0, &parameters);
    save_domain(parameters.save_path, &parameters);
    // printf("Initial data saved!\n");

    // Solve PDE
    // printf("Solving PDE...\n");
    log_message(&parameters, "Solving PDE...");
    solve_PDE(y_0, p, &parameters);
    // printf("PDE solved!\n");

    // Free memory
    // printf("Freeing memory...\n");
    log_message(&parameters, "Freeing memory...");
    free(u);
    free(v);
    free(w);
    free(T);
    free(p);
    free(Y);
    free_parameters(&parameters);
    // printf("Memory freed!\n");

    return 0;
}
