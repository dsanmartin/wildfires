#include "../../include/cu/main.cuh"

int main(int argc, char *argv[]) { 
    int device_id, number_of_SMs;
    char *parameters_file_path;
    double *u, *v, *w, *T, *Y, *p, *p_device, *y_0, *y_0_device;
    int size, size_Y;
    size_t threads_per_block, number_of_blocks;

    // Get parameters input file path
    if (argc == 2) {
        parameters_file_path = argv[1];
    } else {
        printf("Usage: %s <parameters_file_path>\n", argv[0]);
        return 1;
    }

    // Get the number of Streaming Multiprocessors
    cudaGetDevice(&device_id);
    cudaDeviceGetAttribute(&number_of_SMs, cudaDevAttrMultiProcessorCount, device_id);

    // Define kernel parameteres
    threads_per_block = THREADS;
    number_of_blocks = 32 * number_of_SMs;

    // Read parameters file
    Parameters parameters = read_parameters_file(parameters_file_path);
    log_parameters(parameters);

    // Set number of blocks and threads in parameters structure
    parameters.threads_per_block = threads_per_block;
    parameters.number_of_blocks = number_of_blocks;

    // Allocate memory for u, v, w, T, Y, p, y_0
    size = parameters.Nx * parameters.Ny * parameters.Nz;
    size_Y = parameters.Nx * parameters.Ny * parameters.Nz_Y_max;
    u = (double *) malloc(size * sizeof(double));
    v = (double *) malloc(size * sizeof(double));
    w = (double *) malloc(size * sizeof(double));
    T = (double *) malloc(size * sizeof(double));
    p = (double *) malloc(size * sizeof(double));
    Y = (double *) malloc(size_Y * sizeof(double));
    y_0 = (double *) malloc((4 * size + size_Y) * sizeof(double));
    // Allocate memory in device
    cudaMalloc(&y_0_device, (4 * size + size_Y) * sizeof(double));
    cudaMalloc(&p_device, size * sizeof(double));

    // Initialize u, v, w, T, Y, p
    log_message(parameters, "Initial conditions...");
    initial_conditions(u, v, w, T, Y, p, parameters);

    // Create y_0 using T and Y
    log_message(parameters, "Creating y_0...");
    create_y_0(u, v, w, T, Y, y_0, parameters);

    // Save initial data
    log_message(parameters, "Saving initial data...");
    save_data(y_0, p, 0, parameters);
    save_domain(parameters);
    save_topography(parameters);

    // Copy data to device
    cudaMemcpy(y_0_device, y_0, (4 * size + size_Y) * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(p_device, p, size * sizeof(double), cudaMemcpyHostToDevice);

    // Solve PDE
    log_message(parameters, "Solving PDE...");
    solve_PDE(y_0_device, p_device, parameters);

    // Free memory
    log_message(parameters, "Freeing memory...");
    free(u);
    free(v);
    free(w);
    free(T);
    free(p);
    free(Y);
    free(y_0);
    cudaFree(y_0_device);
    cudaFree(p_device);
    free_parameters(parameters);
    printf("Memory freed!\n");

    return 0;
}
