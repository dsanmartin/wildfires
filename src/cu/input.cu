#include "../../include/cu/input.cuh"

void initial_conditions_input(double *u, double *v, double *w, double *T, double *Y, double *p, Parameters parameters) {
    // Fill the arrays with initial conditions. Data is stored in binary files with column major order.
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int Nz_Y_max = parameters.Nz_Y_max;
    char filename_u[FILENAME_SIZE];
    char filename_v[FILENAME_SIZE];
    char filename_w[FILENAME_SIZE];
    char filename_T[FILENAME_SIZE];
    char filename_Y[FILENAME_SIZE];
    char filename_p[FILENAME_SIZE];
    FILE *input_u, *input_v, *input_w, *input_T, *input_Y, *input_p;
    sprintf(filename_u, "%su.bin", parameters.input_path);
    sprintf(filename_v, "%sv.bin", parameters.input_path);
    sprintf(filename_w, "%sw.bin", parameters.input_path);
    sprintf(filename_T, "%sT.bin", parameters.input_path);
    sprintf(filename_Y, "%sY.bin", parameters.input_path);
    sprintf(filename_p, "%sp.bin", parameters.input_path);
    input_u = fopen(filename_u, "rb");
    input_v = fopen(filename_v, "rb");
    input_w = fopen(filename_w, "rb");
    input_T = fopen(filename_T, "rb");
    input_Y = fopen(filename_Y, "rb");
    input_p = fopen(filename_p, "rb");
    if (input_u == NULL || input_v == NULL || input_w == NULL || input_T == NULL || input_Y == NULL || input_p == NULL) {
        // fprintf(stderr, "Error opening initial conditions files.\n");
        // Log message
        log_message(parameters, "Error opening initial conditions files.\n");
        exit(EXIT_FAILURE);
    }
    // Read data from files
    size_t read_u = fread(u, sizeof(double), Nx * Ny * Nz, input_u);
    size_t read_v = fread(v, sizeof(double), Nx * Ny * Nz, input_v);
    size_t read_w = fread(w, sizeof(double), Nx * Ny * Nz, input_w);
    size_t read_T = fread(T, sizeof(double), Nx * Ny * Nz, input_T);
    size_t read_Y = fread(Y, sizeof(double), Nx * Ny * Nz_Y_max, input_Y);
    size_t read_p = fread(p, sizeof(double), Nx * Ny * Nz, input_p);

    if (read_u != Nx * Ny * Nz ||
        read_v != Nx * Ny * Nz ||
        read_w != Nx * Ny * Nz ||
        read_T != Nx * Ny * Nz ||
        read_Y != Nx * Ny * Nz_Y_max ||
        read_p != Nx * Ny * Nz) {
        // fprintf(stderr, "Error reading initial conditions: file size mismatch.\n");
        log_message(parameters, "Error reading initial conditions: file size mismatch.\n");
        exit(EXIT_FAILURE);
    }
    // Close files
    fclose(input_u);
    fclose(input_v);
    fclose(input_w);
    fclose(input_T);
    fclose(input_Y);
    fclose(input_p);
}