#include "../include/output.h"

void save_scalar(char *path, double *x, double *y, double *z, double *f, int Nx, int Ny, int Nz) {
    FILE *output = fopen(path, "w");
    fprintf(output, "x, y, z, f\n");
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                fprintf(output, "%lf, %lf, %lf, %lf\n", x[i], y[j], z[k], f[idx(i, j, k, Nx, Ny, Nz)]);
            }
        }
    }
    fclose(output);
}

void save_data(char *save_path, double *data, int n, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->k_Y_h + 1;
    double *x = parameters->x;
    double *y = parameters->y;
    double *z = parameters->z;
    double T, Y;
    char *filename_T = (char *) malloc(100 * sizeof(char));
    char *filename_Y = (char *) malloc(100 * sizeof(char));
    FILE *output_T, *output_Y;
    int T_index = 0;
    int Y_index = Nx * Ny * Nz;
    sprintf(filename_T, "%sT.csv.%d", save_path, n);
    sprintf(filename_Y, "%sY.csv.%d", save_path, n); 
    output_T = fopen(filename_T, "w");
    output_Y = fopen(filename_Y, "w");
    fprintf(output_T, "x, y, z, T\n");
    fprintf(output_Y, "x, y, z, Y\n");
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                T = data[T_index + idx(i, j, k, Nx, Ny, Nz)];
                fprintf(output_T, "%lf, %lf, %lf, %lf\n", x[i], y[j], z[k], T);
                if (k < Nz_Y) {
                    Y = data[Y_index + idx(i, j, k, Nx, Ny, Nz_Y)];
                    fprintf(output_Y, "%lf, %lf, %lf, %lf\n", x[i], y[j], z[k], Y);
                }
            }
        }
    }
    free(filename_T);
    free(filename_Y);
    fclose(output_T);
    fclose(output_Y);
}
