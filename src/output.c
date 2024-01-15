#include "../include/output.h"

void save_scalar(char *path, double *x, double *y, double *z, double *f, int Nx, int Ny, int Nz) {
    FILE *output = fopen(path, "w");
    fprintf(output, "x, y, z, f\n");
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                fprintf(output, "%lf, %lf, %lf, %lf\n", x[i], y[j], z[k], f[IDX(i, j, k, Nx, Ny, Nz)]);
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
    double u, v, w, T, Y;
    char *filename_u = (char *) malloc(100 * sizeof(char));
    char *filename_v = (char *) malloc(100 * sizeof(char));
    char *filename_w = (char *) malloc(100 * sizeof(char));
    char *filename_T = (char *) malloc(100 * sizeof(char));
    char *filename_Y = (char *) malloc(100 * sizeof(char));
    FILE *output_u, *output_v, *output_w, *output_T, *output_Y;
    int u_index = 0;
    int v_index = Nx * Ny * Nz;
    int w_index = 2 * Nx * Ny * Nz;
    int T_index = 3 * Nx * Ny * Nz;
    int Y_index = 4 * Nx * Ny * Nz;
    sprintf(filename_u, "%su.csv.%d", save_path, n);
    sprintf(filename_v, "%sv.csv.%d", save_path, n);
    sprintf(filename_w, "%sw.csv.%d", save_path, n);
    sprintf(filename_T, "%sT.csv.%d", save_path, n);
    sprintf(filename_Y, "%sY.csv.%d", save_path, n);
    output_u = fopen(filename_u, "w");
    output_v = fopen(filename_v, "w");
    output_w = fopen(filename_w, "w"); 
    output_T = fopen(filename_T, "w");
    output_Y = fopen(filename_Y, "w");
    fprintf(output_u, "x, y, z, u\n");
    fprintf(output_v, "x, y, z, v\n");
    fprintf(output_w, "x, y, z, w\n");
    fprintf(output_T, "x, y, z, T\n");
    fprintf(output_Y, "x, y, z, Y\n");
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                u = data[u_index + IDX(i, j, k, Nx, Ny, Nz)];
                v = data[v_index + IDX(i, j, k, Nx, Ny, Nz)];
                w = data[w_index + IDX(i, j, k, Nx, Ny, Nz)];
                T = data[T_index + IDX(i, j, k, Nx, Ny, Nz)];
                fprintf(output_u, "%lf, %lf, %lf, %lf\n", x[i], y[j], z[k], u);
                fprintf(output_v, "%lf, %lf, %lf, %lf\n", x[i], y[j], z[k], v);
                fprintf(output_w, "%lf, %lf, %lf, %lf\n", x[i], y[j], z[k], w);
                fprintf(output_T, "%lf, %lf, %lf, %lf\n", x[i], y[j], z[k], T);
                if (k < Nz_Y) {
                    Y = data[Y_index + IDX(i, j, k, Nx, Ny, Nz_Y)];
                    fprintf(output_Y, "%lf, %lf, %lf, %lf\n", x[i], y[j], z[k], Y);
                }
            }
        }
    }
    free(filename_u);
    free(filename_v);
    free(filename_w);
    free(filename_T);
    free(filename_Y);
    fclose(output_u);
    fclose(output_v);
    fclose(output_w);
    fclose(output_T);
    fclose(output_Y);
}
