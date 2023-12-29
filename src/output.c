#include "../include/output.h"

void save_data(char *path, double *x, double *y, double *z, double *f, int Nx, int Ny, int Nz) {
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
