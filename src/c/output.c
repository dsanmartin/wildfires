#include "../../include/c/output.h"

void save_data(char *save_path, double *data, double *p, int n, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nz_Y = parameters->k_Y_h + 1;
    char *filename_u = (char *) malloc(100 * sizeof(char));
    char *filename_v = (char *) malloc(100 * sizeof(char));
    char *filename_w = (char *) malloc(100 * sizeof(char));
    char *filename_T = (char *) malloc(100 * sizeof(char));
    char *filename_Y = (char *) malloc(100 * sizeof(char));
    char *filename_p = (char *) malloc(100 * sizeof(char));
    FILE *output_u, *output_v, *output_w, *output_T, *output_Y, *output_p;
    int u_index = parameters->field_indexes.u;
    int v_index = parameters->field_indexes.v;
    int w_index = parameters->field_indexes.w;
    int T_index = parameters->field_indexes.T;
    int Y_index = parameters->field_indexes.Y;
    sprintf(filename_u, "%su.bin.%d", save_path, n);
    sprintf(filename_v, "%sv.bin.%d", save_path, n);
    sprintf(filename_w, "%sw.bin.%d", save_path, n);
    sprintf(filename_T, "%sT.bin.%d", save_path, n);
    sprintf(filename_Y, "%sY.bin.%d", save_path, n);
    sprintf(filename_p, "%sp.bin.%d", save_path, n);
    output_u = fopen(filename_u, "wb");
    output_v = fopen(filename_v, "wb");
    output_w = fopen(filename_w, "wb"); 
    output_T = fopen(filename_T, "wb");
    output_Y = fopen(filename_Y, "wb");
    output_p = fopen(filename_p, "wb");
    fwrite(data + u_index, sizeof(double), Nx * Ny * Nz, output_u);
    fwrite(data + v_index, sizeof(double), Nx * Ny * Nz, output_v);
    fwrite(data + w_index, sizeof(double), Nx * Ny * Nz, output_w);
    fwrite(data + T_index, sizeof(double), Nx * Ny * Nz, output_T);
    fwrite(p, sizeof(double), Nx * Ny * Nz, output_p);
    fwrite(data + Y_index, sizeof(double), Nx * Ny * Nz_Y, output_Y);
    free(filename_u);
    free(filename_v);
    free(filename_w);
    free(filename_T);
    free(filename_Y);
    free(filename_p);
    fclose(output_u);
    fclose(output_v);
    fclose(output_w);
    fclose(output_T);
    fclose(output_Y);
    fclose(output_p);
}

void save_domain(char *save_path, Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    int Nt = parameters->Nt;
    int NT = parameters->NT;
    double *x = parameters->x;
    double *y = parameters->y;
    double *z = parameters->z;
    double *t = parameters->t;
    char *filename_x = (char *) malloc(100 * sizeof(char));
    char *filename_y = (char *) malloc(100 * sizeof(char));
    char *filename_z = (char *) malloc(100 * sizeof(char));
    char *filename_t = (char *) malloc(100 * sizeof(char));
    FILE *output_x, *output_y, *output_z, *output_t;
    sprintf(filename_x, "%sx.bin", save_path);
    sprintf(filename_y, "%sy.bin", save_path);
    sprintf(filename_z, "%sz.bin", save_path);
    sprintf(filename_t, "%st.bin", save_path);
    output_x = fopen(filename_x, "wb");
    output_y = fopen(filename_y, "wb");
    output_z = fopen(filename_z, "wb");
    output_t = fopen(filename_t, "wb");
    fwrite(x, sizeof(double), Nx, output_x);
    fwrite(y, sizeof(double), Ny, output_y);
    fwrite(z, sizeof(double), Nz, output_z);
    for (int n = 0; n < Nt; n++) {
        if (n % NT == 0 || n == Nt - 1) {
            fwrite(&t[n], sizeof(double), 1, output_t);
        }
    }
    free(filename_x);
    free(filename_y);
    free(filename_z);
    free(filename_t);
    fclose(output_x);
    fclose(output_y);
    fclose(output_z);
    fclose(output_t);
}