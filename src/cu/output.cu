#include "../../include/cu/output.cuh"

void save_data(double *data, double *p, int n, Parameters parameters) {
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
    FILE *output_u, *output_v, *output_w, *output_T, *output_Y, *output_p;
    int u_index = parameters.field_indexes.u;
    int v_index = parameters.field_indexes.v;
    int w_index = parameters.field_indexes.w;
    int T_index = parameters.field_indexes.T;
    int Y_index = parameters.field_indexes.Y;
    sprintf(filename_u, "%su.bin.%d", parameters.save_path, n);
    sprintf(filename_v, "%sv.bin.%d", parameters.save_path, n);
    sprintf(filename_w, "%sw.bin.%d", parameters.save_path, n);
    sprintf(filename_T, "%sT.bin.%d", parameters.save_path, n);
    sprintf(filename_Y, "%sY.bin.%d", parameters.save_path, n);
    sprintf(filename_p, "%sp.bin.%d", parameters.save_path, n);
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
    fwrite(data + Y_index, sizeof(double), Nx * Ny * Nz_Y_max, output_Y);
    fclose(output_u);
    fclose(output_v);
    fclose(output_w);
    fclose(output_T);
    fclose(output_Y);
    fclose(output_p);
}

void save_domain(Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    int Nt = parameters.Nt;
    int NT = parameters.NT;
    double *x = parameters.x;
    double *y = parameters.y;
    double *z = parameters.z;
    double *t = parameters.t;
    char filename_x[FILENAME_SIZE];
    char filename_y[FILENAME_SIZE];
    char filename_z[FILENAME_SIZE];
    char filename_t[FILENAME_SIZE];
    FILE *output_x, *output_y, *output_z, *output_t;
    sprintf(filename_x, "%sx.bin", parameters.save_path);
    sprintf(filename_y, "%sy.bin", parameters.save_path);
    sprintf(filename_z, "%sz.bin", parameters.save_path);
    sprintf(filename_t, "%st.bin", parameters.save_path);
    output_x = fopen(filename_x, "wb");
    output_y = fopen(filename_y, "wb");
    output_z = fopen(filename_z, "wb");
    output_t = fopen(filename_t, "wb");
    fwrite(x, sizeof(double), Nx, output_x);
    fwrite(y, sizeof(double), Ny, output_y);
    fwrite(z, sizeof(double), Nz, output_z);
    for (int n = 0; n <= Nt; n++) {
        if (n % NT == 0 || n == Nt) {
            fwrite(&t[n], sizeof(double), 1, output_t);
        }
    }
    fclose(output_x);
    fclose(output_y);
    fclose(output_z);
    fclose(output_t);
}

void save_topography(Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    // int *Nz_Y = parameters->Nz_Y;
    double *topography = parameters.topography;    
    char filename_topography[FILENAME_SIZE];
    // char filename_Y_mask[FILENAME_SIZE];
    FILE *output_topography;//, *output_Y_mask;
    sprintf(filename_topography, "%stopography.bin", parameters.save_path);
    // sprintf(filename_Y_mask, "%sY_mask.bin", parameters->save_path);
    output_topography = fopen(filename_topography, "wb");
    // output_Y_mask = fopen(filename_Y_mask, "wb");
    fwrite(topography, sizeof(double), Nx * Ny, output_topography);
    // fwrite(Nz_Y, sizeof(int), Nx * Ny, output_Y_mask);
    fclose(output_topography);
    // fclose(output_Y_mask);
}