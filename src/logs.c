#include "../include/logs.h"

void log_parameters(Parameters *parameters, int to_file) {
    FILE *output;
    if (to_file)
        output = fopen("data/output/parameters.txt", "w");
    else
        output = stdout;
    fprintf(output, "Domain:\n");
    fprintf(output, "  [%lf, %lf] x [%lf, %lf] x [%lf, %lf] x [%lf, %lf]\n",
        parameters->x_min, parameters->x_max, parameters->y_min, parameters->y_max,
        parameters->z_min, parameters->z_max, parameters->t_min, parameters->t_max);
    fprintf(output, "Grid size:\n");
    fprintf(output, "  Nx: %d, Ny: %d, Nz: %d, Nt: %d\n", parameters->Nx, parameters->Ny, parameters->Nz, parameters->Nt);
    fprintf(output, "  dx: %lf, dy: %lf, dz: %lf, dt: %lf\n", parameters->dx, parameters->dy, parameters->dz, parameters->dt);
    fprintf(output, "mu: %lf\n", parameters->mu);
    if (to_file)
        fclose(output);
}