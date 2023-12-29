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
    fprintf(output, "  NT: %d\n", parameters->NT);
    fprintf(output, "Fluid parameters:\n");
    fprintf(output, "  mu: %lf, alpha: %lf\n", parameters->mu, parameters->alpha);
    fprintf(output, "Temperature initial condition:\n");
    fprintf(output, "  T_inf: %lf, T_hot: %lf\n", parameters->T_inf, parameters->T_hot);
    fprintf(output, "  T0_shape: %s\n", parameters->T0_shape);
    fprintf(output, "  T0_x_start: %lf, T0_x_end: %lf, T0_x_center: %lf, T0_length: %lf\n", 
        parameters->T0_x_start, parameters->T0_x_end, parameters->T0_x_center, parameters->T0_length);
    fprintf(output, "  T0_y_start: %lf, T0_y_end: %lf, T0_y_center: %lf, T0_width: %lf\n", 
        parameters->T0_y_start, parameters->T0_y_end, parameters->T0_y_center, parameters->T0_width);
    fprintf(output, "  T0_z_start: %lf, T0_z_end: %lf, T0_z_center: %lf, T0_height: %lf\n",
        parameters->T0_z_start, parameters->T0_z_end, parameters->T0_z_center, parameters->T0_height);
    if (to_file)
        fclose(output);
}