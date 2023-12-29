#ifndef STRUCTURES_H
#define STRUCTURES_H

typedef struct _parameters {
    int Nx;
    int Ny;
    int Nz;
    int Nt;
    int NT;
    double x_min;
    double x_max;
    double y_min;
    double y_max;
    double z_min;
    double z_max;
    double t_min;
    double t_max;
    double dx;
    double dy;
    double dz;
    double dt;
    double *x;
    double *y;
    double *z;
    double *t;
    double c_p;
    double c_v;
    double rho;
    double mu;
    double k;
    double alpha;
    double delta;
    double nu;
    double T_inf;
    double g;
    double C_s;
    double Pr;
    double C_D;
    double a_v;
    double T_pc;
    double H_R;
    double A;
    double E_A;
    double T_a;
    double h;
    double Y_D;
    double Y_f;
    double T_hot;
    double Y_h;
    int k_Y_h;
    double T_min;
    double T_max;
    double Y_min;
    double Y_max;
    char *U0_type;
    double u_z0;
    double d;
    double u_ast;
    double kappa;
    double u_r;
    double z_r;
    double alpha_u;
    char *T0_shape;
    double T0_x_start;
    double T0_x_end;
    double T0_x_center;
    double T0_length;
    double T0_y_start;
    double T0_y_end;
    double T0_y_center;
    double T0_width;
    double T0_z_start;
    double T0_z_end;
    double T0_z_center;
    double T0_height;
    char *Z_shape;
    double hill_center_x;
    double hill_center_y;
    double hill_length;
    double hill_width;
    double hill_height;
} Parameters;

#endif