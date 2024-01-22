#ifndef STRUCTURES_H
#define STRUCTURES_H

typedef struct _turbulence_indexes {
    int ux;
    int uy;
    int uz;
    int vx;
    int vy;
    int vz;
    int wx;
    int wy;
    int wz;
    int Tx;
    int Ty;
    int Tz;
    int uxx;
    int uyy;
    int uzz;
    int vxx;
    int vyy;
    int vzz;
    int wxx;
    int wyy;
    int wzz;
    int Txx;
    int Tyy;
    int Tzz;
    // int uyx;
    // int uzx;
    // int uxy;
    // int uzy;
    // int uxz;
    // int uyz;
    // int vyx;
    // int vzx;
    // int vxy;
    // int vzy;
    // int vxz;
    // int vyz;
    // int wyx;
    // int wzx;
    // int wxy;
    // int wzy;
    // int wxz;
    // int wyz;
    int fw;
    int fwx;
    int fwy;
    int fwz;
} TurbulenceIndexes;

typedef struct _pressure_indexes  {
    int ux;
    int vy;
    int wz;
} PressureIndexes;

typedef struct _field_indexes {
    int u;
    int v;
    int w;
    int T;
    int Y;
} FieldIndexes;

typedef struct _parameters {
    /* Mesh nodes */
    int Nx; // Number of grid points in x
    int Ny; // Number of grid points in y
    int Nz; // Number of grid points in z
    int Nt; // Number of grid points in t
    int NT; // Number of time steps to be saved
    /* Domain */
    double x_min; // Minimum x value
    double x_max; // Maximum x value
    double y_min; // Minimum y value
    double y_max; // Maximum y value
    double z_min; // Minimum z value
    double z_max; // Maximum z value
    double t_min; // Minimum t value
    double t_max; // Maximum t value
    double dx; // Grid spacing in x
    double dy; // Grid spacing in y
    double dz; // Grid spacing in z
    double dt; // Grid spacing in t
    double *x; // x domain
    double *y; // y domain
    double *z; // z domain
    double *t; // t domain
    double *r; // Indexes in frequency space (x direction)
    double *s; // Indexes in frequency space (y direction)
    double *kx; // Domain in frequency space (x direction)
    double *ky; // Domain in frequency space (y direction)
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
    char U0_type[32];
    double u_z0;
    double d;
    double u_ast;
    double kappa;
    double u_r;
    double z_r;
    double alpha_u;
    char T0_shape[32];
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
    char Z_shape[32];
    double hill_center_x;
    double hill_center_y;
    double hill_length;
    double hill_width;
    double hill_height;
    FieldIndexes field_indexes;
    TurbulenceIndexes turbulence_indexes;
    PressureIndexes pressure_indexes;
} Parameters;

#endif