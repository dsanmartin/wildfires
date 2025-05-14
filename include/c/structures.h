/**
 * @file structures.h
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Structures used in the wildfire simulation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <stdio.h>

/**
 * @brief Structure to store the indexes of the fields.
 * 
 */
typedef struct _field_indexes {
    int u;
    int v;
    int w;
    int T;
    int Y;
} FieldIndexes;

/**
 * @brief Structure to store the indexes of the gradients for pressure field.
 * 
 */
typedef struct _pressure_indexes  {
    int ux;
    int vy;
    int wz;
} PressureIndexes;

/**
 * @brief Structure to store the indexes of the turbulence fields.
 * 
 */
typedef struct _turbulence_indexes {
    // int ux;
    // int uy;
    // int uz;
    // int vx;
    // int vy;
    // int vz;
    // int wx;
    // int wy;
    // int wz;
    int rho;
    int Tx;
    int Ty;
    int Tz;
    // int uxx;
    // int uyy;
    // int uzz;
    // int vxx;
    // int vyy;
    // int vzz;
    // int wxx;
    // int wyy;
    // int wzz;
    int Txx;
    int Tyy;
    int Tzz;
    // int fw;
    int mu_sgs;
    int S_11;
    int S_12;
    int S_13;
    int S_21;
    int S_22;
    int S_23;
    int S_31;
    int S_32;
    int S_33;
    // int fwx;
    // int fwy;
    // int fwz;
} TurbulenceIndexes;

/**
 * @brief Structure to store the log files.
 * 
 */
typedef struct _log_files {
    FILE *parameters;
    FILE *log;
} LogFiles;

// typedef struct _dead_nodes {
//     int i; int j; int k;
// } DeadNodes;

// typedef struct _cut_nodes {
//     int i; int j; int k;
// } CutNodes;

/**
 * @brief Structure to store the parameters of the simulation.
 * 
 */
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
    double rho_inf;    
    double mu;
    double k;
    double alpha;
    double delta;
    double nu;
    double T_inf;
    double g;
    double C_s;
    double Pr;
    double C_d;
    double alpha_s;
    double sigma_s;
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
    int Nz_Y_max; // Maximum number of grid points in z for Y variable
    int *Nz_Y; // Number of grid points in z for Y variable
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
    char topo_shape[32]; // Topography shape (simple_hill or flat)
    double *topography; // Topography field
    double hill_center_x; 
    double hill_center_y; 
    double hill_length; 
    double hill_width; 
    double hill_height; 
    double t_source; 
    // IBM dead nodes and cut nodes
    // int n_dead_nodes;
    // int n_cut_nodes;
    // double *dead_nodes;
    // double *cut_nodes;
    int *cut_nodes;
    double *z_ibm;
    double p_top;
    double u_dead_nodes;
    double v_dead_nodes;
    double w_dead_nodes;
    double T_dead_nodes;
    double Y_dead_nodes;
    // Pressure solver parameters
    int variable_density;
    double pressure_solver_tol;
    int pressure_solver_iter;
    int pressure_solver_log;
    char method[32];
    char save_path[64];
    char sim_id[32];
    int n_threads;
    int z_domain;
    double k_nu_grid;
    int k_transition;
    double z_transition;
    size_t threads_per_block;
    size_t number_of_blocks;
    FieldIndexes field_indexes;
    TurbulenceIndexes turbulence_indexes;
    PressureIndexes pressure_indexes;
    // DeadNodes *dead_nodes;
    // CutNodes *cut_nodes;
    LogFiles log_files;
} Parameters;

#endif