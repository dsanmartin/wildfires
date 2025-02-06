#include "../../include/c/ibm.h"

void flat_terrain(Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            parameters->topography[IDX(i, j, 0, Nx, Ny, 1)] = 0.0;
        }
    }
}

void simple_hill(Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    double *x = parameters->x;
    double *y = parameters->y;
    double hill_center_x = parameters->hill_center_x;
    double hill_center_y = parameters->hill_center_y;
    double hill_length = parameters->hill_length;
    double hill_width = parameters->hill_width;
    double hill_height = parameters->hill_height;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            parameters->topography[IDX(i, j, 0, Nx, Ny, 1)] = hill_height * gaussian(x[i], y[j], 0, hill_center_x, hill_center_y,  0, hill_length, hill_width, 1);
        }
    }
}

void cube(Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    double *x = parameters->x;
    double *y = parameters->y;
    double *z = parameters->z;
    double cube_center_x = 100;
    double cube_center_y = 100;
    double cube_center_z = 0;
    double cube_length = 20;
    double cube_width = 20;
    double cube_height = 3;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                // Modify cut nodes considering a cube in the domain
                if (cube_center_x - cube_length/2 <= x[i] && x[i] <= cube_center_x + cube_length/2 &&
                    cube_center_y - cube_width/2 <= y[j] && y[j] <= cube_center_y + cube_width/2 &&
                    cube_center_z <= z[k] && z[k] <= cube_center_z + cube_height) {
                    parameters->cut_nodes[IDX(i, j, 0, Nx, Ny, 1)] = k;
                }
            }
        }
    }
}

void ibm_parameters(Parameters *parameters) {
    int Nx = parameters->Nx;
    int Ny = parameters->Ny;
    int Nz = parameters->Nz;
    double *z = parameters->z;
    // double dz = parameters->dz;
    double Y_h = parameters->Y_h;
    double topo;
    double dz;
    int Nz_Y_max = 0;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            // Get topography data
            topo = parameters->topography[IDX(i, j, 0, Nx, Ny, 1)];                     
            // Initialize the cut nodes and Nz_Y
            parameters->cut_nodes[IDX(i, j, 0, Nx, Ny, 1)] = 0;
            parameters->Nz_Y[IDX(i, j, 0, Nx, Ny, 1)] = 0;
            for (int k = 0; k < Nz; k++) {      
                if (k < Nz - 1)
                    dz = z[k + 1] - z[k];
                else
                    dz = z[k] - z[k - 1];          
                // Find the cut nodes indexes
                if ((topo - dz) <= z[k] && z[k] <= topo) {
                    parameters->cut_nodes[IDX(i, j, 0, Nx, Ny, 1)] = k;
                }         
                // Distance from the z coordinate to the topography. Negative does not matter
                parameters->z_ibm[IDX(i, j, k, Nx, Ny, Nz)] = z[k] - topo;       
                // Compute Nz for Y variable to avoid storing a lot of zeros
                if (z[k] <= topo + Y_h) { // If the node is inside the hill + fuel height
                    parameters->Nz_Y[IDX(i, j, 0, Nx, Ny, 1)] = k + 1; // Store the index of the last node plus one to include the last node
                }
                // printf("Nz_Y[%d, %d] = %d\n", i, j, parameters->Nz_Y[IDX(i, j, 0, Nx, Ny, 1)]);                
            }                   
            Nz_Y_max = MAX(Nz_Y_max, parameters->Nz_Y[IDX(i, j, 0, Nx, Ny, 1)]); // Find the maximum number of nodes for Y variable  
        }
    }
    parameters->Nz_Y_max = Nz_Y_max;
}

/*
void flat_terrain(Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            parameters.topography[IDX(i, j, 0, Nx, Ny, 1)] = 0.0;
        }
    }
}

void simple_hill(Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    double *x = parameters.x;
    double *y = parameters.y;
    double hill_center_x = parameters.hill_center_x;
    double hill_center_y = parameters.hill_center_y;
    double hill_length = parameters.hill_length;
    double hill_width = parameters.hill_width;
    double hill_height = parameters.hill_height;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            parameters.topography[IDX(i, j, 0, Nx, Ny, 1)] = hill_height * gaussian(x[i], y[j], 0, hill_center_x, hill_center_y,  0, hill_length, hill_width, 1);
        }
    }
}

void ibm_parameters(Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    double *z = parameters.z;
    double dz = parameters.dz;
    double Y_h = parameters.Y_h;
    double topo;
    int Nz_Y_max = 0;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            // Get topography data
            topo = parameters.topography[IDX(i, j, 0, Nx, Ny, 1)];                     
            // Initialize the cut nodes and Nz_Y
            parameters.cut_nodes[IDX(i, j, 0, Nx, Ny, 1)] = 0;
            parameters.Nz_Y[IDX(i, j, 0, Nx, Ny, 1)] = 0;
            for (int k = 0; k < Nz; k++) {                
                // Find the cut nodes indexes
                if ((topo - dz) <= z[k] && z[k] <= topo) {
                    parameters.cut_nodes[IDX(i, j, 0, Nx, Ny, 1)] = k;
                }         
                // Distance from the z coordinate to the topography. Negative does not matter
                parameters.z_ibm[IDX(i, j, k, Nx, Ny, Nz)] = z[k] - topo;       
                // Compute Nz for Y variable to avoid storing a lot of zeros
                if (z[k] <= topo + Y_h) { // If the node is inside the hill + fuel height
                    parameters.Nz_Y[IDX(i, j, 0, Nx, Ny, 1)] = k + 1; // Store the index of the last node plus one to include the last node
                }
                // printf("Nz_Y[%d, %d] = %d\n", i, j, parameters.Nz_Y[IDX(i, j, 0, Nx, Ny, 1)]);                
            }                   
            Nz_Y_max = MAX(Nz_Y_max, parameters.Nz_Y[IDX(i, j, 0, Nx, Ny, 1)]); // Find the maximum number of nodes for Y variable  
        }
    }
    parameters.Nz_Y_max = Nz_Y_max;
    printf("Nz_Y_max = %d\n", parameters.Nz_Y_max);
}
*/

// void topography(Parameters parameters) {
//     if (strcmp(parameters.topo_shape, "simple hill") == 0)
//         simple_hill(parameters);
//     else if (strcmp(parameters.topo_shape, "flat") == 0)
//         flat_terrain(parameters);
//     ibm_parameters(parameters);
// }

/*
void ibm_nodes(Parameters parameters) {
    int Nx = parameters.Nx;
    int Ny = parameters.Ny;
    int Nz = parameters.Nz;
    double dz = parameters.dz;
    double topo;
    double *z = parameters.z;
    double *topography = parameters.topography;
    int max_k = 0, n_dead_nodes = 0;
    parameters.n_cut_nodes = Nx * Ny;
    // Create the cut nodes
    CutNodes cut_nodes[parameters.n_cut_nodes];
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                topo = topography[IDX(i, j, k, Nx, Ny, Nz)];
                // Check if the node is a dead node
                if (z[k] < (topo - dz)) {
                    n_dead_nodes++;
                    max_k = MAX(max_k, k);
                }
                // Check if the node is a cut node
                else {
                    if ((topo - dz) <= z[k] && z[k] <= topo) {
                        cut_nodes[parameters.n_dead_nodes].i = i;
                        cut_nodes[parameters.n_dead_nodes].j = j;
                        cut_nodes[parameters.n_dead_nodes].k = k;
                    }   
                }             
            }
        }
    }
    parameters.cut_nodes = cut_nodes;
    // Create the dead nodes
    parameters.n_dead_nodes = n_dead_nodes;
    DeadNodes dead_nodes[parameters.n_dead_nodes];
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < max_k; k++) {
                topo = topography[IDX(i, j, k, Nx, Ny, Nz)];
                if (z[k] < (topo - dz)) {
                    dead_nodes[parameters.n_dead_nodes].i = i;
                    dead_nodes[parameters.n_dead_nodes].j = j;
                    dead_nodes[parameters.n_dead_nodes].k = k;
                }
            }
        }
    }
}

*/
