#ifndef UTILS_H
#define UTILS_H

#define idxR(i, j, k, Nx, Ny, Nz) (k) + (Nz) * ((j) + (Ny) * (i)) // Indexing macro row-major
#define idxC(i, j, k, Nx, Ny, Nz) (i) + (Nx) * ((j) + (Ny) * (k)) // Indexing macro column-major 
#define IDX(i, j, k, Nx, Ny, Nz) idxR(i, j, k, Nx, Ny, Nz) // Default indexing macro
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

void axpy(double *y, double *x, double a, int size);
void copy(double *destination, double *source, int size);

#endif