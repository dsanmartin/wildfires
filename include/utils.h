#ifndef UTILS_H
#define UTILS_H

#define M_PI 3.14159265358979323846
#define idxR(i, j, k, Nx, Ny, Nz) (k) + (Nz) * ((j) + (Ny) * (i)) // Indexing macro row-major
#define idxC(i, j, k, Nx, Ny, Nz) (i) + (Nx) * ((j) + (Ny) * (k)) // Indexing macro column-major 
#define IDX(i, j, k, Nx, Ny, Nz) idxR(i, j, k, Nx, Ny, Nz) // Default indexing macro
#define FFTWIDX(i, j, k, Nx, Ny, Nz) (j) + (Ny) * (i) + (Nx) * (Ny) * (k) // Indexing macro for FFTW
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#include <stdio.h>
#include <time.h>

void axpy(double *y, double *x, double a, int size);
void caxpy(double *c, double *x, double *y, double a, int size);
void copy(double *destination, double *source, int size);
void copy_slice(double *destination, double *source, int Nx, int Ny, int Nz, int slice);
void fft_freq(double *f, int N, double d);
void generate_current_datetime_string(char *datetime_str, size_t max_size);

#endif