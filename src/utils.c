#include "../include/utils.h"

void axpy(double *y, double *x, double a, int size) {
    for (int i = 0; i < size; i++) {
        y[i] += a * x[i];
    }
}

void caxpy(double *c, double *x, double *y, double a, int size) {
    for (int i = 0; i < size; i++) {
        c[i] = y[i] + a * x[i];
    }
}

void copy(double *destination, double *source, int size) {
    for (int i = 0; i < size; i++) {
        destination[i] = source[i];
    }
}

void copy_slice(double *destination, double *source, int Nx, int Ny, int Nz, int slice) {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            destination[IDX(i, j, 0, Nx, Ny, 1)] = source[IDX(i, j, slice, Nx, Ny, Nz)];
        }
    }
}



/*
DFT sample frequencies
f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd
*/
void fft_freq(double *f, int N, double d) {
    int n;
    if (N % 2 == 0) {
        n = N / 2;
        for (int i = 0; i < n; i++) {
            f[i] = i / (N * d);
        }
        for (int i = n; i < N; i++) {
            f[i] = (i - N) / (N * d);
        }
    } else {
        n = (N - 1) / 2;
        for (int i = 0; i <= n; i++) {
            f[i] = i / (N * d);
        }
        for (int i = n + 1; i < N; i++) {
            f[i] = (i - N) / (N * d);
        }
    }
}


void generate_current_datetime_string(char *datetime_str, size_t max_size) {
    time_t t = time(NULL);
    struct tm *tm = localtime(&t);
    if (tm != NULL) {
        // Format: YYYYMMDDHHMMSS
        strftime(datetime_str, max_size, "%Y%m%d%H%M%S", tm);
    }
}