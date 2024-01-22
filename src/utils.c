#include "../include/utils.h"

void axpy(double *y, double *x, double a, int size) {
    for (int i = 0; i < size; i++) {
        y[i] += a * x[i];
    }
}

void copy(double *destination, double *source, int size) {
    for (int i = 0; i < size; i++) {
        destination[i] = source[i];
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