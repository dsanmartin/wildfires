#include "../include/poisson.h"

void thomas_algorithm(double complex *a, double complex *b, double complex *c, double complex *d, double complex *x, double complex *l, double complex *u, double complex *y, int N) {
    // Create l and u
    l[0] = 0;
    u[0] = b[0];
    for (int i = 1; i < N; i++) {
        l[i] = a[i] / u[i - 1];
        u[i] = b[i] - l[i] * c[i - 1];
    }

    // Solve Ly = d
    y[0] = d[0];
    for (int i = 1; i < N; i++) {
        y[i] = d[i] - l[i] * y[i - 1];
    }

    // Solve Ux = y
    x[N - 1] = y[N - 1] / u[N - 1];
    for (int i = N - 2; i >= 0; i--) {
        x[i] = (y[i] - c[i] * x[i + 1]) / u[i];
    }

}