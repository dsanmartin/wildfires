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