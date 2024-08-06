/**
 * @file utils.cu
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Utility functions for the wildfire simulation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "../../include/cu/utils.cuh"

__global__ 
void caxpy(double *c, double *x, double *y, double a, int size) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;
    for (int i = idx; i < size; i += stride) {
        c[i] = y[i] + a * x[i];
    }
}

__global__
void copy(double *destination, double *source, int size) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = gridDim.x * blockDim.x;
    for (int i = idx; i < size; i += stride) {
        destination[i] = source[i];
    }
}

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

void format_seconds(double seconds, char *formatted_time) {
    int hours = seconds / 3600;
    int minutes = ((int) seconds % 3600) / 60;
    int secs = ((int) seconds % 3600) % 60;
    sprintf(formatted_time, "%02d:%02d:%02d", hours, minutes, secs);
}