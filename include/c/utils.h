/**
 * @file utils.h
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Utility functions for the wildfire simulation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef UTILS_H
#define UTILS_H

#define idxR(i, j, k, Nx, Ny, Nz) (k) + (Nz) * ((j) + (Ny) * (i)) // Indexing macro row-major
#define idxC(i, j, k, Nx, Ny, Nz) (i) + (Nx) * ((j) + (Ny) * (k)) // Indexing macro column-major 
#define IDX(i, j, k, Nx, Ny, Nz) idxR(i, j, k, Nx, Ny, Nz) // Default indexing macro
#define FFTWIDX(i, j, k, Nx, Ny, Nz) (j) + (Ny) * (i) + (Nx) * (Ny) * (k) // Indexing macro for FFTW
// #define FFTWIDX(i, j, k, Nx, Ny, Nz) (k) + (Nz) * ((j) + (Ny) * (i))
#define MAX(x, y) (((x) > (y)) ? (x) : (y)) // Maximum of two numbers
#define MIN(x, y) (((x) < (y)) ? (x) : (y)) // Minimum of two numbers

#include <stdio.h>
#include <time.h>
#include "constants.h"

/**
 * Performs the operation c = a * x + y, where c, x, and y are arrays of doubles.
 *
 * @param c     The output array where the result will be stored.
 * @param x     The input array x.
 * @param y     The input array y.
 * @param a     The scalar value a.
 * @param size  The size of the arrays c, x, and y.
 */
void caxpy(double *c, double *x, double *y, double a, int size);

/**
 * @brief Copies the elements from the source array to the destination array.
 *
 * This function copies the elements from the source array to the destination array.
 *
 * @param destination Pointer to the destination array.
 * @param source Pointer to the source array.
 * @param size The number of elements to copy.
 */
void copy(double *destination, double *source, int size);

/**
 * @brief Calculates the frequencies of the Fast Fourier Transform (FFT).
 *
 * This function calculates the frequencies of the Fast Fourier Transform (FFT)
 * for a given number of points and a specified frequency resolution.
 * f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
 * f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd
 *
 * @param f Pointer to the array where the calculated frequencies will be stored.
 * @param N The number of points for which the frequencies will be calculated.
 * @param d The frequency resolution.
 */
void fft_freq(double *f, int N, double d);

/**
 * @brief Generates a string representation of the current date and time.
 *
 * This function generates a string representation of the current date and time
 * and stores it in the provided buffer.
 *
 * @param datetime_str The buffer to store the generated date and time string.
 * @param max_size The maximum size of the buffer.
 */
void generate_current_datetime_string(char *datetime_str, size_t max_size);

/**
 * @brief Formats the given number of seconds into a human-friendly time format.
 *
 * This function takes the number of seconds as input and formats it into a string
 * representing a human-readable time format. The formatted time is stored in the
 * provided `formatted_time` buffer.
 *
 * @param seconds The number of seconds to format.
 * @param formatted_time The buffer to store the formatted time.
 */
void format_seconds(double seconds, char *formatted_time);

#endif