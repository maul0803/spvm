#include "spmv.h"
#include <stdio.h>
#include <stdlib.h>
typedef struct {
    unsigned int *row_offsets;
    unsigned int *colulmn_indices;
    double *values;
    unsigned int n;
} CSR;

CSR convert_dense_to_CSR(const unsigned int n, const double mat[]){
    CSR csr;
    csr.n = n;
    // Count the number of values in the matrix to create the arrays
    unsigned int number_of_not_zeros = 0;
    for (unsigned int i = 0; i < n * n; i++) {
        if (mat[i] != 0.0) {
            number_of_not_zeros++;
        }
    }
    csr.row_offsets = (unsigned int *)malloc((n + 1) * sizeof(unsigned int)); // n+1 pour le dernier offset
    csr.colulmn_indices = (unsigned int *)malloc(number_of_not_zeros * sizeof(unsigned int));
    csr.values = (double *)malloc(number_of_not_zeros * sizeof(double));

    for (unsigned int i = 0; i < n; i++){
        for 
    }
}
int my_dense(const unsigned int n, const double mat[], double vec[], double result[])
{
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
            result[i] += mat[i * n + j] * vec[j];
        }
    }
    return 0;
}
