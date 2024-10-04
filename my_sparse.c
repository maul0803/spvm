#include "spmv.h"
#include "my_sparse.h"
#include <stdio.h>
#include <stdlib.h>

CSR convert_dense_to_CSR(const unsigned int n, const double mat[]){
    CSR csr;
    // Count the number of values in the matrix to create the arrays
    unsigned int size_row_offsets = n+1;
    unsigned int size_column_indices_values = 0;
    for (unsigned int i = 0; i < n * n; i++) {
        if (mat[i] != 0.0) {
            size_column_indices_values++;
        }
    }
    csr.row_offsets = (unsigned int *)malloc((size_row_offsets) * sizeof(unsigned int)); // n+1 for the last offset
    csr.column_indices = (unsigned int *)malloc(size_column_indices_values * sizeof(unsigned int));
    csr.values = (double *)malloc(size_column_indices_values * sizeof(double));
    csr.size_row_offsets = size_row_offsets;
    csr.size_column_indices_values = size_column_indices_values;
    unsigned int buffer = 0;
    for (unsigned int i = 0; i < n; i++){
        csr.row_offsets[i] = buffer;
        for (unsigned int j = 0; j < n ; j++){
            if (mat[i * n + j] !=0){
                csr.column_indices[buffer] = j;
                csr.values[buffer] = mat[i * n + j];
                buffer ++;
            }
        }
    }
    csr.row_offsets[n] = buffer;
    return csr;
}

double* convert_CSR_to_dense(const CSR csr, unsigned int n) {
    double *mat = (double *) calloc(n * n, sizeof(double));
    // Fill the matrix with 0.0.
    for (unsigned  int i = 0; i < n; i ++){
        for (unsigned int j = 0; j < n; j++){
            mat[i * n + j] = 0.0;
        }
    }
    // Go through all rows
    for (unsigned int i = 0; i < n; i++) {
        // Go through all columns of each row
        for (unsigned int j = csr.row_offsets[i]; j < csr.row_offsets[i + 1]; j++) {
            mat[i * n + csr.column_indices[j]] = csr.values[j];
        }
    }
    return mat;
}

void free_CSR(CSR *csr) {
    free(csr->row_offsets);
    free(csr->column_indices);
    free(csr->values);
}

int my_sparse(CSR *csr, double vec[], double result[]) {
    // Go through all rows
    for (unsigned int i = 0; i < csr->size_row_offsets - 1; i++) {
        result[i] = 0.0;
        // Go through all columns of each row
        for (unsigned int j = csr->row_offsets[i]; j < csr->row_offsets[i + 1]; j++) {
            result[i] += csr->values[j] * vec[csr->column_indices[j]]; // Select the column shortcut
        }
    }
    return 0;
}
