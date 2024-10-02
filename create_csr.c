#include <stdio.h>
#include <stdlib.h>
#include "timer.h"
#include "spmv.h"
#include "my_sparse.h"
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
    csr.row_offsets = (unsigned int *)malloc((size_row_offsets) * sizeof(unsigned int)); // n+1 pour le dernier offset
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

void print_CSR(const CSR *csr, const unsigned int n) {
    printf("Row Offsets (size %u):\n", csr->size_row_offsets);
    for (unsigned int i = 0; i < csr->size_row_offsets; i++) {
        printf("%u ", csr->row_offsets[i]);
    }
    printf("\n\n");

    printf("Column Indices (size %u):\n", csr->size_column_indices_values);
    for (unsigned int i = 0; i < csr->size_column_indices_values; i++) {
        printf("%u ", csr->column_indices[i]);
    }
    printf("\n\n");

    printf("Values (size %u):\n", csr->size_column_indices_values);
    for (unsigned int i = 0; i < csr->size_column_indices_values; i++) {
        printf("%.2f ", csr->values[i]);
    }
    printf("\n\n");
}

int main() {
    unsigned int n = 4;
    //https://www.researchgate.net/figure/A-sparse-matrix-and-its-CSR-format_fig1_273788746
    double mat[] = {
            1.0, 0.0, 2.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            1.0, 0.0, 2.0, 3.0,
            0.0, 1.0, 0.0, 2.0
    };

    printf("Dense Matrix (%u x %u):\n", n, n);
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
            printf("%.2f\t", mat[i * n + j]);
        }
        printf("\n");
    }
    printf("\n");
    CSR csr = convert_dense_to_CSR(n, mat);
    print_CSR(&csr, n);
    free_CSR(&csr);
    return 0;
}