#include <stdio.h>
#include <stdlib.h>
#include "timer.h"
#include "spmv.h"
#include "my_sparse.h"
#include <stdbool.h>
#include "spmv.h"
#include "my_sparse.h"
#include <stdio.h>
#include <stdlib.h>
#include "timer.h"
#include "spmv.h"
#include <stdbool.h>
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
    bool is_first_element = true;
    for (unsigned int i = 0; i < n; i++){
        for (unsigned int j = 0; j < n ; j++){
            if (mat[i * n + j] !=0){
                if (is_first_element){// The first element of the row
                    csr.row_offsets[i] = buffer;
                    is_first_element = false;
                }
                csr.column_indices[buffer] = j;
                csr.values[buffer] = mat[i * n + j];
                buffer ++;
            }
        }
        if (is_first_element){//The row is empty
            csr.row_offsets[i] = csr.row_offsets[i-1];
        }
        is_first_element = true;
    }
    csr.row_offsets[n] = buffer;
    return csr;
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

int my_sparse(CSR *csr, double vec[], double result[]) {
    // Go through all rows
    unsigned int index_column_start = 0;
    unsigned int buffer = 1;
    unsigned int index_column_end = csr->row_offsets[buffer];
    //  Initialisation
    while (index_column_start == index_column_end){
        buffer++;
        index_column_end = csr->row_offsets[buffer];
    }
    for (unsigned int j = index_column_start; j < index_column_end; j++) {
        result[0] += csr->values[j] * vec[csr->column_indices[j]]; // Select the column shortcut
    }
    for (unsigned int i = 1; i < csr->size_row_offsets - 1; i++) {
        result[i] = 0.0;
        printf("row_offsets[i-1]: %.2u | row_offsets[i]: %.2u | i: %d\n", csr->row_offsets[i-1], csr->row_offsets[i], i);
        if (csr->row_offsets[i-1]!=csr->row_offsets[i]){
            buffer = i+1;
            index_column_start = csr->row_offsets[i];
            index_column_end = csr->row_offsets[buffer];
            // Go through all columns of each row
            while (index_column_start == index_column_end){
                buffer++;
                index_column_end = csr->row_offsets[buffer];
            }
            printf("index_column_start: %d, index_column_end: %d\n", index_column_start, index_column_end);
            for (unsigned int j = index_column_start; j < index_column_end; j++) {
                result[i] += csr->values[j] * vec[csr->column_indices[j]]; // Select the column shortcut
            }
        }
    }
    return 0;
}

int main() {
    unsigned int n = 5;
    //https://www.researchgate.net/figure/A-sparse-matrix-and-its-CSR-format_fig1_273788746 (deprecated)
    //https://docs.nvidia.com/nvpl/_static/sparse/storage_format/sparse_matrix.html
    double mat_[] = {
            0.0, 0.3, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            4.0, 5.0, 0.0, 0.0,
            0.0, 6.0, 7.0, 8.0
    };
    double mat[] = {
            1.0, 0.0, 2.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            4.0, 5.0, 0.0, 0.0, 0.0,
            0.0, 6.0, 7.0, 8.0, 0.0
    };
    double vec[] = {1.0, 1.0, 1.0, 1.0};
    double result[] = {0.0, 0.0, 0.0, 0.0};
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

    my_sparse(&csr, vec, result);
    for (unsigned int i = 0; i < n; i++) {
        printf("%.2f\n", result[i]);
    }
    free_CSR(&csr);
    return 0;
}


