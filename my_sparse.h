//
// Created by cytech on 02/10/24.
//

#ifndef SPVM_MY_SPARSE_H
#define SPVM_MY_SPARSE_H

typedef struct {
    unsigned int *row_offsets;
    unsigned int *column_indices;
    double *values;
    unsigned int size_row_offsets;
    unsigned int size_column_indices_values;
} CSR;

CSR convert_dense_to_CSR(unsigned int n, const double mat[]);
double* convert_CSR_to_dense(CSR csr, unsigned int n);
void free_CSR(CSR *csr);
int my_sparse(CSR *csr, double vec[], double result[]);

#endif //SPVM_MY_SPARSE_H
