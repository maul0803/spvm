#include "spmv.h"

#include <stdio.h>

typedef struct {
    unsigned int *row_indices; // Indices des lignes
    unsigned int *col_indices; // Indices des colonnes
    double *values;            // Valeurs non nulles
    unsigned int nnz;          // Nombre d'éléments non nuls
} SparseMatrix;

int my_sparse(const SparseMatrix *mat, const double vec[], double result[]) {
    // Initialiser le tableau result à 0
    for (unsigned int i = 0; i < mat->nnz; i++) {
        result[i] = 0.0;
    }

    // Multiplication matrice creuse-vecteur
    for (unsigned int i = 0; i < mat->nnz; i++) {
        unsigned int row = mat->row_indices[i];
        result[row] += mat->values[i] * vec[mat->col_indices[i]];
    }

    return 0; // Retourner 0 pour indiquer que l'opération s'est bien passée
}
