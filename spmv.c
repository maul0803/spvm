#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_cblas.h>      // CBLAS in GSL (the GNU Scientific Library)
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_vector.h>
#include "timer.h"
#include "spmv.h"
#include "my_sparse.h"
// Added libraries
#include <gsl/gsl_blas.h>
#include <gsl/gsl_spblas.h>
#define DEFAULT_SIZE 1024
#define DEFAULT_DENSITY 25.0

unsigned int populate_sparse_matrix(double mat[], unsigned int n, double density, unsigned int seed)
{
  unsigned int nnz = 0;

  srand(seed);

  for (unsigned int i = 0; i < n * n; i++) {
    if ((rand() % 100) / 100.0 < density) {
      // Get a pseudorandom value between -9.99 e 9.99
      mat[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
      nnz++;
    } else {
      mat[i] = 0;
    }
  }

  return nnz;
}

unsigned int populate_vector(double vec[], unsigned int size, unsigned int seed)
{
  srand(seed);

  for (unsigned int i = 0; i < size; i++) {
    vec[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
  }

  return size;
}

int is_nearly_equal(double x, double y)
{
  const double epsilon = 1e-5 /* some small number */;
  return fabs(x - y) <= epsilon * fabs(x);
  // see Knuth section 4.2.2 pages 217-218
}

unsigned int check_result(double ref[], double result[], unsigned int size)
{
  for(unsigned int i = 0; i < size; i++) {
    if (!is_nearly_equal(ref[i], result[i]))
      return 0;
  }

  return 1;
}

int main(int argc, char *argv[])
{
  int size;        // number of rows and cols (size x size matrix)
  double density;  // aprox. ratio of non-zero values

  if (argc < 2) {
    size = DEFAULT_SIZE;
    density = DEFAULT_DENSITY;
  } else if (argc < 3) {
    size = atoi(argv[1]);
    density = DEFAULT_DENSITY;
  } else {
    size = atoi(argv[1]);
    density = atoi(argv[2]);
  }

  double *mat, *vec, *refsol, *mysol;

  mat = (double *) malloc(size * size * sizeof(double));
  vec = (double *) malloc(size * sizeof(double));
  refsol = (double *) malloc(size * sizeof(double));
  mysol = (double *) malloc(size * sizeof(double));

  unsigned int nnz = populate_sparse_matrix(mat, size, density, 1);
  populate_vector(vec, size, 2);

  printf("Matriz size: %d x %d (%d elements)\n", size, size, size*size);
  printf("%d non-zero elements (%.2lf%%)\n\n", nnz, (double) nnz / (size*size) * 100.0);
  //
  // Dense computation using CBLAS (eg. GSL's CBLAS implementation)
  //
  printf("\nDense computation\n----------------\n");

  timeinfo start, now;
  timestamp(&start);

  cblas_dgemv(CblasRowMajor, CblasNoTrans, size, size, 1.0, mat, size, vec, 1, 0.0, refsol, 1);

  timestamp(&now);
  printf("Time taken by CBLAS dense computation: %ld ms\n", diff_milli(&start, &now));

  //
  // Using your own dense implementation
  //
  printf("\nMy Own Dense computation\n----------------\n");
  timestamp(&start);

  my_dense(size, mat, vec, mysol);

  timestamp(&now);
  printf("Time taken by my dense matrix-vector product: %ld ms\n", diff_milli(&start, &now));

  if (check_result(refsol, mysol, size) == 1)
    printf("Result my dense is ok!\n");
  else
    printf("Result my dense is wrong!\n");
  //
  // Let's try now SpMV: Sparse Matrix - Dense Vector computation
  //
  printf("\nSparse computation\n----------------\n");
  // Convert mat to a sparse format: CSR
  // Use the gsl_spmatrix struct as datatype
  // Sparse computation using GSL's sparse algebra functions
  //
  // Initialisation of the vector
  gsl_vector *gsl_vec = gsl_vector_calloc(size);
  for (int i = 0; i < size; i++){
      double val = vec[i];
      gsl_vector_set(gsl_vec, i, val);
  }
  // Initialisation of the dense matrix
  gsl_matrix *gsl_dense = gsl_matrix_alloc(size, size);
  for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
          double val = mat[i * size + j];
          gsl_matrix_set(gsl_dense, i, j, val);
      }
  }
  // Convert dense to COO format
  gsl_spmatrix *gsl_coo = gsl_spmatrix_alloc(size, size);
  gsl_spmatrix_d2sp(gsl_coo, gsl_dense);
  // Convert COO to CSR format
  gsl_spmatrix *gsl_csr = gsl_spmatrix_crs(gsl_coo);
  // Multiply the matrix and the vector
  gsl_vector *gsl_mysol = gsl_vector_calloc(size);
  // Results
  timestamp(&start);
  gsl_spblas_dgemv(CblasNoTrans, 1, gsl_csr, gsl_vec, 0, gsl_mysol);
  timestamp(&now);
  printf("Time taken by SpMV sparse matrix - vector product: %ld ms\n", diff_milli(&start, &now));

  // Initialisation of the solution
  for (int i = 0; i < size; i++){
      double val = gsl_vector_get(gsl_mysol, i);
      mysol[i] = val;
  }
  if (check_result(refsol, mysol, size) == 1)
      printf("Result SpMV sparse implementation is ok!\n");
  else
      printf("Result SpMV spare implementation is wrong!\n");


  // Your own sparse implementation
  //
  printf("\nMy Own Spare computation\n----------------\n");

  CSR per_mat_csr = convert_dense_to_CSR(size, mat);
  // Compare times (and computation correctness!)
  timestamp(&start);
  my_sparse(&per_mat_csr, vec, mysol);
  timestamp(&now);
  printf("Time taken by my own sparse matrix - vector product: %ld ms\n", diff_milli(&start, &now));
  if (check_result(refsol, mysol, size) == 1)
      printf("Result own spare implementation is ok!\n");
  else
      printf("Result own spare implementation is wrong!\n");
  // Free resources
  free(mat);
  free(vec);
  free(refsol);
  free(mysol);
  return 0;
}
