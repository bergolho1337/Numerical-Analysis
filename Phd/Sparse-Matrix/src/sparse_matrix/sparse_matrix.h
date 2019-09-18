#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

struct csr_matrix
{
    uint32_t num_rows;          // Number of lines
    uint32_t num_cols;          // Number of columns
    uint32_t num_elements;      // Number of non-zero elements
    uint32_t *rows;             // Lines indexes
    uint32_t *cols;             // Columns indexes
    double *data;               // Non-zero elements
};

struct csr_matrix* new_csr_matrix ();
void read_csr_sparse_matrix_from_mtx_file (struct csr_matrix *csr, const char mtx_file[]);
void csr_matrix_vector_multiplication(struct csr_matrix *csr);
void print_csr_sparse_matrix (struct csr_matrix *csr);


#endif