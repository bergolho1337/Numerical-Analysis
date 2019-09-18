#include "sparse_matrix.h"

struct csr_matrix* new_csr_matrix ()
{
    struct csr_matrix *result = (struct csr_matrix*)malloc(sizeof(struct csr_matrix));

    result->num_elements = 0;
    result->num_rows = 0;
    result->num_cols = 0;
    result->rows = NULL;
    result->cols = NULL;
    result->data = NULL;

    return result;
}

void read_csr_sparse_matrix_from_mtx_file (struct csr_matrix *csr, const char mtx_file[])
{

    FILE *file = fopen(mtx_file,"r");
    if (!file)
    {
        fprintf(stderr,"[-] ERROR! Cannot open file '%s'\n",mtx_file);
        exit(EXIT_FAILURE);
    }

    char str[50];
    uint32_t num_rows, num_cols, num_elements;
    uint32_t row, col;
    double value;

    fgets(str, sizeof(str), file);
    fscanf(file,"%u %u %u",&num_rows,&num_cols,&num_elements);

    csr->num_elements = num_elements;
    csr->num_rows = num_rows;
    csr->num_cols = num_cols;
    
    // Alocate memory
    csr->data = (double*)calloc(num_elements,sizeof(double));
    csr->rows = (uint32_t*)calloc(num_elements,sizeof(double));
    csr->cols = (uint32_t*)calloc(num_elements,sizeof(double));

    uint32_t counter = 0;
    uint32_t cur_row = -1;
    while (fscanf(file,"%u %u %lf",&row,&col,&value) != EOF)
    {
        row--; col--;

        csr->cols[counter] = col;
        csr->data[counter] = value;
        if (row != cur_row)
        {
            csr->rows[row] = counter;
            cur_row++;
        }
        counter++;
    }
    csr->rows[cur_row+1] = counter;
    fclose(file);

}

void csr_matrix_vector_multiplication(struct csr_matrix *csr)
{
    double x[csr->num_cols];
    double y[csr->num_rows];

    // Unit vector for testing : Sum of the elements from the line
    for (uint32_t i = 0; i < csr->num_cols; i++)
        x[i] = 1.0;
    
    for (uint32_t i = 0; i < csr->num_rows; i++)
    {
        double sum = 0.0;
        uint32_t start_index = csr->rows[i];
        uint32_t end_index = csr->rows[i+1];
        for (uint32_t j = start_index; j < end_index; j++)
            sum += csr->data[j] * x[csr->cols[j]];
        y[i] = sum;
    }

    for (uint32_t i = 0; i < csr->num_rows; i++)
        printf("%g\n",y[i]);
}

void print_csr_sparse_matrix (struct csr_matrix *csr)
{
    uint32_t num_elements = csr->num_elements;

    printf("data\n");
    for (uint32_t i = 0; i < num_elements; i++)
    {
        printf("%g ",csr->data[i]);
    }
    printf("\n");

    printf("cols\n");
    for (uint32_t i = 0; i < num_elements; i++)
    {
        printf("%u ",csr->cols[i]);
    }
    printf("\n");

    printf("rows\n");
    for (uint32_t i = 0; i < num_elements; i++)
    {
        printf("%u ",csr->rows[i]);
    }
    printf("\n");

    printf("\nFormatted output\n");
    for (uint32_t i = 0; i < csr->num_rows; i++) 
    {
        uint32_t start_row = csr->rows[i];
        uint32_t end_row = csr->rows[i+1];
        double col_values[csr->num_cols];
        memset(col_values,0.0,sizeof(double)*csr->num_cols);

        for (uint32_t j = start_row; j < end_row; j++)
        {
            uint32_t col = csr->cols[j];
            double value = csr->data[j];
            
            col_values[col] = value;
        }

        for (uint32_t j = 0; j < csr->num_cols; j++)
            printf("%g ",col_values[j]);
        printf("\n");
    }  
}