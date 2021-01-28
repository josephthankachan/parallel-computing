#ifndef _MATRIX_H_
#define _MATRIX_H_

typedef struct matrix_st {
	unsigned rows;
	unsigned cols;
	double* data;
} matrix;

void mat_init(matrix* mat, unsigned rows, unsigned cols);

double* get_mat_ele(matrix *mat, unsigned row, unsigned col);
void set_mat_ele(matrix *mat, unsigned row, unsigned col, double val);
void mat_free(matrix* mat);
void mat_make_identity(matrix* mat);
void mat_make_transpose(matrix* src_mat, matrix* dest_mat);
void mat_multiply_row(matrix* l, matrix* r, matrix* out, unsigned row, unsigned col);
void mat_copy(matrix* src_mat, matrix* dest_mat);
void mat_copy_from_array(matrix* dest_mat, double* input_data);
void mat_copy_to_array(matrix* src_mat, double* output_data);
void print_mat(matrix* mat);
void mat_copy_from_epiphany(void *dev, unsigned row, unsigned col, unsigned eMatrix, matrix* dest_mat);
void mat_copy_to_epiphany(void *dev, unsigned row, unsigned col, matrix* src_mat, unsigned eMatrix);
#endif  /* _MATRIX_H_ */
