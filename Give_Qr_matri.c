#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Mat.h"
#include "e-hal.h"
#include "Comm.h"
#include <string.h>


/* Givens Rotation program*/

int find_para(double* x, double* y, double* c, double* s)
{
    if((*y) != 0.0) {
        double distance;
        distance = sqrt((*x)*(*x) + (*y)*(*y));
        *s = (*y) / distance;
        *c = (*x) / distance;
        return 0;
        }
    else {
        return 1;
    }
}

void givens_qrde_1(matrix *input, matrix *R)
{
    double c = 0.0, s = 0.0;
    matrix r_temp, G;
    mat_init(&r_temp, R->rows, R->cols);
    mat_init(&G, input->rows, input->cols);
    mat_copy(input , R);
    for(unsigned col = 0; col < input->cols - 1; col++) {
        for(unsigned row = col + 1; row < input->rows; row++) {
            if(!find_para(get_mat_ele(R, col, col), get_mat_ele(R, row, col), &c, &s)) {
                mat_make_identity(&G);
                set_mat_ele(&G, col, col, c);
                set_mat_ele(&G, col, row, s);
                set_mat_ele(&G, row, col, -s);
                set_mat_ele(&G, row, row, c);

                mat_copy(R, &r_temp);
                for(unsigned i= 0; i < R->rows; i++) {
                    mat_multiply_row(&G, &r_temp, R, i, col);
                    mat_multiply_row(&G, &r_temp, R, i, row);

                }
            }
        }
    }
    mat_free(&r_temp);
    mat_free(&G);
}

void epi_givens_qrde_single_core_2(void *dev, shm_t *shm, void *e_mem, matrix *input, matrix *R)
{
    double c = 0.0, s = 0.0;
    unsigned clr = 0;
    unsigned set = 0;
    unsigned ready = 0;
    unsigned go = 0;
    matrix r_temp, G;
    mat_init(&r_temp, R->rows, R->cols);
    mat_init(&G, input->rows, input->cols);
    mat_copy(input , R);
    for(unsigned col = 0; col < input->cols - 1; col++) {
        for(unsigned row = col + 1; row < input->rows; row++) {
            if(!find_para(get_mat_ele(R, col, col), get_mat_ele(R, row, col), &c, &s)) {
                mat_make_identity(&G);
                set_mat_ele(&G, col, col, c);
                set_mat_ele(&G, col, row, s);
                set_mat_ele(&G, row, col, -s);
                set_mat_ele(&G, row, row, c);

                mat_copy(R, &r_temp);
                mat_copy_to_epiphany(dev, 0, 0, &G, Amat);
                mat_copy_to_epiphany(dev, 0, 0, &r_temp, Bmat);
                mat_copy_to_epiphany(dev, 0, 0, R, Cmat);
                e_write(dev, 0, 0, Rcomp, &row, sizeof(row));
                e_write(dev, 0, 0, Ccomp, &col, sizeof(col));
                shm->go[1] = 1;
                
		if(E_ERR == e_write(e_mem, 0, 0, (off_t)0, shm, sizeof(shm_t)))
		            printf("Can't write to memory!\n");

                while(!ready) {
                    e_read(dev, 0, 0, READY, &ready, sizeof(ready));
                }

                ready = 0;
                e_write(dev, 0, 0, READY, &ready, sizeof(ready));
                mat_copy_from_epiphany(dev, 0, 0, Cmat, R);

            }
        }
    }
    shm->active[1] = 0;
    if(E_ERR == e_write(e_mem, 0, 0, (off_t)0, shm, sizeof(shm_t)))
		printf("Can't write to memory!\n");
    printf("matrix R is :\n");
    printf("___________________\n");
    printf("\n");
    print_mat(R);
    mat_free(&r_temp);
    mat_free(&G);
}

void epi_givens_qrde_two_core_3(void *dev, shm_t *shm, void *e_mem, matrix *input, matrix *R)
{
    double c = 0.0, s = 0.0;
    unsigned clr = 0;
    unsigned set = 0;
    unsigned ready1 = 0, ready2 = 0;
    unsigned go = 0;
    matrix r_temp, G;
    mat_init(&r_temp, R->rows, R->cols);
    mat_init(&G, input->rows, input->cols);
    mat_copy(input , R);
    for(unsigned col = 0; col < input->cols - 1; col++) {
        for(unsigned row = col + 1; row < input->rows; row++) {
            if(!find_para(get_mat_ele(R, col, col), get_mat_ele(R, row, col), &c, &s)) {
                mat_make_identity(&G);
                set_mat_ele(&G, col, col, c);
                set_mat_ele(&G, col, row, s);
                set_mat_ele(&G, row, col, -s);
                set_mat_ele(&G, row, row, c);

                mat_copy(R, &r_temp);
                mat_copy_to_epiphany(dev, 0, 0, &G, Amat);
                mat_copy_to_epiphany(dev, 0, 0, &r_temp, Bmat);
                mat_copy_to_epiphany(dev, 0, 1, &G, Amat);
                mat_copy_to_epiphany(dev, 0, 1, &r_temp, Bmat);
                mat_copy_to_epiphany(dev, 0, 0, R, Cmat);
                e_write(dev, 0, 0, Rcomp, &row, sizeof(row));
                e_write(dev, 0, 0, Ccomp, &col, sizeof(col));
                e_write(dev, 0, 1, Rcomp, &row, sizeof(row));
                e_write(dev, 0, 1, Ccomp, &col, sizeof(col));
                shm->go[1] = 1;
                shm->go[2] = 1;
                if(E_ERR == e_write(e_mem, 0, 0, (off_t)0, shm, sizeof(shm_t)))
		            printf("Can't write to memory!\n");

                while(!(ready1&&ready2)) {
                    e_read(dev, 0, 0, READY, &ready1, sizeof(ready1));
                    e_read(dev, 0, 1, READY, &ready2, sizeof(ready2));
                }

                ready1 = 0;
                e_write(dev, 0, 0, READY, &ready1, sizeof(ready1));
                e_write(dev, 0, 1, READY, &ready1, sizeof(ready1));
            }
            mat_copy_from_epiphany(dev, 0, 0, Cmat, R);
        }
    }
    
    shm->active[1] = 0;
    shm->active[2] = 0;
    if(E_ERR == e_write(e_mem, 0, 0, (off_t)0, shm, sizeof(shm_t)))
		printf("Can't write to memory!\n");
    printf("matrix R is :\n");
    printf("___________________\n");
    printf("\n");
    print_mat(R);
    mat_free(&r_temp);
    mat_free(&G);
}

void qr_deompostion_on_four_core(void *dev, shm_t *shm, void *e_mem, matrix *input, matrix *R)
{
    double c = 0.0, s = 0.0;
    unsigned clr = 0;
    unsigned set = 0;
    unsigned ready1 = 0, ready2 = 0, ready3 = 0, ready4 = 0;
    unsigned go = 0;
    matrix r_temp, G;
    mat_init(&r_temp, R->rows, R->cols);
    mat_init(&G, input->rows, input->cols);
    mat_copy(input , R);
    for(unsigned col = 0; col < input->cols - 1; col++) {
        for(unsigned row = col + 1; row < input->rows; row++) {
            if(!find_para(get_mat_ele(R, col, col), get_mat_ele(R, row, col), &c, &s)) {
                mat_make_identity(&G);
                set_mat_ele(&G, col, col, c);
                set_mat_ele(&G, col, row, s);
                set_mat_ele(&G, row, col, -s);
                set_mat_ele(&G, row, row, c);

                mat_copy(R, &r_temp);
                mat_copy_to_epiphany(dev, 0, 0, &G, Amat);
                mat_copy_to_epiphany(dev, 0, 0, &r_temp, Bmat);
                mat_copy_to_epiphany(dev, 0, 1, &G, Amat);
                mat_copy_to_epiphany(dev, 0, 1, &r_temp, Bmat);
                mat_copy_to_epiphany(dev, 0, 2, &G, Amat);
                mat_copy_to_epiphany(dev, 0, 2, &r_temp, Bmat);
                mat_copy_to_epiphany(dev, 0, 3, &G, Amat);
                mat_copy_to_epiphany(dev, 0, 3, &r_temp, Bmat);
                mat_copy_to_epiphany(dev, 0, 0, R, Cmat);
                e_write(dev, 0, 0, Rcomp, &row, sizeof(row));
                e_write(dev, 0, 0, Ccomp, &col, sizeof(col));
                e_write(dev, 0, 1, Rcomp, &row, sizeof(row));
                e_write(dev, 0, 1, Ccomp, &col, sizeof(col));
                e_write(dev, 0, 2, Rcomp, &row, sizeof(row));
                e_write(dev, 0, 2, Ccomp, &col, sizeof(col));
                e_write(dev, 0, 3, Rcomp, &row, sizeof(row));
                e_write(dev, 0, 3, Ccomp, &col, sizeof(col));
                shm->go[0] = 1;
                shm->go[1] = 1;
                shm->go[2] = 1;
                shm->go[3] = 1;
                if(E_ERR == e_write(e_mem, 0, 0, (off_t)0, shm, sizeof(shm_t)))
		            printf("Can't write to memory!\n");

                while(!(ready1&&ready2&&ready3&&ready4)) {
                    e_read(dev, 0, 0, READY, &ready1, sizeof(ready1));
                    e_read(dev, 0, 1, READY, &ready2, sizeof(ready2));
                    e_read(dev, 0, 2, READY, &ready3, sizeof(ready3));
                    e_read(dev, 0, 3, READY, &ready4, sizeof(ready4));
                    //printf("Epiphany core data not ready \n ");
                }

                ready1 = 0;
                ready2 = 0;
                ready3 = 0;
                ready4 = 0;
                e_write(dev, 0, 0, READY, &ready1, sizeof(ready1));
                e_write(dev, 0, 1, READY, &ready1, sizeof(ready2));
                e_write(dev, 0, 2, READY, &ready1, sizeof(ready3));
                e_write(dev, 0, 3, READY, &ready1, sizeof(ready4));

                //printf("%x \n", (off_t) Cmat + col*SZ1);
            }
            // printf("Column %d \n", col);
            mat_copy_from_epiphany(dev, 0, 0, Cmat, R);
        }
    }
    
    /* idle core setting */
    shm->active[0] = 0;
    shm->active[1] = 0;
    shm->active[2] = 0;
    shm->active[3] = 0;
    if(E_ERR == e_write(e_mem, 0, 0, (off_t)0, shm, sizeof(shm_t)))
		printf("ERROR WHILE WRITING TO MEMMMORY!\n");
    print_mat(R);
    mat_free(&r_temp);
    mat_free(&G);
    
}
  
/* ****************************************************** */
/* Matrix program */


size_t SZ = INPUT_ROWS*INPUT_COLS*sizeof(double);
size_t SZ1 = INPUT_ROWS*sizeof(double);

void mat_init(matrix *mat, unsigned rows, unsigned cols)
{
    if(NULL == (mat->data = malloc(rows * cols * sizeof(double)))){
        printf("Couldn't allocate memory for matrix \n");
        exit (1);
    }

    mat->rows = rows;
    mat->cols = cols;    
}


double* get_mat_ele(matrix *mat, unsigned row, unsigned col)
{
    return (mat->data + (col * mat->rows + row));
}

void set_mat_ele(matrix *mat, unsigned row, unsigned col, double val)
{
    mat->data[(col * mat->rows) + row] = (double) val;
}

void mat_free(matrix *mat)
{
    free(mat->data);

    mat->rows = 0;
    mat->cols = 0;

}

void mat_make_identity(matrix *mat)
{
    for(unsigned r = 0; r < mat->rows; r++) {
        for(unsigned c = 0; c < mat->cols; c++) {
            if(c == r)
                set_mat_ele(mat, r, c,  (double) 1);
            else
                set_mat_ele(mat, r, c, (double) 0);
        }
    }
}

void mat_make_transpose(matrix *src_mat, matrix *dest_mat)
{
    for(unsigned r =0; r < src_mat->rows; r++ )  {
        for(unsigned c = 0; c < src_mat->cols; c++)  {
            set_mat_ele(dest_mat, r, c, *(get_mat_ele(src_mat, c, r)));
        }
    }
}

void mat_multiply_row(matrix *l, matrix *r, matrix *out, unsigned row, unsigned col)
{
    double val;
    val = 0.0;
    for(unsigned i =0; i < l->rows; i++) {
        val += (*get_mat_ele(l, row, i))*(*get_mat_ele(r, i, col));
    }
    set_mat_ele(out, row, col, val);
}

void mat_copy(matrix *src_mat, matrix *dest_mat)
{
    memcpy(dest_mat->data, src_mat->data, src_mat->rows*src_mat->cols*sizeof(double));
}

void mat_copy_from_array(matrix *dest_mat, double *input_data)
{
    memcpy(dest_mat->data, input_data, dest_mat->rows*dest_mat->cols*sizeof(double));
}

void mat_copy_to_array(matrix *src_mat, double *output_data)
{
    memcpy(output_data, src_mat->data, src_mat->rows*src_mat->cols*sizeof(double));
}

void print_mat(matrix *mat)
{
    printf("[ \n");
    for(int i=0; i < mat->rows; i++) {
        for(int j = 0; j < mat->cols; j++) {
            printf("%lf, ", *(get_mat_ele(mat, i, j)));
        }
        printf("\n");
    }
    printf("] \n");
}


/* Methods to copy matrix data to and from epiphany core */
void mat_copy_from_epiphany(void *dev, unsigned row, unsigned col,unsigned eMatrix, matrix *dest_mat)
{
    e_read(dev, row, col, (off_t) eMatrix, dest_mat->data, SZ);
}

void mat_copy_to_epiphany(void *dev, unsigned row, unsigned col, matrix *src_mat, unsigned eMatrix)
{
    e_write(dev, row, col, (off_t) eMatrix, src_mat->data, SZ);
}























































