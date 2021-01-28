
#ifndef _GIVENS_QR_H
#define _GIVENS_QR_H
/* function declarations */
int find_para(double* x, double* y, double* c, double* s);
void givens_qrde_1(matrix* input, matrix* R);
void epi_givens_qrde_single_core_2(void *dev, shm_t *shm, void *e_mem, matrix *input, matrix *R);
void epi_givens_qrde_two_core_3(void *dev, shm_t *shm, void *e_mem, matrix *input, matrix *R);
void qr_deompostion_on_four_core(void *dev, shm_t *shm, void *e_mem, matrix *input, matrix *R);
#endif
