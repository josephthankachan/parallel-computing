


/* Sequential QR decomposition program  */
/* header files for operation*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Comm.h"
#include "Mat.h"
#include "Give_Qr.h"

#include <e-loader.h>
#include "e-hal.h"

#define FAIL(...) { fprintf(stderr, __VA_ARGS__); exit(1); }
#define SHM_OFFSET 0x01000000

/* Give Input matrix */

double input_data[256] = {
 5.124575, 7.898948, 7.946213, 3.592056, 6.081988, 7.166724, 7.788865, 9.798719, 4.258254, 1.054377, 7.652938, 3.552302, 2.828004, 1.762015, 2.638748, 8.069129,
 2.522771, 1.012173, 8.776333, 6.716550, 4.709401, 1.210326, 1.280568, 3.150382, 9.152067, 2.220713, 4.314055, 2.801475, 0.355802, 4.751442, 2.452192, 8.104572,
 8.535061, 1.283802, 0.308171, 1.881512, 6.372757, 4.749593, 5.185006, 8.183881, 0.145697, 9.572844, 4.968377, 2.484793, 3.972274, 0.722540, 4.236591, 8.985795,
 3.747190, 1.811495, 6.834885, 9.154308, 5.894849, 0.924535, 4.470197, 0.745557, 2.168096, 6.291060, 6.742811, 0.373962, 4.722142, 6.529717, 9.824009, 7.356618,
 6.461459, 9.131959, 7.072122, 0.483971, 1.670396, 8.949241, 4.018287, 3.546607, 8.672633, 0.835284, 1.438314, 5.034810, 6.390177, 3.856355, 2.282590, 2.083286,
 9.260374, 6.927009, 4.018142, 0.385422, 1.099953, 3.724873, 9.097155, 9.176734, 6.069993, 8.241838, 9.264992, 1.061532, 3.459950, 2.737167, 2.078934, 3.678259,
 6.428624, 0.599322, 6.768367, 8.416201, 2.400816, 6.837505, 0.179495, 4.902166, 7.376048, 9.049492, 7.510052, 1.206717, 5.326238, 2.406306, 4.195184, 9.811339,
 8.330168, 8.243491, 5.997973, 8.827311, 5.079985, 6.432596, 3.651846, 9.114349, 1.086232, 7.299483, 4.097865, 8.284714, 0.099244, 9.042177, 1.308067, 3.394185,
 6.530458, 2.984305, 6.574254, 5.616988, 9.094666, 5.588702, 5.533086, 7.406141, 4.161861, 0.238815, 9.648564, 6.822000, 0.370360, 4.103983, 8.835959, 5.350417,
 6.909482, 7.122869, 4.538856, 0.413352, 0.563682, 8.213211, 3.641620, 3.405092, 1.848284, 4.093572, 8.746218, 8.198015, 4.431138, 8.115096, 3.752868, 9.309067,
 4.522346, 4.953636, 6.954018, 4.096635, 1.076255, 2.290696, 4.178736, 1.015898, 5.956050, 5.301290, 3.462898, 9.802859, 9.976361, 8.454410, 2.433835, 3.138024,
 6.622229, 8.626863, 4.287552, 5.543465, 5.243287, 7.461806, 2.883427, 3.698015, 4.793037, 1.159368, 1.629081, 2.267589, 0.671444, 9.046292, 4.284649, 2.955473,
 4.077998, 8.902470, 5.945147, 3.556053, 8.246637, 4.732425, 0.047543, 5.104174, 6.272343, 6.639598, 9.528616, 1.437121, 6.459073, 0.140771, 9.660500, 8.103308,
 7.141063, 7.821795, 4.805562, 5.816616, 8.541925, 0.016668, 2.015325, 4.703988, 5.374593, 0.995246, 5.829943, 5.670019, 8.307987, 6.370590, 9.468715, 1.161778,
 7.457226, 0.099300, 7.410219, 0.187684, 1.784756, 1.066501, 8.859654, 6.758221, 9.879946, 3.690862, 4.245235, 3.524249, 2.964446, 4.373680, 2.311288, 2.500651,
 4.752189, 0.592439, 2.946092, 0.001511, 6.350062, 6.068312, 2.986027, 0.368879, 3.349494, 5.183024, 7.386629, 0.150104, 9.933439, 1.343407, 0.192133, 6.214067
};

double *output_data;

e_platform_t platform;
e_epiphany_t dev;
e_mem_t emem;

double sequential_qrde_parallella(void)
{

    matrix input_matrix, R;
    clock_t start, end;
    double seconds;

    /* matrix initialisation then values of matrix copying */

    mat_init(&input_matrix, INPUT_ROWS, INPUT_COLS);
    mat_init(&R, INPUT_ROWS, INPUT_COLS);
    
    mat_copy_from_array(&input_matrix, input_data);
    printf("The input matrix : \n");
    printf("_____________________\n");
    printf("\n");
    print_mat(&input_matrix);
    
    start = clock();

   /* calling QR decomposition */
  
    givens_qrde_1(&input_matrix, &R);
    end = clock();
    seconds = (double)(end - start)/CLOCKS_PER_SEC;
   /* coping data to the output matrix */

    output_data = malloc(INPUT_ROWS*INPUT_COLS*sizeof(double));

    mat_copy_to_array(&R, output_data);
    printf("The  R matrix :\n");
    printf("____________________\n ");
    printf("\n ");
    print_mat(&R);

    
    printf("Time taken for computation is %lf  us \n", seconds*1000000 );
    /*finaly free the matrix */

    mat_free(&input_matrix);
    mat_free(&R);
    free(output_data);
    return seconds*1000000;
}

unsigned long long sequential_qrde_epiphany(void)
{
    unsigned long long us = 0;
    shm_t shm;
    unsigned clr = 0;
    unsigned set = 1;
    unsigned corenum = 1;
    unsigned done = 0, active = 0;
    int i, j;
    /*matrix initalisation*/
    matrix input_matrix, R;
    mat_init(&input_matrix, INPUT_ROWS, INPUT_COLS);
    mat_init(&R, INPUT_ROWS, INPUT_COLS);
   /* matrix copy from array and print the input matrix*/
    mat_copy_from_array(&input_matrix, input_data);
    printf("The input matrix is :\n");
    printf("_________________________\n");
    printf("\n");
    print_mat(&input_matrix);

	if(E_OK != e_init(NULL))
		FAIL("Can't init!\n");

	e_reset_system();

	e_get_platform_info(&platform);
	
	if(E_OK != e_open(&dev, 0, 0, ROW, COL))
		FAIL("Can't open workgroup!\n");
    
	if(E_OK != e_alloc(&emem, SHM_OFFSET, sizeof(shm_t)))
		FAIL("Can't alloc!\n");

	if(E_OK != e_load_group("e_seq.elf", &dev, 0, 0, ROW, COL, E_FALSE))
	  	FAIL("Can't load! \n")
    
    for(i = 0; i < ROW; i++) {
        for(j = 0; j < COL; j++) {
            e_write(&dev, i, j, ACTIVE_ADDR, &clr, sizeof(clr));
            e_write(&dev, i, j, DONE_ADDR, &clr, sizeof(clr));
            e_write(&dev, i, j, CORENUM_ADDR, &corenum, sizeof(corenum));
            shm.go[corenum] = 0;
            corenum++;
        }
    }
    e_write(&dev, 0, 0, ACTIVE_ADDR, &set, sizeof(set));

    if(E_ERR == e_write(&emem, 0, 0, (off_t)0, &shm, sizeof(shm_t)))
		FAIL("Can't write to memory!\n");

    e_start_group(&dev);

    epi_givens_qrde_single_core_2(&dev, &shm, &emem, &input_matrix, &R);
    
    while(!done) {
        e_read(&dev, 0, 0, DONE_ADDR, &done, sizeof(done));
    }
    
    printf("Sequential computation on a single core completed \n");

	if(E_ERR == e_read(&emem, 0, 0, (off_t)0, &shm,	sizeof(shm_t)))
	 	FAIL("Can't poll!\n");
	
    us = shm.total_us[1];
    if(E_OK != e_free(&emem)) FAIL("Can't free!\n");
	if(E_OK != e_close(&dev)) FAIL("Can't close!\n");
	if(E_OK != e_finalize())  FAIL("Can't finalize!\n");
    return us;
}

unsigned long long parallel_qrde_epiphany(void)
{
    unsigned long long us = 0;
    shm_t shm;
    unsigned clr = 0;
    unsigned set = 1;
    unsigned corenum = 1;
    unsigned done = 0, active = 0;
    int i, j;
    matrix input_matrix, R;
    /*matrix initialisation*/
    mat_init(&input_matrix, INPUT_ROWS, INPUT_COLS);
    mat_init(&R, INPUT_ROWS, INPUT_COLS);
    /* copy matrix from array and printing input matrix*/   
    mat_copy_from_array(&input_matrix, input_data);
    printf("The input matrix is \n");
    printf("_________________________\n ");
    printf("\n");
    print_mat(&input_matrix);

	if(E_OK != e_init(NULL))
		FAIL("Can't init!\n");
	e_reset_system();

	e_get_platform_info(&platform);
	
	if(E_OK != e_open(&dev, 0, 0, ROW, COL))
		FAIL("Can't open workgroup!\n");
    
	if(E_OK != e_alloc(&emem, SHM_OFFSET, sizeof(shm_t)))
		FAIL("Can't alloc!\n");

	if(E_OK != e_load_group("e_two_core.elf", &dev, 0, 0, ROW, COL, E_FALSE))
	  	FAIL("Can't load! \n")
    
    	for(i = 0; i < ROW; i++) {
        	for(j = 0; j < COL; j++) {
            		e_write(&dev, i, j, ACTIVE_ADDR, &clr, sizeof(clr));
            		e_write(&dev, i, j, DONE_ADDR, &clr, sizeof(clr));
            		e_write(&dev, i, j, CORENUM_ADDR, &corenum, sizeof(corenum));
            		shm.go[corenum] = 0;
            		corenum++;
        	}
    	}
    shm.active[1] = 1;
    shm.active[2] = 1;


    if(E_ERR == e_write(&emem, 0, 0, (off_t)0, &shm, sizeof(shm_t)))
		FAIL("Can't write to memory!\n");

    e_start_group(&dev);

    epi_givens_qrde_two_core_3(&dev, &shm, &emem, &input_matrix, &R);
    
    while(!done) {
        e_read(&dev, 0, 0, DONE_ADDR, &done, sizeof(done));
    }
    printf("Sequential computation on two cores completed \n");


	if(E_ERR == e_read(&emem, 0, 0, (off_t)0, &shm,	sizeof(shm_t)))
	 	FAIL("Can't poll!\n");
	
    us = (shm.total_us[1] + shm.total_us[2])/2;
    
    if(E_OK != e_free(&emem)) FAIL("Can't free!\n");
	if(E_OK != e_close(&dev)) FAIL("Can't close!\n");
	if(E_OK != e_finalize())  FAIL("Can't finalize!\n");
    return us;
}

unsigned long long parallel_qr_decmp_4_core(void)
{
    unsigned long long us = 0;
    shm_t shm;
    unsigned clr = 0;
    unsigned set = 1;
    unsigned corenum = 0;
    unsigned done1 = 0, done2 = 0,done3 = 0, done4 = 0, active = 0;
    int i, j;
    matrix input_matrix, R, Q;

    /* Init the matrices */
    mat_init(&input_matrix, INPUT_ROWS, INPUT_COLS);
    mat_init(&R, INPUT_ROWS, INPUT_COLS);
    
    /* Copy data to the input matrix */
    mat_copy_from_array(&input_matrix, input_data);
    printf("Input matrix : \n");
    print_mat(&input_matrix);

	if(E_OK != e_init(NULL))
		FAIL("ERROR WHILE INITIALIZING THE BOARD!\n");
	/* reset epiphany chip */
	e_reset_system();

	e_get_platform_info(&platform);
	
	if(E_OK != e_open(&dev, 0, 0, ROW, COL))
		FAIL("ERROR WHILE OPENING THE WORKGROUP!\n");
    
    /* initialize and allocate shared memory buffer */
	if(E_OK != e_alloc(&emem, SHM_OFFSET, sizeof(shm_t)))
		FAIL("ERROR WHILE ALLOCATING THE SHARED MEMMORY! \n");

	if(E_OK != e_load_group("epiphany_four_core_util.elf", &dev, 0, 0, ROW, COL, E_FALSE))
	  	FAIL("ERROR WHILE LOADING THE GROUPS \n")
    
    for(i = 0; i < ROW; i++) {
        for(j = 0; j < COL; j++) {
            e_write(&dev, i, j, ACTIVE_ADDR, &clr, sizeof(clr));
            e_write(&dev, i, j, DONE_ADDR, &clr, sizeof(clr));
            e_write(&dev, i, j, CORENUM_ADDR, &corenum, sizeof(corenum));
            shm.go[corenum] = 0;
            corenum++;
        }
    }
    shm.active[0] = 1;
    shm.active[1] = 1;
    shm.active[2] = 1;
    shm.active[3] = 1;


    if(E_ERR == e_write(&emem, 0, 0, (off_t)0, &shm, sizeof(shm_t)))
		FAIL("ERROR WHILE WRITING TO MEMORY!\n");

    /* Start the program on epiphany core */
    e_start_group(&dev);

    qr_deompostion_on_four_core(&dev, &shm, &emem, &input_matrix, &R);
    
   
    while(!done1 || !done2 || !done3 || !done4) {
        e_read(&dev, 0, 0, DONE_ADDR, &done1, sizeof(done1));
        e_read(&dev, 0, 1, DONE_ADDR, &done2, sizeof(done2));
        e_read(&dev, 0, 2, DONE_ADDR, &done3, sizeof(done3));
        e_read(&dev, 0, 3, DONE_ADDR, &done4, sizeof(done4));
    }
    printf("Sequential computation on four cores completed \n");

	if(E_ERR == e_read(&emem, 0, 0, (off_t)0, &shm,	sizeof(shm_t)))
	 	FAIL("ERROR WHILE POLLING!\n");
	
    us = (shm.total_us[0] + shm.total_us[1]+ shm.total_us[2] + shm.total_us[3])/4;
    printf("core1 %5.2llu core2  %5.2llu core3  %5.2llu core4  %5.2llu us ",shm.total_us[0],shm.total_us[1],shm.total_us[2],shm.total_us[3]);
    
    /* Free epiphany resources */
    if(E_OK != e_free(&emem)) FAIL("ERROR WHILE FREEING!\n");
	if(E_OK != e_close(&dev)) FAIL("ERROR WHILE CLOSING!\n");
	if(E_OK != e_finalize())  FAIL("ERROR WHILE FINALIZING!\n");
    return us;
}



























	
































    











	



    

    






