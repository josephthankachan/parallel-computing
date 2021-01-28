/* code by JOSEPH_THANKACHAN */


#include <stdio.h>
#include <stdlib.h>
#include "Sequence_main.h"

int main(int argc, char** argv)
{
    double sequence_parallella;
    unsigned long long sequence_epiphany;
    unsigned long long parallel_epiphany_two_core;
    unsigned long long prllel_epiphany_four_core;
    
    /* Parallella core sequential execution */
    sequence_parallella = sequential_qrde_parallella();
    printf("Parallella core Sequential execytion is : \n");
    printf("%lf us \n", sequence_parallella);
    printf("************************************************************************************************************************************************************************\n");

    /* Single epiphany core sequential execution */
    sequence_epiphany = sequential_qrde_epiphany();
    printf("Epiphany core sequential execution is : \n");
    printf("%llu us \n", sequence_epiphany);
    printf("************************************************************************************************************************************************************************\n ");

    /* Two epiphany cores parallel exection*/
    parallel_epiphany_two_core = parallel_qrde_epiphany();
    printf("Parallel execution on 2 cores is : \n");
    printf("%llu us \n", parallel_epiphany_two_core);
    printf("************************************************************************************************************************************************************************\n");


/* Parallel execution on four epiphany cores */
    printf("************************************************************************************************************************************************************************\n");
    printf("QR decomposition on epiphany four cores\n");
    prllel_epiphany_four_core = parallel_qr_decmp_4_core();
    printf("Parallel execution on two cores \n");
    printf("%llu us \n", prllel_epiphany_four_core);
    printf("************************************************************************************************************************************************************************\n");

}

