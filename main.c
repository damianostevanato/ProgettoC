#include <stdio.h>
#include "bmp.h"
#include <math.h>
#include <string.h>
#include "ip_lib.h"

/*
Compilazione:
gcc -c ip_lib.c -o ip_lib.o -lm -Wall
gcc ip_lib.o bmp.o main.c -o main -lm -Wall
Compilazione per verifica memory leaks:
gcc -c ip_lib.c -o ip_lib.o -lm -Wall -g -O1
gcc ip_lib.o bmp.o main.c -o main -lm -Wall -g -O1
valgrind -v --leak-check=full ./main
*/
int main(){
    /*dichiarazione veriabili*/
    ip_mat *mat,*mat2,*concat,*sum,*sub,*scalare,*s_scalare;

    /*operazioni su variabili*/
    mat=ip_mat_create(4,4,3,1.0);
    mat2=ip_mat_create(4,4,3,2.0);
    sum=ip_mat_sum(mat,mat2);
    sub=ip_mat_sub(mat,mat2);
    concat= ip_mat_concat(mat,mat2,1);
    scalare = ip_mat_mul_scalar(mat,5);
    s_scalare = ip_mat_add_scalar(mat,5);

    /*stampa delle variabili*/
    printf("mat1 \n");
    ip_mat_show(mat);
    printf("mat2 \n");
    ip_mat_show(mat2);
    printf("somma mat1 + mat2 \n");
    ip_mat_show(sum);
    printf("sottrazione mat1-mat2 \n");
    ip_mat_show(sub);
    printf("concatenazione mat1 mat2 \n");
    ip_mat_show(concat);
    printf("moltiplicazione di mat per uno scalare c=5 \n");
    ip_mat_show(scalare);
    printf("somma di mat ad uno scalare c=5 \n");
    ip_mat_show(s_scalare);

    /*free della memoria*/
    ip_mat_free(mat);
    ip_mat_free(mat2);
    ip_mat_free(sum);
    ip_mat_free(sub);
    ip_mat_free(concat);
    ip_mat_free(scalare);
    ip_mat_free(s_scalare);
    return 0;
}