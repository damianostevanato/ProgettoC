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
    ip_mat *mat,*mat2,*concat,*scalare,*s_scalare;
    mat=ip_mat_create(4,4,3,1.0);
    mat2=ip_mat_create(4,4,3,2.0);
    concat=ip_mat_sub(mat,mat2);
    scalare = ip_mat_mul_scalar(mat,5);
    s_scalare = ip_mat_add_scalar(mat,5);
    ip_mat_show(scalare);
    ip_mat_show(s_scalare);
    ip_mat_free(mat);
    ip_mat_free(s_scalare);
    ip_mat_free(mat2);
    ip_mat_free(concat);
    ip_mat_free(scalare);
    return 0;
}