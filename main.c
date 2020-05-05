#include <stdio.h>
#include "bmp.h"
#include <math.h>
#include <string.h>
#include "ip_lib.h"

int main(){
    ip_mat *mat,*mat2,*concat,*somma;
    mat=ip_mat_create(4,4,3,1.0);
    mat2=ip_mat_create(4,4,3,2.0);
    concat=ip_mat_concat(mat,mat2,0);
    somma = ip_mat_sum(mat,mat2);
    ip_mat_show(mat);
    printf("----------------");
    ip_mat_show(mat2);
    printf("----------------");
    ip_mat_show(somma);
    printf("----------------");
    ip_mat_free(mat);
    ip_mat_free(mat2);
    ip_mat_free(somma);
    ip_mat_free(concat);
    return 0;
}