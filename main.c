#include <stdio.h>
#include "bmp.h"
#include <math.h>
#include <string.h>
#include "ip_lib.h"

int main(){
    ip_mat *mat,*mat2,*concat,*scalare;
    mat=ip_mat_create(4,4,3,1.0);
    mat2=ip_mat_create(4,4,3,2.0);
    concat=ip_mat_sub(mat,mat2);
    scalare = ip_mat_mul_scalar(mat,5);

    ip_mat_show(scalare);
    ip_mat_free(mat);
    ip_mat_free(mat2);
    ip_mat_free(concat);
    ip_mat_free(scalare);
    return 0;
}