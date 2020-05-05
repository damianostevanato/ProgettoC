#include <stdio.h>
#include "bmp.h"
#include <math.h>
#include <string.h>
#include "ip_lib.h"

int main(){
    ip_mat *mat,*mat2,*concat;
    mat=ip_mat_create(4,4,3,1.0);
    mat2=ip_mat_create(4,4,1,2.0);
    concat=ip_mat_concat(mat,mat2,2);
    ip_mat_show(concat);
    ip_mat_free(mat);
    ip_mat_free(mat2);
    ip_mat_free(concat);
    return 0;
}