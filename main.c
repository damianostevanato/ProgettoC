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
    /*---------------PROVA PARTE 1---------------*/
    /*dichiarazione veriabili*/
    /*ip_mat *mat,*mat2,*concat,*sum,*sub,*scalare,*s_scalare,*mean;*/

    /*operazioni su variabili*/
    /*mat=ip_mat_create(4,4,3,1.0);
    mat2=ip_mat_create(4,4,3,2.0);
    sum=ip_mat_sum(mat,mat2);
    sub=ip_mat_sub(mat,mat2);
    concat= ip_mat_concat(mat,mat2,1);
    scalare = ip_mat_mul_scalar(mat,5);
    s_scalare = ip_mat_add_scalar(mat,5);
    mean=ip_mat_mean(mat,mat2);*/

    /*stampa delle variabili*/
    /*printf("mat1 \n");
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
    printf("mean\n:");
    ip_mat_show(mean);*/

    /*free della memoria*/
    /*ip_mat_free(mat);
    ip_mat_free(s_scalare);
    ip_mat_free(mat2);
    ip_mat_free(sum);
    ip_mat_free(sub);
    ip_mat_free(concat);
    ip_mat_free(scalare);
    ip_mat_free(mean);*/
    /*-----------PROVA IMMAGINE-----------------*/
    /*char *file="flower.bmp";
    ip_mat *t=NULL;  
    ip_mat *subset=NULL;  
    Bitmap *img=NULL;
    img=bm_load(file);
    t=bitmap_to_ip_mat(img);
    subset=ip_mat_subset(t,0,6,0,6);
    ip_mat_show(subset);
    ip_mat_free(t);
    bm_free(img);*/

    /*----------------------PROVA SECONDA PARTE----------------------*/
    char *flower = "flower2.bmp";
    char *mongolfiera= "mongolfiere.bmp";
    char *c_flower = "corrupted_flower.bmp";
    char *b_flower = "brighten_flower.bmp";
    char *bl_flower = "blend_flower.bmp";
    char *g_flower = "gray_flower.bmp";

    ip_mat *mongolfiera_mat,*flower_mat,*corrupted_flower,*bright_flower,*blend,*gray_scale;
    Bitmap *flower_bmp,*mongolfiera_bmp;

    flower_bmp = bm_load(flower);
    flower_mat = bitmap_to_ip_mat(flower_bmp);
    bm_free(flower_bmp);
    mongolfiera_bmp = bm_load(mongolfiera);
    mongolfiera_mat = bitmap_to_ip_mat(mongolfiera_bmp);
    bm_free(mongolfiera_bmp);

    corrupted_flower = ip_mat_corrupt(flower_mat,255);
    clamp(corrupted_flower,0,255);
    bright_flower = ip_mat_brighten(flower_mat,100);
    clamp(bright_flower,0,255);
    gray_scale = ip_mat_to_gray_scale(flower_mat);
    blend = ip_mat_blend(flower_mat,mongolfiera_mat,0.5);
  
    flower_bmp = ip_mat_to_bitmap(corrupted_flower);
    bm_save(flower_bmp,c_flower);
    bm_free(flower_bmp);
    flower_bmp = ip_mat_to_bitmap(bright_flower);
    bm_save(flower_bmp,b_flower);
    bm_free(flower_bmp);
    flower_bmp = ip_mat_to_bitmap(gray_scale);
    bm_save(flower_bmp,g_flower);
    bm_free(flower_bmp);
    flower_bmp = ip_mat_to_bitmap(blend);
    bm_save(flower_bmp,bl_flower);
    bm_free(flower_bmp);

    ip_mat_free(mongolfiera_mat);
    ip_mat_free(flower_mat);
    ip_mat_free(corrupted_flower);
    ip_mat_free(bright_flower);
    ip_mat_free(blend);
    ip_mat_free(gray_scale);
   
    return 0;
}