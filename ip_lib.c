/*
 Created by Sebastiano Vascon on 23/03/20.
*/

/*Comands to compile and link
 *gcc -c ip_lib.c -o ip_lib.o
 *gcc ip_lib.o bmp.o test.c -otest -lm
 */
#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"


/* Inizializza una ip_mat con dimensioni h w e k. Ogni elemento è inizializzato a v.
 * Inoltre crea un vettore di stats per contenere le statische sui singoli canali.
 * */
ip_mat * ip_mat_create(unsigned int h, unsigned int w,unsigned  int k, float v){
    ip_mat *t;
    int i,j,l;

    /*Allocazione della struttura ip_mat*/
    t=(ip_mat*)malloc(sizeof(ip_mat));

    t->w=w;
    t->h=h;
    t->k=k;

    /*creazione della matrice 3D di valori float tramite un array frastagliato a 3 dimensioni*/
    t->data=(float***)malloc(h*sizeof(float**));
    for(i=0;i<h;i++){
        t->data[i]=(float**)malloc(w*sizeof(float*));
        for(j=0;j<w;j++){
            t->data[i][j]=(float*)malloc(k*sizeof(float));
        }
    }

    /*inizializzazione della matrice*/
    for(i=0;i<h;i++){
        for(j=0;j<w;j++){
            for(l=0;l<k;l++){
                set_val(t,i,j,l,v);
            }
        }
    }

    /*creazione e inizializzazione della struttura stat*/
    t->stat=(stats*)malloc(sizeof(stats)*k);

    return t;
}

/* Libera la memoria (data, stat e la struttura) */
void ip_mat_free(ip_mat *a){
    int i,j;
    /*Libero la memoria allocata da stat*/
    free(a->stat);
    for(i=0;i<a->h;i++){
        for(j=0;j<a->w;j++){
            free(a->data[i][j]);
        }
        free(a->data[i]);
    }
    free(a->data);
    free(a);
}

/* Calcola il valore minimo, il massimo e la media per ogni canale
 * e li salva dentro la struttura ip_mat stats
 * */
void compute_stats(ip_mat * t){
    int i,j,l;
    for(l=0;l<t->k;l++){
        t->stat[l].max=t->data[0][0][l];
        t->stat[l].min=t->data[0][0][l];
        float acc=0;
        for(i=0;i<t->h;i++){
            for(j=0;j<t->w;j++){
                if(t->stat[l].max<t->data[i][j][l]){
                    t->stat[l].max=t->data[i][j][l];
                }
                if(t->stat[l].min>t->data[i][j][l]){
                    t->stat[l].min=t->data[i][j][l];
                }
                acc+=t->data[i][j][l];
            }
        }
        t->stat[l].mean=acc/((t->w)*(t->h));
    }
}

/* Inizializza una ip_mat con dimensioni w h e k.
 * Ogni elemento è generato da una gaussiana con media mean e varianza var */
void ip_mat_init_random(ip_mat * t, float mean, float var){
    int i,j,l;
    for(l=0;l<t->k;l++){
        for(i=0;i<t->h;i++){
            for(j=0;j<t->w;j++){
                float gaus=get_normal_random()*var+mean;
                set_val(t,i,j,l,gaus);
            }
        }
    }
    compute_stats(t);
}

/* Crea una copia di una ip_mat e lo restituisce in output */
ip_mat * ip_mat_copy(ip_mat * in){
    ip_mat *out;
    int i,j,l;
    float val;
    out=ip_mat_create(in->h,in->w,in->k,0.0);
    for(i=0;i<in->h;i++){
        for(j=0;j<in->w;j++){
            for(l=0;l<in->k;l++){
                val=get_val(in,i,j,l);
                set_val(out,i,j,l,val);
            }
        }
    }
    compute_stats(out);
    return out;
}

ip_mat* copy_concat(ip_mat *a, ip_mat *b, int dim){
    ip_mat *out=NULL;
    int i,j,l;
    float vala,valb;
    switch(dim){
        case 0:
        out=ip_mat_create(a->h+b->h,a->w,a->k,0.0);
        for(i=0;i<a->h;i++){
            for(j=0;j<a->w;j++){
                for(l=0;l<a->k;l++){
                    vala=get_val(a,i,j,l);
                    set_val(out,i,j,l,vala);
                }
            }
        }
        for(i=a->h;i<b->h+a->h;i++){
            for(j=0;j<b->w;j++){
                for(l=0;l<b->k;l++){
                    valb=get_val(b,i-a->h,j,l);
                    set_val(out,i,j,l,valb);
                }
            }
        }
        break;

        case 1:
        out=ip_mat_create(a->h,a->w+b->w,a->k,0.0);
        for(i=0;i<a->h;i++){
            for(j=0;j<a->w;j++){
                for(l=0;l<a->k;l++){
                    vala=get_val(a,i,j,l);
                    set_val(out,i,j,l,vala);
                }
            }
        }
        for(i=0;i<b->h;i++){
            for(j=a->w;j<b->w+a->w;j++){
                for(l=0;l<b->k;l++){
                    valb=get_val(b,i,j-a->w,l);
                    set_val(out,i,j,l,valb);
                }
            }
        }
        break;

        case 2:
        out=ip_mat_create(a->h,a->w,a->k+b->k,0.0);
        for(i=0;i<a->h;i++){
            for(j=0;j<a->w;j++){
                for(l=0;l<a->k;l++){
                    vala=get_val(a,i,j,l);
                    set_val(out,i,j,l,vala);
                }
            }
        }
        for(i=0;i<b->h;i++){
            for(j=0;j<b->w;j++){
                for(l=a->k;l<b->k+a->k;l++){
                    valb=get_val(b,i,j,l-a->k);
                    set_val(out,i,j,l,valb);
                }
            }
        }
        break;
    }
    return out;
}

/*restituisce una sottomatrice nxmxk di una struttura dati ip_mat*/
ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end){

    if(row_end > t->h || col_end > t->w || row_start > t->h ||col_start > t->w ){
        printf("errore subset!");
        exit(2);
    }else{
        ip_mat *subset;
        subset = ip_mat_create(row_end-row_start,col_end-col_start,t->k,0.0);
        int i,j,l,val;
        for(i=row_start;i<row_end;i++){
            for(j=col_start;j<col_end;j++){
                for(l=0;l<t->k;l++){
                    val = get_val(t,i,j,l);
                    set_val(subset,i-row_start,j-col_start,l,val);
                }
            }
        }
        compute_stats(subset);
        return subset;
    }
}


/* Concatena due ip_mat su una certa dimensione.
 * Ad esempio:
 * ip_mat_concat(ip_mat * a, ip_mat * b, 0);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h + b.h
 *      out.w = a.w = b.w
 *      out.k = a.k = b.k
 *
 * ip_mat_concat(ip_mat * a, ip_mat * b, 1);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h = b.h
 *      out.w = a.w + b.w
 *      out.k = a.k = b.k
 *
 * ip_mat_concat(ip_mat * a, ip_mat * b, 2);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h = b.h
 *      out.w = a.w = b.w
 *      out.k = a.k + b.k
 * */
ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione){
    ip_mat *out=NULL;
    switch(dimensione){
        case 0:
        if(a->w!=b->w||a->k!=b->k){
            printf("Errore ip_mat_concat!!\n");
            exit(3);
        }
        else{
            out=copy_concat(a,b,dimensione);
        }
        break;

        case 1:
        if(a->h!=b->h || a->k!=b->k){
            printf("Errore ip_mat_concat!!\n");
            exit(4);
        }
        else{
            out=copy_concat(a,b,dimensione);
        }
        break;

        case 2:
        if(a->w!=b->w || a->h!=b->h){
            printf("Errore ip_mat_concat!!\n");
            exit(5);
        }
        else{
            out=copy_concat(a,b,dimensione);
        }
        
        default:
        printf("Errore ip_mat_concat!!\n");
        exit(7);
        break;
    }
    return out;
}

void ip_mat_show(ip_mat * t){
    unsigned int r,l,c;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->w,t->h,t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(r=0;r<t->h;r++) {
            for (c = 0; c < t->w; c++) {
                printf("%f ", get_val(t,r,c,l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_show_stats(ip_mat * t){
    unsigned int k;

    compute_stats(t);

    for(k=0;k<t->k;k++){
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}

ip_mat * bitmap_to_ip_mat(Bitmap * img){
    unsigned int i=0,j=0;

    unsigned char R,G,B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }

    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t){

    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            bm_set_pixel(b, j,i, (unsigned char) get_val(t,i,j,0),
                    (unsigned char) get_val(t,i,j,1),
                    (unsigned char) get_val(t,i,j,2));
        }
    }
    return b;
}


float get_val(ip_mat * a, unsigned int r,unsigned int c,unsigned int k){
    if(r<a->h && c<a->w && k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[r][c][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w && k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(){
    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    return cos(2*PI*y2)*sqrt(-2.*log(y1));

}