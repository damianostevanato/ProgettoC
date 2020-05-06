/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"

/*--------------------PARTE PRIMA---------------------------*/

/* Inizializza una ip_mat con dimensioni h w e k. Ogni elemento è inizializzato a v.
 * Inoltre crea un vettore di stats per contenere le statische sui singoli canali.
 * */
ip_mat * ip_mat_create(unsigned int h, unsigned int w,unsigned  int k, float v){
    ip_mat *t;
    unsigned int i,j,l;

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
    unsigned int i,j;
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
    unsigned int i,j,l;
    float acc;
    for(l=0;l<t->k;l++){
        t->stat[l].max=t->data[0][0][l];
        t->stat[l].min=t->data[0][0][l];
        acc=0;
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
    unsigned int i,j,l;
    for(l=0;l<t->k;l++){
        for(i=0;i<t->h;i++){
            for(j=0;j<t->w;j++){
                float gaus=get_normal_random()*var+mean;
                set_val(t,i,j,l,gaus);
            }
        }
    }
}

/* Crea una copia di una ip_mat e lo restituisce in output */
ip_mat * ip_mat_copy(ip_mat * in){
    ip_mat *out;
    unsigned int i,j,l;
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
    return out;
}

/*restituisce una sottomatrice nxmxk di una struttura dati ip_mat*/
ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end){

    if(row_end > t->h || col_end > t->w || row_start > t->h ||col_start > t->w ){
        printf("errore subset!");
        exit(2);
    }else{
        ip_mat *subset;
        unsigned int i,j,l;
        float val;
        subset = ip_mat_create(row_end-row_start,col_end-col_start,t->k,0.0);
        for(i=row_start;i<row_end;i++){
            for(j=col_start;j<col_end;j++){
                for(l=0;l<t->k;l++){
                    val = get_val(t,i,j,l);
                    set_val(subset,i-row_start,j-col_start,l,val);
                }
            }
        }
        return subset;
    }
}


ip_mat* copy_concat(ip_mat *a, ip_mat *b, int dim){
    ip_mat *out=NULL;
    unsigned int i,j,l;
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
        break;
        
        default:
        printf("Errore ip_mat_concat!!\n");
        exit(6);
        break;
    }
    return out;
}
/* esegue la somma dei valori di ogni canale tra 2 strutture di tipo ip_mat*/
ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b){
    if( (a->h == b->h) && (a->w == b->w) && (a->k == b->k)){
        ip_mat *sum;
        unsigned int i,j,l;
        float somma;
        sum = ip_mat_create(a->h,a->w,a->k,0.0);
        for(i=0;i<sum->h;i++){
            for(j=0;j<sum->w;j++){
                for(l=0;l<sum->k;l++){
                    somma = 0.0;
                    somma = get_val(a,i,j,l)+get_val(b,i,j,l);
                    set_val(sum,i,j,l,somma);
                } 
            }
        }
        compute_stats(sum);
        return sum;
    }else{
        printf("Errore sum!");
        exit(7);
    }
}

/* Esegue la sottrazione di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output.
 * */
ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b){
    if( (a->h == b->h) && (a->w == b->w) && (a->k == b->k)){
        ip_mat *dif;
        unsigned int i,j,l;
        float differenza; 
        dif = ip_mat_create(a->h,a->w,a->k,0.0);
        for(i=0;i<dif->h;i++){
            for(j=0;j<dif->w;j++){
                for(l=0;l<dif->k;l++){
                    differenza = 0.0;
                    differenza = get_val(a,i,j,l)-get_val(b,i,j,l);
                    set_val(dif,i,j,l,differenza);
                } 
            }
        }
        return dif;
    }else{
        printf("Errore sub!");
        exit(8);
    }
}
/*moltiplica ogni elemento di ogni canale per uno scalare C*/
ip_mat * ip_mat_mul_scalar(ip_mat *a, float c){
    ip_mat *x;
    unsigned int i,j,l;
    float molt;
    x = ip_mat_create(a->h,a->w,a->k,0.0);
    for(i=0;i<x->h;i++){
        for(j=0;j<x->w;j++){
            for(l=0;l<x->k;l++){
                molt = get_val(a,i,j,l)*c;
                set_val(x,i,j,l,molt);
            }
        }
    }
    return x;
}
/*somma ad ogni elemento di ogni canale uno scalare C*/
ip_mat *  ip_mat_add_scalar(ip_mat *a, float c){
    ip_mat *x;
    unsigned int i,j,l;
    x = ip_mat_create(a->h,a->w,a->k,0.0);
    for(i=0;i<x->h;i++){
        for(j=0;j<x->w;j++){
            for(l=0;l<x->k;l++){
                set_val(x,i,j,l,(get_val(a,i,j,l)+c));
            }
        }
    }
    return x;
}

/* Calcola la media di due ip_mat a e b e la restituisce in output.*/
ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b){
    ip_mat *out,*somma;
    out=ip_mat_create(a->h,a->w,a->k,0.0);
    somma=ip_mat_create(a->h,a->w,a->k,0.0);
    somma=ip_mat_sum(a,b);
    out=ip_mat_mul_scalar(somma,0.5);
    ip_mat_free(somma);
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

/*-----------------------------PARTE SECONDA---------------------------*/

/*Calcola la media di tutti i pixel su ogni singolo canale tramite compute_stats
e ritorna una nuova immagine che su ogni pixel ha come valore la media appena calcolata 
per ogni corrispettivo canale
*/
ip_mat * ip_mat_to_gray_scale(ip_mat * in){
    ip_mat *gray_scale;
    unsigned int i,j,l;
    float somma,media;
    gray_scale= ip_mat_create(in->h,in->w,in->k,0.0);
    for(i=0;i<in->h;i++){
        for(j=0;j<in->w;j++){
            somma =0.0;
            for(l=0;l<in->k;l++){
                somma = somma+get_val(in,i,j,l);
            }
            media = somma/3.0;
            set_val(gray_scale,i,j,0,media);
            set_val(gray_scale,i,j,1,media);
            set_val(gray_scale,i,j,2,media);
        }
    }
    return gray_scale;
}
/*Aumenta la luminosita di un immagine sommando un certo valore a tutti i pixel
ritorna la nuova immagine con luminosita aumentata
utilizza la funzione ip_mat_add_scalar
*/
ip_mat * ip_mat_brighten(ip_mat * a, float bright){
    return ip_mat_add_scalar(a,bright);
}
/*esegue la sovrapposizione di 2 immagini delle stesse dimensioni hxwxk
formula : out = alpha*A+(1-alpha)*B
nb: alpha dovrebbe essere compreso tra 0-1
nb: a e b devono essere uguali?
*/
ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha){
    ip_mat *a_scalar_alpha,*b_scalar_alpha,*out;
    a_scalar_alpha = ip_mat_mul_scalar(a,alpha);
    b_scalar_alpha = ip_mat_mul_scalar(b,(1-alpha));
    out = ip_mat_sum(a_scalar_alpha,b_scalar_alpha);
    ip_mat_free(a_scalar_alpha);
    ip_mat_free(b_scalar_alpha);
    return out;
}
/*crea una matrice di numeri gaussiano casuali, moltiplica ogni singolo canale di ogni singolo
pixel della matrice gaussiana per amount, il risultato viene sommato alla matrice immagine data
in input
utilizza ip_mat_create,ip_mat_mul_scalar, ip_mat_sum
*/
ip_mat * ip_mat_corrupt(ip_mat * a, float amount){
    ip_mat *random_gaussian,*gaussian_amount,*out;
    random_gaussian = ip_mat_create(a->h,a->w,a->k,0.0);
    ip_mat_init_random(random_gaussian,0,1);
    gaussian_amount = ip_mat_mul_scalar(random_gaussian,amount);
    out = ip_mat_sum(a,gaussian_amount);
    ip_mat_free(random_gaussian);
    ip_mat_free(gaussian_amount);
    return out;
}
/*---------------------------PARTE TERZA----------------------------*/
/*
Controlla se i valori di ogni singolo canale di ogni singolo pixel sono compresi tra low e high
se val>high -> val = high 
se val<low -> val = low
*/
void clamp(ip_mat * t, float low, float high){
    unsigned int i,j,l;
    for(i=0;i<t->h;i++){
        for(j=0;j<t->w;j++){
            for(l=0;l<t->k;l++){
                if(get_val(t,i,j,l)>high){
                    set_val(t,i,j,l,high);
                }else if(get_val(t,i,j,l)<low){
                    set_val(t,i,j,l,low);
                }
            }
        }
    }
}