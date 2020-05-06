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

/*restituisce una sottomatrice nxmxk di una struttura dati ip_mat*/
ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end){

    if(row_end > t->h || col_end > t->w || row_start > t->h ||col_start > t->w ){
        //printf("errore subset! %d %d %d %d | %d,%d",row_start,row_end,col_start,col_end,t->h,t->w);
        exit(2);
    }else{
        ip_mat *subset;
        subset = ip_mat_create(row_end-row_start,col_end-col_start,t->k,0.0);
        int i,j,l;
        float val;
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
        sum = ip_mat_create(a->h,a->w,a->k,0.0);
        int i,j,l;
        float somma;
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
        dif = ip_mat_create(a->h,a->w,a->k,0.0);
        int i,j,l;
        float differenza;
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
    x = ip_mat_create(a->h,a->w,a->k,0.0);
    int i,j,l;
    float molt;
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
    x = ip_mat_create(a->h,a->w,a->k,0.0);
    int i,j,l;
    float somma;
    for(i=0;i<x->h;i++){
        for(j=0;j<x->w;j++){
            for(l=0;l<x->k;l++){
                somma = get_val(a,i,j,l)+c;
                set_val(x,i,j,l,somma);
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


/*---------------------------PARTE TERZA----------------------------*/

float prod_mat(ip_mat *a, ip_mat *k,int layer){
    if(a->h!=k->h || a->w!=k->w){
        printf("Errore prod_mat\n");
        exit(8);
    }
    else{
        float acc=0.;
        int i,j;
        if(k->k==1){
            for(i=0;i<a->h;i++){
                for(j=0;j<a->w;j++){
                    acc=acc+(a->data[i][j][layer]*k->data[i][j][0]);
                }
            }
        }
        else{
            if(a->k==k->k){
                for(i=0;i<a->h;i++){
                    for(j=0;j<a->w;j++){
                        acc=acc+(a->data[i][j][layer]*k->data[i][j][layer]);
                    }
                }
            }
            else{
                printf("Errore prod_mat\n");
                exit(9);
            }
        }
        return acc;
    }
}

/* Effettua la convoluzione di un ip_mat "a" con un ip_mat "f".
 * La funzione restituisce un ip_mat delle stesse dimensioni di "a".
 * */
ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f){
    int pad_h,pad_w,i,j,l;
    ip_mat *padded=NULL,*out=NULL;
    pad_h=(f->h-1)/2;
    pad_w=(f->w-1)/2;
    out=ip_mat_create(a->h,a->w,a->k,0.0);
    padded=ip_mat_padding(a,pad_h,pad_w);    
    for(i=0;i<out->h;i++){
        for(j=0;j<out->w;j++){
            for(l=0;l<out->k;l++){
                float val=prod_mat(ip_mat_subset(padded,i,i+f->h,j,j+f->w),f,l);
                set_val(out,i,j,l,val);
            }
        }
    }
    clamp(out,0.,255.);
    return out;
}

/* Aggiunge un padding all'immagine. Il padding verticale è pad_h mentre quello
 * orizzontale è pad_w.
 * L'output sarà un'immagine di dimensioni:
 *      out.h = a.h + 2*pad_h;
 *      out.w = a.w + 2*pad_w;
 *      out.k = a.k
 * con valori nulli sui bordi corrispondenti al padding e l'immagine "a" riportata
 * nel centro
 * */
ip_mat * ip_mat_padding(ip_mat * a, int pad_h, int pad_w){
    ip_mat *out;
    int i,j,l;
    out=ip_mat_create(a->h+(pad_h*2),a->w+(pad_w*2),a->k,0.0);
    for(i=pad_h;i<a->h+pad_h;i++){
        for(j=pad_w;j<a->w+pad_w;j++){
            for(l=0;l<a->k;l++){
                set_val(out,i,j,l,get_val(a,i-pad_h,j-pad_w,l));
            }
        }
    }
    return out;
}

/* Crea un filtro di sharpening */
ip_mat * create_sharpen_filter(){
    ip_mat *filter=ip_mat_create(3,3,1,0.);
    set_val(filter,0,0,0,0.);
    set_val(filter,0,1,0,-1.);
    set_val(filter,0,2,0,0.);
    set_val(filter,1,0,0,-1.);
    set_val(filter,1,1,0,5.);
    set_val(filter,1,2,0,-1.);
    set_val(filter,2,0,0,0.);
    set_val(filter,2,1,0,-1.);
    set_val(filter,2,2,0,0.);
    return filter;
}

/* Crea un filtro per rilevare i bordi */
ip_mat * create_edge_filter(){
    ip_mat *filter=ip_mat_create(3,3,1,0.);
    set_val(filter,0,0,0,-1.);
    set_val(filter,0,1,0,-1.);
    set_val(filter,0,2,0,-1.);
    set_val(filter,1,0,0,-1.);
    set_val(filter,1,1,0,8.);
    set_val(filter,1,2,0,-1.);
    set_val(filter,2,0,0,-1.);
    set_val(filter,2,1,0,-1.);
    set_val(filter,2,2,0,-1.);
    return filter;
}

/* Crea un filtro per aggiungere profondità */
ip_mat * create_emboss_filter(){
    ip_mat *filter=ip_mat_create(3,3,1,0.);
    set_val(filter,0,0,0,-2.);
    set_val(filter,0,1,0,-1.);
    set_val(filter,0,2,0,0.);
    set_val(filter,1,0,0,-1.);
    set_val(filter,1,1,0,1.);
    set_val(filter,1,2,0,1.);
    set_val(filter,2,0,0,0.);
    set_val(filter,2,1,0,1.);
    set_val(filter,2,2,0,2.);
    return filter;
}

/* Crea un filtro medio per la rimozione del rumore */
ip_mat * create_average_filter(int w, int h, int k){
        if(w%2==0||h%2==0){
        printf("Errore convoluzione!!");
        exit(9);
    }
    else{
        ip_mat *filter=ip_mat_create(h,w,k,0.);
        float val=1./(w*h);
        int i,j,l;
        for(i=0;i<filter->h;i++){
            for(j=0;j<filter->w;j++){
                for(l=0;l<filter->k;l++){
                    set_val(filter,i,j,l,val);
                }
            }
        }
        return filter;
    }

}

/* Crea un filtro gaussiano per la rimozione del rumore */
ip_mat * create_gaussian_filter(int w, int h, int k, float sigma){
    ip_mat *filter=ip_mat_create(h,w,k,0.0);
    int cx,cy;
    cx=(w-1)/2;
    cy=(h-1)/2;
    int i,j,l;
    float x,y,sum=0.;
    for(l=0;l<filter->k;l++){
        sum=0.;
        for(i=0;i<filter->h;i++){
            for(j=0;j<filter->w;j++){
                    x=i-cx;
                    y=j-cy;
                    float val=(1/(2*PI*(sigma*sigma)))*pow(E_NEPERO,(-(((x*x)+(y*y))/(2*sigma*sigma))));
                    set_val(filter,i,j,l,val);
                    sum+=val;
            }
        }
        for(i=0;i<filter->h;i++){
            for(j=0;j<filter->w;j++){
                float val=get_val(filter,i,j,l)/sum;
                set_val(filter,i,j,l,val);
            }
        }
    }
    return filter;
}

/* Nell'operazione di clamping i valori <low si convertono in low e i valori >high in high.*/
void clamp(ip_mat * t, float low, float high){
    int i,j,l;
    for(l=0;l<t->k;l++){
        for(i=0;i<t->h;i++){
            for(j=0;j<t->w;j++){
                if(get_val(t,i,j,l)<0.)
                    set_val(t,i,j,l,0.);
                if(get_val(t,i,j,l)>255.)
                    set_val(t,i,j,l,255.);
            }
        }
    }
}

/* Effettua una riscalatura dei dati tale che i valori siano in [0,new_max].
 * Utilizzate il metodo compute_stat per ricavarvi il min, max per ogni canale.
 *
 * I valori sono scalati tramite la formula valore-min/(max - min)
 *
 * Si considera ogni indice della terza dimensione indipendente, quindi l'operazione
 * di scalatura va ripetuta per ogni "fetta" della matrice 3D.
 * Successivamente moltiplichiamo per new_max gli elementi della matrice in modo da ottenere un range
 * di valori in [0,new_max].
 * */
void rescale(ip_mat * t, float new_max){
    compute_stats(t);
    int i,j,l;
    for(l=0;l<t->k;l++){
        for(i=0;i<t->h;i++){
            for(j=0;j<t->w;j++){
                set_val(t,i,j,l,(get_val(t,i,j,l)-t->stat[l].min)/(t->stat[l].max-t->stat[l].min));
                set_val(t,i,j,l,get_val(t,i,j,l)*new_max);
            }
        }
    }    
}