/*
 Created by Sebastiano Vascon on 23/03/20.
*/
#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"
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
    compute_stats(t);
    return t;
}

/* Libera la memoria (data, stat e la struttura) */
void ip_mat_free(ip_mat *a){
    unsigned int i,j;
    if(a!=NULL){
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
                float gaus=get_normal_random(mean,var);
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
    compute_stats(out);
    return out;
}

/* Restituisce una sotto-matrice, ovvero la porzione individuata da:
 * t->data[row_start...row_end][col_start...col_end][0...k]
 * La terza dimensione la riportiamo per intero, stiamo in sostanza prendendo un sottoinsieme
 * delle righe e delle colonne.
 * */
ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end){

    if(row_end > t->h || col_end > t->w || row_start > t->h ||col_start > t->w ){
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
    compute_stats(out);
    return out;
}

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
    compute_stats(out);
    return out;
}




/**** PARTE 1: OPERAZIONI MATEMATICHE FRA IP_MAT ****/
/* Esegue la somma di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output. */
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
        compute_stats(dif);
        return dif;
    }else{
        printf("Errore sub!");
        exit(8);
    }
}

/* Moltiplica un ip_mat per uno scalare c. Si moltiplica c per tutti gli elementi di "a"
 * e si salva il risultato in un nuovo tensore in output. */
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
    compute_stats(x);
    return x;
}

/* Aggiunge ad un ip_mat uno scalare c e lo restituisce in un nuovo tensore in output. */
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
    compute_stats(x);
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
    compute_stats(out);
    return out;
}


/**** PARTE 2: SEMPLICI OPERAZIONI SU IMMAGINI ****/

/* Operazione di corruzione con rumore gaussiano:
 * Aggiunge del rumore gaussiano all'immagine, il rumore viene enfatizzato
 * per mezzo della variabile amount.
 * out = a + gauss_noise*amount
 * */

ip_mat * ip_mat_corrupt(ip_mat * a, float amount){
    unsigned int i,j,l;
    ip_mat *out;
    out = ip_mat_create(a->h,a->w,a->k,0.0);
    for(l=0;l<out->k;l++){
        for(i=0;i<out->h;i++){
            for(j=0;j<out->w;j++){
                float gauss = get_normal_random(0,2*0.4);
                set_val(out,i,j,l,get_val(a,i,j,l)+gauss*amount);
                printf("%f\n",gauss);
            }
        }
    }
    clamp(out,0,255);
    compute_stats(out);
    return out;
}

/* Converte un'immagine RGB ad una immagine a scala di grigio.
 * Quest'operazione viene fatta calcolando la media per ogni pixel sui 3 canali
 * e creando una nuova immagine avente per valore di un pixel su ogni canale la media appena calcolata.
 * Avremo quindi che tutti i canali saranno uguali.
 * */
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
    compute_stats(gray_scale);
    return gray_scale;
}

/* Operazione di brightening: aumenta la luminosità dell'immagine
 * aggiunge ad ogni pixel un certo valore*/
ip_mat * ip_mat_brighten(ip_mat * a, float bright){
    ip_mat *out;
    out=ip_mat_add_scalar(a,bright);
    compute_stats(out);
    return out;
}

/* Effettua la fusione (combinazione convessa) di due immagini */
ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha){
    ip_mat *a_scalar_alpha,*b_scalar_alpha,*out;
    a_scalar_alpha = ip_mat_mul_scalar(a,alpha);
    b_scalar_alpha = ip_mat_mul_scalar(b,(1-alpha));
    out = ip_mat_sum(a_scalar_alpha,b_scalar_alpha);
    ip_mat_free(a_scalar_alpha);
    ip_mat_free(b_scalar_alpha);
    compute_stats(out);
    return out;
}


/**** PARTE 3: CONVOLUZIONE E FILTRI *****/

/*moltiplica un layer di una matrice per un layer del filtro*/
float prod_mat(ip_mat *a, ip_mat *k,unsigned int layer){
    if(a->h!=k->h || a->w!=k->w){
        printf("Errore prod_mat\n");
        exit(8);
    }
    else{
        float acc=0.;
        unsigned int i,j;
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
    int pad_h,pad_w;
    unsigned int i,j,l;
    ip_mat *padded=NULL,*out=NULL,*sub=NULL;
    pad_h=(f->h-1)/2;
    pad_w=(f->w-1)/2;
    out=ip_mat_create(a->h,a->w,a->k,0.0);
    padded=ip_mat_padding(a,pad_h,pad_w);    
    for(i=0;i<out->h;i++){
        for(j=0;j<out->w;j++){
            for(l=0;l<out->k;l++){
                sub = ip_mat_subset(padded,i,i+f->h,j,j+f->w);
                float val=prod_mat(sub,f,l);
                ip_mat_free(sub);
                set_val(out,i,j,l,val);
            }
        }
    }
    ip_mat_free(padded);
    compute_stats(out);
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
ip_mat * ip_mat_padding(ip_mat * a, unsigned int pad_h, unsigned int pad_w){
    ip_mat *out;
    unsigned int i,j,l;
    out=ip_mat_create(a->h+(pad_h*2),a->w+(pad_w*2),a->k,0.0);
    for(i=pad_h;i<a->h+pad_h;i++){
        for(j=pad_w;j<a->w+pad_w;j++){
            for(l=0;l<a->k;l++){
                set_val(out,i,j,l,get_val(a,i-pad_h,j-pad_w,l));
            }
        }
    }
    compute_stats(out);
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
    compute_stats(filter);
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
    compute_stats(filter);
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
ip_mat * create_average_filter(unsigned int w, unsigned int h, unsigned int k){
        if(w%2==0||h%2==0){
        printf("Errore convoluzione!!");
        exit(9);
    }
    else{
        float val=1./(w*h);
        unsigned int i,j,l;
        ip_mat *filter=ip_mat_create(h,w,k,0.);
        for(i=0;i<filter->h;i++){
            for(j=0;j<filter->w;j++){
                for(l=0;l<filter->k;l++){
                    set_val(filter,i,j,l,val);
                }
            }
        }
        compute_stats(filter);
        return filter;
    }

}

/* Crea un filtro gaussiano per la rimozione del rumore */
ip_mat * create_gaussian_filter(unsigned int w, unsigned int h, unsigned int k, float sigma){
    int cx,cy;
    unsigned int i,j,l;
    float x=0,y=0,sum=0.;
    ip_mat *filter=ip_mat_create(h,w,k,0.0);
    cx=(w-1)/2;
    cy=(h-1)/2;
    for(l=0;l<filter->k;l++){
        sum=0.;
        for(i=0;i<filter->h;i++){
            for(j=0;j<filter->w;j++){
                    float val=(1/(2*PI*(sigma*sigma)))*exp(-(((x*x)+(y*y))/(2*sigma*sigma)));
                    x=i-cx;
                    y=j-cy;
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
    compute_stats(filter);
    return filter;
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
    unsigned int i,j,l;
    compute_stats(t);
    for(l=0;l<t->k;l++){
        for(i=0;i<t->h;i++){
            for(j=0;j<t->w;j++){
                set_val(t,i,j,l,(get_val(t,i,j,l)-t->stat[l].min)/(t->stat[l].max-t->stat[l].min));
                set_val(t,i,j,l,get_val(t,i,j,l)*new_max);
            }
        }
    }
}

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

void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0;i<t->h;i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", get_val(t,i,j,l));
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

    compute_stats(out);

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

float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(float media, float std){

    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float num = cos(2*PI*y2)*sqrt(-2.*log(y1));

    return media + num*std;
}
