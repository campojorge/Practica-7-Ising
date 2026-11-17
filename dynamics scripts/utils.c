#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "Random_utils.h"
#include <math.h>


void vecinos(int sitio, Parametros p, int v[]){

    int x,y;

    y = sitio / p.L;
    x = sitio % p.L;

    /* Condiciones periódicas */
    v[0] = ((y + p.L - 1) % p.L) * p.L + x; //Arriba
    v[1] = ((y + 1) % p.L) * p.L + x; //Abajo
    v[2] = y * p.L + (x + p.L - 1) % p.L; //Izquierda
    v[3] = y * p.L + (x + 1) % p.L; //Derecha
}

double delta_energia(Parametros p, int red[], int sitio){

    int v[4];
    int arriba,abajo,derecha,izquierda,suma_vecinos,i;
    double delta;

    vecinos(sitio,p,v);
    suma_vecinos=0;

    for(i=0;i<4;i++){
        suma_vecinos+=red[v[i]];
    }

    delta=2.0 * red[sitio] * (p.J * suma_vecinos + p.H);
    return delta;
}

double delta_energia2(Parametros p, int clase){
    double delta;

    if(clase<5){
        if(clase == 0){
            delta=2.0  * (p.J *4  + p.H);
        }
        else if(clase == 1){
            delta=2.0  * (p.J *2  + p.H);
        }
        else if(clase == 2){
            delta=2.0  * (p.J *(0)  + p.H);
        }
        else if(clase == 3){
            delta=2.0  * (p.J *(-2)  + p.H);
        }
        else if(clase == 4){
            delta=2.0  * (p.J *(-4)  + p.H);
        }
    }
    else{
        if(clase == 5){
            delta=-2.0  * (p.J *4  + p.H);
        }
        else if(clase == 6){
            delta=-2.0  * (p.J *2  + p.H);
        }
        else if(clase == 7){
            delta=-2.0  * (p.J *(0)  + p.H);
        }
        else if(clase == 8){
            delta=-2.0  * (p.J *(-2)  + p.H);
        }
        else if(clase == 9){
            delta=-2.0  * (p.J *(-4)  + p.H);
        }
    }

    return delta;
}

   

double probabilidad_glauber(double delta_E, Parametros p){
    double prob;
    //Tratemos primero los casos especiales
    //1-Si ambas configuraicones tienen la misma energia, entonces la probabilidad de inversion es 1/2
    if (delta_E == 0.0) {
        prob=0.5;
    }
    //2-Para temperatura 0 un cambio que disminuye la energia se acepta siempre y uno que la aumenta se rechaza siempre
    if (p.kBT == 0.0) {
        if(delta_E<0.0){
            prob=1.0;
        }
        else{
            prob=0.0;
        }
    }

    //Caso general
   
    if (delta_E > 0.0 && p.kBT>0.0) {
        prob=1/(1+exp(delta_E/p.kBT));
    }

    return prob;
}

void intento_glauber(Parametros p, int red[]){

    int sitio;
    double delta_E,prob;

    //Elegimos un espin cualquiera
    sitio = (int)(Random() * p.L * p.L);


    //Calculamos la diferencia de energía y con ello la probabilidad
    delta_E = delta_energia(p, red);
    prob = probabilidad_glauber(delta_E, p);

    if (Random() < prob) {
        red[sitio] = -red[sitio];
    }
}

double magnetizacion(Parametros p,int red[]){
    int i;
    double M;
    M=0;
    for(i=0;i<(p.L*p.L);i++){
        M+=red[i];
    }

    return M;
}

int clase_sitio(Parametros p,int red[],int sitio){
    int v[4];
    int n_arriba,i,clase;

    n_arriba=0;

    vecinos(sitio, p, v);

    for (i = 0; i < 4; i++) {
        if (red[v[i]] == 1) {
            n_arriba++;
        }
    }

    /*
    Clases 0..4:
        espín +1 con 4,3,2,1,0 vecinos arriba.

    Clases 5..9:
        espín -1 con 4,3,2,1,0 vecinos arriba.
    */

    if(red[sitio] == 1){
        clase=n_arriba;

    }
    else if(red[sitio] == -1){
        clase=(5-n_arriba)+4;
    }

    return clase;
}


/*  MANEJO DE LAS LISTAS */

void insertar(Nfold *nf, int sitio, int c){
    int pos;
    pos = nf->cantidad[c];

    nf->lista[c][pos] = sitio;
    nf->clase[sitio] = c;
    nf->posicion[sitio] = pos;
    nf->cantidad[c]++;
}

void eliminar(Nfold *nf, int sitio){
    int c,pos,ultimo;

    c = nf->clase[sitio];
    pos = nf->posicion[sitio];

    /*
    Rellenamos el hueco con el último elemento
    de la misma clase.
    */
    nf->cantidad[c]--;

    ultimo = nf->lista[c][nf->cantidad[c]];

    nf->lista[c][pos] = ultimo;
    nf->posicion[ultimo] = pos;
}


/*INICIALIZACIÓN*/

void inicializar_nfold(Parametros p,Nfold *nf,int red[]){

    int c,sitio,s,clase;
    double delta_E;


    // Número total de espines
    nf->N = p.L * p.L;

    // Preparar las diez clases
    for (c = 0; c < 10; c++){
        nf->cantidad[c] = 0;
        delta_E=delta_energia2(p,c);
        nf->prob[c] = probabilidad_glauber(delta_E, p);
    }

    // Introducir cada sitio de la red en su clase
    for (sitio = 0; sitio < nf->N; sitio++) {
        clase = clase_sitio(p, red, sitio);
        insertar(nf, sitio, clase);
    }
}



void paso_nfold(const Parametros p,int red[],Nfold *nf,double (*uniforme)(void),uint64_t *intentos,int *delta_M){
    int c,c_elegida,i,j,ind,sitio,n_afectados,repetido;
    int v[4];
    int afectados[5];
    double peso[10];
    double W,Q,m,objetivo,acumulada;

    W=0.0;
    *intentos = 0;
    *delta_M = 0;

    /* Peso de cada clase: n_c * p_c */
    for (c = 0; c < 10; c++) {
        peso[c] = nf->cantidad[c] * nf->prob[c];
        W += peso[c];
    }


    /* Probabilidad total de aceptar un intento */
    Q = W / nf->N;


    /*
    Tiempo de espera en intentos:

        m = floor(log(r) / log(1-Q)) + 1

    Se toma r = 1-u, con u en [0,1).
    */

    if (Q >= 1.0) {
        m = 1.0;
    } 
    else {
        m = floor(log1p(-uniforme()) / log1p(-Q)) + 1.0;
    }


    /* Elegir la clase */
    objetivo = uniforme() * W;
    acumulada = 0.0;

    c_elegida = -1;

    for (c = 0; c < 10; c++) {
        if (peso[c] > 0.0) {
            c_elegida = c;
        }
        acumulada += peso[c];
        if (objetivo < acumulada) {
            break;
        }
    }

    /* Elegir un espín dentro de esa clase */
    ind = (int)(uniforme() * nf->cantidad[c_elegida]);

    sitio = nf->lista[c_elegida][ind];

    n_afectados = 1;
    afectados[0] = sitio;

    vecinos(sitio, p->L, v);

    for (i = 0; i < 4; i++) {
        repetido = 0;

        for (int j = 0; j < n_afectados; ++j) {
            if (afectados[j] == v[i]) {
                repetido = 1;
            }
        }

        /* Para L=2 algunos vecinos aparecen repetidos */
        if (!repetido) {
            afectados[n_afectados] = v[i];
            n_afectados++;
        }
    }

    /* Sacar los sitios afectados de sus clases */
    for (i = 0; i < n_afectados; i++) {
        eliminar(nf, afectados[i]);
    }

    /* Girar el espín seleccionado */


    red[sitio] = -red[sitio];

    *delta_M = 2 * red[sitio];

    /* Introducir los sitios en sus nuevas clases */
    for (int i = 0; i < n_afectados; ++i) {
        int c = clase_sitio(p, red, afectados[i]);

        insertar(nf, afectados[i], c);
    }

    *intentos = (uint64_t)m;
}
