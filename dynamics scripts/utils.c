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
   
    if ( p.kBT>0.0) {
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
    delta_E = delta_energia(p, red,sitio);
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

    
    /*Clases 0..4:
        espín +1 con 4,3,2,1,0 vecinos arriba.

    Clases 5..9:
        espín -1 con 4,3,2,1,0 vecinos arriba.*/
    

    if(red[sitio] == 1){
        clase=4 - n_arriba;
    }
    else if(red[sitio] == -1){
        clase=(5-n_arriba)+4;
    }

    return clase;
}


//  MANEJO DE LAS LISTAS 

void insertar(Nfold *nf, int sitio, int c){
    int pos;
    pos = nf->NNCLS[c];

    nf->LOC[c][pos] = sitio;
    nf->LOOK[sitio] = c;
    nf->LOOC[sitio] = pos;
    nf->NNCLS[c]++;
}

void eliminar(Nfold *nf, int sitio){
    int c,pos,ultimo;

    c = nf->LOOK[sitio];
    pos = nf->LOOC[sitio];

    /*
    Rellenamos el hueco con el último elemento
    de la misma clase.*/
    
    nf->NNCLS[c]--;

    ultimo = nf->LOC[c][nf->NNCLS[c]];

    nf->LOC[c][pos] = ultimo;
    nf->LOOC[ultimo] = pos;
}


/*INICIALIZACIÓN*/

void inicializar_nfold(Parametros p,Nfold *nf,int red[]){

    int c,sitio,s,clase;
    double delta_E;


    // Número total de espines
    

    // Preparar las diez clases
    for (c = 0; c < 10; c++){
        nf->NNCLS[c] = 0;
        delta_E=delta_energia2(p,c);
        nf->prob[c] = probabilidad_glauber(delta_E, p);
    }

    // Introducir cada sitio de la red en su clase
    for (sitio = 0; sitio < nf->N; sitio++) {
        clase = clase_sitio(p, red, sitio);
        insertar(nf, sitio, clase);
    }
}



void paso_nfold(Parametros p,int red[],Nfold *nf,long long *intentos,int *delta_M){
    
    
    int c,c_elegida,i,j,ind,sitio,n_afectados,repetido;
    int v[4];
    int afectados[5];
    double peso[10];
    double m,objetivo,acumulada;
    double W[10],Q[10];
    int clase_f[5];

    for (i = 0; i<10;i++)
    W[i] = 0;

    *intentos = 0;
    *delta_M = 0;

    // Peso de cada clase: n_c * p_c 
    for(i = 0;i<10; i++)
    {
      for (c = 0; c <= i; c++) 
    {
        peso[c] = nf->NNCLS[c] * nf->prob[c];
        W[i] += peso[c];
    }  
    }
    
    // Probabilidad total de aceptar un intento 
    for(i=0;i<10;i++)
    Q[i] = W[i] / nf->N;


    // Tiempo de espera en intentos:

        //m = floor(log(Random()) / log(1-Q)) + 1

    //Se toma r = 1-u, con u en [0,1).
    

    if (Q[9] >= 1.0) {
        m = 1.0;
    } 
    else {
        m = floor(log(1 - Random()) / log1p(-Q[9])) + 1.0;
    }


    //Elegir la clase 
    objetivo = Random() * Q[9];
    acumulada = 0.0;

    c_elegida = -1;

    for(i = 0;i<10;i++)
    {
        if(objetivo < Q[i])
        {
        c_elegida = i;
        break;
        }
    }

    //Elegir un espín dentro de esa clase 
    ind =(int)(Random() * nf->NNCLS[c_elegida]);

    sitio = nf->LOC[c_elegida][ind];

    afectados[0] = sitio;

    vecinos(sitio, p, v);

    for (i = 0; i < 4; i++) 
    {
        afectados[i+1] = v[i];
    }

    //Sacar los sitios afectados de sus clases 
    for (i = 0; i < 5; i++) {
        eliminar(nf, afectados[i]);
    } 

    // Introducir los sitios en sus nuevas clases 
    for (i = 0; i < 5; i++) {
        clase_f[i] = clase_sitio(p, red, afectados[i]);
    }

    act_clase(sitio,nf,clase_f);
    for(i=0; i<5; i++)
    {
        insertar(nf,afectados[i],clase_f[i]);
    }

    red[sitio] = -red[sitio];
    *delta_M = 2 * red[sitio];
    *intentos = m;
}

void act_clase(int sitio, Nfold *nf, int *clase_f)
{
    int i;

    if ( clase_f[0] < 5) 
    {
        clase_f[0]+= 5;
        for(i = 1;i<5;i++)
        {
            clase_f[i] += 1;
        } 
    }
    else 
    {
        clase_f[0] -= 5;
        for(i =1;i<5;i++)
        {
            clase_f[i] -=1;
        }
    }
    
}

double calcular_correlacion(int red[],int red_inicial[],int N) {
    int i;
    double suma;
    suma = 0.0;

    for (i = 0; i < N; i++) {
        suma += red_inicial[i] * red[i];
    }

    return suma / N;
}