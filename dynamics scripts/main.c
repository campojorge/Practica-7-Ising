#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "Random_utils.h"
#include <time.h>



/*int main(){

    FILE *f;

    Parametros p;
    p.L=8;
    p.J=1.0;
    p.H=-0.3;
    p.kBT=1.8;

    int N,n_sim,i,j,k,M,intento;
    int red[p.L*p.L];
    char filename1[256];

    N=(int)p.L*p.L;
    n_sim=1;
    ini_ran(123456);

    snprintf(filename1, sizeof(filename1),"../data_results/N=%d_J=%lf_H=%lf_kBT=%lf.csv",N ,p.J,p.H,p.kBT);
    f=fopen(filename1, "w");

    if (f == NULL) {
        printf("Error al abrir el fichero");
        return 1;
    }

    fprintf(f, "# simulacion t M\n");

    for(k=0;k<n_sim;k++){

        intento=0;

        //Definimos la configuracion inicial: Todos los espines hacia arriba
        M=N;
        for(i=0;i<(p.L*p.L);i++){
                red[i]=1;
        }
        

        while(M!=0){
            intento_glauber(p,red);
            M=magnetizacion(p,red);
            intento++;
            fprintf(f, "%d %.12g %d \n",k,(double)intento/N,M);
            
        }
    }
    fclose(f);

    return 0;
}*/

static Nfold nf_global;      // fuera de main


int main(){

    FILE *f;

    Parametros p;
    Nfold *nf;
    long long intentos;
    int delta_M;

    p.L=8;
    p.J=1.0;
    p.H=-0.3;
    p.kBT=1.8;

    int N,n_sim,i,j,k,t;
    int red[p.L*p.L];
    char filename1[256];

    
    N=(int)p.L*p.L;
    n_sim=1;
    ini_ran(123456);
    int M;
    nf->N = N;

    snprintf(filename1, sizeof(filename1),"../data_results/N-FW_N=%d_J=%lf_H=%lf_kBT=%lf.csv",N ,p.J,p.H,p.kBT);
    f=fopen(filename1, "w");

    if (f == NULL) {
        printf("Error al abrir el fichero");
        return 1;
    }

    fprintf(f, "# simulacion t M\n");

    for(k=0;k<n_sim;k++){

        t=0;
        //Definimos la configuracion inicial: Todos los espines hacia arriba
        M=N;
        for(i=0;i<(p.L*p.L);i++){
                red[i]=1;
        }
        
        inicializar_nfold(p, nf, red);
        while(M>0){
            paso_nfold(p, red,nf,&intentos,&delta_M);
            M += delta_M;
            t+=intentos;
            fprintf(f, "%d %.12g %d \n",k,(double)intentos/N,M);
            
        }
    }
    fclose(f);

    return 0;
}


