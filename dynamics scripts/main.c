#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "Random_utils.h"
#include <time.h>



int main(){

    FILE *f;
    FILE *g;

    Parametros p;

    int N,n_sim,i,k,M,ite;
    char filename1[256],filename2[256];
    double t_CPU,intervalo,corr,t_fin;
    p.L=8;
    p.J=1.0;
    p.H=-0.3;
    p.kBT=50;
    n_sim=1000;
    intervalo=0.1;
    t_fin=10.0;

    int red[p.L*p.L],red_inicial[p.L*p.L];

    N=(int)p.L*p.L;
    ini_ran(123456);

    snprintf(filename1, sizeof(filename1),"../data_results/N=%d_J=%lf_H=%lf_kBT=%lf.csv",N ,p.J,p.H,p.kBT);
    snprintf(filename2, sizeof(filename2),"../data_results/correlacion_N=%d_J=%lf_H=%lf_kBT=%lf.csv",N ,p.J,p.H,p.kBT);
    f=fopen(filename1, "w");
    g=fopen(filename2, "w");

    if (f == NULL || g == NULL) {
        printf("Error al abrir el fichero");
        return 1;
    }

    fprintf(f, "# simulacion iteraciones t_CPU MCSS\n");
    fprintf(g, "# simulacion MCSS C\n");

    for(k=0;k<n_sim;k++){
        ite=0;
        clock_t inicio = clock();
        //Definimos la configuracion inicial: Todos los espines hacia arriba
        M=N;
        for(i=0;i<(p.L*p.L);i++){
                red[i]=red_inicial[i]=1;
        }
        
        corr=calcular_correlacion(red,red_inicial,N);
        fprintf(g, "%d %.12g %.10f\n",k,(double)ite/N,corr);
        while(M>0){
            intento_glauber(p,red);
            M=magnetizacion(p,red);
            ite++;
            corr=calcular_correlacion(red,red_inicial,N);

            if((ite/N)<t_fin && ite%(int)(intervalo*N) == 0){
                fprintf(g, "%d %.12g %.10f\n",k,(double)ite/N,corr);
            }
        }
        clock_t fin=clock();
        t_CPU = (double)(fin - inicio) / CLOCKS_PER_SEC;
        fprintf(f, "%d %d %.6f %.12g\n",k,ite,t_CPU,(double)ite/N);
    }
    fclose(f);
    fclose(g);

    return 0;
}


/*

int main(){
    FILE *g;
    FILE *f;

    Parametros p;
    Nfold nf;
    long long intentos,t;
    int N,n_sim,i,k,delta_M,M,ite;
    char filename1[256],filename2[256];
    double t_CPU,intervalo,t_fin,corr;
    

    p.L=8;
    p.J=1.0;
    p.H=-0.3;
    p.kBT=50;
    n_sim=1000;
    intervalo=0.1;
    t_fin=10.0;

    int red[p.L*p.L],red_inicial[p.L*p.L];
  
    N=(int)p.L*p.L;
    nf.N = N;
    ini_ran(123456);


    snprintf(filename1, sizeof(filename1),"../data_results/N-FW_N=%d_J=%lf_H=%lf_kBT=%lf.csv",N ,p.J,p.H,p.kBT);
    snprintf(filename2, sizeof(filename2),"../data_results/correlacion_N-FW_N=%d_J=%lf_H=%lf_kBT=%lf.csv",N ,p.J,p.H,p.kBT);
    f=fopen(filename1, "w");
    g=fopen(filename2, "w");

    if (f == NULL || g == NULL) {
        printf("Error al abrir el fichero");
        return 1;
    }

    fprintf(f, "# simulacion iteraciones t_CPU MCSS\n");
    fprintf(g, "# simulacion MCSS C\n");

    for(k=0;k<n_sim;k++){
        ite=0;
        t=0;
        clock_t inicio = clock();
        //Definimos la configuracion inicial: Todos los espines hacia arriba
        M = N;
        for(i=0;i<(p.L*p.L);i++){
                red[i]=red_inicial[i]=1;
        }
        
        inicializar_nfold(p, &nf, red);
        corr=calcular_correlacion(red,red_inicial,N);
        fprintf(g, "%d %.12g %.10f\n",k,(double)ite/N,corr);

        while(M>0){
            paso_nfold(p,red,&nf,&intentos,&delta_M);
            M += delta_M;
            
            t+=intentos;
            ite++;

            corr=calcular_correlacion(red,red_inicial,N);

            if((ite/N)<t_fin && ite%(int)(intervalo*N) == 0){
                fprintf(g, "%d %.12g %.10f\n",k,(double)ite/N,corr);
            }
        }
        clock_t fin=clock();
        t_CPU = (double)(fin - inicio) / CLOCKS_PER_SEC;
        fprintf(f, "%d %d %.6f %.12g\n",k,ite,t_CPU,(double)t/N);
    }
    fclose(f);
    fclose(g);

    return 0;
}

*/