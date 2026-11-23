#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "Random_utils.h"
#include <time.h>

/*

int main(){

    FILE *g;

    Parametros p;

    int N,n_sim,i,k,M,ite;
    char filename2[256];
    double intervalo,corr,t_fin;
    p.L=8;
    p.J=0;
    p.H=0;
    p.kBT=1.8;
    n_sim=10000;
    intervalo=0.1;
    t_fin=10.0;

    int red[p.L*p.L],red_inicial[p.L*p.L];

    N=(int)p.L*p.L;
    ini_ran(123456);

    snprintf(filename2, sizeof(filename2),"../data_results/correlacion_N=%d_J=%lf_H=%lf_kBT=%lf.csv",N ,p.J,p.H,p.kBT);
    g=fopen(filename2, "w");

    if (g == NULL) {
        printf("Error al abrir el fichero");
        return 1;
    }

    fprintf(g, "# simulacion MCSS C\n");

    for(k=0;k<n_sim;k++){
        ite=0;
        //Definimos la configuracion inicial: Todos los espines hacia arriba
        M=N;
        for(i=0;i<(p.L*p.L);i++){
                red[i]=red_inicial[i]=1;
        }
        
        corr=calcular_correlacion(red,red_inicial,N);
        fprintf(g, "%d %.12g %.10f\n",k,(double)ite/N,corr);
        while((double)(ite/N)<t_fin){
            intento_glauber(p,red);
            M=magnetizacion(p,red);
            ite++;
            corr=calcular_correlacion(red,red_inicial,N);
            if(ite%(int)(intervalo*N) == 0){
                fprintf(g, "%d %.12g %.10f\n",k,(double)ite/N,corr);
            }
        }
    }
    fclose(g);
    return 0;
}

*/

int main(){
    FILE *g;

    Parametros p;
    Nfold nf;
    long long intentos,t;
    int N,n_sim,i,k,delta_M,M,ite,medida;
    char filename2[256];
    double intervalo,t_fin,corr,corr_antes,t_medida;

    p.L=8;
    p.J=1.0;
    p.H=-0.3;
    p.kBT=50;
    n_sim=10000;
    intervalo=0.1;
    t_fin=10.0;

    int red[p.L*p.L],red_inicial[p.L*p.L];

    N=p.L*p.L;
    nf.N=N;
    ini_ran(123456);

    snprintf(filename2, sizeof(filename2),"../data_results/correlacion_N-FW_N=%d_J=%lf_H=%lf_kBT=%lf.csv",N,p.J,p.H,p.kBT);

    g=fopen(filename2,"w");

    if(g==NULL){
        printf("Error al abrir el fichero");
        return 1;
    }

    fprintf(g,"# simulacion MCSS C\n");

    for(k=0;k<n_sim;k++){
        ite=0;
        t=0;
        medida=1;

        M=N;
        for(i=0;i<N;i++){
            red[i]=red_inicial[i]=1;
        }

        inicializar_nfold(p,&nf,red);

        corr=calcular_correlacion(red,red_inicial,N);
        fprintf(g,"%d %.12g %.10f\n",k,0.0,corr);

        while((double)t/N<t_fin){

            // Guardamos la correlación antes del cambio.
            corr_antes=corr;

            paso_nfold(p,red,&nf,&intentos,&delta_M);

            M+=delta_M;
            t+=intentos;
            ite++;

            corr=calcular_correlacion(red,red_inicial,N);

            // Mismos tiempos de medida que en Glauber.
            t_medida=ceil(medida*intervalo*N-1e-9)/N;

            // Guardar todas las medidas atravesadas por el salto.
            while(t_medida<=(double)t/N && t_medida<=t_fin){

                if(t_medida<(double)t/N){
                    fprintf(g,"%d %.12g %.10f\n",k,t_medida,corr_antes);
                }
                else{
                    fprintf(g,"%d %.12g %.10f\n",k,t_medida,corr);
                }

                medida++;
                t_medida=ceil(medida*intervalo*N-1e-9)/N;
            }
        }
    }

    fclose(g);
    return 0;
}