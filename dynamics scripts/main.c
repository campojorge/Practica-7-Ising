#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "Random_utils.h"
#include <time.h>

/*

int main(){

    FILE *f;

    Parametros p;

    int N,n_sim,i,k,M,ite;
    char filename1[256];
    double t_CPU;
    p.L=8;
    p.J=1.0;
    p.H=-0.3;
    p.kBT=1.8;
    n_sim=50000;

    int red[p.L*p.L];

    N=(int)p.L*p.L;
    ini_ran(123456);

    snprintf(filename1, sizeof(filename1),"../data_results/N=%d_J=%lf_H=%lf_kBT=%lf.csv",N ,p.J,p.H,p.kBT);
    f=fopen(filename1, "w");

    if (f == NULL) {
        printf("Error al abrir el fichero");
        return 1;
    }

    fprintf(f, "# simulacion iteraciones t_CPU MCSS\n");

    for(k=0;k<n_sim;k++){
        ite=0;
        clock_t inicio = clock();
        //Definimos la configuracion inicial: Todos los espines hacia arriba
        M=N;
        for(i=0;i<(p.L*p.L);i++){
                red[i]=1;
        }
        
        while(M>0){
            intento_glauber(p,red);
            M=magnetizacion(p,red);
            ite++;
        }
        clock_t fin=clock();
        t_CPU = (double)(fin - inicio) / CLOCKS_PER_SEC;
        fprintf(f, "%d %d %.6f %.12g\n",k,ite,t_CPU,(double)ite/N);
    }
    fclose(f);

    return 0;
}

*/

/*
int main(){
    FILE *f;

    Parametros p;
    Nfold nf;
    long long intentos,t;
    int N,n_sim,i,k,delta_M,M,ite;
    char filename1[256];
    double t_CPU;
    

    p.L=8;
    p.J=1.0;
    p.H=-0.3;
    p.kBT=1.8;
    n_sim=50000;


    int red[p.L*p.L];
  
    N=(int)p.L*p.L;
    nf.N = N;
    ini_ran(123456);


    snprintf(filename1, sizeof(filename1),"../data_results/N-FW_N=%d_J=%lf_H=%lf_kBT=%lf.csv",N ,p.J,p.H,p.kBT);
    f=fopen(filename1, "w");

    if (f == NULL) {
        printf("Error al abrir el fichero");
        return 1;
    }

    fprintf(f, "# simulacion iteraciones t_CPU MCSS\n");

    for(k=0;k<n_sim;k++){
        ite=0;
        t=0;
        clock_t inicio = clock();
        //Definimos la configuracion inicial: Todos los espines hacia arriba
        M = N;
        for(i=0;i<(p.L*p.L);i++){
                red[i]=1;
        }
        
        inicializar_nfold(p, &nf, red);

        while(M>0){
            paso_nfold(p,red,&nf,&intentos,&delta_M);
            M += delta_M;
            
            t+=intentos;
            ite++;
        }
        clock_t fin=clock();
        t_CPU = (double)(fin - inicio) / CLOCKS_PER_SEC;
        fprintf(f, "%d %d %.6f %.12g\n",k,ite,t_CPU,(double)t/N);
    }
    fclose(f);

    return 0;
}
*/
int main()
{

    FILE *file_M;
    FILE *file_E;

    Parametros p;
    Nfold nf;
    long long intentos,t;
    int N,n_sim,i,k,delta_M,M,ite,medida;
    char filename1[256];
    char filename2[256];    
    double intervalo,t_fin,corr,corr_antes,t_medida;

    p.L=8;
    p.J=1.0;
    p.H=-0.3;
    p.kBT=1.3615; //T = 0.6Tc
    n_sim=10000;
    intervalo=0.1;
    t_fin=10.0;

    int red[p.L*p.L],red_inicial[p.L*p.L];

    N=p.L*p.L;
    nf.N=N;
    ini_ran(123456);

    snprintf(filename1, sizeof(filename1),"../data_results/Test2_Glauber_M(t)_N=%d_J=%lf_H=%lf_kBT=%lf.csv",N,p.J,p.H,p.kBT);
    snprintf(filename2, sizeof(filename2),"../data_results/Test2_Glauber_E(t)_N=%d_J=%lf_H=%lf_kBT=%lf.csv",N,p.J,p.H,p.kBT);

    file_M=fopen(filename1,"w");
    file_E=fopen(filename2,"w");

    if(file_M==NULL){
        printf("Error al abrir el fichero");
        return 1;
    }

    //1. GENERAMOS LA CONFIGURACIÓN CON TODO ESPINES HACIA ARRIBA.

    //int i;
    for (i=0;i<N;i++)
    {
        red[i]=red_inicial[i]=1;
    }

    int Magnet[(int)t_fin];
    int Energia[(int)t_fin];
    int E = 0;
    //2. GENERAMOS UN BUCLE CON 5 ITERACCIONES PARA CADA VALOR DE LA TEMPERATURA A BARRER. 

        // TENEMOS QUE HACER 50000 SIMULACIONES PARA CADA TEMPERATURA, Y GUARDAR EL PROMEDIO DE <M(t)>, <E(t)> Y <H(t)> EN UN FICHERO. 
        for(n_sim = 0;n_sim<50000;n_sim++)
        {
            //4. HACEMOS UN BUCLE PARA REALIZAR EL ESCAPE DEL SISTEMA HASTA QUE SE ALCANCE EL TIEMPO FINAL.
            t=0;
            while((double)t/N<t_fin) // t_fin = 10 MCSS
            {
                intento_glauber(p,red);
                M=magnetizacion(p,red);
                E = energy(p,red);
                t+=1;
                Magnet[i] = M;
                Energia[i] = E;
            }
        fprintf(file_M,"%d %d %d %d %d %d %d %d %d %d\n",Magnet[0],Magnet[1],Magnet[2],Magnet[3],Magnet[4],Magnet[5],Magnet[6],Magnet[7],Magnet[8],Magnet[9]);
        fprintf(file_E,"%d %d %d %d %d %d %d %d %d %d\n",Energia[0],Energia[1],Energia[2],Energia[3],Energia[4],Energia[5],Energia[6],Energia[7],Energia[8],Energia[9]);
    
        }
 }