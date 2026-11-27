#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "Random_utils.h"
#include <time.h>
const double KB = 1.380649e-23;

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

    FILE *file;

    Parametros p;
    Nfold nf;
    long long intentos,t;
    int N,n_sim,i,k,delta_M,M,ite,medida;
    char filename[256];    
    double intervalo,t_fin,corr,corr_antes,t_medida;

    p.L=8;
    p.J=1.0;
    p.H=-0.3;
    p.kBT=50; //T = 0.6Tc
    n_sim=10000;
    intervalo=0.1;
    t_fin=10000.0;

    int red[p.L*p.L],red_inicial[p.L*p.L];

    N=p.L*p.L;
    nf.N=N;
    ini_ran(123456);

    snprintf(filename, sizeof(filename),"../data_results/Test2_Glauber_N=%d_J=%lf_H=%lf_kBT=%lf.csv",N,p.J,p.H,p.kBT);

    file=fopen(filename,"w");
    fprintf(file,"# t MCSS M_med M_var E_med E_var Cv ErrM ErrE ErrCv\n");

    if(file==NULL){
        printf("Error al abrir el fichero");
        return 1;
    }

    //1. GENERAMOS LA CONFIGURACIÓN CON TODO ESPINES HACIA ARRIBA.

    //int i;
    for (i=0;i<N;i++)
    {
        red[i]=red_inicial[i]=1;
    }
    int SM[(int)t_fin],SE[(int)t_fin],SM2[(int)t_fin],SE2[(int)t_fin]; 
    int E;
    double Cv, ErrM, ErrE, ErrCv;
    M=magnetizacion(p,red);
    E = energy(p,red);

    // INICIALIZAMOS TODAS EN 0

    for (i=0;i<(int)t_fin;i++)
    {
    SM[i] = 0;
    SE[i] = 0;
    SM2[i] = 0;
    SE2[i] = 0;
    }

    // TENEMOS QUE HACER 50000 SIMULACIONES PARA CADA TEMPERATURA, Y GUARDAR EL PROMEDIO DE <M(t)>, <E(t)> Y <H(t)> EN UN FICHERO. 
   /* for(n_sim = 0;n_sim<50000;n_sim++)
    {
        //4. HACEMOS UN BUCLE PARA REALIZAR EL ESCAPE DEL SISTEMA HASTA QUE SE ALCANCE EL TIEMPO FINAL.
        t=0;
        while((double)t/N<t_fin) // t_fin = 10000 MCSS
        {
            intento_glauber(p,red);
               
            if(t%N == 0)
            {
                M=magnetizacion(p,red);
                E = energy(p,red);
                SM[(int)t/N] += M;
                SE[(int)t/N] += E;
                SM2[(int)t/N] += M*M;
                SE2[(int)t/N] += E*E;
            }
                t+=1;
        }
        // HAY QUE HACER UN RESET EN LA RED (INICIALIZARLA SPIN UP)
        for (i=0;i<N;i++)
            {
                red[i]=1;
            }
    }

    for (i=0;i<(int)t_fin;i++)
    {
        SM[i] /= 50000;
        SE[i] /= 50000;
        SM2[i] = SM2[i]/50000 - SM[i]*SM[i];
        SE2[i] = SE2[i]/50000 - SE[i]*SE[i];
        Cv = (SE2[i]-SE[i]*SE[i])/(p.L*p.L*p.kBT*p.kBT/KB);
        ErrM = sqrt(SM2[i]/50000);
        ErrE = sqrt(SE2[i]/50000);
        ErrCv = sqrt((SE2[i]/50000 - Cv*Cv)/50000);
        fprintf(file,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n",(double)i,SM[i],SM2[i],SE[i],SE2[i],Cv,ErrM,ErrE,ErrCv);
    }*/
 }