#include "utils.h"
#include <math.h>

double delta_energia(Parametros p, int red[p.L][p.L], Posicion pos){

    int arriba,abajo,derecha,izquierda,suma_vecinos;
    double delta;

    //Buscamos los cuatro vecinos del espin seleccionado(Condiciones de contorno periódicas)
    if (pos.y == 0) {
        abajo = p.L - 1;
    } else {
        abajo = pos.y - 1;
    }

    if (pos.y == p.L - 1) {
        arriba = 0;
    } else {
        arriba = pos.y + 1;
    }

    if (pos.x == 0) {
        izquierda = p.L - 1;
    } else {
        izquierda = pos.x - 1;
    }

    if (pos.x == p.L - 1) {
        derecha = 0;
    } else {
        derecha = pos.x + 1;
    }

    suma_vecinos = red[arriba][pos.x] + red[abajo][pos.x]+ red[pos.y][izquierda] + red[pos.y][derecha];

    delta=2.0 * red[pos.y][pos.x] * (p.J * suma_vecinos + p.H);
    return delta
}
   

double probabilidad_glauber(double delta_E, double kBT){
    double prob;
    //Tratemos primero los casos especiales
    //1-Si ambas configuraicones tienen la misma energia, entonces la probabilidad de inversion es 1/2
    if (delta_E == 0.0) {
        prob=0.5;
    }
    //2-Para temperatura 0 un cambio que disminuye la energia se acepta siempre y uno que la aumenta se rechaza siempre
    if (kBT == 0.0) {
        if(delta_E<0.0){
            prob=1.0;
        }
        else{
            prob=0.0;
        }
    }

    //Caso general
   
    if (delta_E > 0.0 && kBT>0.0) {
        prob=1/(1+exp(delta_E/kBT));
    }
    
    return prob;
}

int intento_glauber(const Parametros *p, int red[p->L][p->L],
                    double (*uniforme)(void))
{
    Posicion pos;
    pos.y = (int)(uniforme() * p->L);
    pos.x = (int)(uniforme() * p->L);

    double delta_E = delta_energia(p, red, pos);
    double probabilidad = probabilidad_glauber(delta_E, p->kBT);

    if (uniforme() < probabilidad) {
        red[pos.y][pos.x] = -red[pos.y][pos.x];
        return 1;
    }
    return 0;
}

size_t paso_glauber(const Parametros *p, int red[p->L][p->L],
                    double (*uniforme)(void))
{
    size_t aceptados = 0;
    size_t n_espines = (size_t)p->L * p->L;

    for (size_t intento = 0; intento < n_espines; ++intento) {
        aceptados += intento_glauber(p, red, uniforme);
    }
    return aceptados;
}