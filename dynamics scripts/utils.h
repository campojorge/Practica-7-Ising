#include <stddef.h>
#define PRACTICA7_UTILS_H


typedef struct {
    int L;       // Lado de la red
    double J;    // Acoplamiento entre primeros vecinos
    double H;    // Campo externo 
    double kB;  // Constate de boltsmann
} Parametros;

typedef struct {
    int x;  // Columna
    int y;  // Fila
} Posicion;


double delta_energia(Parametros p, int red[p.L][p.L], Posicion pos);
double probabilidad_glauber(double delta_E, double kBT);


int intento_glauber(const Parametros *p, int red[p->L][p->L],
                    double (*uniforme)(void));

/* Realiza L*L intentos con reemplazo y devuelve el número aceptado. */
size_t paso_glauber(const Parametros *p, int red[p->L][p->L],
                    double (*uniforme)(void));

