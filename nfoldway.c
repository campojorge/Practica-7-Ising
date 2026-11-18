/*
 * Algoritmo n-fold way para el modelo de Ising 2D (Sec. 5.2 del paper).
 * Escape del estado metaestable: desde todos los espines UP hasta M = 0.
 *
 * Requiere random.c (Parisi-Rapuano) en el mismo directorio, con NormRANu
 * definido en el. Ese fichero debe poder incluirse tal cual (con sus propios
 * #include si los necesita, p. ej. stdlib.h para srand).
 *
 * Compilar:  gcc -std=c99 -O2 -Wall -Wextra -o nfold nfold_pr.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "random.c"

#define NCLASES 10

/* Si ya tienes este tipo definido en otro fichero, elimina esta definicion */
typedef struct {
    int L;
    double J;
    double H;
    double kBT;
} Parametros;

/*
 * Estado del n-fold way: solo se guardan las clases de los espines.
 *   nncls[i]  : n_i, numero de espines de la clase i       (NNCLS)
 *   loc       : posiciones de los espines de cada clase    (LOC, 10 x N)
 *   look      : indice de cada espin dentro de su clase    (LOOK)
 *   looc      : clase de cada espin                        (LOOC)
 * Se cumple loc[looc[i]*N + look[i]] = i.
 *
 * Clases (indices 0..9; en la Tabla 1 del paper son 1..10):
 *   0..4 : espin UP   con 4,3,2,1,0 vecinos UP
 *   5..9 : espin DOWN con 4,3,2,1,0 vecinos UP
 */
typedef struct {
    int N;
    int nncls[NCLASES];
    int *loc;
    int *look;
    int *looc;
} EstadoNFold;

/* ------------------------------------------------------------------ */
/* Probabilidades de inversion p_i de cada clase                       */
/* ------------------------------------------------------------------ */
void calcular_probabilidades(Parametros p, double prob[NCLASES]){

    int c, u, suma_vecinos;
    double s, menos_Eold, menos_Enew;

    for (c = 0; c < NCLASES; c++) {

        /* Orientacion del espin y numero de vecinos UP (Tabla 1) */
        if (c < 5) {
            s = 1.0;
        } else {
            s = -1.0;
        }
        u = 4 - (c % 5);
        suma_vecinos = 2 * u - 4;

        /* Columnas -E_old y -E_new de la Tabla 1 */
        menos_Eold = s * (p.J * suma_vecinos + p.H);
        menos_Enew = -menos_Eold;

        /* p_i = exp(-Enew/kT) / [exp(-Eold/kT) + exp(-Enew/kT)]  (Ec. 2),
           escrito como 1/(1+exp(x)) para evitar desbordamientos */
        prob[c] = 1.0 / (1.0 + exp((menos_Eold - menos_Enew) / p.kBT));
    }
}

/* ------------------------------------------------------------------ */
/* Vecinos de un espin (condiciones de contorno periodicas)            */
/* El espin i esta en la fila y = i / L y la columna x = i % L          */
/* ------------------------------------------------------------------ */
void obtener_vecinos(Parametros p, int i, int vecinos[4]){

    int x, y, arriba, abajo, izquierda, derecha;

    y = i / p.L;
    x = i % p.L;

    if (y == 0) {
        abajo = p.L - 1;
    } else {
        abajo = y - 1;
    }

    if (y == p.L - 1) {
        arriba = 0;
    } else {
        arriba = y + 1;
    }

    if (x == 0) {
        izquierda = p.L - 1;
    } else {
        izquierda = x - 1;
    }

    if (x == p.L - 1) {
        derecha = 0;
    } else {
        derecha = x + 1;
    }

    vecinos[0] = arriba * p.L + x;
    vecinos[1] = abajo * p.L + x;
    vecinos[2] = y * p.L + izquierda;
    vecinos[3] = y * p.L + derecha;
}

/* ------------------------------------------------------------------ */
/* Gestion del estado                                                  */
/* ------------------------------------------------------------------ */
int estado_crear(Parametros p, EstadoNFold *e){

    e->N = p.L * p.L;
    e->loc = (int *)malloc(sizeof(int) * (size_t)e->N * NCLASES);
    e->look = (int *)malloc(sizeof(int) * (size_t)e->N);
    e->looc = (int *)malloc(sizeof(int) * (size_t)e->N);

    if (e->loc == NULL || e->look == NULL || e->looc == NULL) {
        return -1;
    }
    return 0;
}

void estado_liberar(EstadoNFold *e){
    free(e->loc);
    free(e->look);
    free(e->looc);
}

/* Todos los espines UP: todos en la clase 0 (n_1 = N y n_i = 0 para i > 1) */
void estado_inicial(EstadoNFold *e){

    int i;

    for (i = 0; i < NCLASES; i++) {
        e->nncls[i] = 0;
    }
    e->nncls[0] = e->N;

    for (i = 0; i < e->N; i++) {
        e->loc[i] = i;
        e->look[i] = i;
        e->looc[i] = 0;
    }
}

/* Cambia un espin de clase actualizando NNCLS, LOC, LOOK y LOOC */
void mover_spin(EstadoNFold *e, int spin, int clase_nueva){

    int N, clase_vieja, pos_vieja, ultima_pos, ultimo_spin, pos_nueva;

    N = e->N;
    clase_vieja = e->looc[spin];
    pos_vieja = e->look[spin];

    /* Se saca de la clase vieja: el ultimo espin ocupa su hueco */
    ultima_pos = e->nncls[clase_vieja] - 1;
    ultimo_spin = e->loc[clase_vieja * N + ultima_pos];
    e->loc[clase_vieja * N + pos_vieja] = ultimo_spin;
    e->look[ultimo_spin] = pos_vieja;
    e->nncls[clase_vieja]--;

    /* Se anade al final de la clase nueva */
    pos_nueva = e->nncls[clase_nueva];
    e->loc[clase_nueva * N + pos_nueva] = spin;
    e->looc[spin] = clase_nueva;
    e->look[spin] = pos_nueva;
    e->nncls[clase_nueva]++;
}

/* ------------------------------------------------------------------ */
/* Una ejecucion del n-fold way hasta M = 0. Devuelve el tiempo en MCSS */
/* ------------------------------------------------------------------ */
double nfold_escape(Parametros p, EstadoNFold *e, const double prob[NCLASES]){

    int N, M, j, k, v, nb, spin, sentido;
    int vecinos[4];
    long long t;                    /* tiempo en mcs */
    double Q[NCLASES];
    double acumulado, Q10, objetivo, mm;
    double r_tilde, r_bar, r;

    N = e->N;
    M = N;
    t = 0;
    estado_inicial(e);

    while (M > 0) {

        /* Tres numeros aleatorios: r_tilde, r_bar y r */
        r_tilde = Random();
        do {
            r_bar = Random();
        } while (r_bar <= 0.0);     /* evita ln(0) */
        r = Random();

        /* Los 10 valores Q_i (Ec. 42) */
        acumulado = 0.0;
        for (j = 0; j < NCLASES; j++) {
            acumulado += (double)e->nncls[j] * prob[j] / N;
            Q[j] = acumulado;
        }
        Q10 = Q[NCLASES - 1];

        if (!(Q10 > 0.0)) {
            fprintf(stderr, "Error: Q10 = 0 (underflow); temperatura demasiado baja.\n");
            exit(EXIT_FAILURE);
        }

        /* Tiempo para salir de la configuracion actual (Ec. 44) */
        if (Q10 >= 1.0) {
            mm = 1.0;
        } else {
            mm = floor(log(r_bar) / log1p(-Q10)) + 1.0;
        }
        if (!(mm < 1.0e18)) {
            fprintf(stderr, "Error: tiempo de escape fuera de rango.\n");
            exit(EXIT_FAILURE);
        }

        /* Clase k que cumple Q_{k-1} <= r_tilde*Q10 < Q_k (Ec. 45) */
        objetivo = r_tilde * Q10;
        k = 0;
        for (j = 0; j < NCLASES; j++) {
            if (e->nncls[j] > 0 && prob[j] > 0.0) {
                k = j;
                if (objetivo < Q[j]) {
                    break;
                }
            }
        }

        /* Con r se elige uno de los n_k espines de la clase k */
        j = (int)(r * e->nncls[k]);
        if (j >= e->nncls[k]) {
            j = e->nncls[k] - 1;
        }
        spin = e->loc[k * N + j];

        /* Se invierte el espin: clase k -> k+5 (UP) o k-5 (DOWN).
           Sus 4 vecinos cambian de clase en +1 (UP) o -1 (DOWN). */
        if (k < 5) {
            sentido = +1;
            M -= 2;
        } else {
            sentido = -1;
            M += 2;
        }
        mover_spin(e, spin, k + 5 * sentido);

        obtener_vecinos(p, spin, vecinos);
        for (v = 0; v < 4; v++) {
            nb = vecinos[v];
            mover_spin(e, nb, e->looc[nb] + sentido);
        }

        /* Se avanza el tiempo en m intentos */
        t += (long long)mm;
    }

    return (double)t / N;           /* MCSS = mcs / N */
}

int main(void){

    Parametros p;
    EstadoNFold estado;
    double prob[NCLASES];
    double tiempo, suma, suma2, media, desv;
    int K, i, semilla;

    /* Parametros de la simulacion */
    p.L = 32;                       /* se requiere L >= 3 */
    p.J = 1.0;
    p.H = -3.0;
    p.kBT = 0.5;
    K = 1000;                       /* numero de escapes */
    semilla = 42;

    if (p.L < 3 || p.kBT <= 0.0) {
        fprintf(stderr, "Error: se requiere L >= 3 y kBT > 0.\n");
        return EXIT_FAILURE;
    }

    ini_ran(semilla);
    calcular_probabilidades(p, prob);

    if (estado_crear(p, &estado) != 0) {
        fprintf(stderr, "Error: sin memoria.\n");
        estado_liberar(&estado);
        return EXIT_FAILURE;
    }

    suma = 0.0;
    suma2 = 0.0;
    for (i = 0; i < K; i++) {
        tiempo = nfold_escape(p, &estado, prob);
        suma += tiempo;
        suma2 += tiempo * tiempo;
    }

    media = suma / K;
    desv = sqrt(suma2 / K - media * media);
    printf("Tiempo medio de escape (K=%d): %.6g MCSS (desv. est. = %.6g)\n",
           K, media, desv);

    estado_liberar(&estado);
    return EXIT_SUCCESS;
}