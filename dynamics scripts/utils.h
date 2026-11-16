
typedef struct {
    int L;       // Lado de la red
    double J;    // Acoplamiento entre primeros vecinos
    double H;    // Campo externo 
    double kBT;  // Energia termica 
} Parametros;

typedef struct {
    int x;  // Columna
    int y;  // Fila
} Posicion;


double delta_energia(Parametros p, int red[p.L][p.L], Posicion pos);
double probabilidad_glauber(double delta_E, Parametros p);
void intento_glauber(Parametros p, int red[p.L][p.L]);
void paso_glauber(Parametros p, int red[p.L][p.L]);
double magnetizacion(Parametros p,int red[p.L][p.L]);

