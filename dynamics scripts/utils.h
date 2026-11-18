#define L_MAX 128
#define N_MAX (L_MAX * L_MAX)

typedef struct {
    int L;       // Lado de la red
    double J;    // Acoplamiento entre primeros vecinos
    double H;    // Campo externo 
    double kBT;  // Energia termica 
} Parametros;

typedef struct {
    int N;

    /* Número de espines y probabilidad de cada clase */
    int NNCLS[10];
    double prob[10];

    /* Arrays normales: no necesitan malloc */
    int LOC[10][N_MAX];
    int LOOK[N_MAX];
    int LOOC[N_MAX];

} Nfold;


void vecinos(int sitio, Parametros p, int v[]);
double delta_energia(Parametros p, int red[], int sitio);
double delta_energia2(Parametros p, int clase);
double probabilidad_glauber(double delta_E, Parametros p);
void intento_glauber(Parametros p, int red[]);
double magnetizacion(Parametros p,int red[]);

 // FUNCIONES N-FOLD WAY

int clase_sitio(Parametros p,int red[],int sitio);
void insertar(Nfold *nf, int sitio, int c);
void eliminar(Nfold *nf, int sitio);
void inicializar_nfold(Parametros p,Nfold *nf, int red[]);
void act_clase(int sitio, Nfold *nf, int *clase_f);
void paso_nfold(Parametros p,int red[],Nfold *nf,long long *intentos,int *delta_M);


