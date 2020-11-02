/* Programma di simulazione di dinamica molecolare */


# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <math.h>

# define Rcut 2.5
# define EMPTY -1
# define Natoms 864
# define DIM 3
# define dt 0.001
# define PI 3.1415926535
# define NCLMAX 10000 /* Maximum number of linked-list cells */

double L;


void Initialization (double **, double **, double **, double, double, int *);
void Rescale (double **, double);
void Evolution (double **, double **, double **, double *);
void HalfStep (double **, double **);
void Acceleration (double **, double **, double *);
void Observables (double **, double *, double *, double, int, FILE *);
double instantaneousTemperature (double **);
double Gaussian(void);

int main () {
    
    int stepCount, stepLimit, stepPrint, stepScale, seed, i, d;
    double **v, **r, **a, kin, pot, mass, T, rho;
    
    FILE *fp;
    
    seed = time(0);
    srand48(seed);
    
    fp = fopen("drift.dat","w+");
    
    /*printf ("\n Inserire il tempo computazionale massimo(int 10) : ");
     scanf ("%d", &stepLimit);
     
     printf ("\n Inserire la massa delle particelle (double 1.0): ");
     scanf ("%lf", &mass);
     
     printf ("\n Inserire la temperatura iniziale del sistema in unità ridotte (double 1.0): ");
     scanf ("%lf", &T);
     
     printf ("\n Inserire la densita' di particelle in unità ridotte (double 0.8): ");
     scanf ("%lf", &rho);*/
    
    stepLimit = 1000;
    mass = 1.0;
    T = 0.1;
    rho = 0.2;
    
    stepPrint = (int) (stepLimit / 1000);
    stepScale = (int) (stepLimit / 5);
    
    r = ( double ** ) malloc ( Natoms * sizeof ( double * ) );
    v = ( double ** ) malloc ( Natoms * sizeof ( double * ) );
    a = ( double ** ) malloc ( Natoms * sizeof ( double * ) );
    
    for (i=0;i<Natoms;i++){ 
        
        r[i] = ( double * ) malloc ( DIM * sizeof ( double ) );
        v[i] = ( double * ) malloc ( DIM * sizeof ( double ) );
        a[i] = ( double * ) malloc ( DIM * sizeof ( double ) );
    }
    
    if ((a == NULL) || (r == NULL) || (v == NULL ))
        printf("\n si è verificato un errore nella allocazione della memoria");
    
    Initialization (r, v, a, T, rho, &seed);
    
    for (stepCount=1; stepCount<=stepLimit; stepCount++){
        
        Evolution (r, v, a, &pot);
               
        if (stepCount % stepPrint == 0)
            Observables(v, &kin, &pot, mass, stepCount,fp);
        
        /*if (stepCount % stepScale == 0)
         Rescale(v, T);*/
        
        //fprintf (fp,"%d\t%lf\n",stepCount, instantaneousTemperature(v));
    }
    
    
    for (i=0;i<Natoms;i++){
        
        free (r[i]);
        free (v[i]);
        free (a[i]);
    }
    
    free (r);
    free (v);
    free (a);
    
    fclose(fp);
    
}

/* Questa funzione, chiamata una sola volta dal pogramma inizializza gli array di posizione e velocità riscalando il momento totale a zero */

void Initialization (double **r, double **v, double **a, double T, double rho, int *seed){
    
    int i, d, k, M=1, nX, nY, nZ;
    double sumv[DIM] = {0.}, vSum[DIM]={0.0}, e[DIM];
    double firstCell[4][3] = {
        {0.25, 0.25, 0.25},
        {0.75, 0.75, 0.25},
        {0.75, 0.25, 0.75},
        {0.25, 0.75, 0.75}};
    
    
    L = pow(Natoms / rho, 1.0/3);
    
    while (4 * M * M * M < Natoms) M++; // M^3 è così il numero di celle unitarie    
    
    double l = L / M; // a è la lunghezza di un lato di questi cubi
    
    int n = 0;
    
    
    for (nX=0; nX<M;nX++)
        
        for (nY=0; nY<M;nY++)
            
            for (nZ=0; nZ<M;nZ++)
                
                for (k=0; k<4; k++)
                    
                    if (n<Natoms){
                        
                        r[n][0] = (nX + firstCell[k][0]) * l;
                        r[n][1] = (nY + firstCell[k][1]) * l;
                        r[n][2] = (nZ + firstCell[k][2]) * l;
                        
                        n++;
                    }
    
    
    
    for (i=0;i<Natoms;i++)
        for (d=0;d<DIM;d++){
            
            v[i][d] = Gaussian();
            vSum[d] += v[i][d];
        }
    
    for (d=0;d<DIM;d++) vSum[d] /= Natoms;
    
    for (i=0;i<Natoms;i++)
        for (d=0;d<DIM;d++){
            
            v[i][d] -= vSum[d]; /* momento totale == NULL */
            a[i][d] = 0.; /* atomi inizialment fermi */
        }
    
    Rescale (v, T);
}

void Rescale (double **v, double T){
    
    int i, d;
    double vSum2=0.0, fs;
    
    for (i=0;i<Natoms;i++)
        for(d=0;d<DIM;d++)
            
            vSum2 += v[i][d] * v[i][d];
    
    fs = sqrt (3 * T * Natoms / vSum2); /* fattore di scala per tener conto della temperatura iniziale del sistema */
    
    for (i=0;i<Natoms;i++)
        for(d=0;d<DIM;d++)
            
            v[i][d] *= fs;
    
}



void Evolution (double **r, double **v, double **a, double *pot){
    
    int i, d;
    
    HalfStep (v, a);
    
    for (i = 0; i < Natoms; i++)
        for (d=0; d<DIM; d++)
            
            r[i][d] += v[i][d] * dt;
    
    Acceleration(r, a, pot);
    HalfStep(v, a);
}

/* Mezzo passo per il velocity Verlet */

void HalfStep (double **v, double **a){
    
    int i, d;
    
    for (i=0; i<Natoms; i++)
        for (d=0; d<DIM; d++)  v[i][d] += 0.5 * dt * a[i][d];
}
/****************************************************************************************/

/* La funzione Acceleration è il cuore del programma, la scatola viene divisa in celle con lunghezza maggiore del raggio di cutoff, le celle vengolo linkate con delle liste e in questo modo data una particella di head si passa alle successive e se ne considerano le interazioni con le 26 celle adiacenti e con se stessa. Si calcola l'energia potenziale del sistema */
/****************************************************************************************/

void Acceleration(double **r, double **a, double *pot){
    
    
    int i, j, d;
    double d2, d2inv, d6inv, f, U, ecut, dr[DIM];
    
    
    ecut = 4 * ((1 / ( pow (Rcut, 12))) - (1 / pow(Rcut,6)));
    U = 0.0;

    for (i=0;i<Natoms;i++)
        for (d=0;d<DIM;d++)
            
            a[i][d] = 0.0;

    
    for (i=0;i<Natoms;i++)
        for (j=i+1;j<Natoms;j++){
                            
                for (d2=0.0, d=0; d<DIM; d++){
                    
                    dr[d] = r[i][d] - r[j][d];
                    dr[d] -= floor((dr[d] / L) + 0.5) * L;
                    
                    d2 += dr[d] * dr[d];
                }
                
                if (d2 < Rcut * Rcut){
                    
                    d2inv = 1. / d2;
                    d6inv = d2inv * d2inv * d2inv;
                    f = 24. * d2inv * d6inv * (2. * d6inv -1.);
                    
                    for (d=0; d<DIM; d++){
                        
                        a[i][d] += f * dr[d];
                        a[j][d] -= f * dr[d];
                    }
                    
                    U += 4. * d6inv * (d6inv - 1.) - ecut;
                                        
                } /* endif d2 */
        }
     
    *pot = U;

}

/* In Observables vengono calcolati e stampati gli osservabili del sistema */

void Observables (double **v, double *kin, double *pot, double mass, int stepCount, FILE *fp){
    
    double kinE=0., ETOT;
    int i, d;
    
    for (i=0; i<Natoms; i++)
        for (d=0; d<DIM; d++)
            
            kinE += v[i][d] * v[i][d];
    
    kinE *= 0.5 * mass;
    
    *kin = kinE;
    
    ETOT = *pot + *kin;
    
    //printf ("TSTEP = %d\tPOT = %lf\tKIN = %lf\tTOT = %lf\n", stepCount, *pot, *kin, ETOT);
    fprintf (fp, "%.15lf\n",ETOT);
    
    
    
}
/***************************************************************************************/
/* the function generates normal distributed number with Box-Muller algorithm */

double Gaussian (void){
    
    double x, y, s;
    
    x = ((double)lrand48()/RAND_MAX);
    y = ((double)lrand48()/RAND_MAX);
    s = sqrt (-2.0 * log(x)) * cos (2.0 * PI * y);
    
    return s;
    
}

/******************************************************************************/


double instantaneousTemperature (double **v) {
    
    double sum = 0;
    int i, d;
    
    for (i = 0; i < Natoms; i++)
        for (d = 0; d < DIM; d++)
            
            sum += v[i][d] * v[i][d];
    
    return sum / (3 * Natoms);
}

/******************************************************************************/

