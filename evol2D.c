/***********************************************************************
 *  pack.c                                                              *
 *  programma per la soluzione dell'equazione di Schroedinger           *
 *  dipendente dal tempo per un pacchetto gaussiano                     *
 *  compilare con cc -o pack pack.c -lm                                 *
 *  uscita per grafico 3D gnuplot nei files                             *
 *  repsi.dat, impsi.dat pro.dat                                        *
 *  rispettivamente per la parte reale, immaginaria ed il modulo quadro *
 *  delle funzione d'onda                                               *
 ***********************************************************************/

#include <stdio.h>
#include <math.h>

#define PI     3.1415926
#define NX     100
#define NY     100
#define x0    -150.           // Posizione iniziale del pacchetto
#define y0    -150.
#define w      1.            // Larghezza del pacchetto
#define kx0    2.            // Velocita' di gruppo del pacchetto
#define ky0    2.
#define xmin  -200.           // Valore minimo di x
#define xmax   200.           // Valore massimo di x
#define ymin   -200.
#define ymax   200.
#define Nt     5000        // Numero di punti temporali
#define dt     0.1        // intervallo temporale
#define stampa 100          // intervallo di stampa su file

double potential (int, int);

int main(){
    
    int i, j, t, k=0;
    double x, y, dx, dx2, dy, dy2, alfa;
    double aux, auy, b, R[NX][NY], I[NX][NY], Rnew[NX][NY], Inew[NX][NY], P1, P2;
    FILE *pro, *repsi, *impsi;
    
    
    /* PARAMETRI */
    
    dx = (double) (xmax - xmin) / NX; // intervallo x
    dy = (double) (ymax - ymin) / NY; // intervallo y
    //dx=dy=0.5;
    dx2 = dx * dx;
    dy2 = dy * dy;
    
    repsi = fopen("repsi.dat","w");
    impsi = fopen("impsi.dat","w");
    pro = fopen("mod.dat","w");
    
    /*** VALORI INIZIALI ***/
    
    for (i=0;i<NX;i++) R[i][0] = R[i][NY-1] = I[i][0] = I[i][NY-1] = 0.;
    for (j=0;j<NY;j++) R[0][j] = R[NX-1][j] = I[0][j] = I[NX-1][j] = 0.;
    
    for (i=1;i<NX-1;i++){
    for (j=1; j<NY-1; j++) {
        
        x = xmin + i * dx;
        y = ymin + j * dy;
        aux = x - x0;
        auy = y - y0;
        b = pow(2. * PI * w * w, -0.5);
        
        // Valore iniziale R[t=0]
        R[i][j] = b * cos(kx0 * aux + ky0 * auy) * exp(-(aux * aux + auy * auy) / (4. * w * w));
        
        // Valore iniziale I[t=dt/2]
        I[i][j] = b * sin(kx0 * (aux - 0.25 * kx0 * dt) + ky0 * (auy - 0.25 * ky0 * dt)) * exp(-((aux - 0.5 * kx0 * dt) * (aux - 0.5 * kx0 * dt) + (auy - 0.5 * ky0 * dt) * (auy - 0.5 * ky0 * dt)) / (4. * w * w));
        
    }
    }
    
    
    for (i=0;i<NX;i++){
    for (j=0;j<NY;j++) {
    
        fprintf (repsi, "%f\n", R[i][j]);
        fprintf (impsi, "%f\n", I[i][j]);
        fprintf (pro, "%lf\t%lf\t%e\n", (i) * dx, j * dy, R[i][j]*R[i][j]+I[i][j]*I[i][j]);
        }
    }
    
    /*** LOOP TEMPORALE ***/
    
    alfa = dt / (2 * dx2);
    
    for (t=0;t<Nt;t++) {
        
        for (i=1;i<NX-1;i++){
        for (j=1; j<NY-1; j++) {
            //Rnew[i][j] = R[i][j] - (0.5 * dt) * ((I[i+1][j] - 2. * I[i][j] + I[i-1][j]) / dx2 + (I[i][j+1] - 2. * I[i][j] + I[i][j-1]) / dy2) + potential(i, j) * I[i][j] * dt;
            Rnew[i][j] = R[i][j] + 2. * ((4. * alfa + dt * potential (i, j)) * I[i][j] - alfa * (I[i+1][j] + I[i-1][j] + I[i][j+1] + I[i][j-1]));
            R[i][j] = Rnew[i][j];
        }}
        
        for (i=1;i<NX-1;i++){
        for (j=1;j<NY-1;j++) {
            //Inew[i][j] = I[i][j] + (0.5 * dt) * ((R[i+1][j] - 2. * R[i][j] + R[i-1][j]) / dx2 + (R[i][j+1] - 2. * R[i][j] + R[i][j-1]) / dy2) - potential(i, j) * R[i][j] * dt;
            Inew[i][j] = I[i][j] - 2. * ((4. * alfa + dt * potential (i, j)) * R[i][j] + alfa * (R[i+1][j] + R[i-1][j] + R[i][j+1] + R[i][j-1]));
            I[i][j] = Inew[i][j];
        }}
        
        if(t%stampa == 0) {
            for (i=1;i<NX-1;i++){
            for (j=1;j<NY-1;j++) {
                
                fprintf(repsi,"%f\n",R[i][j]);
                fprintf(impsi,"%f\n",I[i][j]);
                
                P1 = Rnew[i][j] * R[i][j] + I[i][j] * I[i][j];
                //P2 = Rnew[i][j] * Rnew[i][j] + Inew[i][j] * I[i][j];
                
                fprintf(pro,"%lf\t%lf\t%e\n", (i + k * NX - 1) * dx, j * dy, P1);
                //fprintf(pro,"%lf\t%lf\t%e\n", (i + 0.5 + k * NX - 1) * dx, j * dy, P2);
            }}
            fprintf(repsi,"\n");
            fprintf(impsi,"\n");
            fprintf(pro,"\n\n");
        
            //k++;
        }
        
    }
    
    fclose(repsi);
    fclose(impsi);
    fclose(pro);
}

double potential (int i, int j){
    
    double r2;
    r2=(NX*NX+NY*NY)/(20.*20.);
    
    if(((i-NX/2)*(i-NX/2)+(j-NY/2)*(j-NY/2))<r2){ //potenziale circolare di 1/20 del lato della mesh
        return 0.;}
    else{return 0.;}
}