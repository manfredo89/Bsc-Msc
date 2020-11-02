#include <stdio.h>
#include <math.h>

#define PI     3.1415926
#define NX     1000
#define x0    -50.           // Posizione iniziale del pacchetto
#define w      10.            // Larghezza del pacchetto
#define k0     5.            // Velocita' di gruppo del pacchetto
#define xmin  -200.           // Valore minimo di x
#define xmax   200.           // Valore massimo di x
#define Nt     40000         // Numero di punti temporali
#define dt     0.01         // intervallo temporale
#define stampa 1000         // intervallo di stampa su file

double potential (int);

int main(){
    
    int i,j,k=0;
    double x,dx,dx2;
    double aux,bux,R[NX]={0.},I[NX]={0.},Rnew[NX]={0.},Inew[NX]={0.},P1,P2;
    FILE *mod,*repsi,*impsi;
    
    /* PARAMETRI */
    
    dx  =(xmax-xmin)/NX; // intervallo x
    dx2 =dx*dx;
    
    repsi=fopen("repsi.dat","w");
    impsi=fopen("impsi.dat","w");
    mod=fopen("mod.dat","w");
    
    /*** VALORI INIZIALI ***/
    
    fprintf(repsi,"%lf\t%lf\n", 0., R[0]);
    fprintf(impsi,"%lf\t%lf\n", dx/2, I[0]);
    
    for(j=1; j<NX-1; j++) {
        
        x=xmin+j*dx;
        aux=x-x0;
        bux=pow(2.*PI*w*w,-0.25);
        
        // Valore iniziale R[t=0]
        R[j]=bux*cos(k0*aux)*exp(-aux*aux/(4.*w*w));
        
        // Valore iniziale I[t=dt/2]
        I[j]=bux*sin(k0*(aux-0.5*k0*dt))*exp(-(aux-0.5*k0*dt)*(aux-0.5*k0*dt)/(4.*w*w));
        
        fprintf(repsi,"%lf\t%lf\n", j*dx, R[j]);
        fprintf(impsi,"%lf\t%lf\n", (2*j+1)/2*dx, I[j]);
    }
    
    fprintf(repsi,"%lf\t%lf\n", 998*dx, R[NX-1]);
    fprintf(impsi,"%lf\t%lf\n", (1999/2)*dx, I[NX-1]);
    
    /*** LOOP TEMPORALE ***/
    
    for(i=0; i<Nt; i++) {
        
        for(j=2; j<NX; j++) {
            Rnew[1]= R[1] - 0.5*(I[2]-2.*I[1]+I[0])/dx2*dt;
            R[j-1]=Rnew[j-1];
            Rnew[j]=R[j]-0.5*(I[j+1]-2.*I[j]+I[j-1])/dx2*dt + potential(j) * I[j] * dt;
        }
        
        for(j=2; j<NX; j++) {
            Inew[1]= I[1] + 0.5*(R[2]-2.*R[1]+R[0])/dx2*dt;
            I[j-1]=Inew[j-1];
            Inew[j]=I[j]+0.5*(R[j+1]-2.*R[j]+R[j-1])/dx2*dt - potential(j) * R[j] * dt;
        }
        
        if(i%stampa == 0) {
            
            
            for(j=1; j<NX-1; j++) {
                
                fprintf(repsi,"%lf\t%lf\n", k*(xmax-xmin)+j*dx, R[j]);
                fprintf(impsi,"%lf\t%lf\n", k*(xmax-xmin)+dx*(2*j+1)/2, I[j]);
                
                
                if (Rnew[j]*R[j]<0.) Rnew[j]=copysign(Rnew[j],R[j]);
                if (Inew[j]*I[j]<0.) Inew[j]=copysign(Inew[j],I[j]);
                
                P1=Rnew[j]*R[j]+I[j]*I[j];
                P2=R[j]*R[j]+Inew[j]*I[j];
                
                fprintf(mod,"%lf\t%lf\t%lf\n", k*(xmax-xmin)+j*dx, sqrt(P2), potential(j));
                fprintf(mod,"%lf\t%lf\t%lf\n", k*(xmax-xmin)+dx*(2*j+1)/2, sqrt(P1), potential(j));
            }
            
            k++;
        }
        
    }
    
    fclose(repsi);
    fclose(impsi);
    fclose(mod);
}

double potential (int i){
    
    if ((i<499)||(i>501)) return 0.0;
    else return 0.0;
}
