/***********************************************************************
*  fourier.c                                                           
*  programma per la soluzione dell'equazione di Fourier                
*  compilare con cc -lm fourier.c                                          
*  uscita per grafico 3D gnuplot nel file fourier.dat                   
*  comando: splot 'fourier.dat'          
* 
*  l'equazione del calore di Fourier
* 
* 		rho c dT(x,t)/dt = lambda grad^2 T(x,t)
* 
*  viene discretizzata nel tempo con un metodo di Eulero
* 
* 		T(x,t+dt) = T(x,t) + const grad^2 T(x,t) dt
* 		const = lambda / (rho c)
* 
*  con la discretizzazione spaziale del Laplaciano 
* 
* 		grad^2 T = [T(x+dx) -2 T(x) + T(x-dx)] / (dx*dx) 
* 
*  bada bene: nel programma specifico, dx vale 1.    
*                                         
***********************************************************************/

#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

#define Ngriglia 101  /* dimensioni della griglia spaziale */
#define Ntmax 30000   /* iterazioni temporali */
#define lambda 0.12	  /* conducibilita' termica */
#define calsp 0.113	  /* calore specifico */
#define dens 7.8	  /* densita' */

int main(){ 
  	int i,j;
  	double cost; // costante fisica di integrazione
  	double T[Ngriglia][2]; // temperatura     
  	FILE *uscita;

	/* apri il file di uscita*/
  	uscita = fopen("fourier2.dat","w");
  	
  	/* condizione iniziale: a t=0 tutti i punti sono a 100 C */
  	for(i=0; i<Ngriglia; i++) T[i][0]=100.0 * sin ((double)i * M_PI / Ngriglia);
  	/* eccetto agli estremi che sono a 0 C */
  	for(j=0; j<2; j++) T[0][j] = T[Ngriglia-1][j] = 0.;
	
	/* costanti fisiche */
  	cost=lambda/(calsp*dens);
  	
  	/* loop temporale */	     
  	for(i=1; i<=Ntmax; i++){ 
  		
  		/* T[x][0] contiene la soluzione al tempo t,
  		 * T[x][1] conterra` la soluzione al tempo t+dt;
  		 * il loop implementa la discretizzazione di
  		 * T(x,t+dt) =  T(x,t) + const grad^2 T dt       */
    	for(j=1; j<(Ngriglia-1); j++)              
      		T[j][1] = T[j][0] 
      				+ cost*(T[j+1][0]+T[j-1][0]-2.0*T[j][0]);
      	/* bada bene: nella discretizzazione del Laplaciano
      	 * abbiamo assunto dx = 1 !!!                       */
    
    	/* salva ogni 1000 passi */
    	if((i%1000==0) || (i==1)){              
    		for(j=0 ; j<Ngriglia; j++) // formato 3D per gnuplot 
    			fprintf(uscita, "%f\n", T[j][1]); 
      		fprintf(uscita, "\n"); // linea vuota per gnuplot 
    	}
    	
    	/* copia la soluzione al tempo t+dt in T[x][0] */
    	for(j=0; j<Ngriglia; j++) T[j][0]=T[j][1];  
    
    }
  	
  	fprintf(stderr,"Dati in fourier.dat\n");
  	
  	fclose(uscita);

	return 0;
}
