/* Programma di generazione di numeri pseudocasuali distribuiti in modo gaussiano 
 con l'algoritmo Box-Muller */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#define Natoms 864

double Gaussian();

main (){

    FILE *fp;
    int i, seed = time(0);
    srand48(seed);
	
    fp = fopen ("gaussian.dat","w+");
    
    for (i=0;i<Natoms;i++){        
       
        fprintf (fp,"%lf\t%lf\t%lf\n",Gaussian(), Gaussian(), Gaussian());
	}
    	
  fclose(fp);

}

double Gaussian (void){
    
    double x, y, s;
    
    x = ((double)lrand48()/RAND_MAX);
    y = ((double)lrand48()/RAND_MAX);
    s = sqrt (-2.0 * log(x)) * cos (2.0 * M_PI * y);
    
    return s;
    
}
