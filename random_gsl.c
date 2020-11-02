/* Programma di simulazione di dinamica molecolare */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define Natoms 3864


main () {
    
    double p[Natoms];
    int i, seed = time(0);
    srand(seed);
    FILE *fp;
    fp = fopen ("gaussian.dat","w+");
    
    gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
    
    for (i=0;i<Natoms;i++){
        
        p[i] = gsl_ran_gaussian(r,1.0);
    }
    for (i=0;i<Natoms;i++) fprintf (fp,"%d\t%d\t%lf \n",i,i,p[i]);

    fclose (fp);

}
