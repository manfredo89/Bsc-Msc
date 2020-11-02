#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

main (){
    
    FILE *fp;
    int i;
    double val[2500], tot=0., med, de2;
    
    fp = fopen ("energy1.dat","r");
    
    for (i=0;i<2500;i++){
        fscanf(fp,"%lf",&val[i]);
        tot += val[i];
    }
    
    med = tot / 2500.;
    
    for (i=0;i<2500;i++){
        
        de2 = (val[i] - med) * (val[i] - med) / 2500.;
    }
 
    printf ("varianza %.10lf \n", sqrt(de2));
    fclose (fp);
}