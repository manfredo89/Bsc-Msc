# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <math.h>

void main (){
    
    FILE *fp;
    double val[1000], sigma2=0.;
    int i;
    
    fp = fopen ("drift.dat","r");
    
    for (i=0;i<1000;i++){
    
        fscanf(fp,"%lf",&val[i]);
    }
    
    for (i=1;i<1000;i++){
        sigma2 += (val[0] - val[i]) * (val[0] - val[i]);
    }
    printf ("varianza = %.10lf\n", sqrt(sigma2/1000));
    fclose (fp);

}