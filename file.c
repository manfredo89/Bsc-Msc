/* Questo programma apre due file last e data ognuno con una colonna i dati
 e copia il contenuto dei due file su un terzo file last disponendo i dati 
 in due colonne */ 

#include <stdio.h>
#include <stdlib.h>

int main (){
    
    FILE *fp0,*fp1,*fp2;
    double d0,d1;
    int N0,N1;
    int i;
    
    fp0=fopen("time.dat","r");
    fp1=fopen("cudatime.dat","r");
    fp2=fopen("speed.dat","w");
    
    for(i=0;i<200;i++){
    
        fscanf(fp0,"%d\t",&N0);
        fscanf(fp0,"%lf",&d0);
        
		fscanf(fp1,"%d\t",&N1);
        fscanf(fp1,"%lf",&d1);
    
    if (N0 == N1) fprintf(fp2,"%d\t%lf\n",N0,d0/d1);
    else printf ("errore, controllare \n");
    
    }
    
    fclose(fp0);
    fclose(fp1);
    fclose(fp2);
}