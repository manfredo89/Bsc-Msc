#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define n 10000 //numero di punti array spaziale
#define V0 -64. /*Profondità della buca in eV*/
#define Estep 3 //incremento di energia in eV per la ricerca degli autovalori (per una ricerca piu accurata è necessario infittire)

double cond (double x,int i,double E0);
double k (int i,double E0);
double Psi (double *array, double a,double dx,double E0);

main(){
           double a,dx,c,E0,DE,fE,fdE,E1;
           int count,i,flag,j;
           double *array,check,Erip,A;
           float B;

           
           FILE *uscita;
           
           array=(double *)malloc(n*sizeof(double));
           uscita=fopen("numerov.dat","w");
           
           a=0.195;  /*semilarghezza della buca in nm*/
           dx=(4*a)/n; /*dimensione dell' unità della mesh spaziale*/
           E0=-64.0;  /*energia di partenza*/
           DE=0.00001;
           flag=0;
           Erip=0;
           
           j=1;
           do{
           
           count=1;
           do{
           fE=Psi(array,a,dx,E0);  /*calcolo il valore del punto di controllo (3/4 a)*/
           fdE=(Psi(array,a,dx,(E0+DE))-Psi(array,a,dx,E0))/DE;  /*calcolo la sua derivata prima rispetto all'energia per newton raphson */
           E1=((E0-fE/fdE)); /*metodo di newton-raphson*/
           check=(array[n/4*3-1]-array[n/4*3])*(array[n/4*3+1]+array[n/4*3+2]); /*la funzione check mi assicura che le due funzioni abbiano la derivata spaziare dello stesso segno*/                   
            //printf("e1 is:%lf\n",E1);
           if(fabs((E1-E0)/E1)<0.001 && check>0)           /*soglia di accettazione*/
           { count=20; flag=1;
           printf("E1 is:%f\n",E1); 
           } 

           
           E0=E1;                           /*aggiornamento del valore dell'energia calcolato*/      
           count=count+1; 
           }while(count<20);
           
           fE=Psi(array,a,dx,E1);
           if(flag==1 && Erip!=E1){
           A=0;/*reset valore integrale*/
           /*integrale tramite metodo dei trapezi per la normalizzazione*/
           for(i=0;i<n-1;i++){A=A+(array[i]+array[i+1])*(array[i]+array[i+1]);}//somma modulo quadro
           A=sqrt(A*0.5*dx);//radice quadrata e moltiplicazione per il passo reticolare
           /*stampa*/                        
           for(i=0;i<n;i++){
             array[i]=array[i]/A;//normalizzazione
             fprintf(uscita,"%lf %lf\n",i*dx-2*a,array[i]); 
             }
             fprintf(uscita,"\n");
             }
           Erip=E1; /*evita di stampare 2 volte la stessa energia su file*/
           E0=V0+j*Estep;//incremento di # eV l'energia per la ricerca di un nuovo autovalore
           //printf("E0 is:%f\n",E0);
           flag=0;
           j+=1;
           }while(E0<0);
           
printf("numerov.dat\n"); 
fclose(uscita);
           
  

}

/*Algoritmo di numerov*/
double Psi (double *array, double a,double dx,double E0){
       int i;
       double c;
    /*calcolo della funzione d'onda in una regione esterna alla buca in questo caso in x=2a ed in x=-2a*/   
   array[0]=cond(-2*a,0,E0); 
   array[1]=cond(-2*a+dx,0,E0);
  /* array[0]=-cond(-2*a,0,E0);     per soluzioni dispari
   array[1]=-cond(-2*a+dx,0,E0);*/
   array[n-1]=cond(-2*a,n-1,E0);
   array[n-2]=cond(-2*a+dx,n-1,E0);
   
   /*Soluzione destra che a partire dalle condizioni fissate agli estremi a -2a implementa la soluzione fino a 3/4 a*/
    for(i=1;i<(n/4.*3);i++){
             array[i+1]=(2*array[i]-array[i-1]-(5./6.)*k(i,E0)*array[i]*(dx*dx)-k(i-1,E0)*array[i-1]*(((dx*dx))/12.))/(1+((dx*dx)*k(i+1,E0))/12.);
             }
    /*Soluzione sinistra che a partire dalle condizioni fissate agli estremi a 2a implementa la soluzione fino a 3/4 a*/
   for(i=n-2;i>(n/4*3+1);i--){
             array[i-1]=(2*array[i]-array[i+1]-(5./6.)*k(i,E0)*array[i]*(dx*dx)-k(i+1,E0)*array[i+1]*(((dx*dx))/12.))/(1+((dx*dx)*k(i-1,E0))/12.);
             }

c=(array[n/4*3]-array[n/4*3+1])/array[n/4*3+1];
/*C è il punto di controllo, se è zero, la condizione di continuità della derivata prima è soddisfatta*/

return c;
}

/*funzione che determina il valore della funzione d'onda agli estremi della mesh*/
double cond (double x,int i, double E0){
       double P,t;
       t=sqrt(-2*13.1234*E0);
       P=exp(t*x);
       
       return P;
}   
/*funzione K(X)*/
double k (int i, double E){
       double s;
              if((i<n/4-1) || i>(3*n/4-1)){
                       s=2*13.1234*E;
                       }else{
                   s=2*13.1234*(-V0+E);
                   }
           return s; 
}
      
