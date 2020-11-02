


/*Soluzione numerica dell'equazione del calore tramite metodo di Eulero per un sistema a 2 strati, gli strati esterni possiedono le stesse proprietà fisiche,
mentre lo strato centrale ha una conducibilità termica circa la metà*/
#include <stdlib.h>
#include <stdio.h>
#define Ngriglia 101 /* dimensioni di ciascuna delle griglie spaziali, un intervallo spazioale corrisponde ad un cm*/

#define Ntmax 10000000 /* iterazioni temporali, 1 intervallo corrisponde a 10^-4 s*/
#define lambda1 1.4 /* conducibilita' materiale 1 termica W/mK*/
#define calsp1 0.88 /* calore specifico materiale 1 J KG/K^-1*/
#define dens1 1400 /* densita' materiale 1 Kg/m^3 */
#define lambda2 0.8 /* conducibilita' materiale 2 termica W/mK*/
#define calsp2 0.88 /*calore specifico materiale 2 J KG/K^-1*/
#define dens2 1400 /* densita' materiale 2*/
int main(){
int i,j;
double cost1,cost2; // costanti fisiche
double Parete1[Ngriglia][2],Parete2[Ngriglia][2],Parete3[Ngriglia][2]; //i tre diversi array delle temperature delle pareti, 
                                                                         //ogni passo discretocorrisponde ad 1 cm.
FILE *uscita;
/* apri il file di uscita*/
uscita = fopen("fourier.dat","w");

/* condizioni iniziali per i 3 diversi strati in °C*/
for(i=0; i<(Ngriglia); i++) Parete1[i][0]=100.;
for(i=0; i<(Ngriglia); i++) Parete2[i][0]=0.;
for(i=0; i<(Ngriglia); i++) Parete3[i][0]=0.;

/*Condizioni a contorno per le pareti esterne*/
Parete1[0][1]=100.;

Parete3[Ngriglia-1][1]=0;

/* costanti fisiche del materiale 1 e del materiale 2 */
cost1=lambda1/(calsp1*dens1);
cost2=lambda2/(calsp2*dens2);


/*Condizioni di interfaccia iniziali*/

Parete1[Ngriglia-1][0]=Parete2[0][0]=(cost1*Parete1[Ngriglia-2][0]+cost2*Parete2[1][0])/(cost1+cost2);
Parete2[Ngriglia-1][0]=Parete3[0][0]=(cost2*Parete2[Ngriglia-2][0]+cost1*Parete3[1][0])/(cost2+cost1);

/* loop temporale */

for(i=1; i<=Ntmax; i++){

/*metodo di Eulero per la parete 1*/
      for(j=1; j<(Ngriglia-1); j++){
              Parete1[j][1] = Parete1[j][0]+ cost1*(Parete1[j+1][0]+Parete1[j-1][0]-2.0*Parete1[j][0]);
              }
/*metodo di Eulero per la parete 2*/
      for(j=1; j<(Ngriglia-1); j++){
              Parete2[j][1] = Parete2[j][0]+ cost2*(Parete2[j+1][0]+Parete2[j-1][0]-2.0*Parete2[j][0]);
              }
/*metodo di Eulero per la parete 3*/
      for(j=1; j<(Ngriglia-1); j++){
              Parete3[j][1] = Parete3[j][0]+ cost1*(Parete3[j+1][0]+Parete3[j-1][0]-2.0*Parete3[j][0]);
              }

/*aggiornamento condizioni di interfaccia iniziali*/
Parete1[Ngriglia-1][1]=Parete2[0][1]=(cost1*Parete1[Ngriglia-2][1]+cost2*Parete2[1][1])/(cost1+cost2);
Parete2[Ngriglia-1][1]=Parete3[0][1]=(cost2*Parete2[Ngriglia-2][1]+cost1*Parete3[1][1])/(cost2+cost1);

/* salva ogni 200000 passi, che corrispondono a 20 s */
if((i%200000==0) || (i==1)){
         for(j=0 ; j<Ngriglia; j++)
                fprintf(uscita, "%f\n", Parete1[j][1]);
         for(j=0 ; j<Ngriglia; j++)
                fprintf(uscita, "%f\n", Parete2[j][1]);
         for(j=0 ; j<Ngriglia; j++)
                fprintf(uscita, "%f\n", Parete3[j][1]);
         fprintf(uscita, "\n"); 
}
/* copia la soluzione al tempo per ciascuna parete a t+dt  */
for(j=0; j<Ngriglia; j++){ Parete1[j][0]=Parete1[j][1]; Parete3[j][0]=Parete3[j][1];Parete2[j][0]=Parete2[j][1];}

}
fprintf(stderr,"Dati in fourier.dat\n");
fclose(uscita);
return 0;
}
