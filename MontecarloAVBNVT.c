#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define AVB 1 //porre a 0.5 per disattivare l'AVB

#define N_ATOMS 200  //numero di particelle del sistema


#define TEMPERATURE 1   //temperatura del sistema
#define LBOX 12          //lato del box
#define M 8          //numero di celle per direzione
#define KBOLTZ 1     //costante di Boltzmann

/*informazioni potenziale Square Well*/
#define DELTA 0.25        //parte attrattiva
#define H -1.              //profondità della buca
#define SIGMA 1             //parte attrattiva

#define DeltaR 0.1    //traslazione massima della particella

double randN(double inizio, double fine);   //genera numeri casuali tra inizio e fine
void disposizione(double coords[][3]);             //dispone le particelle in un reticolo cubico per la configurazione iniziale
void MAPS(int MAP[]);                         //realizza una mappa delle celle vicine con condizioni periodiche a contorno
void Links(int HEAD[], int LIST[], double coords[][3]);        //realizza 
double ENERGYPARTICLE(int N,double coords[][3], int LIST[],int MAP[],int HEAD[],double RX,double RY,double RZ,int *flag); //calcola l'energia di singola particella
double moveMC(int HEAD[],int LIST[],double coords[][3],int MAP[],int *flag,int *mossa_accettata, int *mossa_rifiutata); //funzione del Montecarlo NVT
double AVBbonding(int HEAD[],int LIST[],double coords[][3],int MAP[],int *flag,int *mossa_accettata, int *mossa_rifiutata);    //funzioni per l'AVB
double AVBunbonding(int HEAD[],int LIST[],double coords[][3],int MAP[],int *flag,int *mossa_accettataAVB, int *mossa_rifiutataAVB);
int RANDOMneighbors(int N,double coords[][3], int LIST[],int MAP[],int HEAD[],double RX,double RY,double RZ);


main(void){
    printf("MonteCarlo AVB\n");
    double SystemEnergy,Lato,Energy,a,N;
    int naccept, n,i,MAP[26*M*M*M],randomNumber;
    int HEAD[M*M*M], LIST[N_ATOMS],j,acc,rig,*mossa_accettata,*mossa_rifiutata,*flag,controllo,accAVB,rigAVB,*mossa_accettataAVB,*mossa_rifiutataAVB;
    double coords[N_ATOMS][3];
    FILE *uscita;
    /* apri il file di uscita*/
    uscita = fopen("AVBMC.dat","w");
    
    mossa_accettata=&acc;
    mossa_rifiutata=&rig;
    mossa_accettataAVB=&accAVB;
    mossa_rifiutataAVB=&rigAVB;
    flag=&controllo;
    
    *mossa_accettata=0;
    *mossa_rifiutata=0;
    *mossa_accettataAVB=0;
    *mossa_rifiutataAVB=0;
    controllo=0;
    Lato=LBOX;
    naccept=0;
    srand(time(0));
    disposizione(coords);  //configurazione iniziale delle particelle in un reticolo cubico
    MAPS(MAP);       //linkaggio celle con boundary condition
    Links(HEAD,LIST,coords);  //creazione cell list
    Energy=0;
    for(i=0;i<N_ATOMS;i++){           //calcolo l'energia iniziale chiamando la funzione di energia di singola particella per ogni atomo
    Energy+=ENERGYPARTICLE(i,coords,LIST,MAP,HEAD,coords[i][0],coords[i][1],coords[i][2],flag);
    }
    SystemEnergy=0.5*Energy; //energia del sistema
    printf("energy is:%f\n",SystemEnergy);
      
     for(j=0;j<2000;j++){                                   
    for(i=0;i<2000;i++){ 
                        N=randN(0,AVB);
                                  if(N<=0.5){
                                      a=moveMC(HEAD,LIST,coords,MAP,flag,mossa_accettata,mossa_rifiutata);
                                      SystemEnergy+=a;
                                     }
                                  if(N>0.5 && N<=0.75){
                                      SystemEnergy+=AVBbonding(HEAD,LIST,coords,MAP,flag,mossa_accettataAVB,mossa_rifiutataAVB);
                                     }
                                      if(N>0.75){
                                      SystemEnergy+=AVBunbonding(HEAD,LIST,coords,MAP,flag,mossa_accettataAVB,mossa_rifiutataAVB);
                                     }
            //if(j==0)fprintf(uscita, "%d %f %d %d %d %d\n",i+1000*j,SystemEnergy/N_ATOMS,*mossa_accettata,*mossa_rifiutata,*mossa_accettataAVB,*mossa_rifiutataAVB);
                         }
                          printf("e1 is:%f\n",SystemEnergy/(N_ATOMS+1));
                          
                    fprintf(uscita, "%d\t%f\n",i+2000*j,SystemEnergy/(N_ATOMS+1));    
            }
                           

fclose(uscita);
}


/*generazione numeri casuali*/
double randN(double inizio,double fine){
       return (fine-inizio)*rand()/RAND_MAX+inizio;
}

/*disposizione iniziale*/
void disposizione(double coords[][3]){
      int n3=2,i;
      double ix,iy,iz,L,L2;
      L=LBOX;
      L2=LBOX/2.;
      while ((n3*n3*n3)<N_ATOMS) n3++; //calcola il lato ottimale del reticolo cubico per n atomi

      ix=iy=iz=0;
  /* Assign particle positions */
  for (i=0;i<N_ATOMS;i++) {
    coords[i][0] = ((double)ix+0.5)*L/n3;  //L/n3 è il passo reticolare
    coords[i][1] = ((double)iy+0.5)*L/n3;
    coords[i][2] = ((double)iz+0.5)*L/n3;
    ix++;
    if (ix==n3) {
      ix=0;
      iy++;
      if (iy==n3) {
	iy=0;
	iz++;
      }
    }
  }
}
/*data la terna di coordinate X,Y,Z dei cubi che compongono il box, ICELL restituisce l'indice di cella*/
int ICELL(int IX,int IY,int IZ){
  return  (IX+3*M)%M+((IY+3*M)%M)*M+((IZ+3*M)%M)*M*M;
}
/*l'array map è una lista 26*M*M*M che contiene tutte le celle vicine di una cella*/
void MAPS (int MAP[]){
     int IX,IY,IZ,IMAP;
     for(IX=0;IX<M;IX++){
             for(IY=0;IY<M;IY++){
                    for(IZ=0;IZ<M;IZ++){   
                         IMAP=26*(ICELL(IX,IY,IZ));  //valore della cella centrale è l'indice che punta alla prima locazione di memoria della lista MAP
                                                       // che contiene i primi vicini con condizioni periodiche a contorno
                         
                         MAP[IMAP + 0] = ICELL(IX-1,IY-1,IZ);
                         MAP[IMAP + 1] = ICELL(IX-1,IY,IZ);
                         MAP[IMAP + 2]= ICELL (IX-1,IY+1,IZ);
                         MAP[IMAP + 3] = ICELL (IX+1,IY-1,IZ);
                         MAP[IMAP + 4]= ICELL (IX+1,IY,IZ);
                         MAP[IMAP + 5]= ICELL (IX+1,IY+1,IZ);
                         MAP[IMAP + 6]= ICELL (IX,IY-1,IZ);
                         MAP[IMAP + 7]= ICELL (IX,IY+1,IZ);
                         MAP[IMAP + 8] = ICELL(IX-1,IY-1,IZ-1);
                         MAP[IMAP + 9] = ICELL(IX-1,IY,IZ-1);
                         MAP[IMAP + 10]= ICELL (IX-1,IY+1,IZ-1);
                         MAP[IMAP + 11] = ICELL (IX+1,IY-1,IZ-1);
                         MAP[IMAP + 12]= ICELL (IX+1,IY,IZ-1);
                         MAP[IMAP + 13]= ICELL (IX+1,IY+1,IZ-1);
                         MAP[IMAP + 14]= ICELL (IX,IY-1,IZ-1);
                         MAP[IMAP + 15]= ICELL (IX,IY+1,IZ-1);
                         MAP[IMAP + 16] = ICELL(IX-1,IY-1,IZ+1);
                         MAP[IMAP + 17] = ICELL(IX-1,IY,IZ+1);
                         MAP[IMAP + 18]= ICELL (IX-1,IY+1,IZ+1);
                         MAP[IMAP + 19] = ICELL (IX+1,IY-1,IZ+1);
                         MAP[IMAP + 20]= ICELL (IX+1,IY,IZ+1);
                         MAP[IMAP + 21]= ICELL (IX+1,IY+1,IZ+1);
                         MAP[IMAP + 22]= ICELL (IX,IY-1,IZ+1);
                         MAP[IMAP + 23]= ICELL (IX,IY+1,IZ+1);
                         MAP[IMAP + 24]= ICELL ( IX, IY, IZ +1);
                         MAP[IMAP + 25]= ICELL ( IX, IY, IZ -1);
                    }
             }
     }
    
}

/*funzione che realizza una linked list delle particelle*/
void Links(int HEAD[], int LIST[], double coords[][3]){
     int I,a,b,c,ICELL;
     double CELLI,L2,L;
     L2=LBOX/2.;
     L=LBOX;
     CELLI=M;
     for(ICELL=0;ICELL<(M*M*M);ICELL++){
                                        HEAD[ICELL]=-1;
                                        }
     for(I=0;I<N_ATOMS;I++){
                               a=(int)(((coords[I][0])/(L))*CELLI);
                               b=(int)(((coords[I][1])/(L))*CELLI)*M;
                               c=(int)(((coords[I][2])/(L))*CELLI)*M*M;
                               ICELL=a+b+c; 
                               LIST[I]=HEAD[ICELL];
                               HEAD[ICELL]=I;
                               }
}
/*data ma particella con indice N e le sue coordinate, questa funzione restituisce l'energia di singola particella*/
double ENERGYPARTICLE(int N,double coords[][3], int LIST[],int MAP[],int HEAD[],double RX,double RY,double RZ,int *flag){
       double L2,RXIJ,RYIJ,RZIJ,distanceij,Ep,L,CELLI,epsilon,epsilonHC;
       int ICELL,a,b,c,I,D;
       int JCELL,JCELLO,NABOR,J;
       L2=LBOX/2.;
       L=LBOX;
       CELLI=M;
       epsilon=(SIGMA+DELTA);
       epsilonHC=SIGMA;
       a=(int)(((RX)/(L))*CELLI);
       b=(int)(((RY)/(L))*CELLI)*M;
       c=(int)(((RZ)/(L))*CELLI)*M*M;
       ICELL=a+b+c;                            /*dalle coordinate della particella N questa funzione trova l'indice della cella in cui si trova*/
       Ep=0;
       *flag=0; //funzione di controllo sovrapposizione inizializzata a zero
           I = HEAD[ICELL];                              /*Estra la particella capolista e la confronta*/
                  if((I>-1)&&(I!=N)){
                  RXIJ=coords[I][0]-RX;                  /*N può essere sia il capolista di quella lista, sia essere nell'array List*/
                  RYIJ=coords[I][1]-RY;
                  RZIJ=coords[I][2]-RZ;
                  RXIJ-=rint(RXIJ/L)*L;
                  RYIJ-=rint(RYIJ/L)*L;                      //minimum image convention per la distanza
                  RZIJ-=rint(RZIJ/L)*L;
                  distanceij=sqrt(RXIJ*RXIJ+RYIJ*RYIJ+RZIJ*RZIJ);
                  
                  if((distanceij<epsilon) && (distanceij!=0)){
                                          if(distanceij<=epsilonHC){(*flag)+=1;}
                                          Ep=H+Ep;}}
                                          
                                          
           if(I>-1){
           J=LIST[I];
           while(J>-1){                                   /*confronta N con le altre particelle presenti nella lista*/
           if(J!=N){
           RXIJ=coords[J][0]-RX;
           RYIJ=coords[J][1]-RY;
           RZIJ=coords[J][2]-RZ;
           RXIJ-=rint(RXIJ/L)*L;
           RYIJ-=rint(RYIJ/L)*L;
           RZIJ-=rint(RZIJ/L)*L;
           distanceij=sqrt(RXIJ*RXIJ+RYIJ*RYIJ+RZIJ*RZIJ);
           
           if((distanceij<epsilon)&& (distanceij!=0)){
                                     if(distanceij<=epsilonHC){(*flag)+=1;}
                                          Ep=Ep+H;}}
           
           J=LIST[J];
           }
           }
           JCELLO=26*(ICELL);                                       //Confronto delle distanze con le liste delle celle vicine
           for(NABOR=0;NABOR<26;NABOR++){
               JCELL = MAP[JCELLO + NABOR];
               I = HEAD[JCELL];                              
                if((I>-1)&&(I!=N)){
                  RXIJ=coords[I][0]-RX;
                  RYIJ=coords[I][1]-RY;
                  RZIJ=coords[I][2]-RZ;
                  RXIJ-=rint(RXIJ/L)*L;
                  RYIJ-=rint(RYIJ/L)*L;
                  RZIJ-=rint(RZIJ/L)*L;
                  distanceij=sqrt(RXIJ*RXIJ+RYIJ*RYIJ+RZIJ*RZIJ);
                 
                  if((distanceij<epsilon) && (distanceij!=0)){
                                          if(distanceij<=epsilonHC){(*flag)+=1;}
                                          Ep=H+Ep;}}
           if(I>-1){
           J=LIST[I];
           while(J>-1){
           if(J!=N){
           RXIJ=coords[J][0]-RX;
           RYIJ=coords[J][1]-RY;
           RZIJ=coords[J][2]-RZ;
           RXIJ-=rint(RXIJ/L)*L;
           RYIJ-=rint(RYIJ/L)*L;
           RZIJ-=rint(RZIJ/L)*L;
           distanceij=sqrt(RXIJ*RXIJ+RYIJ*RYIJ+RZIJ*RZIJ);
           
           if((distanceij<(epsilon))&& (distanceij!=0)){
                                          if(distanceij<=epsilonHC){(*flag)+=1;}
                                          Ep=H+Ep;}}
           
           J=LIST[J];
           }
           }
           }
          
return Ep;                                   //restituisce il valore dell'energia
}



double moveMC(int HEAD[],int LIST[],double coords[][3],int MAP[],int *flag,int *mossa_accettata, int *mossa_rifiutata){
       int I,a,b,c,F,ICELL,JCELL,J;
       double NEWX,NEWY,NEWZ,CELLI;
       double OLDX,OLDY,OLDZ,Beta;
       double Eold,Enew,DE,N,distance;
       double ACCTRS,L2,L;
       Beta=1./(TEMPERATURE);
       CELLI=M;
       L2=(LBOX)/2.;
       L=LBOX;
       N=N_ATOMS;
       
       I=(int)(randN(0,N));              //scelta random della particella
       OLDX=coords[I][0];                              
       OLDY=coords[I][1];
       OLDZ=coords[I][2];
       NEWX=coords[I][0]+2*DeltaR*((randN(0,1)-0.5));             //nuove coordinate, spostamento massimo di DeltaR
       NEWY=coords[I][1]+2*DeltaR*((randN(0,1)-0.5));
       NEWZ=coords[I][2]+2*DeltaR*((randN(0,1)-0.5));
       
      NEWX-=floor(NEWX/L)*L;                   //minimum image convention
      NEWY-=floor(NEWY/L)*L;
      NEWZ-=floor(NEWZ/L)*L;
      
       

       
       /*estrazione indice di cella di partenza*/
       a=(int)(((coords[I][0])/(L))*CELLI);
       b=(int)(((coords[I][1])/(L))*CELLI)*M;
       c=(int)(((coords[I][2])/(L))*CELLI)*M*M;
       ICELL=a+b+c;
       
       /*estrazione indice cella di arrivo*/
       a=(int)(((NEWX)/(L))*CELLI);
       b=(int)(((NEWY)/(L))*CELLI)*M;
       c=(int)(((NEWZ)/(L))*CELLI)*M*M;
       JCELL=a+b+c;
       
       
       Eold=ENERGYPARTICLE(I,coords,LIST,MAP,HEAD,OLDX,OLDY,OLDZ,flag);             //energia della particella estratta nella vecchia configirazione
       
       Enew=ENERGYPARTICLE(I,coords,LIST,MAP,HEAD,NEWX,NEWY,NEWZ,flag);               //energia della particella estratta nella nuova configirazione
      
       DE=Enew-Eold;
       ACCTRS=exp(-Beta*DE);                 //calcola soglia di accettazione
       
       N=randN(0,1);                     //estrazione numero random e confronto con la soglia di accettazione
       if(N<ACCTRS && *flag==0){
                    if(ICELL==JCELL){        //se la cella di partenza e arrivo sono le stesse aggiona le coordinate
                    coords[I][0]=NEWX;
                    coords[I][1]=NEWY;
                    coords[I][2]=NEWZ;}else{
                                           if(I==HEAD[ICELL]){                  //se la particella mossa era capolista,fai capolista la particella successiva
                                                 HEAD[ICELL]=LIST[I];
                                            }
                                           else if(I!=HEAD[ICELL]){            //se la particella mossa non era capolista, bisogna cancellarla dalla lista delle 
                                                 J=HEAD[ICELL];                      //particelle in quella cella
                                                 while(J>-1){
                                                 if(LIST[J]==I){
                                                      LIST[J]=LIST[I];
                                                      }      
                                                 J = LIST[J];}
                                            }

                                           LIST[I]=HEAD[JCELL];                //fai comunque capolista la particella mossa
                                           HEAD[JCELL]=I;}
                                           coords[I][0]=NEWX;                    //aggiornamento coordinate
                                           coords[I][1]=NEWY;
                                           coords[I][2]=NEWZ;
                                           (*mossa_accettata)+=1;
                                           
                                           
                   }else{DE=0;
                   (*mossa_rifiutata)+=1;
                   }         //rifiuto mossa
 
return DE;
}           

double AVBbonding(int HEAD[],int LIST[],double coords[][3],int MAP[],int *flag,int *mossa_accettataAVB, int *mossa_rifiutataAVB){
       int a,b,c;
       int I,J,Ni,ICELL,JCELLnew,JCELLold,K,N;
       double POT,Va,V,A,B;
       double r,dr,R,S,T,DE,Beta;
       double RXI,RYI,RZI,RXJ,RYJ,RZJ,RXIJ,RYIJ,RZIJ,CELLI,ACCTRS,distance;
       double RXJnew,RYJnew,RZJnew;
       double EJold,EJnew,distanceij,L,L2;
       double rho,phi,theta;
       
          POT=H;
          V=LBOX*LBOX*LBOX;
          Va=4*M_PI*SIGMA*SIGMA*SIGMA/3.;
          r=SIGMA;
          N=N_ATOMS-1;
          dr=DELTA;
          CELLI=M;
          L2=(LBOX)/2.;
          L=LBOX;
          
                                /*scelta di due particelle I e J con I!=J, e J non deve essere contenuta nelvolume di interazione di I*/
          do{
          do{
          I=(int)(randN(0,N));
          RXI=coords[I][0];
          RYI=coords[I][1];
          RZI=coords[I][2];
          
          J=(int)(randN(0,N));
          RXJ=coords[J][0];
          RYJ=coords[J][1];
          RZJ=coords[J][2];
       
          RXIJ-=rint(RXIJ/L)*L;
          RYIJ-=rint(RYIJ/L)*L;
          RZIJ-=rint(RZIJ/L)*L;
          distanceij=sqrt(RXIJ*RXIJ+RYIJ*RYIJ+RZIJ*RZIJ);}while(distanceij<((SIGMA+DELTA))&& I==J);
          
          
          //calolo il numeri di primi vicini dividendo l'energia per il valore del porenziale, non deve essere maggiore 
          Ni=(int)(ENERGYPARTICLE(I,coords,LIST,MAP,HEAD,RXI,RYI,RZI,flag)/POT);}while(Ni>11); //di 11 ossia del numeri di primi vicini -1 per una sfera
         
         
          EJold=ENERGYPARTICLE(J,coords,LIST,MAP,HEAD,RXJ,RYJ,RZJ,flag);//energia della particella j nella vecchia configurazione
          
          
              rho=r+randN(0,1)*dr;             //disposizione casuale di J nel volume di interesse di I, la scelta random è eseguita in coordinate polari
              theta=randN(0,1)*2*M_PI;         //in una 'buccia' intorno alla particella I
              phi=randN(0,1)*2*M_PI;
              RXJnew=RXI+rho*sin(theta)*cos(phi);
              RYJnew=RYI+rho*sin(theta)*sin(phi);
              RZJnew=RZI+rho*cos(theta);
              
              RXJnew-=floor(RXJnew/L)*L;
              RYJnew-=floor(RYJnew/L)*L;
              RZJnew-=floor(RZJnew/L)*L;
       
       
          
          
          
          EJnew=ENERGYPARTICLE(J,coords,LIST,MAP,HEAD,RXJnew,RYJnew,RZJnew,flag); //energia nuova configurazione
          DE=EJnew-EJold;
          A=Va/((Ni+1)*(4*M_PI*V-Va));
          Beta=1./(KBOLTZ*TEMPERATURE);
          B=exp(-Beta*DE);
          ACCTRS=(N-Ni-1)*A*B;                     //calcolo soglia di accettazione per l'AVB bonding
         
           /*estrazione indice di cella di partenza*/
          a=(int)(((RXI)/(L))*CELLI);
          b=(int)(((RYI)/(L))*CELLI)*M;
          c=(int)(((RZI)/(L))*CELLI)*M*M;
          ICELL=a+b+c;
       
          /*estrazione indice cella di arrivo*/
          a=(int)(((RXJ)/(L))*CELLI);
          b=(int)(((RYJ)/(L))*CELLI)*M;
          c=(int)(((RZJ)/(L))*CELLI)*M*M;
          JCELLold=a+b+c;
          
           /*estrazione indice di cella di partenza*/
          a=(int)(((RXJnew)/(L))*CELLI);
          b=(int)(((RYJnew)/(L))*CELLI)*M;
          c=(int)(((RZJnew)/(L))*CELLI)*M*M;
          JCELLnew=a+b+c;
          T=randN(0,1);
          
          
          if(T<ACCTRS && *flag==0){                                          //aggiornamento Cell list come visto precedentemente
                                           if(J==HEAD[JCELLold]){ 
                                                HEAD[JCELLold]=LIST[J];
                                                LIST[J]=-1;
                                                }
                                           else if(J!=HEAD[JCELLold]){
                                                 K=HEAD[JCELLold];
                                                 while(K>-1){
                                                      if(LIST[K]==J){
                                                      LIST[K]=LIST[J];
                                                      LIST[J]=-1;
                                                       }      
                                                 K = LIST[K];}
                                            }
                
                                           LIST[J]=HEAD[JCELLnew];
                                           HEAD[JCELLnew]=J;
                                           coords[J][0]=RXJnew;
                                           coords[J][1]=RYJnew;
                                           coords[J][2]=RZJnew;
                                           (*mossa_accettataAVB)+=1;
                                           
                                           }else{DE=0;
                                           (*mossa_rifiutataAVB)+=1;}
                                           
return DE;
}

double AVBunbonding(int HEAD[],int LIST[],double coords[][3],int MAP[],int *flag,int *mossa_accettataAVB, int *mossa_rifiutataAVB){
    int a,b,c;
       int I,J,Ni,ICELL,JCELLnew,JCELLold,K,N;
       double POT,Va,V,A,B;
       double S,T,DE,Beta;
       double RXI,RYI,RZI,RXJ,RYJ,RZJ,RXIJ,RYIJ,RZIJ,CELLI,ACCTRS,distance;
       double RXJnew,RYJnew,RZJnew;
       double EJold,EJnew,distanceij,L,L2;
       
          POT=H;
          V=LBOX*LBOX*LBOX;
          Va=4*M_PI*SIGMA*SIGMA*SIGMA/3.;
          N=N_ATOMS-1;
          CELLI=M;
          L=LBOX;
          L2=LBOX/2.;
    do{
    I=(int)(randN(0,N));            //scegli la particella I
          RXI=coords[I][0];
          RYI=coords[I][1];
          RZI=coords[I][2];
          
          a=(int)(((RXI)/(L))*CELLI);
          b=(int)(((RYI)/(L))*CELLI)*M;
          c=(int)(((RZI)/(L))*CELLI)*M*M;
          ICELL=a+b+c;
          if(LIST[I]==-1) return 0; //la particella non è legata
          
    J=RANDOMneighbors(I,coords,LIST,MAP,HEAD,RXI,RYI,RZI);}while(J==-1);  //scegli un suo vicino nel volume di interessa tramite la funzione RANDOMneighbors
    
    RXJ=coords[J][0];  //memorizza vecchie coordinate
    RYJ=coords[J][1];
    RZJ=coords[J][2];
    
    
    
    EJold=ENERGYPARTICLE(J,coords,LIST,MAP,HEAD,RXJ,RYJ,RZJ,flag);
    Ni=(int)(EJold/H); //numero di particelle nel volume di interazione di J
    
    do{
    RXJnew=randN(0,1)*L;       //scegli a caso le nuove nel box
    RYJnew=randN(0,1)*L;
    RZJnew=randN(0,1)*L;
    RXJnew-=floor(RXJnew/L)*L;
    RYJnew-=floor(RYJnew/L)*L;
    RZJnew-=floor(RZJnew/L)*L;
    
    RXIJ=RXJnew-RXI;
    RYIJ=RYJnew-RYI;
    RZIJ=RZJnew-RZI;
    
    RXIJ-=rint(RXIJ/L)*L;
    RYIJ-=rint(RYIJ/L)*L;
    RZIJ-=rint(RZIJ/L)*L;
    distanceij=sqrt(RXIJ*RXIJ+RYIJ*RYIJ+RZIJ*RZIJ);}while(distanceij<((SIGMA+DELTA))); //assicura che le nuove coordinate siano fuori dal volume di interazione di I
    

    EJnew=ENERGYPARTICLE(J,coords,LIST,MAP,HEAD,RXJnew,RYJnew,RZJnew,flag); //nuova enegia della particella j
     
    
    DE=EJnew-EJold;
          
          
          A=(Ni/((N+1-Ni)*Va))*(4*M_PI*V-Va);  //soglia di accettazione mossa
          Beta=1./(KBOLTZ*TEMPERATURE);
          B=exp(-Beta*DE);
          ACCTRS=A*B;
          
       
          /*estrazione indice cella di arrivo*/
          a=(int)(((RXJ)/(L))*CELLI);
          b=(int)(((RYJ)/(L))*CELLI)*M;
          c=(int)(((RZJ)/(L))*CELLI)*M*M;
          JCELLold=a+b+c;
          
           /*estrazione indice di cella di partenza*/
          a=(int)(((RXJnew)/(L))*CELLI);
          b=(int)(((RYJnew)/(L))*CELLI)*M;
          c=(int)(((RZJnew)/(L))*CELLI)*M*M;
          JCELLnew=a+b+c;
          
          T=randN(0,1);
          if(T<ACCTRS && *flag==0){                            //se accetto aggiorno la cell list
                                           if(J==HEAD[JCELLold]){ 
                                                    HEAD[JCELLold]=LIST[J];
                                                    LIST[J]=-1;}
                                            else if(J!=HEAD[JCELLold]){
                                                    K=HEAD[JCELLold];
                                                    while(K>-1){
                                                        if(LIST[K]==J){
                                                        LIST[K]=LIST[J];
                                                        LIST[J]=-1;
                                                         }     
                                                        K = LIST[K];}
                                            }
               
                                           LIST[J]=HEAD[JCELLnew];
                                           HEAD[JCELLnew]=J;
                                           coords[J][0]=RXJnew;
                                           coords[J][1]=RYJnew;
                                           coords[J][2]=RZJnew;
                                           (*mossa_accettataAVB)+=1;
                                           }else{DE=0;(*mossa_rifiutataAVB)+=1;}
                                           
return DE;
    
    
    
          
}          
  /*scelta di un vicino di N a caso*/        
int RANDOMneighbors(int N,double coords[][3], int LIST[],int MAP[],int HEAD[],double RX,double RY,double RZ){
       double L2,RXIJ,RYIJ,RZIJ,distanceij,L,CELLI,epsilon,D;
       int ICELL,a,b,c,I,i;
       int JCELL,JCELLO,NABOR,J,vicini[30];
       L2=LBOX/2.;
       L=LBOX;
       CELLI=M;
       epsilon=(SIGMA+DELTA);
       for(i=0;i<20;i++) vicini[i]=-1;
       
       a=(int)(((RX)/(L))*CELLI);             //dalle coordinate estraggo indice di cella
       b=(int)(((RY)/(L))*CELLI)*M;
       c=(int)(((RZ)/(L))*CELLI)*M*M;
       ICELL=a+b+c;
       i=0;
       
           I = HEAD[ICELL];                              
           if((I>-1)&&(I!=N)){
                  RXIJ=coords[I][0]-coords[N][0];  //confronto con il capolista
                  RYIJ=coords[I][1]-coords[N][1];
                  RZIJ=coords[I][2]-coords[N][2];
                  RXIJ-=rint(RXIJ/L)*L;
                  RYIJ-=rint(RYIJ/L)*L;
                  RZIJ-=rint(RZIJ/L)*L;
                  distanceij=sqrt(RXIJ*RXIJ+RYIJ*RYIJ+RZIJ*RZIJ);
                  if(distanceij<epsilon && distanceij!=0){vicini[i]=I;i+=1;}}  //memorizzo la particella nell'array se le particelle sono nello stesso volume
                                                     //di interazione
                                          
           if(I>-1){
           J=LIST[I];
           while(J>-1){                                 //stessa cosa con le particella della lista e delle particelle delle liste delle altre 
           if(J!=N){                                                       //celle vicine
           RXIJ=coords[J][0]-coords[N][0];
           RYIJ=coords[J][1]-coords[N][1];
           RZIJ=coords[J][2]-coords[N][2];
           RXIJ-=rint(RXIJ/L)*L;
           RYIJ-=rint(RYIJ/L)*L;
           RZIJ-=rint(RZIJ/L)*L;
           distanceij=sqrt(RXIJ*RXIJ+RYIJ*RYIJ+RZIJ*RZIJ);
           if(distanceij<epsilon && distanceij!=0){vicini[i]=J; i+=1;}}
                                    
           
           J=LIST[J];
           }
           }
           JCELLO=16*(ICELL);
           for(NABOR=0;NABOR<26;NABOR++){
               JCELL = MAP[JCELLO + NABOR];
               I = HEAD[JCELL];                              
                 if((I>-1)&&(I!=N)){
                  RXIJ=coords[I][0]-coords[N][0];
                  RYIJ=coords[I][1]-coords[N][1];
                  RZIJ=coords[I][2]-coords[N][2];
                  RXIJ-=rint(RXIJ/L)*L;
                  RYIJ-=rint(RYIJ/L)*L;
                  RZIJ-=rint(RZIJ/L)*L;
                  distanceij=sqrt(RXIJ*RXIJ+RYIJ*RYIJ+RZIJ*RZIJ);
                  if(distanceij<epsilon  && distanceij!=0){vicini[i]=I; i+=1;}}
                                          
           if(I>-1){
           J=LIST[I];
           while(J>-1){
           if(J!=N){
           RXIJ=coords[J][0]-coords[N][0];
           RYIJ=coords[J][1]-coords[N][1];
           RZIJ=coords[J][2]-coords[N][2];
           RXIJ-=rint(RXIJ/L)*L;
           RYIJ-=rint(RYIJ/L)*L;
           RZIJ-=rint(RZIJ/L)*L;
           distanceij=sqrt(RXIJ*RXIJ+RYIJ*RYIJ+RZIJ*RZIJ);
           if(distanceij<epsilon && distanceij!=0){vicini[i]=J;i+=1;}}
                                  
           J=LIST[J];
           }
           }
           }
          
           a=(int)(int)randN(0,1)*(i-1);  //scelgo a caso tra 0 e i, l'indice dell'array che contiene l'indice della particella nel volume di I
           
           J=vicini[a];
           
return J;
} 





