/* Questo programma consente di calcolare gli autovalori energetici di un sistema tipo buca di potenziale
 finita mediante l'algoritmo di Numerov. Compilare con gcc -lm -o numerov Numerov.c 
 Viene generato un file "numerov.dat" con le autofunzioni trovate, è possibile visualizzarle con gnuplot
 con il comando plot "numerov.dat" index i con i numero dell'autovalore (ad esempio per lo stato
 fondamentale digitare plot "numerov.dat" index 0 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define mesh 10000        /* impostare una mesh con modulo4 = 0 */
#define a 0.195          /* larghezza di un quarto di intervallo */
#define mh2 13.1234      /* massa / (costante di planck / 2PI)^2 in unit‡ ridotte */
#define cl 1e-9          /* intervallo di confidenza per la ricerca degli autivalori */
#define epass 0.01        /* divisione dell'intevallo energetico */
#define mc 1200          /* matching point */

    double wfunction (double, int, double, FILE *);

    double wfunction (double e, int parity, double POT, FILE *fp){

    int i;
    double v[mesh] = {0.}, k[mesh];
    double dx, ddx12, xmax, norm=0.;
    double yl[mc], yr[mesh - mc], y[mesh], control;

    xmax = 4*a;                 /* lunghezza totale intervallo */
    dx = xmax / mesh;           /* passo di integrazione */
    ddx12 = dx * dx / 12.;      /* costante per l'algoritmo di Numerov */


        for (i = (mesh/4); i < (3*mesh)/4; ++i) v[i] = POT;       /* potenziale v[i] per ogni posizione */

        for (i = 0; i < mesh; ++i) k[i] = 2. * (e - v[i])*mh2;    /* funzione k */

        yl[0] = exp(sqrt(-2*e*mh2)*(-2*a));                 /* condizioni iniziali: la funzione d'onda */
        yl[1] = exp(sqrt(-2*e*mh2)*(-2*a + dx));            /* tende esponenzialmente a zero agli estremi */
        yr[mesh - mc - 1] = parity*exp(sqrt(-2*e*mh2)*(-2*a));
        yr[mesh - mc - 2] = parity*exp(sqrt(-2*e*mh2)*(-2*a + dx));

        /* calcolo della funzione d'onda left */
        for (i=1; i < (mc - 1); i++){
        yl[i+1] = (2*yl[i] - yl[i-1] - ddx12*(k[i-1]*yl[i-1] + 10*yl[i]*k[i]))/(1+ddx12*k[i+1]);
            fprintf (fp, "\n%lf\t%lf", i*dx, yl[i+1]);
}
        
        /* calcolo della funzione d'onda right */
        for (i=(mesh - mc - 2); i>0; i--){
        yr[i-1] = (2*yr[i] - yr[i+1] - ddx12*(k[mc+1+i]*yr[i+1] + 10*k[mc+i]*yr[i]))/(1+ddx12*k[mc+i-1]);
            fprintf (fp, "\n%lf\t%lf", (i+mc)*dx, yr[i-1]);
        }
        fprintf (fp, "\n\n");
        
        //printf ("e=%lf yl[mc-1]= %lf yr[0]= %lf diff= %lf\n",e,yl[mc-1],yr[0], yl[mc-1]-yr[0]);

        control = (yl[mc-1]-yr[0]);//(((yl[mc-1]-yl[mc-2])/yl[mc-1])-((yr[1]-yr[0])/y[0]));///fabs(((yl[mc-1]-yl[mc-2])/yl[mc-1])+((yr[1]-yr[0])/y[0]));

        
        /* normalizzazione */
        for (i=0; i<mc; i++) norm += fabs(yl[i]);
        for (i=0; i<(mesh-mc); i++) norm += fabs(yr[i]);
        
        //printf ("-----------------e= %.10lf--------------------------------------------------\n", e);
        
        for (i=0; i<mesh; i++) {              /* inserisco le due funzioni left e right in uno stesso array e plotto*/
            
            if (i<mc) y[i] = yl[i]/norm;
            else y[i] = yr[i - mc]/norm;
            //if ((i==mc-2)||(i==mc-1)||(i==mc)||(i==mc+1)) printf ("y[%d] = %.15lf \n", i, y[i]);
        }

        /*if (fabs(control)<=cl){
            for (i=0;i<mesh;i++)
                fprintf (fp, "\n%lf\t%lf", i*dx, y[i]);

            fprintf (fp, "\n\n");}*/

        return control;

    }


    int main(){

    int i=0, parity;
    double e, emin, emax, val, POT, contder, contder2;
        
        FILE *fp;
        
        fp = fopen ("numerov.dat","w+");
        
        printf ("\n");
        printf ("Programma per la ricerca degli autovalori di una buca di potenziale finita\n\n");
        printf ("Inserire la profondita' della buca di potenziale (-64.0) Vo = ");
        scanf ("%lf", &POT);
        
        e = POT + epass;    /* parto da un valore dell'energia simile a quello del potenziale */

        while (fabs(e) > epass){

            emin = e;
            emax = emin + epass;

            parity = 2 * ((i + 1) % 2) - 1;

            contder = wfunction (emin, parity, POT, fp);
            contder2 = wfunction (emax, parity, POT, fp);
            //printf ("contder= %lf contder2= %lf product= %lf\n",contder,contder2, contder2*contder);

            if (contder*contder2<0){

                val = emax; /* conservo il limite dell'intervallo */

                while (fabs(contder)>cl){

                    e = (emin + emax)/2.;        /* bisezione */
                    contder = wfunction (e, parity, POT, fp);
                    contder2 = wfunction (emax, parity, POT, fp);

                    if (contder*contder2<0.) emin = e;
                    else emax = e;

                }
                
                printf ("l'autovalore corrispondente al livello energetico %i e' %lf\n ", i, e);
                i++;        /* il primo ciclo corrisponde al ground state, quelli successivi a stati eccitati */
                e = val + epass;    /* ricomincio a suddividere l'intervallo energetico da dove avevo iniziato */
        
            }

            else e+=epass;    /* se le derivate agli estremi dell'intervallo hanno segno concorde passo oltre */

        }
        
        fclose (fp);

    }
