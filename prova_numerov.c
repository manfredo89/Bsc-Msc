#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define mesh 2000        /* impostare una mesh con modulo4 = 0 */
#define a 0.195          /* larghezza di un quarto di intervallo */
#define mh2 13.1234      /* massa / (costante di planck / 2PI)^2 in unità ridotte */
#define cl 1E-7         /* intervallo di confidenza per la ricerca degli autivalori */
#define epass 0.1        /* divisione dell'intevallo energetico */
#define POT -64

    struct doublepoint {
        double x1;
        double x2;
    };

    struct doublepoint *wfunction (double, int, FILE *);

    struct doublepoint *wfunction (double e, int parity, FILE *fp){

    int i;
    double v[mesh] = {0.}, k[mesh];
    double dx, ddx12, xmax, norm=0.;
    double yl[(mesh/2)], yr[(mesh/2)], y[mesh];
    static struct doublepoint contder;

    xmax = 4*a;                 /* lunghezza totale intervallo */
    dx = xmax / mesh;           /* passo di integrazione */
    ddx12 = dx * dx / 12.;      /* costante per l'algoritmo di Numerov */


        for (i = (mesh/4); i < (3*mesh)/4; ++i) v[i] = POT;       /* potenziale v[i] per ogni posizione */

        for (i = 0; i < mesh; ++i) k[i] = 2. * (e - v[i])*mh2;    /* funzione k */

        yl[0] = exp(sqrt(-2*e*mh2)*(-2*a));                 /* condizioni iniziali: la funzione d'onda */
        yl[1] = exp(sqrt(-2*e*mh2)*(-2*a + dx));            /* tende esponenzialmente a zero agli estremi */
        yr[(mesh/2)-1] = parity*exp(sqrt(-2*e*mh2)*(-2*a));
        yr[(mesh/2)-2] = parity*exp(sqrt(-2*e*mh2)*(-2*a + dx));

        /* calcolo della funzione d'onda left */
        for (i=1; i < (mesh/2)-1; i++)
        yl[i+1] = (2*yl[i] - yl[i-1] - ddx12*(k[i-1]*yl[i-1] + 10*yl[i]*k[i]))/(1+ddx12*k[i+1]);
        
        /* calcolo della funzione d'onda right */
        for (i=(mesh/2)-2; i>=0; i--)
        yr[i-1] = (2*yr[i] - yr[i+1] - ddx12*(k[(mesh/2)+1+i]*yr[i+1] + 10*k[(mesh/2)+i]*yr[i]))/(1+ddx12*k[(mesh/2)+i-1]);

        /* normalizzazione */
        for (i=0; i<(mesh/2); i++) norm += fabs(yl[i]) + fabs(yr[i]);
        
        for (i=0; i<mesh; i++) {              /* inserisco le due funzioni left e right in uno stesso array e plotto*/

            if (i<mesh/2) y[i] = yl[i]/norm;
            else y[i] = yr[i - (mesh/2)]/norm;
        }

        contder.x1 = copysign (yl[(mesh/2)-1] - yr[0], yl[(mesh/2)-1]);    /*  */
        contder.x2 = copysign(((yl[(mesh/2)-1] - yl[(mesh/2)-2])/dx)-((yr[1]-yr[0])/dx),   /*  */
                              ((yl[(mesh/2)-1] - yl[(mesh/2)-2])/dx));

        /* la struct contiene: 1) differenza tra left e right --> continuità
                               2) differenza tra derivate left e right --> derivata continua
                                  il segno della derivata left viene passato per usare il metodo di bisezione*/
        
        for (i=0;i<mesh;i++)
                fprintf (fp, "\n%lf\t%lf", i*dx, y[i]);
                fprintf (fp, "\n\n");

        return &contder;

    }


    int main(){

    int parity;
    double e;
    struct doublepoint contder;
        
        FILE *fp;
        
        fp = fopen ("provanumerov.dat","w+");
        
        e = -55.0;
        parity = -1;
        
        contder = *wfunction (e, parity, fp);


    
        fclose (fp);

    }
