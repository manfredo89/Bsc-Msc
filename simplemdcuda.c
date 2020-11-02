/* Programma di simulazione di dinamica molecolare */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#define Rcut 2.5 // cutoff distance
#define Natoms 108 // 4*(x^3) --> 32, 108, 256, 500, 864, ...
#define DIM 3
#define dt 0.005
#define rho 0.2 // number density
#define eCut -0.01631689114 // 4 * ((1 / ( pow (Rcut, 12))) - (1 / pow(Rcut,6))) // cutoff energy
#define mass 1.0
#define PI 3.1415926535

typedef struct {
    
    double x, y, z; // coordinates
    double w; // free parameter
    
    } M_double4;


/****************************************************************************************/
/* first half kick for Verlet integration */

void HalfStep (M_double4 *v, M_double4 *a){
    
    int i;
    
    for (i=0;i<Natoms;i++){
        
        v[i].x += 0.5 * dt * a[i].x;
        v[i].y += 0.5 * dt * a[i].y;
        v[i].z += 0.5 * dt * a[i].z;
        v[i].w = 0.5 * mass * ((v[i].x * v[i].x) + (v[i].y * v[i].y) + (v[i].z * v[i].z)); // kinetic energy
    }
}
/****************************************************************************************/
/* positions updates */

void Position (M_double4 *r, M_double4 *v, M_double4 *a){
    
    int i;
    
     for (i=0;i<Natoms;i++){
         
        v[i].x += 0.5 * dt * a[i].x;
        v[i].y += 0.5 * dt * a[i].y;
        v[i].z += 0.5 * dt * a[i].z;
        
        r[i].x += v[i].x * dt;
        r[i].y += v[i].y * dt;
        r[i].z += v[i].z * dt;
    }
}
/****************************************************************************************/
/* the function generates normal distributed number with Box-Muller algorithm */

double Gaussian (void){
    
    double x, y, s;
    
    x = ((double)lrand48()/RAND_MAX);
    y = ((double)lrand48()/RAND_MAX);
    s = sqrt (-2.0 * log(x)) * cos (2.0 * PI * y);
    
    return s;
    
}

/****************************************************************************************/
/* calculus of forces with Lennard Jones potential and cutoff energy (see Molecular Simulation by D. Frenkel) */

void Acceleration(M_double4 *r, M_double4 *a){
    
    int i, j, axis;
    double d2, d2inv, d6inv, f, dr[DIM];
    double L = pow(Natoms / rho, 1.0/3);

     for (i=0;i<Natoms;i++){
    
        a[i].x = 0.0;
        a[i].y = 0.0;
        a[i].z = 0.0;
        a[i].w = 0.0;
     }
       
       for (i=0;i<Natoms;i++){ 
        for (j=i+1;j<Natoms;j++){
            
            dr[0] = r[i].x - r[j].x;
		    dr[1] = r[i].y - r[j].y;
		    dr[2] = r[i].z - r[j].z;
            
		    dr[0] -= floor((dr[0] / L) + 0.5) * L;
		    dr[1] -= floor((dr[1] / L) + 0.5) * L;
		    dr[2] -= floor((dr[2] / L) + 0.5) * L;
            
               for (d2=0.0, axis=0; axis<DIM; axis++) d2 += dr[axis] * dr[axis];
                
            if (d2 < Rcut * Rcut){
                    
                d2inv = 1.0 / d2;
                d6inv = d2inv * d2inv * d2inv;
                f = 24.0 * d2inv * d6inv * (2.0 * d6inv - 1.0);
                        
                        a[i].x += f * dr[0];
                        a[i].y += f * dr[1];
                        a[i].z += f * dr[2];
                        a[j].x -= f * dr[0];
                        a[j].y -= f * dr[1];
                        a[j].z -= f * dr[2];
                
                a[i].w += 4.0 * d6inv * (d6inv - 1.0) - eCut;
                
                  } /* endif d2 */
        }
            } /* endfor j */
}
/****************************************************************************************/

void Rescale (M_double4 *v, double T){
    
    int i;
    double vSum2=0.0, fs;
    
    for (i=0;i<Natoms;i++) vSum2 += (v[i].x * v[i].x) + (v[i].y * v[i].y) + (v[i].z * v[i].z);
    
    fs = sqrt (3.0 * T * Natoms / vSum2); /* scaling factor to set the temperature */
    
    for (i=0;i<Natoms;i++){
        
        v[i].x *= fs;
        v[i].y *= fs;
        v[i].z *= fs;
    }
}


/****************************************************************************************/
/* The function initializes the r (fcc lattice), v, a arrays and sets c.o.m. speed to zero */

void Initialization (M_double4 *r, M_double4 *v, M_double4 *a, double T, int *seed){
    
    int i, axis, k, M=1, nX, nY, nZ;
    double vSum[DIM]={0.0};
    double L = pow(Natoms / rho, 1.0/3);
    double firstCell[4][3] = {
        {0.25, 0.25, 0.25},
        {0.75, 0.75, 0.25},
        {0.75, 0.25, 0.75},
        {0.25, 0.75, 0.75}};
    
    while (4 * M * M * M < Natoms) M++; // M^3 will be the nuber of boxes to contain all the Natoms
    
    double l = L / M; // is the single box lenght
    
    int n = 0;
    
    
    for (nX=0; nX<M;nX++)
        
        for (nY=0; nY<M;nY++)
            
            for (nZ=0; nZ<M;nZ++)
                
                for (k=0; k<4; k++)
                    
                    if (n<Natoms){
                        
                        r[n].x = (nX + firstCell[k][0]) * l;
                        r[n].y = (nY + firstCell[k][1]) * l;
                        r[n].z = (nZ + firstCell[k][2]) * l;
                        
                        n++;
                    }
    
    
    for (i=0;i<Natoms;i++){
        
        v[i].x = Gaussian();
        v[i].y = Gaussian();
        v[i].z = Gaussian();
        
        vSum[0] += v[i].x;
        vSum[1] += v[i].y;
        vSum[2] += v[i].z;
        
    }
    
    for (axis=0;axis<DIM;axis++) vSum[axis] /= Natoms;
    
    for (i=0;i<Natoms;i++){
        
        v[i].x -= vSum[0]; /* total momentum = NULL */
        v[i].y -= vSum[1];
        v[i].z -= vSum[2];
        
        a[i].x = 0.0;
        a[i].y = 0.0;
        a[i].z = 0.0;
    }
    
    Rescale (v, T);
}
/****************************************************************************************/

void Evolution (M_double4 *r, M_double4 *v, M_double4 *a){
    
    Position(r, v, a);
    Acceleration(r, a);
    HalfStep(v, a);
    
}

/****************************************************************************************/

int main () {
    
    int stepCount, stepLimit, stepPrint, seed, i;
    double T, mtime;
    double POT, KIN;
    long seconds, useconds;
    M_double4 *r, *v, *a;
    struct timeval start, end;
    FILE *fp;
    
    gettimeofday(&start, NULL);
    fp = fopen ("energy1.dat","w+");
    seed = time(0);
    srand48(seed);
    
    stepLimit = 50000;
    stepPrint = 20;
    T = 0.001;
    
    /* allocating memory on the host */
    r = ( M_double4 * ) malloc ( Natoms * sizeof ( M_double4 ) );
    v = ( M_double4 * ) malloc ( Natoms * sizeof ( M_double4 ) );
    a = ( M_double4 * ) malloc ( Natoms * sizeof ( M_double4 ) );
    
    Initialization (r, v, a, T, &seed);
    
    /* main cycle */
    for (stepCount=1; stepCount<=stepLimit; stepCount++){
        
        if ((stepCount%stepPrint) == 0){
            
            POT = 0.0;
            KIN = 0.0;
            
            for (i=0;i<Natoms;i++){
                
                POT += a[i].w; // potential energy
                KIN += v[i].w; // kinetic energy
            }
            
           //printf ("TSTEP = %d\tPOT = %lf\tKIN = %lf\tTOT = %lf \n", stepCount, POT, KIN, POT + KIN);
            //fprintf (fp, "%lf\t%lf\t%lf\t%d\n", POT, KIN, POT + KIN, stepCount);
            fprintf (fp, "%.10lf\n",POT+KIN);
            
            }
        
        Evolution (r, v, a);
    }
    
    gettimeofday(&end, NULL);
	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = (double)seconds + (double)useconds/1000000.0;
	
	printf ("execution time is %lf seconds with dt = %lf \n",mtime, dt);

    fclose(fp);
    
    free (r); // freeing host memory
    free (v);
    free (a);
    
}
