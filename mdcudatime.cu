/* Programma di simulazione di dinamica molecolare */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <cuda.h>

#define Rcut 2.5f // cutoff distance
#define DIM 3
#define dt 0.0005f
#define eCut -0.01631689114f // 4 * ((1 / ( pow (Rcut, 12))) - (1 / pow(Rcut,6))) // cutoff energy
#define mass 1.0f
#define PI 3.1415926535f
#define NUM_THREAD 64  // Number of threads per block
#define NUM_BLOCK (int) ceil (Natoms/(float)NUM_THREAD)  // Numb of thread blocks

typedef struct {
    
    float x, y, z; // coordinates
    float w; // free parameter
    
    } M_double4;


/****************************************************************************************/
/* first half kick for Verlet integration */

__global__ void HalfStep (M_double4 *v, M_double4 *a, int Natoms){
    
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    
    if (i < Natoms){
        
        v[i].x += 0.5f * dt * a[i].x;
        v[i].y += 0.5f * dt * a[i].y;
        v[i].z += 0.5f * dt * a[i].z;
        v[i].w = 0.5f * mass * ((v[i].x * v[i].x) + (v[i].y * v[i].y) + (v[i].z * v[i].z)); // kinetic energy
    }
}
/****************************************************************************************/
/* positions updates */

__global__ void Position (M_double4 *r, M_double4 *v, M_double4 *a, int Natoms){
    
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    
    if (i < Natoms){
        
        v[i].x += 0.5f * dt * a[i].x;
        v[i].y += 0.5f * dt * a[i].y;
        v[i].z += 0.5f * dt * a[i].z;
        
        r[i].x += v[i].x * dt;
        r[i].y += v[i].y * dt;
        r[i].z += v[i].z * dt;
        
        }
    
}
/****************************************************************************************/
/* the function generates normal distributed number with Box-Muller algorithm */

float Gaussian (void){
    
    float x, y, s;
    
    x = ((float)lrand48()/RAND_MAX);
    y = ((float)lrand48()/RAND_MAX);
    s = sqrt (-2.0f * log(x)) * cos (2.0f * PI * y);
    
    return s;
    
}

/****************************************************************************************/
/* calculus of forces with Lennard Jones potential and cutoff energy (see Molecular Simulation by D. Frenkel) */

__global__ void Acceleration(M_double4 *r, M_double4 *a, int Natoms, float rho){
    
    
    int i, j, axis;
    float d2, d2inv, d6inv, f, dr[DIM];
    float L = powf(Natoms / rho, 1.0f/3);
    
    i = blockIdx.x*blockDim.x+threadIdx.x;
     
     if (i < Natoms){
    
        a[i].x = 0.0f;
        a[i].y = 0.0f;
        a[i].z = 0.0f;
        a[i].w = 0.0f;
       
        for (j=0;j<Natoms;j++){
            if (i == j) continue;
            
            	    dr[0] = r[i].x - r[j].x;
		    dr[1] = r[i].y - r[j].y;
		    dr[2] = r[i].z - r[j].z;
            
		    dr[0] -= floorf((dr[0] / L) + 0.5f) * L;
		    dr[1] -= floorf((dr[1] / L) + 0.5f) * L;
		    dr[2] -= floorf((dr[2] / L) + 0.5f) * L;
                
                for (d2=0.0f, axis=0; axis<DIM; axis++) d2 += dr[axis] * dr[axis];
            
            if (d2 < Rcut * Rcut){
                    
                d2inv = 1.0f / d2;
                d6inv = d2inv * d2inv * d2inv;
                f = 24.0f * d2inv * d6inv * (2.0f * d6inv -1.0f);
                        
                        a[i].x += f * dr[0];
                        a[i].y += f * dr[1];
                        a[i].z += f * dr[2];
                        
                        a[i].w += 4.0f * d6inv * (d6inv - 1.0f) - eCut;
               
                } /* endif d2 */
            } /* endfor j */
	} /* endfor i */
}

/****************************************************************************************/

void Rescale (M_double4 *v, double T, int Natoms){
    
    int i;
    float vSum2=0.0f, fs;
    
    for (i=0;i<Natoms;i++) vSum2 += (v[i].x * v[i].x) + (v[i].y * v[i].y) + (v[i].z * v[i].z);
    
    fs = sqrt (3.0f * T * Natoms / vSum2); /* scaling factor to set the temperature */
    
    for (i=0;i<Natoms;i++){
        
        v[i].x *= fs;
        v[i].y *= fs;
        v[i].z *= fs;
    }
}

/****************************************************************************************/
/* The function initializes the r (fcc lattice), v, a arrays and sets c.o.m. speed to zero */

void Initialization (M_double4 *r, M_double4 *v, M_double4 *a, float T, int *seed, int Natoms, float rho){
    
    int i, axis, k, M=1, nX, nY, nZ;
    float vSum[DIM]={0.0f};
    float L = pow(Natoms / rho, 1.0f/3); // total box lenght
    float firstCell[4][3] = {
        {0.25f, 0.25f, 0.25f},
        {0.75f, 0.75f, 0.25f},
        {0.75f, 0.25f, 0.75f},
        {0.25f, 0.75f, 0.75f}};
    
    while (4 * M * M * M < Natoms) M++; // M^3 will be the nuber of boxes to contain all the Natoms
    
    float l = L / M; // is the single box lenght
    
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
        
        a[i].x = 0.0f;
        a[i].y = 0.0f;
        a[i].z = 0.0f;
    }
    
    Rescale (v, T, Natoms);

}
/****************************************************************************************/

void Evolution (M_double4 *r, M_double4 *v, M_double4 *a, int Natoms, float rho){
    
    dim3 dimGrid (NUM_BLOCK, 1, 1);
    dim3 dimBlock (NUM_THREAD, 1, 1);
    
    Position <<<dimGrid, dimBlock>>>(r, v, a, Natoms);
    Acceleration <<<dimGrid, dimBlock>>>(r, a, Natoms, rho);
    HalfStep <<<dimGrid, dimBlock>>>(v, a, Natoms);
    
}
/****************************************************************************************/

int main () {
    
    int stepCount, stepLimit, seed, Natoms;
    float T, mtime;
    float rho;
    long seconds, useconds;
    M_double4 *h_r, *h_v, *h_a;
    M_double4 *d_r, *d_v, *d_a;
    struct timeval start, end;
    cudaEvent_t gpu_start, gpu_stop;
    float gpu_runtime;
    FILE *fp;
    
    seed = time(0);
    srand48(seed);
  
    fp = fopen("cuda_time.dat","w+");
    
    stepLimit = 100;
    T = 0.5f;
    rho = 0.2f;
    
    while (rho<1.0f){
    
        for (Natoms=10;Natoms<=1500;Natoms+=100){
    
            size_t size = NUM_BLOCK * NUM_THREAD * sizeof ( M_double4 );
    
        
            /* allocating memory on the host */
            h_r = ( M_double4 * ) malloc ( size );
            h_v = ( M_double4 * ) malloc ( size );
            h_a = ( M_double4 * ) malloc ( size );
    
            /* allocating memory on the device */
            cudaMalloc (&d_r, size);
            cudaMalloc (&d_v, size);
            cudaMalloc (&d_a, size);
    
            Initialization (h_r, h_v, h_a, T, &seed, Natoms, rho);
     
            /* copying data from host to device */
            cudaMemcpy (d_r, h_r, size, cudaMemcpyHostToDevice);
            cudaMemcpy (d_v, h_v, size, cudaMemcpyHostToDevice);
            cudaMemcpy (d_a, h_a, size, cudaMemcpyHostToDevice);
        
            gettimeofday (&start, NULL);
            
            cudaEventCreate (&gpu_start);
            cudaEventCreate (&gpu_stop);
            cudaEventRecord (gpu_start, 0);
        
            /* main cycle */
            for (stepCount=0; stepCount<=stepLimit; stepCount++){
           
                Evolution (d_r, d_v, d_a, Natoms, rho);
            }
            
            cudaEventRecord (gpu_stop, 0);
            cudaEventSynchronize (gpu_stop);
            cudaEventElapsedTime (&gpu_runtime, gpu_start, gpu_stop);
    
            gettimeofday (&end, NULL);
            seconds  = end.tv_sec  - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
            mtime = (float)seconds + (float)useconds/1000000.0f;
        
            fprintf (fp, "%d\t%.5lf\n", Natoms, gpu_runtime /*1000 * mtime*/ / (stepLimit + 1));
    
            free (h_r); // freeing host memory
            free (h_v);
            free (h_a);
        
            cudaFree (d_r); // freeing device memory
            cudaFree (d_v);
            cudaFree (d_a);

        }
    
        fprintf (fp,"\n\n");
        rho += 0.5f;
    }
	
        printf ("\nComputational times in file cuda_time.dat\n", mtime);

        fclose(fp);
    
        
}
/****************************************************************************************/
