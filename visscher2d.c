/********************************************************************
 * Solution for time dependent Schrodinger equation *
 * for a two dimensional gaussian wave packet entering a *
 * double slit (Young experiment) *
 * slits have 5 units width separated 3 units *
 * initial condition gaussian wave packet *
 * ps(x,0)=exp(i*k0*x) exp(-0.5(x-x0)^2/2sig0^2)) *
 * *exp(-0.5(y-y0)^2/2sig0^2)) *
 * centered at x0=0.0 and y0=-7.0 *
 * see Visscher Computers in Physics 5, 596 nov/dic, (1991) *
 * Askar and Cakmak J.Chem Phys 68,N06,2794-2798,(1978) *
 * psr real part wf. defined at t=0,dt,2st,... *
 * psi imag part wf. defined at t=0.dt, 1.5dt,2.5dt,... *
 * ps(i,1): old iteration in time *
 * ps(i,2): new iteration in time *
 *******************************************************************/

#include "stdio.h"
#include "math.h"

#define D 500

int main()
{
    double psr[91][91][2], psi[91][91][2], v[91][91], p2[91][91];
    int nit; void initial(double psr[][91][2], double psi[][91][2]);
    void potential(double v[][91]);
    void solution(double v[][91], double psr[][91][2], double psi[][91][2], double p2[][91], int nit);
    
    /* input a positive integer which is proportional to the time
     you want to see the position of the wave packet.	 */
    /*printf(" Enter a positive integer from 1(initial time)\n");*/
    /*printf("to 1800 to get wave packet position at that time:\n");*/
    scanf("%d", &nit);
    
    /* initializes the constant values and the wave packet */
    initial(psr,psi);
    
    /* two-dimensional slit potential multiplied by dt is set */
    potential(v);
    
    /* time dependent Schrodinger equation is solved */
    solution(v,psr,psi,p2,nit);
}


void initial(double psr[][91][2], double psi[][91][2])
{
    double a1, a2, a4, dtx, yy, y, x, excr, exci;
    double time, sig0, dx, dx2, k0x, k0y, dy, x0, y0, dt, xx;
    int i, j, nnx, nny;
    
    time=0.;
    sig0=1.2;
    dx=0.2;
    dx2=dx*dx;
    k0x=0.;
    k0y=2.5;
    dy=dx;
    dt=0.0025;
    dtx=dt/dx2;
    x0=0.;
    y0=-7.;
    nnx=89;
    nny=89;
    yy=-9.;
    x=0.;
    
    for(j=0;j<=nny+1;j++){
        xx=-9.;
        for(i=0;i<=nnx+1;i++){
            y=k0x*xx+k0y*yy;
            excr=cos(y);
            a1=((xx-x0)/sig0)*((xx-x0)/sig0);
            a2=((yy-y0)/sig0)*((yy-y0)/sig0);
            a4=exp(-0.5*(a1+a2));
            psr[i][j][0]=a4*excr;
            xx=xx+dx;
        }
        yy=yy+dy;
    }
    time=time+dt/2.;
    yy=-9.;
    for(j=0;j<=nny+1;j++){
        xx=-9.;
        for(i=0;i<=nnx+1;i++){
            y=k0x*xx+k0y*yy;
            exci=sin(y);
            a1=((xx-x0)/sig0)*((xx-x0)/sig0);
            a2=((yy-y0)/sig0)*((yy-y0)/sig0);
            a4=exp(-0.5*(a1+a2));
            psi[i][j][0]=a4*exci;
            xx=xx+dx;
        }
        yy=yy+dy;
    }
    
}

void potential(double v[][91])
{
    /* sets slit is a potential with aperture multiplied by dt */
    double xx, yy, dx, dy, dtx, dt;
    int i, j, nnx, nny;
    dx=0.2;
    dy=dx;
    dt=0.0025;
    nnx=89;
    nny=89;
    /* v[i]*dt */
    yy=-9.;
    for(j=0;j<=nny+1;j++){
        xx=-9.;
        for(i=0;i<=nnx+1;i++){
            if (j==35 && (i<=47 && i>=45 || i<=38 || i>=54 )){
                v[i][j]=170.*dt;
            }
            else {
                v[i][j]=0.;
            }
            xx=xx+dx;
        }
        yy=yy+dy;
    }
}
void solution(double v[][91], double psr[][91][2], double psi[][91][2], double p2[][91], int nit)
{
    /* solves time dependent Schrodinger equation */
    double a1, a2, dx2, dtx, dx, dt, xx, yy, dy, aa, sol;
    int n, i, j, nnx, nny;
    char s[]="run1.0001";
    FILE *out;
    
    dx=0.2;
    dy=dx;
    dt=0.0025;
    dx2=dx*dx;
    dtx=dt/dx2;
    nnx=89;
    nny=89;
    /* next time steps for real and imaginary parts of wave packet
     nit (the input) gives the last time step and the probability
     is plotted (with the potential multiplied by a factor) */
    for(n=1;n<=nit;n++){
        /* compute real part of wave packet and probability */
        for(j=1;j<=nny;j++){
            for(i=1;i<=nnx;i++){
                a2=v[i][j]*psi[i][j][0]+2.0*dtx*psi[i][j][0];
                a1=psi[i+1][j][0]+psi[i-1][j][0]+psi[i][j+1][0]+psi[i][j-1][0];
                psr[i][j][1]=psr[i][j][0]-dtx*a1+2.*a2;
                p2[i][j]=psr[i][j][0]*psr[i][j][1]+psi[i][j][0]*psi[i][j][0];
            }
            /* at x edges derivative is zero */
            psr[0][j][1]=psr[1][j][1];
            psr[nny+1][j][1]=psr[nny][j][1];
        }
        
        /* imaginary part of wave packet is next */
        for(j=0;j<=nny+1;j++){
            for(i=1;i<=nnx;i++){
            a2=v[i][j]*psr[i][j][1]+2.0*dtx*psr[i][j][1];
            a1=psr[i+1][j][1]+psr[i-1][j][1]+psr[i][j+1][1]+psr[i][j-1][1];
            psi[i][j][1]=psi[i][j][0]+dtx*a1-2.*a2 ;
        }
        /* at x edges derivative is zero */
            psi[0][j][1]=psi[1][j][1];
            psi[nny+1][j][1]=psi[nny][j][1];
            }
        /* new iterations are now the old ones, recycle */
            for(j=0;j<=nny+1;j++){
                for(i=0;i<=nnx+1;i++){
                    psi[i][j][0]=psi[i][j][1];
                    psr[i][j][0]=psr[i][j][1]; 
                }
            }
            if(n<10){
                s[8]= n+48;
            }
            if(n<100 && n>9){
                s[7]=(n/10)+48;
                s[8]= (n%10)+48;
            }
            if(n<1000 && n>99){
                s[6]=(n/100)+48;
                s[7]= (((n%100))/10)+48;
                s[8]= (n%10)+48;
            }
            out=fopen(s,"wb");
            for(i=1;i<=nny;i=i+1){
                for(j=1;j<=nny;j=j+1){
                    sol = p2[i][j]+v[i][j]/340.;
                    if(sol<1.2e-10)sol=1.2e-10;
                    fprintf(out, "%e\n", sol);
                }
                fprintf(out, "\n");
            }
            fclose(out);
            
            }
            
            }
