/* Esempio 5:
 
 M(RT)^2 per un fluido di Lennard Jones
 
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h> /* necessario per la funzione rand() */


const int npmax = 500;
/* const double KB=8.617e-5;  energie in eV */
const double KB=1.;  /* Energie in eps_lj, temperature in eps_lj/KB, lunghezze in sigma_lj */

double delta,rho,sigma_lj,eps_lj,l,li,l2;
double Pi,T,beta,acc,vbox;
int npart;


struct particle
{
    double x,y,z;
};

struct table
{
    int point[200];
};





double total_potential(),pair_pot(),vlj(),vtail(),total_virial();
struct table grnow(),addtable();

struct table addtable(struct table table1, struct table table2)
{
    int i;
    struct table temp;
    for(i=0;i<200;i++)
    {
        temp.point[i]=table1.point[i]+table2.point[i];
    };
    return temp;
}

main()
{
    struct particle walker[npmax];    /* walker */
    struct table gr;
    double etot,vout;
    double accum,dr,l22;
    double dv,vave,vave2,vsum,vsum2,vnow; /* stimatori */
    double dp,pave,pave2,psum,psum2,pnow;
    double vbulk;
    double gr_norm,r,norm,rout,rin;
    float delta_in, sigma_lj_in, eps_lj_in,rho_in,T_in;
    int nstep1,nstep2;
    int i;
    FILE *out;
    out=fopen("gr.out","w+");
    Pi=4.*atan(1.);
    
    printf("Fluido di Lennard-Jones in 3D: \n");
    
    printf("Lennard-Jones sigma: ");
    scanf("%f",&sigma_lj_in);
    sigma_lj=sigma_lj_in;
    printf("Lennard-Jones epsilon: ");
    scanf("%f",&eps_lj_in);
    eps_lj=eps_lj_in;
    printf("Numero di particelle: ");
    scanf("%i",&npart);
    printf("Densita:");
    scanf("%f",&rho_in);
    rho=rho_in;
    printf("Temperatura:");
    scanf("%f",&T_in);
    T=T_in;
    printf("Numero di passi di equilibratura: ");
    scanf("%i",&nstep1);
    printf("Numero di passi: ");
    scanf("%i",&nstep2);
    printf("Lunghezza dello step di Metropolis: ");
    scanf("%f",&delta_in);
    delta=delta_in;
    
    
    l=pow((double)npart/rho,1./3.);
    li=1./l;
    l2=0.5*l;
    
    printf ("Lato della scatola di simulazione %10.8f \n",l);
    beta = 1./(KB*T);
    printf ("beta: %10.5f \n",beta);
    vbox=0.;
    l22=l2*l2-1.e-10;
    vbox=vlj(l22);
    vout=vtail(l);
    vbulk=0.5*vbox*(Pi*npart/6.-1);
    printf ("Correzione di coda: %10.5e\n",vout);
    printf ("Potenziale a L/2: %10.5e %10.5e\n",vbox,vlj(l2*l2));
    printf ("Correzione interna: %10.5e\n",vbulk);
    
    srand(45);
    get_initial_positions(walker);
    vsum=0.;
    vsum2=0.;
    vnow=0.;
    accum=0.;
    for(i=0;i<nstep1;i++)    /* loop principale */
    {
        advance(walker);
        if(i%100==0) printf("Step: %6i    Acceptance:  %10.8f \n",\
                            i+1,accum/(i+1));
        vnow=total_potential(walker);
        accum=accum+acc/(double)npart;
        vsum=vsum+vnow;
        vsum2=vsum2+vnow*vnow;
    };
    printf("\n Fine della fase di equilibratura \n");
    vsum=0.;
    vsum2=0;
    vnow=0.;
    accum=0.;
    psum=0.;
    psum2=0.;
    for(i=0;i<201;i++) gr.point[i]=0;
    for(i=0;i<nstep2;i++)
    {
        advance(walker);
        vnow=total_potential(walker)/(double)npart+(vout+vbulk);
        pnow=rho/beta+total_virial(walker)/(l*l*l);
        gr=addtable(gr,grnow(walker));
        accum=accum+acc/(double)npart;
        vsum=vsum+vnow;
        vsum2=vsum2+vnow*vnow;
        psum=psum+pnow;
        psum2=psum2+pnow*pnow;
        if((i+1)%100==0) printf(\
                                "Step: %6i    Acceptance:  %10.5f  Potential (now: %10.8f ave: %10.8f )\n",\
                                i+1,accum/(i+1), vnow,\
                                vsum/(i+1));
    };
    vave=vsum/(double)nstep2;
    vave2=vsum2/(double)nstep2;
    pave=psum/(double)nstep2;
    pave2=psum2/(double)nstep2;
    dv = sqrt(fabs(vave*vave-vave2)/(double)nstep2);  /* calcola l' errore */
    dp = sqrt(fabs(pave*pave-pave2)/(double)nstep2);
    printf("Energia potenziale media: %10.8f +-  %10.8f\n",\
           vave,dv);
    etot=1.5*KB*T+vave;
    printf("Energia totale media: %10.8f +-  %10.8f\n",etot,dv);
    printf("Pressione media: %10.8f +-  %10.8f\n",\
           pave,dp);
    
    printf("Rapporto di accettazione: %10.8f \n",accum/(double)nstep2);
    for(i=0;i<200;i++)
    {
        dr=l*0.5/201.;
        r=(i+0.5)*dr;
        rin=i*dr;
        rout=(i+1)*dr;
        norm=4./3.*Pi*(rout*rout*rout-rin*rin*rin)*rho;
        gr_norm=(double)gr.point[i]/((double)nstep2*norm*npart);
        fprintf(out,"%10.5f %10.5f\n",r,gr_norm);
    };
}

get_initial_positions(struct particle walker[npmax])


{
    struct particle trial_position;
    
    int i,j;
    double rr,dx,dy,dz;
    int overlap;
    i=0;
    while(i<npart)
    {
        trial_position.x=l*(0.5-(double)rand()/(double)RAND_MAX);
        trial_position.y=l*(0.5-(double)rand()/(double)RAND_MAX);
        trial_position.z=l*(0.5-(double)rand()/(double)RAND_MAX);
        overlap=0;
        for(j=0;j<i;j++)
        {
            dx=trial_position.x-walker[j].x;
            dx=dx-l*rint(dx*li);
            dy=trial_position.y-walker[j].y;
            dy=dy-l*rint(dy*li);
            dz=trial_position.z-walker[j].z;
            dz=dz-l*rint(dz*li);
            rr=dx*dx+dy*dy+dz*dz;
            if(rr<0.5*sigma_lj*sigma_lj) overlap=1;
        };
        if(!overlap)
        {
            walker[i].x=trial_position.x;
            walker[i].y=trial_position.y;
            walker[i].z=trial_position.z;
            i++;
        };
    };
}

advance(struct particle walker[npmax])

{
    struct particle new;
    int i,j;
    double dx,dy,dz,arg;
    double Vold,V,p,csi;
    acc=0.;
    for(j=0;j<npart;j++)
    {
        Vold=pair_pot(j,walker[j].x,walker[j].y,walker[j].z,walker);
        dx=delta*(0.5-(double)rand()/(double)RAND_MAX);
        
        new.x=walker[j].x+dx;
        new.x=new.x-l*rint(new.x*li);
        dy=delta*(0.5-(double)rand()/(double)RAND_MAX);
        new.y=walker[j].y+dy;
        new.y=new.y-l*rint(new.y*li);
        dz=delta*(0.5-(double)rand()/(double)RAND_MAX);
        new.z=walker[j].z+dz;
        new.z=new.z-l*rint(new.z*li);
        V=pair_pot(j,new.x,new.y,new.z,walker);
        arg=V-Vold;
        p=exp(-beta*arg);
        csi=(double)rand()/(double)RAND_MAX;
        
        if(p>csi)
        {
            acc=acc+1.;
            walker[j].x=new.x;
            walker[j].y=new.y;
            walker[j].z=new.z;
        };
        /*      printf("arg: %10.5f p: %10.5f acc: %10.5f \n",arg,p,acc); */
    };
}

double total_potential(struct particle walker[npmax])
{
    int i,j;
    double dx,dy,dz,rr,v;
    v=0.;
    for(i=0;i<npart;i++)
    {
        for(j=0;j<i;j++)
        {
            dx=walker[i].x-walker[j].x;
            dx=dx-l*rint(dx*li);
            dy=walker[i].y-walker[j].y;
            dy=dy-l*rint(dy*li);
            dz=walker[i].z-walker[j].z;
            dz=dz-l*rint(dz*li);
            rr=dx*dx+dy*dy+dz*dz;
            v=v+vlj(rr);
        };
    };
    return v;
}

double pair_pot(int i,double x, double y, double z, struct particle walker[npmax])
{
    int j;
    double dx,dy,dz,rr,v;
    
    v=0.;
    for(j=0;j<i;j++)
	{
        dx=x-walker[j].x;
        dx=dx-l*rint(dx*li);
        dy=y-walker[j].y;
        dy=dy-l*rint(dy*li);
        dz=z-walker[j].z;
        dz=dz-l*rint(dz*li);
        rr=dx*dx+dy*dy+dz*dz;
        v=v+vlj(rr);
	};
    for(j=i+1;j<npart;j++)
	{
        dx=x-walker[j].x;
        dx=dx-l*rint(dx*li);
        dy=y-walker[j].y;
        dy=dy-l*rint(dy*li);
        dz=z-walker[j].z;
        dz=dz-l*rint(dz*li);
        rr=dx*dx+dy*dy+dz*dz;
        v=v+vlj(rr);
	};
    return v;
}

double vlj(double r2)
{
    double r2i,r6i,r12i,v,l22;
    double rcore;
    
    l22=l2*l2;
    if(r2>l22)
    {
        v=0;
        return v;
    };
    
    rcore=0.3*sigma_lj;
    
    if(r2<rcore*rcore) r2=rcore*rcore;
    r2i=sigma_lj*sigma_lj/r2;
    r6i=r2i*r2i*r2i;
    r12i=r6i*r6i;
    v=4.*eps_lj*(r12i-r6i)-vbox;
    return v;
}

double virial(double r2)
{
    double r2i,r6i,r12i,v,l22;
    double rcore;
    
    l22=l2*l2;
    if(r2>l22)
    {
        v=0;
        return v;
    };
    
    rcore=0.3*sigma_lj;
    
    if(r2<rcore*rcore) r2=rcore*rcore;
    r2i=sigma_lj*sigma_lj/r2;
    r6i=r2i*r2i*r2i;
    r12i=r6i*r6i;
    v=4.*eps_lj*(12.*r12i-6.*r6i);
    return v;
}

struct table grnow(struct particle walker[npmax])
{
    int i,j;
    double dx,dy,dz,rr;
    double r,dh;
    int idx;
    struct table gg;
    dh=l2/201.;
    
    for(i=0;i<201;i++) gg.point[i]=0;
    for(i=0;i<npart;i++)
    {
        for(j=0;j<i;j++)
        {
            dx=walker[i].x-walker[j].x;
            dx=dx-l*rint(dx*li);
            dy=walker[i].y-walker[j].y;
            dy=dy-l*rint(dy*li);
            dz=walker[i].z-walker[j].z;
            dz=dz-l*rint(dz*li);
            rr=dx*dx+dy*dy+dz*dz;
            r=sqrt(rr);
            if(r<l2)
            {
                idx=(int)(r/dh);
                gg.point[idx]=gg.point[idx]+2;
            };
        };
    };
    return gg;
}

double vtail(double l)
{
    double aux,aux3,aux9,v;
    
    aux=2.*sigma_lj*li;
    aux3=aux*aux*aux;
    aux9=aux3*aux3*aux3;
    v=8./9.*Pi*eps_lj*rho*(aux9-3.*aux3);
    return v;
}

double total_virial(struct particle walker[npmax])
{
    int i,j;
    double dx,dy,dz,rr,v;
    v=0.;
    for(i=0;i<npart;i++)
    {
        for(j=0;j<i;j++)
        {
            dx=walker[i].x-walker[j].x;
            dx=dx-l*rint(dx*li);
            dy=walker[i].y-walker[j].y;
            dy=dy-l*rint(dy*li);
            dz=walker[i].z-walker[j].z;
            dz=dz-l*rint(dz*li);
            rr=dx*dx+dy*dy+dz*dz;
            v=v+virial(rr)/3.;
        };
    };
    return v;
}