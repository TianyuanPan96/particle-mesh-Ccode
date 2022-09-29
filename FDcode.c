#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define NDIV (1+IMM1/NTAB)

void cholesky (double *mm, int Chainlength, double *cc, int t);
float ran1(long *idum);
float gasdev(long *idum);
long initRan();
float stochastic();

int main (int argc, const char * argv[])
{

//VARIABLES!
	int N, Nbb, Nsc,tmax, k, t, i, j, n,nb, PMFcount;
	double epsilon, kappa, dt, D,L;
    double dx,dy,dz,Fs,rr,r, p, r6, ratio, coeff,c,c2,c3,PMF_force,PMF_energy,PMFave_force,PMFave_energy;
		
//File Input
	FILE *inputfile;
	inputfile = fopen("Input.txt", "r");
	fscanf(inputfile, "epsilon = %lf\n", &epsilon);
	fscanf(inputfile, "kappa = %lf\n", &kappa);
	fscanf(inputfile, "dt = %lf\n", &dt);
	fscanf(inputfile, "tmax = %d\n", &tmax);
	fscanf(inputfile, "Nbb = %d\n", &Nbb);
	fscanf(inputfile, "Nsc = %d\n", &Nsc);
	fscanf(inputfile, "D = %lf\n", &D);
	fclose(inputfile);
	
	N = Nbb*(Nsc+1);
	L = 2.0*(double)Nbb;
	PMFcount = 0;
	PMFave_force = 0.0;
	PMFave_energy = 0.0;
	
	time_t theTime;
	time(&theTime);
	
//open Files
	char* str = malloc(sizeof(char)*30);
	char* str1 = malloc(sizeof(char)*30);
	sprintf(str, "RFD%d_%d_edot%d.xyz", (int)(D*10), Nbb, Nsc);
	sprintf(str1, "PFD%d_%d_edot%d.txt", (int)(D*10), Nbb, Nsc);
	
	FILE *outputfile, *Output;
	outputfile = fopen(str1, "w");
	
	//fprintf(outputfile, ctime(&theTime));
	fprintf(outputfile, "epsilon = %lf\n", epsilon);
	fprintf(outputfile, "kappa = %lf\n", kappa);
	fprintf(outputfile, "dt = %lf\n", dt);
	fprintf(outputfile, "tmax = %d\n", tmax);
	fprintf(outputfile, "Nbb = %d\n", Nbb);
	fprintf(outputfile, "Nsc = %d\n", Nsc);
	fprintf(outputfile, "D = %lf\n", D);
	fprintf(outputfile, "\n");
	fprintf(outputfile, "\n");
	fclose(outputfile);
	
	double *Chainrxa = calloc(N, sizeof(double));
	double *Chainrya = calloc(N, sizeof(double));
	double *Chainrza = calloc(N, sizeof(double));
	double *Chainfxa = calloc(N, sizeof(double));
	double *Chainfya = calloc(N, sizeof(double));
	double *Chainfza = calloc(N, sizeof(double));
	
	double *Chainrxb = calloc(N, sizeof(double));
	double *Chainryb = calloc(N, sizeof(double));
	double *Chainrzb = calloc(N, sizeof(double));
	double *Chainfxb = calloc(N, sizeof(double));
	double *Chainfyb = calloc(N, sizeof(double));
	double *Chainfzb = calloc(N, sizeof(double));
	
    
    //Initialize Random Matrix
	double *RRa = calloc(N*3, sizeof(double));
	double *RRb = calloc(N*3, sizeof(double));
    
    p = sqrt(2.0*dt);
    //Initialize Random Variable
	long *idum = malloc(sizeof(long));
    
    *idum = initRan();
    
    outputfile = fopen(str1, "a");
    fprintf(outputfile, "SEED %ld\n", *idum);
    fclose(outputfile);
    
    
//Initialize Polymer Chain as
	
	for(i = 0; i<Nbb; ++i)
	{
		Chainrxa[i] = 0.0;
		Chainrya[i] = 0.0;
		Chainrza[i] = 1.0+2.0*(double)i;
	}
	for(i = 0; i<(N-Nbb); ++i)
	{
		Chainrxa[i+Nbb] = -2.0-2.0*(double)(i%Nsc);
		Chainrya[i+Nbb] = 0.0;
		Chainrza[i+Nbb] = 1.0+2.0*(double)(i/Nsc);
	}
	
	//Initialize Polymer Chain b
	
	for(i = 0; i<Nbb; ++i)
	{
		Chainrxb[i] = D;
		Chainryb[i] = 0.0;
		Chainrzb[i] = 1.0+2.0*(double)i;
	}
	for(i = 0; i<(N-Nbb); ++i)
	{
		Chainrxb[i+Nbb] = D+2.0+2.0*(double)(i%Nsc);
		Chainryb[i+Nbb] = 0.0;
		Chainrzb[i+Nbb] = 1.0+2.0*(double)(i/Nsc);
	}
	
    
    
    ////Simulation loop
    
    for(t = 0; t<tmax; ++t)
    {
        for(i = 0; i<N; ++i)
        {
            //zero the forces at each time step
            Chainfxa[i] = 0.0; Chainfya[i] = 0.0; Chainfza[i] = 0.0;
			Chainfxb[i] = 0.0; Chainfyb[i] = 0.0; Chainfzb[i] = 0.0;
        }
        
#include "Interactions.h"
        
        //Set random velocities
        for(i = 0; i<N*3; ++i)
        {
            RRa[i] = gasdev(idum);
			RRb[i] = gasdev(idum);
        }
        //Update positions based on forces, flows, random displacements
        for(i = Nbb; i<N; ++i)
        {
            Chainrxa[i] += dt*Chainfxa[i]+p*RRa[3*i];
            Chainrya[i] += dt*Chainfya[i]+p*RRa[3*i+1];
            Chainrza[i] += dt*Chainfza[i]+p*RRa[3*i+2];
			if(Chainrza[i] > L) Chainrza[i]-=L;
			if(Chainrza[i] < 0.0) Chainrza[i]+=L;
			
			Chainrxb[i] += dt*Chainfxb[i]+p*RRb[3*i];
			Chainryb[i] += dt*Chainfyb[i]+p*RRb[3*i+1];
			Chainrzb[i] += dt*Chainfzb[i]+p*RRb[3*i+2];
			if(Chainrzb[i] > L) Chainrzb[i]-=L;
			if(Chainrzb[i] < 0.0) Chainrzb[i]+=L;
        }
        //Every so often, I output an xyz file
        if(t%10000==0 && t>0)
        {
			
			
			if(t>100000)
			{
				PMFave_force += PMF_force/10000.0;
				PMFave_energy += PMF_energy/10000.0;
				PMFcount++;
				printf("%d %lf %lf %lf\n", t, D, PMFave_force/(double)PMFcount, PMFave_energy/(double)PMFcount);
			}
			PMF_force = 0.0;
			PMF_energy = 0.0;
			
            Output = fopen(str, "a");
            fprintf(Output, "%d\n%d %lf\n", 2*N, t, D);
            for(i = 0; i<N; ++i)
            {
                fprintf(Output, "A %lf %lf %lf\n", Chainrxa[i], Chainrya[i], Chainrza[i]);
            }
			for(i = 0; i<N; ++i)
			{
				fprintf(Output, "B %lf %lf %lf\n", Chainrxb[i], Chainryb[i], Chainrzb[i]);
			}
            fclose(Output);
        }
        
    }
	

	
    return 0;
}



float ran1(long *idum)
{
	int j;
	long k;
	static long idum2 = 123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;
	
	if(*idum <= 0)
	{
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2 = (*idum);
		for(j=NTAB+7;j>=0;--j)
		{
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if(*idum<0) *idum+=IM1;
			if(j<NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if(*idum<0) *idum += IM1;
	k=idum2/IQ2;
	if(*idum<0) idum2+= IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if(iy<1) iy += IMM1;
	if((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

float gasdev(long *idum)
{
	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;
	
	if (*idum < 0) iset = 0;
	if (iset == 0)
	{
		do
		{
			v1 = 2.0*ran1(idum)-1.0;
			v2 = 2.0*ran1(idum)-1.0;
			rsq = v1*v1+v2*v2;
		}
		while(rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	}
	else
	{
		iset = 0;
		return gset;
	}
}


void cholesky(double *mm, int Chainlength, double *cc, int t)
{
	int n = Chainlength*3;
	int i,j,k,error;
	double sum;
    error = 0;
	for(i = 0; i<n; ++i)
	{
		for(j = 0; j<n; ++j)
		{
			cc[i+n*j] = 0;
		}
	}
	for(i=0; i<n;++i)
	{
		for(j=i;j<n;++j)
		{
			sum = mm[i+n*j];
			for(k=i-1; k>=0;--k)
			{
				sum -= cc[k+n*i]*cc[k+n*j];
			}
			if(i==j)
			{
				if(sum<=0)
				{
	//				printf("%d, %d, %d, %f\n", t, i,j,sum);
                    error = 1;
				}
				cc[i+n*i]=sqrt(sum);
			}
			else
			{
				if(cc[i+n*i] == 0)
				{
                    error = 1;
				}
				else
				{
					cc[i+n*j] = sum/cc[i+n*i];
				}
			}
		}
	}
    if(error == 1)
    {
        printf("%d\n", t);
        for(i = 0; i<n; ++i)
        {
            for(j = 0; j<n; ++j)
            {
                cc[i+n*j] = 0.0;
                if(i==j) cc[i+n*j] = 1.0;
            }
        }
    }
}
long initRan() {
    time_t seconds;
    time(&seconds); //switched from time due to crashes... unsure why...
    return -1*(unsigned long)(seconds); //Dividing by 100 keeps within correct bounds? not sure why this works
}









