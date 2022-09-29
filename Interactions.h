//
//  Input.h
//  
//
//  Created by Charles Sing on 1/24/14.
//
//
//Spring force
for(i = 0; i<(N-Nbb); ++i)
{
	if(i%Nsc > 0)
	{
		dx = Chainrxa[Nbb+i] - Chainrxa[Nbb+i-1];
		dy = Chainrya[Nbb+i] - Chainrya[Nbb+i-1];
		dz = Chainrza[Nbb+i] - Chainrza[Nbb+i-1];
		if(dz>L/2) dz -= L; if(dz<-L/2) dz += L;
		rr = dx*dx+dy*dy+dz*dz;
		r = sqrt(rr);
		Fs = -kappa*(r-2.0);
		Chainfxa[Nbb+i] += Fs*dx/r;
		Chainfya[Nbb+i] += Fs*dy/r;
		Chainfza[Nbb+i] += Fs*dz/r;
		Chainfxa[Nbb+i-1] += -Fs*dx/r;
		Chainfya[Nbb+i-1] += -Fs*dy/r;
		Chainfza[Nbb+i-1] += -Fs*dz/r;
	}
	else if(i%Nsc==0)
	{
		nb = i/Nsc;
		dx = Chainrxa[Nbb+i] - Chainrxa[nb];
		dy = Chainrya[Nbb+i] - Chainrya[nb];
		dz = Chainrza[Nbb+i] - Chainrza[nb];
		if(dz>L/2) dz -= L; if(dz<-L/2) dz += L;
		rr = dx*dx+dy*dy+dz*dz;
		r = sqrt(rr);
		Fs = -kappa*(r-2.0);
		Chainfxa[Nbb+i] += Fs*dx/r;
		Chainfya[Nbb+i] += Fs*dy/r;
		Chainfza[Nbb+i] += Fs*dz/r;
	}
    
}

for(i = 0; i<(N-Nbb); ++i)
{
	if(i%Nsc > 0)
	{
		dx = Chainrxb[Nbb+i] - Chainrxb[Nbb+i-1];
		dy = Chainryb[Nbb+i] - Chainryb[Nbb+i-1];
		dz = Chainrzb[Nbb+i] - Chainrzb[Nbb+i-1];
		if(dz>L/2) dz -= L; if(dz<-L/2) dz += L;
		rr = dx*dx+dy*dy+dz*dz;
		r = sqrt(rr);
		Fs = -kappa*(r-2.0);
		Chainfxb[Nbb+i] += Fs*dx/r;
		Chainfyb[Nbb+i] += Fs*dy/r;
		Chainfzb[Nbb+i] += Fs*dz/r;
		Chainfxb[Nbb+i-1] += -Fs*dx/r;
		Chainfyb[Nbb+i-1] += -Fs*dy/r;
		Chainfzb[Nbb+i-1] += -Fs*dz/r;
	}
	else if(i%Nsc==0)
	{
		nb = i/Nsc;
		dx = Chainrxb[Nbb+i] - Chainrxb[nb];
		dy = Chainryb[Nbb+i] - Chainryb[nb];
		dz = Chainrzb[Nbb+i] - Chainrzb[nb];
		if(dz>L/2) dz -= L; if(dz<-L/2) dz += L;
		rr = dx*dx+dy*dy+dz*dz;
		r = sqrt(rr);
		Fs = -kappa*(r-2.0);
		Chainfxb[Nbb+i] += Fs*dx/r;
		Chainfyb[Nbb+i] += Fs*dy/r;
		Chainfzb[Nbb+i] += Fs*dz/r;
	}
	
}

//Lennard-Jones interaction
for(i = 0; i<N; ++i)
{
    for(j = i+1; j<N; ++j)
    {
        dx = Chainrxa[i] - Chainrxa[j];
        dy = Chainrya[i] - Chainrya[j];
        dz = Chainrza[i] - Chainrza[j];
		if(dz>L/2) dz -= L; if(dz<-L/2) dz += L;
        rr = dx*dx+dy*dy+dz*dz;
        if(rr<5.04)
        {
            ratio = 4.00/rr;
            r6 = ratio*ratio*ratio;
            coeff = (epsilon/rr)*(48*r6*r6-24*r6);
            Chainfxa[i] += coeff*dx;
            Chainfya[i] += coeff*dy;
            Chainfza[i] += coeff*dz;
            Chainfxa[j] -= coeff*dx;
            Chainfya[j] -= coeff*dy;
            Chainfza[j] -= coeff*dz;
        }
        
    }
}

for(i = 0; i<N; ++i)
{
	for(j = i+1; j<N; ++j)
	{
		dx = Chainrxb[i] - Chainrxb[j];
		dy = Chainryb[i] - Chainryb[j];
		dz = Chainrzb[i] - Chainrzb[j];
		if(dz>L/2) dz -= L; if(dz<-L/2) dz += L;
		rr = dx*dx+dy*dy+dz*dz;
		if(rr<5.04)
		{
			ratio = 4.00/rr;
			r6 = ratio*ratio*ratio;
			coeff = (epsilon/rr)*(48*r6*r6-24*r6);
			Chainfxb[i] += coeff*dx;
			Chainfyb[i] += coeff*dy;
			Chainfzb[i] += coeff*dz;
			Chainfxb[j] -= coeff*dx;
			Chainfyb[j] -= coeff*dy;
			Chainfzb[j] -= coeff*dz;
		}
		
	}
}

for(i = 0; i<N; ++i)
{
	for(j = 0; j<N; ++j)
	{
		dx = Chainrxa[i] - Chainrxb[j];
		dy = Chainrya[i] - Chainryb[j];
		dz = Chainrza[i] - Chainrzb[j];
		if(dz>L/2) dz -= L; if(dz<-L/2) dz += L;
		rr = dx*dx+dy*dy+dz*dz;
		if(rr<5.04)
		{
			ratio = 4.00/rr;
			r6 = ratio*ratio*ratio;
			// coeff = (12*epsilon/rr)*(r6*r6-r6);
			coeff = (epsilon/rr)*(48*r6*r6-24*r6);
			Chainfxa[i] += coeff*dx;
			Chainfya[i] += coeff*dy;
			Chainfza[i] += coeff*dz;
			Chainfxb[j] -= coeff*dx;
			Chainfyb[j] -= coeff*dy;
			Chainfzb[j] -= coeff*dz;
			PMF_force -= coeff*dx;
			PMF_energy += (4*epsilon)*(r6*r6-r6)+epsilon;
		}
		
	}
}
