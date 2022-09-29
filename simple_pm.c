#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0 / IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)
#define NDIV (1 + IMM1 / NTAB)

float ran1(long *idum);
float gasdev(long *idum);
long initRan();
void initialize_polymer(int n, int N, int n_A, double bb, double L, double Lx,
                        double *rx, double *ry, double *rz, int *Particle_type, long *idum);

int main(int argc, const char *argv[])
{
    int natom_perchain, natom_total, natom_A, natom_B, nbb, nsc;
    int nchain, nchain_A, nchain_B;
    int nbonds_perchain, nbonds_total;

    int pm_type;
    double kappa, Nbar_sqrt, chiN, Re, ReReRe, bondlen_sqr;
    int maxsite_x, maxsite_y, maxsite_z, nsites;
    double grid_size;

    double box_x, box_y, box_z;
    double box_volume;
    double TimestepSize, SavingInterval, TotalTimestep, EquilibrationStep;

    // File Input
    FILE *inputfile;
    inputfile = fopen("Input.txt", "r");

    fscanf(inputfile, "nchain_A = %d\n", &nchain_A);
    fscanf(inputfile, "nchain_B = %d\n", &nchain_B);
    fscanf(inputfile, "natom_perchain = %d\n", &natom_perchain);
    fscanf(inputfile, "pm_type = %d\n", &pm_type);      // PM type (0- or 1-order)
    fscanf(inputfile, "Nbar_sqrt = %lf\n", &Nbar_sqrt); // overlap parameter
    fscanf(inputfile, "chiN = %lf\n", &chiN);           // F-H parameter
    fscanf(inputfile, "kappa = %lf\n", &kappa);         // imcompressible parameter

    fscanf(inputfile, "box_size = %lf\n", &box_x);
    fscanf(inputfile, "grid_size = %lf\n", &grid_size);
    fscanf(inputfile, "TimestepSize = %lf\n", &TimestepSize);          // size of time steps
    fscanf(inputfile, "SavingInterval = %d\n", &SavingInterval);       // saving interval
    fscanf(inputfile, "TotalTimestep = %d\n", &TotalTimestep);         // # of steps
    fscanf(inputfile, "EquilibrationStep = %d\n", &EquilibrationStep); // equil. steps

    // fscanf(inputfile, "Nbb = %d\n", &Nbb);
    // fscanf(inputfile, "Nsc = %d\n", &Nsc);
    fclose(inputfile);

    nchain = nchain_A + nchain_B;          // number of B polymer
    natom_A = nchain_A * natom_perchain;   // total number of A atoms
    natom_B = nchain_B * natom_perchain;   // total number of B atoms
    natom_total = nchain * natom_perchain; // total # of molecules(polymers)
    nbonds_perchain = natom_perchain - 1;
    nbonds_total = natom_total - nchain; // total number of bonds

    box_y = box_z = box_x;
    box_volume = box_x * box_y * box_z;
    maxsite_x = maxsite_y = maxsite_z = round(box_x / grid_size);

    double density, density_A, density_B;
    density_A = natom_A / nsites;    // ideal density of A at homogeneous
    density_B = natom_B / nsites;    // ideal density of B at homogeneous
    density = density_A + density_B; // total density of A and B
    ReReRe = natom_perchain * Nbar_sqrt / density;
    Re = cbrt(ReReRe);
    bondlen_sqr = Re * Re / nbonds_perchain;
}

float ran1(long *idum)
{
    // random number between 0 and 1
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0)
    {
        if (-(*idum) < 1)
            *idum = 1;
        else
            *idum = -(*idum);
        idum2 = (*idum);
        for (j = NTAB + 7; j >= 0; --j)
        {
            k = (*idum) / IQ1;
            *idum = IA1 * (*idum - k * IQ1) - k * IR1;
            if (*idum < 0)
                *idum += IM1;
            if (j < NTAB)
                iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ1;
    *idum = IA1 * (*idum - k * IQ1) - k * IR1;
    if (*idum < 0)
        *idum += IM1;
    k = idum2 / IQ2;
    if (*idum < 0)
        idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if (iy < 1)
        iy += IMM1;
    if ((temp = AM * iy) > RNMX)
        return RNMX;
    else
        return temp;
}

float gasdev(long *idum)
{
    // random number with normal distribution N(0,1)
    float ran1(long *idum);
    static int iset = 0;
    static float gset;
    float fac, rsq, v1, v2;

    if (*idum < 0)
        iset = 0;
    if (iset == 0)
    {
        do
        {
            v1 = 2.0 * ran1(idum) - 1.0;
            v2 = 2.0 * ran1(idum) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;
        return v2 * fac;
    }
    else
    {
        iset = 0;
        return gset;
    }
}

long initRan()
{
    time_t seconds;
    time(&seconds);
    return -1 * (unsigned long)(seconds);
}

// initialization of polymer
void initialize_polymer(int n, int N, int n_A, double bb, double L, double Lx,
                        double *rx, double *ry, double *rz, int *Particle_type, long *idum)
{
    double theta, phi;
    int index;
    double b = sqrt(sqrt(bb) / 3.0);
    for (int i = 0; i < n; ++i)
    {
        // Randomizing positions/locations of the first bead(=index 0, 1N, 2N-th...) of each polymer
        rx[i * N] = ran1(idum) * Lx;
        ry[i * N] = ran1(idum) * L;
        rz[i * N] = ran1(idum) * L;
        if (i < n_A)
        {
            Particle_type[i * N] = 0;
        }
        else
        {
            Particle_type[i * N] = 1;
        } // set type A=0, B=1 to distinguish

        // Randomizing positions of the rest beads within the chain spatially
        for (int j = 1; j < N; ++j)
        {
            index = i * N + j;                                     // from index 1(=2nd bead of each chain)
            theta = ran1(idum) * 2.0 * 3.141592;                   // Randomly get theta angle(xy plane)
            phi = acos(2.0 * ran1(idum) - 1.0);                    // randomly get cos() value ranging from -1 to 1 to set 2Pi rotation angle from z axis
            rx[index] = rx[index - 1] + b * cos(theta) * sin(phi); // move this much in x direction relative to previous bead
            ry[index] = ry[index - 1] + b * sin(theta) * sin(phi); // move this much in y direction relative to previous bead
            rz[index] = rz[index - 1] + b * cos(phi);              // move this much in z direction relative to previous bead
            Particle_type[index] = Particle_type[index - 1];       // assign the same type for homopolymer

            // If our new chain goes outside of the periodic box, we wrap around to the other side of the box
            rx[index] -= Lx * floor(rx[index] / Lx);
            ry[index] -= L * floor(ry[index] / L);
            rz[index] -= L * floor(rz[index] / L);

        } // iterate to complete one chain
    }     // done initialization
}