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
void initialize_polymer(int nchain, int natom_perchain, int itype, double bb, double box_size,
                        int start_idx, int mol_start_idx,
                        int *particle_type, int *molecule_id, double **coords, long *idum);
void get_mesh_density(int natom_total, double grid_size, int maxsite_1d, int pm_type,
                      double density_A_ideal, double density_B_ideal,
                      double grid_shift[3], double **coords, int *particle_type, double ****density_grids);

int main(int argc, const char *argv[])
{
    int natom_perchain, natom_total, natom_A, natom_B, nbb, nsc;
    int nchain, nchain_A, nchain_B;
    int nbonds_perchain, nbonds_total;

    int pm_type;
    double kappa, Nbar_sqrt, chiN, Re, ReReRe, bondlen_sqr;
    double grid_size;

    int save_every, tot_mc_cycles, equil_cycles;

    // File Input
    FILE *inputfile;
    inputfile = fopen("Input.txt", "r");

    fscanf(inputfile, "nchain_A = %d\n", &nchain_A);             // # of A polymers
    fscanf(inputfile, "nchain_B = %d\n", &nchain_B);             // # of B polymers
    fscanf(inputfile, "natom_perchain = %d\n", &natom_perchain); // # of atoms per chain
    fscanf(inputfile, "pm_type = %d\n", &pm_type);               // PM type (0- or 1-order)
    fscanf(inputfile, "Nbar_sqrt = %lf\n", &Nbar_sqrt);          // overlap parameter
    fscanf(inputfile, "Re = %lf\n", &Re);                        // End-to-end distance for a polymer
    fscanf(inputfile, "chiN = %lf\n", &chiN);                    // F-H parameter
    fscanf(inputfile, "kappa = %lf\n", &kappa);                  // imcompressible parameter

    fscanf(inputfile, "grid_size = %lf\n", &grid_size);        // mesh size
    fscanf(inputfile, "save_every = %d\n", &save_every);       // saving interval
    fscanf(inputfile, "tot_mc_cycles = %d\n", &tot_mc_cycles); // # of steps (exclude equil. steps)
    fscanf(inputfile, "equil_cycles = %d\n", &equil_cycles);   // # of equil. steps

    // fscanf(inputfile, "Nbb = %d\n", &Nbb);
    // fscanf(inputfile, "Nsc = %d\n", &Nsc);
    fclose(inputfile);

    nchain = nchain_A + nchain_B;          // # of B polymer
    natom_A = nchain_A * natom_perchain;   // total # of A atoms
    natom_B = nchain_B * natom_perchain;   // total # of B atoms
    natom_total = nchain * natom_perchain; // total # of polymers
    nbonds_perchain = natom_perchain - 1;  // # of bonds per chain
    nbonds_total = natom_total - nchain;   // total # of bonds

    double density, density_A_ideal, density_B_ideal;

    double box_size, box_volume;
    int maxsite_1d, nsites;
    ReReRe = Re * Re * Re;
    density = natom_perchain * Nbar_sqrt / ReReRe;
    bondlen_sqr = Re * Re / nbonds_perchain;
    // box_volume = density * natom_total;
    // box_size = cbrt(box_volume);
    // For debug:
    box_size = 10.0;
    box_volume = box_size * box_size * box_size;
    maxsite_1d = round(box_size / grid_size);
    nsites = maxsite_1d * maxsite_1d * maxsite_1d;
    density_A_ideal = natom_A / nsites; // ideal density of A at homogeneous
    density_B_ideal = natom_B / nsites; // ideal density of B at homogeneous

    // Memory allocation
    // type of particle, 1d array (natoms)
    int *particle_type = calloc(natom_total, sizeof(int));
    // molecule ID, 1d array (natoms)
    int *molecule_id = calloc(natom_total, sizeof(int));
    // coordinates of particles, 2d array (natoms, 3)
    double **coords = (double **)calloc(natom_total, sizeof(double *));
    for (int i = 0; i < natom_total; i++)
        coords[i] = (double *)calloc(3, sizeof(double));
    // density profile, 4d array (num_atom_types, num_grid, num_grid, num_grid)
    double ****density_grids = (double ****)calloc(2, sizeof(double ***));
    for (int itype = 0; itype < 2; itype++)
    {
        density_grids[itype] = (double ***)calloc(maxsite_1d, sizeof(double **));
        for (int i = 0; i < maxsite_1d; i++)
        {
            density_grids[itype][i] = (double **)calloc(maxsite_1d, sizeof(double *));
            for (int j = 0; j < maxsite_1d; j++)
            {
                density_grids[itype][i][j] = (double *)calloc(maxsite_1d, sizeof(double));
            }
        }
    }

    // Initialize random number generator
    long *idum = malloc(sizeof(long));
    *idum = initRan();

    // Initialize polymers
    int start_idx = 0;
    int mol_start_idx = 0;
    initialize_polymer(nchain_A, natom_perchain, 0, bondlen_sqr,
                       box_size, start_idx, mol_start_idx, particle_type, molecule_id, coords, idum);
    start_idx = nchain_A * natom_perchain;
    mol_start_idx = nchain_A;
    initialize_polymer(nchain_B, natom_perchain, 1, bondlen_sqr,
                       box_size, start_idx, mol_start_idx, particle_type, molecule_id, coords, idum);

    // Write the initial configuration to a file
    FILE *init_config;
    init_config = fopen("ini.cfg", "w");

    fprintf(init_config, "#Initial configuration\n");
    fprintf(init_config, "%d atoms\n", natom_total);
    fprintf(init_config, "2 atom types\n");
    fprintf(init_config, "%d bonds\n", nbonds_total);
    fprintf(init_config, "1 bond types\n");
    fprintf(init_config, "0.0 %lf xlo xhi\n", box_size);
    fprintf(init_config, "0.0 %lf ylo yhi\n", box_size);
    fprintf(init_config, "0.0 %lf zlo zhi\n", box_size);
    fprintf(init_config, "0 0 0 xy xz yz\n");
    fprintf(init_config, "\n");
    fprintf(init_config, "Atoms # bond\n");
    fprintf(init_config, "\n");

    for (int i = 0; i < natom_total; ++i)
    {
        fprintf(init_config, "%d %d %d %lf %lf %lf\n", i + 1, molecule_id[i], particle_type[i] + 1,
                coords[i][0], coords[i][1], coords[i][2]);
    }
    fclose(init_config);
    // fprintf(outputfile, "SEED %ld\n", *idum); // print out seed, so that it is possible to exactly re-run simulation
    // fclose(outputfile);
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
void initialize_polymer(int nchain, int natom_perchain, int itype, double bb, double box_size,
                        int start_idx, int mol_start_idx,
                        int *particle_type, int *molecule_id, double **coords, long *idum)
{
    double theta, phi;
    int index;
    // double b = sqrt(bb);
    double b = 1;
    for (int i = 0; i < nchain; ++i)
    {
        // Randomizing positions/locations of the first bead(=index 0, 1N, 2N-th...) of each polymer
        coords[i * natom_perchain + start_idx][0] = ran1(idum) * box_size;
        coords[i * natom_perchain + start_idx][1] = ran1(idum) * box_size;
        coords[i * natom_perchain + start_idx][2] = ran1(idum) * box_size;
        particle_type[i * natom_perchain + start_idx] = itype;
        molecule_id[i * natom_perchain + start_idx] = mol_start_idx + i + 1;

        // Randomizing positions of the rest beads within the chain spatially
        for (int j = 1; j < natom_perchain; ++j)
        {
            // from index 1(=2nd bead of each chain)
            index = start_idx + i * natom_perchain + j;
            // Randomly get theta angle(xy plane)
            theta = ran1(idum) * 2.0 * M_PI;
            // randomly get cos() value ranging from -1 to 1 to set 2Pi rotation angle from z axis
            phi = acos(2.0 * ran1(idum) - 1.0);
            // move this much in x direction relative to previous bead
            coords[index][0] = coords[index - 1][0] + b * cos(theta) * sin(phi);
            // move this much in y direction relative to previous bead
            coords[index][1] = coords[index - 1][1] + b * sin(theta) * sin(phi);
            // move this much in z direction relative to previous bead
            coords[index][2] = coords[index - 1][2] + b * cos(phi);
            // assign the atom type
            particle_type[index + start_idx] = itype;
            // assign the molecule ID
            molecule_id[index + start_idx] = mol_start_idx + i + 1;

            // If our new chain goes outside of the periodic box, we wrap around to the other side of the box
            coords[index][0] -= box_size * floor(coords[index][0] / box_size);
            coords[index][1] -= box_size * floor(coords[index][1] / box_size);
            coords[index][2] -= box_size * floor(coords[index][2] / box_size);
            // For debug:
            // printf("%lf\t", coords[index][0]);
            // printf("%lf\t", coords[index][1]);
            // printf("%lf\t", coords[index][2]);
            // printf("\n");

        } // iterate to complete one chain
    }     // done initialization
}

// calculation of density profile of the meshgrids
void get_mesh_density(int natom_total, double grid_size, int maxsite_1d, int pm_type,
                      double density_A_ideal, double density_B_ideal,
                      double grid_shift[3], double **coords, int *particle_type, double ****density_grids)
{
    double rx, ry, rz;
    int xgrid, ygrid, zgrid;
    int atom_type;
    double density_temp[2];

    density_temp[0] = density_A_ideal;
    density_temp[1] = density_B_ideal;
    if (pm_type == 0)
    {
        for (int i = 0; i < natom_total; ++i)
        {
            atom_type = particle_type[i];
            rx = coords[i][0];
            ry = coords[i][1];
            rz = coords[i][2];
            // identify the grid that the particle is in, taking a grid shift into account
            xgrid = rx / grid_size + grid_shift[0];
            xgrid = xgrid % maxsite_1d;
            ygrid = ry / grid_size + grid_shift[1];
            ygrid = ygrid % maxsite_1d;
            zgrid = rz / grid_size + grid_shift[2];
            zgrid = zgrid % maxsite_1d;
            // taking the normalization factor into account
            density_grids[atom_type][xgrid][ygrid][zgrid] += 1 / density_temp[atom_type];
        }
    }
}

void calc_total_particlemesh_energy()
{
}