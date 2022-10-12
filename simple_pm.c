#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
double periodic(double ri, double bl);
double dist_calc(double ra_x, double rb_x, double ra_y, double rb_y, double ra_z, double rb_z, double bl);

void initialize_polymer(int nchain, int natom_perchain, int itype, double bb, double box_size,
                        int start_idx, int mol_start_idx, int *particle_type, int *molecule_id,
                        double **coords, int **bonds, int **atom_bond_list, long *idum);
int metro_crit(double enrg_diff, long *idum);
void get_mesh_density(int natom_total, double grid_size, int maxsite_1d, int pm_type,
                      double density_A_ideal, double density_B_ideal,
                      double grid_shift[3], double **coords, int *particle_type, double ****density_grids);
double calc_total_particlemesh_energy(double ****density_grids, double chi, double kappa,
                                      double density, int maxsite_1d);
double calc_singlemesh_energy(double phiA, double phiB, double chi, double kappa, double density);
double calc_total_bond_energy_harmonic(double **coords, int **bonds, int nbonds_total,
                                       double kspring, double box_size);
double calc_atom_bond_energy_harmonic(double **coords, int **atom_bond_list, int iatom,
                                      double kspring, double box_size);
void mc_atom_displ(int iatom, double displ, int *naccept, double ****density_grids, double density_A_ideal, double density_B_ideal,
                   int *particle_type, double **coords, int **atom_bond_list, double box_size, double grid_size,
                   int maxsite_1d, double chi, double kappa, double density, double kspring, int natom_total,
                   double grid_shift[3], long *idum);

int main(int argc, const char *argv[])
{
    int natom_perchain, natom_total, natom_A, natom_B, nbb, nsc;
    int nchain, nchain_A, nchain_B;
    int nbonds_perchain, nbonds_total;

    int pm_type;
    double kappa, kspring, chi, bondlen, bondlen_sqr;
    double grid_size;

    int save_every, tot_mc_cycles, equil_cycles, grid_shift_every;

    // File Input
    FILE *inputfile;
    inputfile = fopen("Input.txt", "r");

    fscanf(inputfile, "nchain_A = %d\n", &nchain_A);             // # of A polymers
    fscanf(inputfile, "nchain_B = %d\n", &nchain_B);             // # of B polymers
    fscanf(inputfile, "natom_perchain = %d\n", &natom_perchain); // # of atoms per chain
    fscanf(inputfile, "pm_type = %d\n", &pm_type);               // PM type (0- or 1-order)
    fscanf(inputfile, "bondlen = %lf\n", &bondlen);              // bond length           // End-to-end distance for a polymer
    fscanf(inputfile, "chi = %lf\n", &chi);                      // F-H parameter
    fscanf(inputfile, "kappa = %lf\n", &kappa);                  // imcompressible parameter
    fscanf(inputfile, "kspring = %lf\n", &kspring);              // spring constant

    fscanf(inputfile, "grid_size = %lf\n", &grid_size);              // mesh size
    fscanf(inputfile, "save_every = %d\n", &save_every);             // saving interval
    fscanf(inputfile, "tot_mc_cycles = %d\n", &tot_mc_cycles);       // # of steps (exclude equil. steps)
    fscanf(inputfile, "equil_cycles = %d\n", &equil_cycles);         // # of equil. steps
    fscanf(inputfile, "grid_shift_every = %d\n", &grid_shift_every); // frequency of shift grid

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
    bondlen_sqr = bondlen * bondlen;
    // Re = sqrt(nbonds_perchain * bondlen_sqr);
    // ReReRe = Re * Re * Re;
    // density = natom_perchain * Nbar_sqrt / ReReRe;
    density = 1.0;

    // spring constant is 1.5 / b^2
    kspring = kspring / bondlen_sqr;
    box_volume = density * natom_total;
    box_size = cbrt(box_volume);
    maxsite_1d = round(box_size / grid_size);
    nsites = maxsite_1d * maxsite_1d * maxsite_1d;
    density_A_ideal = (double)natom_A / nsites; // ideal density of A at homogeneous
    density_B_ideal = (double)natom_B / nsites; // ideal density of B at homogeneous

    FILE *outputfile;
    outputfile = fopen("output.txt", "w");
    fprintf(outputfile, "nchain_A = %d\n", nchain_A);
    fprintf(outputfile, "nchain_B = %d\n", nchain_B);
    fprintf(outputfile, "natom_perchain = %d\n", natom_perchain);
    fprintf(outputfile, "natom_total = %d\n", natom_total);
    fprintf(outputfile, "nbonds_total = %d\n", nbonds_total);
    fprintf(outputfile, "pm_type = %d\n", pm_type);
    fprintf(outputfile, "density = %lf\n", density);
    fprintf(outputfile, "bondlen = %lf\n", bondlen);
    fprintf(outputfile, "chi = %lf\n", chi);
    fprintf(outputfile, "kappa = %lf\n", kappa);
    fprintf(outputfile, "kspring = %lf\n", kspring);
    fprintf(outputfile, "box_size = %lf\n", box_size);
    fprintf(outputfile, "grid_size = %lf\n", grid_size);
    fprintf(outputfile, "maxsite_1d = %d\n", maxsite_1d);
    fprintf(outputfile, "nsites = %d\n", nsites);
    fprintf(outputfile, "density_A_ideal = %lf\n", density_A_ideal);
    fprintf(outputfile, "density_B_ideal = %lf\n", density_B_ideal);

    fclose(outputfile);

    // Memory allocation
    // type of particle, 1d array (natoms)
    int *particle_type = calloc(natom_total, sizeof(int));
    // molecule ID, 1d array (natoms)
    int *molecule_id = calloc(natom_total, sizeof(int));
    // coordinates of particles, 2d array (natoms, 3)
    double **coords;
    coords = (double **)calloc(natom_total, sizeof(double *));
    for (int i = 0; i < natom_total; ++i)
    {
        coords[i] = (double *)calloc(3, sizeof(double));
    }
    // coords[0][0] = 0;
    // bond table, 2d array (nbonds, 3)
    int **bonds;
    bonds = (int **)calloc(nbonds_total, sizeof(int *));
    for (int i = 0; i < nbonds_total; ++i)
    {
        bonds[i] = (int *)calloc(2, sizeof(int));
    }
    // atom bond list, each entry lists the atoms that the given atom is bonded with
    int **atom_bond_list;
    atom_bond_list = (int **)calloc(natom_total, sizeof(int *));
    for (int i = 0; i < natom_total; ++i)
    {
        // Assuming maximum 4 bonds per atom!!!
        atom_bond_list[i] = (int *)calloc(4, sizeof(int));
    }
    // Then initialize this list to -1 to make sure check can be performed
    for (int i = 0; i < natom_total; ++i)
    {
        // Assuming maximum 4 bonds per atom!!!
        for (int j = 0; j < 4; ++j)
        {
            atom_bond_list[i][j] = -1;
        }
    }
    // density profile, 4d array (num_atom_types, num_grid, num_grid, num_grid)
    double ****density_grids = calloc(2, sizeof(double ***));
    for (int itype = 0; itype < 2; ++itype)
    {
        density_grids[itype] = calloc(maxsite_1d, sizeof(double **));
        for (int i = 0; i < maxsite_1d; ++i)
        {
            density_grids[itype][i] = calloc(maxsite_1d, sizeof(double *));
            for (int j = 0; j < maxsite_1d; ++j)
            {
                density_grids[itype][i][j] = calloc(maxsite_1d, sizeof(double));
            }
        }
    }

    // Initialize random number generator
    long *idum = malloc(sizeof(long));
    *idum = initRan();

    // Initialize polymers
    int start_idx = 0;
    int mol_start_idx = 0;
    initialize_polymer(nchain_A, natom_perchain, 0, bondlen_sqr, box_size, start_idx,
                       mol_start_idx, particle_type, molecule_id, coords, bonds, atom_bond_list, idum);
    start_idx = natom_A;
    mol_start_idx = nchain_A;
    initialize_polymer(nchain_B, natom_perchain, 1, bondlen_sqr, box_size, start_idx,
                       mol_start_idx, particle_type, molecule_id, coords, bonds, atom_bond_list, idum);

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
        fprintf(init_config, "%d %d %d %f %f %f\n", i + 1, molecule_id[i], particle_type[i] + 1,
                coords[i][0], coords[i][1], coords[i][2]);
    }

    fprintf(init_config, "\n");
    fprintf(init_config, "Bonds\n");
    fprintf(init_config, "\n");

    for (int i = 0; i < nbonds_total; ++i)
    {
        fprintf(init_config, "%d 1 %d %d\n", i + 1, bonds[i][0] + 1, bonds[i][1] + 1);
    }
    fclose(init_config);

    // Get the initial density profile, with grid shift = 0
    double grid_shift[3];
    grid_shift[0] = 0;
    grid_shift[1] = 0;
    grid_shift[2] = 0;
    get_mesh_density(natom_total, grid_size, maxsite_1d, pm_type, density_A_ideal, density_B_ideal,
                     grid_shift, coords, particle_type, density_grids);

    // equilibration run
    // Write the thermo info to a file
    double mesh_energy, bond_energy, total_energy;
    FILE *thermo;
    char *thermo_filename = malloc(sizeof(char) * 30);
    sprintf(thermo_filename, "thermo_N%d_A%d_B%d.txt", natom_perchain, nchain_A, nchain_B);
    thermo = fopen(thermo_filename, "w");
    fprintf(thermo, "t\tmesh_energy\tbond_energy\ttotal_energy\tacceptance_rate\n");
    fclose(thermo);

    // Write the density info to a file
    FILE *density_profile;
    char *density_filename = malloc(sizeof(char) * 30);
    sprintf(density_filename, "dens_profile_N%d_A%d_B%d.txt", natom_perchain, nchain_A, nchain_B);
    density_profile = fopen(density_filename, "w");
    fprintf(density_profile, "type\tx\ty\tz\tdensity\n");
    fclose(density_profile);

    thermo = fopen(thermo_filename, "a");
    density_profile = fopen(density_filename, "a");
    int t;
    int naccept = 0;
    int nattempt = 0;
    double acc_rate = 0.0;
    for (t = 0; t <= equil_cycles; ++t)
    {
        // use a new grid_shift every so often
        if (t % grid_shift_every == 0)
        {
            for (int itype = 0; itype < 2; ++itype)
            {
                for (int i = 0; i < maxsite_1d; ++i)
                {
                    for (int j = 0; j < maxsite_1d; ++j)
                    {
                        memset(density_grids[itype][i][j], 0, maxsite_1d * sizeof(double));
                    }
                }
            }
            // move by atmost half the grid size is enough to mimic a random grid discretization
            grid_shift[0] = 0.5 * ran1(idum) * grid_size;
            grid_shift[1] = 0.5 * ran1(idum) * grid_size;
            grid_shift[2] = 0.5 * ran1(idum) * grid_size;
            get_mesh_density(natom_total, grid_size, maxsite_1d, pm_type, density_A_ideal, density_B_ideal,
                             grid_shift, coords, particle_type, density_grids);

            // For debug:
            // for (int itype = 0; itype < 2; ++itype)
            // {
            //     for (int i = 0; i < maxsite_1d; ++i)
            //     {
            //         for (int j = 0; j < maxsite_1d; ++j)
            //         {
            //             for (int k = 0; k < maxsite_1d; ++k)
            //             {
            //                 printf("%f\n", density_grids[itype][i][j][k]);
            //             }
            //         }
            //     }
            // }
        }

        // MC atom displacement move for each particle in the system
        for (int i = 0; i < natom_total; ++i)
        {
            mc_atom_displ(i, bondlen, &naccept, density_grids, density_A_ideal, density_B_ideal,
                          particle_type, coords, atom_bond_list, box_size, grid_size, maxsite_1d, chi, kappa, density,
                          kspring, natom_total, grid_shift, idum);
            nattempt += 1;
        }
        // For debug:
        // mc_atom_displ(0, 0.5 * bondlen, &naccept, density_grids, density_A_ideal, density_B_ideal,
        //               particle_type, coords, atom_bond_list, box_size, grid_size, maxsite_1d, chi, kappa, density,
        //               kspring, natom_total, grid_shift, idum);
        // nattempt += 1;

        // save thermo info every so often
        if (t % save_every == 0)
        {
            acc_rate = (double)naccept / nattempt;
            mesh_energy = calc_total_particlemesh_energy(density_grids, chi, kappa, density, maxsite_1d);
            bond_energy = calc_total_bond_energy_harmonic(coords, bonds, nbonds_total, kspring, box_size);
            total_energy = mesh_energy + bond_energy;
            fprintf(thermo, "%d\t%lf\t%lf\t%lf\t%lf\n", t, mesh_energy, bond_energy, total_energy, acc_rate);
        }
    }

    printf("End of Equilibration\n");
    fprintf(thermo, "\nEnd of Equilibration\n\n");

    // production run
    naccept = 0;
    nattempt = 0;

    // Write the initial configuration to xyz file
    FILE *traj;
    char *traj_name = malloc(sizeof(char) * 30);
    sprintf(traj_name, "traj_N%d_A%d_B%d.xyz", natom_perchain, nchain_A, nchain_B);
    traj = fopen(traj_name, "w");
    // fprintf(traj, "#trajectory file\n");
    fclose(traj);

    traj = fopen(traj_name, "a");

    for (t = 0; t <= tot_mc_cycles; ++t)
    {
        // use a new grid_shift every so often
        if (t % grid_shift_every == 0)
        {
            for (int itype = 0; itype < 2; ++itype)
            {
                for (int i = 0; i < maxsite_1d; ++i)
                {
                    for (int j = 0; j < maxsite_1d; ++j)
                    {
                        memset(density_grids[itype][i][j], 0, maxsite_1d * sizeof(double));
                    }
                }
            }
            grid_shift[0] = 0.5 * ran1(idum) * grid_size;
            grid_shift[1] = 0.5 * ran1(idum) * grid_size;
            grid_shift[2] = 0.5 * ran1(idum) * grid_size;
            get_mesh_density(natom_total, grid_size, maxsite_1d, pm_type, density_A_ideal, density_B_ideal,
                             grid_shift, coords, particle_type, density_grids);
        }

        // For debug:
        // for (int itype = 0; itype < 2; ++itype)
        // {
        //     for (int i = 0; i < maxsite_1d; ++i)
        //     {
        //         for (int j = 0; j < maxsite_1d; ++j)
        //         {
        //             for (int k = 0; k < maxsite_1d; ++k)
        //             {
        //                 printf("%f\n", density_grids[itype][i][j][k]);
        //             }
        //         }
        //     }
        // }

        // MC atom displacement move for each particle in the system
        for (int i = 0; i < natom_total; ++i)
        {
            mc_atom_displ(i, bondlen, &naccept, density_grids, density_A_ideal, density_B_ideal,
                          particle_type, coords, atom_bond_list, box_size, grid_size, maxsite_1d, chi, kappa,
                          density, kspring, natom_total, grid_shift, idum);
            nattempt += 1;
        }
        // For debug:
        // mc_atom_displ(0, 0.5 * bondlen, &naccept, density_grids, density_A_ideal, density_B_ideal,
        //               particle_type, coords, atom_bond_list, box_size, grid_size, maxsite_1d, chi, kappa, density,
        //               kspring, natom_total, grid_shift, idum);
        // nattempt += 1;
        // save trajectory, thermo and density profile every so often
        if (t % save_every == 0)
        {
            printf("Timestep %d finished\n", t);
            // trajectory
            fprintf(traj, "%d\n%d 1.500000\n", natom_total, t);
            for (int i = 0; i < natom_total; ++i)
            {
                fprintf(traj, "%d %lf %lf %lf\n", particle_type[i], coords[i][0], coords[i][1], coords[i][2]);
            }
            // thermo
            acc_rate = (double)naccept / nattempt;
            mesh_energy = calc_total_particlemesh_energy(density_grids, chi, kappa, density, maxsite_1d);
            bond_energy = calc_total_bond_energy_harmonic(coords, bonds, nbonds_total, kspring, box_size);
            total_energy = mesh_energy + bond_energy;
            fprintf(thermo, "%d\t%lf\t%lf\t%lf\t%lf\n", t, mesh_energy, bond_energy, total_energy, acc_rate);
            // density profile
            fprintf(density_profile, "t = %d\n", t);
            for (int itype = 0; itype < 2; ++itype)
            {
                for (int i = 0; i < maxsite_1d; ++i)
                {
                    for (int j = 0; j < maxsite_1d; ++j)
                    {
                        for (int k = 0; k < maxsite_1d; ++k)
                        {
                            fprintf(density_profile, "%d\t%d\t%d\t%d\t%f\n", itype, i, j, k, density_grids[itype][i][j][k]);
                        }
                    }
                }
            }
        }
    }

    fclose(traj);
    fclose(thermo);
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

double periodic(double ri, double bl)
{
    if (ri < 0)
        ri += bl;
    else if (ri > bl)
        ri -= bl;

    return ri;
}

double dist_calc(double ra_x, double rb_x, double ra_y, double rb_y, double ra_z, double rb_z, double bl)
{
    double dpos[3], bl_half;

    bl_half = bl / 2.0;

    dpos[0] = ra_x - rb_x;
    dpos[1] = ra_y - rb_y;
    dpos[2] = ra_z - rb_z;

    if (fabs(dpos[0]) > bl_half)
        dpos[0] = fabs(dpos[0]) - bl;
    if (fabs(dpos[1]) > bl_half)
        dpos[1] = fabs(dpos[1]) - bl;
    if (fabs(dpos[2]) > bl_half)
        dpos[2] = fabs(dpos[2]) - bl;

    return sqrt(dpos[0] * dpos[0] + dpos[1] * dpos[1] + dpos[2] * dpos[2]);
}
// initialization of polymer
void initialize_polymer(int nchain, int natom_perchain, int itype, double bb, double box_size,
                        int start_idx, int mol_start_idx, int *particle_type, int *molecule_id,
                        double **coords, int **bonds, int **atom_bond_list, long *idum)
{
    double theta, phi;
    int index;
    double b = sqrt(bb);
    // double b = 1;
    int ibond;
    ibond = start_idx - mol_start_idx;
    for (int i = 0; i < nchain; ++i)
    {
        // Randomizing positions/locations of the first bead(=index 0, 1N, 2N-th...) of each polymer
        index = i * natom_perchain + start_idx;
        coords[index][0] = ran1(idum) * box_size;
        coords[index][1] = ran1(idum) * box_size;
        coords[index][2] = ran1(idum) * box_size;
        particle_type[index] = itype;
        molecule_id[index] = mol_start_idx + i + 1;
        // For debug:
        // printf("%d\t", index);
        // printf("%d\t", particle_type[index]);
        // printf("%d\t", molecule_id[index]);
        // printf("%lf\t", coords[index][0]);
        // printf("%lf\t", coords[index][1]);
        // printf("%lf\t", coords[index][2]);
        // printf("\n");

        // Randomizing positions of the rest beads within the chain spatially
        for (int j = 1; j < natom_perchain; ++j)
        {
            // from index 1(=2nd bead of each chain)
            index = start_idx + i * natom_perchain + j;
            // assign the atom type
            particle_type[index] = itype;
            // assign the molecule ID
            molecule_id[index] = mol_start_idx + i + 1;
            // create an entry for a bond
            bonds[ibond][0] = index - 1;
            bonds[ibond][1] = index;
            ibond += 1;
            // Add the bond info to both atom[index] and atom[index-1]
            atom_bond_list[index][0] = index - 1;
            atom_bond_list[index - 1][1] = index;
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

            // If our new chain goes outside of the periodic box, we wrap around to the other side of the box
            coords[index][0] = periodic(coords[index][0], box_size);
            coords[index][1] = periodic(coords[index][1], box_size);
            coords[index][2] = periodic(coords[index][2], box_size);
            // For debug:
            // printf("%d\t", index);
            // printf("%d\t", particle_type[index]);
            // printf("%d\t", molecule_id[index]);
            // printf("%lf\t", coords[index][0]);
            // printf("%lf\t", coords[index][1]);
            // printf("%lf\t", coords[index][2]);
            // printf("\n");
            // printf("%d %d\n", bonds[ibond][0], bonds[ibond][1]);

        } // iterate to complete one chain
    }     // done initialization
}

int metro_crit(double enrg_diff, long *idum)
{
    int res = 0;
    double ran, prob;
    if (enrg_diff <= 0)
    {
        res = 1;
    }
    else if ((enrg_diff > 0.0) && (enrg_diff <= 40.0))
    {
        ran = ran1(idum);
        prob = exp(-enrg_diff);
        if (ran <= prob)
        {
            res = 1;
        }
    }

    return res;
}

// calculation of density profile of the meshgrids
void get_mesh_density(int natom_total, double grid_size, int maxsite_1d, int pm_type,
                      double density_A_ideal, double density_B_ideal,
                      double grid_shift[3], double **coords, int *particle_type, double ****density_grids)
{
    double rx, ry, rz;
    int xgrid, ygrid, zgrid;
    int atom_type;

    // TODO: set the density grid to zero within the function
    // memset(density_grids, 0, 2 * maxsite_1d * maxsite_1d * maxsite_1d * sizeof(double));
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
            if (atom_type == 0)
            {
                density_grids[0][xgrid][ygrid][zgrid] += 1 / density_A_ideal;
            }
            else if (atom_type == 1)
            {
                density_grids[1][xgrid][ygrid][zgrid] += 1 / density_B_ideal;
            }
        }
    }
    // For debug:
    // double total_des = 0;
    // printf("Grid shift: %f %f %f\n", grid_shift[0], grid_shift[1], grid_shift[2]);
    // for (int itype = 0; itype < 2; ++itype)
    // {
    //     for (int i = 0; i < maxsite_1d; ++i)
    //     {
    //         for (int j = 0; j < maxsite_1d; ++j)
    //         {
    //             for (int k = 0; k < maxsite_1d; ++k)
    //             {
    //                 printf("%f\n", density_grids[itype][i][j][k]);
    //                 total_des += density_grids[itype][i][j][k];
    //             }
    //         }
    //     }
    // }
    // printf("Total density: %f\n", total_des);
}

double calc_total_particlemesh_energy(double ****density_grids, double chi, double kappa,
                                      double density, int maxsite_1d)
{
    double enrg, enrg_pm;
    double pre_factor, kappaN;
    double phiA, phiB;

    enrg_pm = 0.0;
    for (int i = 0; i < maxsite_1d; ++i)
    {
        for (int j = 0; j < maxsite_1d; ++j)
        {
            for (int k = 0; k < maxsite_1d; ++k)
            {
                phiA = density_grids[0][i][j][k];
                phiB = density_grids[1][i][j][k];
                enrg = chi * phiA * phiB + kappa / 2 * (1 - phiA - phiB) * (1 - phiA - phiB);
                enrg_pm += enrg;
            }
        }
    }
    enrg_pm *= density;
    return enrg_pm;
}

double calc_singlemesh_energy(double phiA, double phiB, double chi, double kappa, double density)
{
    double res;
    res = density * (chi * phiA * phiB + kappa / 2 * (1 - phiA - phiB) * (1 - phiA - phiB));
    return res;
}

double calc_total_bond_energy_harmonic(double **coords, int **bonds, int nbonds_total,
                                       double kspring, double box_size)
{
    double res = 0;
    double rx1, ry1, rz1;
    double rx2, ry2, rz2;
    int i, iatm1, iatm2;
    double dist;
    for (i = 0; i < nbonds_total; ++i)
    {
        iatm1 = bonds[i][0];
        iatm2 = bonds[i][1];
        rx1 = coords[iatm1][0];
        ry1 = coords[iatm1][1];
        rz1 = coords[iatm1][2];
        rx2 = coords[iatm2][0];
        ry2 = coords[iatm2][1];
        rz2 = coords[iatm2][2];
        dist = dist_calc(rx1, rx2, ry1, ry2, rz1, rz2, box_size);

        res += kspring * dist * dist;
    }
    return res;
}

double calc_atom_bond_energy_harmonic(double **coords, int **atom_bond_list, int iatom,
                                      double kspring, double box_size)
{
    double rx1, ry1, rz1;
    double rx2, ry2, rz2;
    double dist;
    int atm_j;

    double res = 0;
    rx1 = coords[iatom][0];
    ry1 = coords[iatom][1];
    rz1 = coords[iatom][2];
    for (int j = 0; j < 4; ++j)
    {
        atm_j = atom_bond_list[iatom][j];
        // check if the index is -1 (meaning nothing recorded),
        // if so, do not calculate
        if (atm_j != -1)
        {
            rx2 = coords[atm_j][0];
            ry2 = coords[atm_j][1];
            rz2 = coords[atm_j][2];
            dist = dist_calc(rx1, rx2, ry1, ry2, rz1, rz2, box_size);

            res += kspring * dist * dist;
        }
    }

    return res;
}

void mc_atom_displ(int iatom, double displ, int *naccept, double ****density_grids, double density_A_ideal, double density_B_ideal,
                   int *particle_type, double **coords, int **atom_bond_list, double box_size, double grid_size,
                   int maxsite_1d, double chi, double kappa, double density, double kspring, int natom_total,
                   double grid_shift[3], long *idum)
{
    double rx, ry, rz;
    int xgrid, ygrid, zgrid;
    int atom_type;
    double kappaN, pre_factor;

    atom_type = particle_type[iatom];
    rx = coords[iatom][0];
    ry = coords[iatom][1];
    rz = coords[iatom][2];
    // identify the grid that the particle is originally in, taking a grid shift into account
    xgrid = rx / grid_size + grid_shift[0];
    xgrid = xgrid % maxsite_1d;
    ygrid = ry / grid_size + grid_shift[1];
    ygrid = ygrid % maxsite_1d;
    zgrid = rz / grid_size + grid_shift[2];
    zgrid = zgrid % maxsite_1d;

    double theta, phi;
    double rx_new, ry_new, rz_new;
    int xgrid_new, ygrid_new, zgrid_new;
    // Randomly get theta angle(xy plane)
    theta = ran1(idum) * 2.0 * M_PI;
    // randomly get cos() value ranging from -1 to 1 to set 2Pi rotation angle from z axis
    phi = acos(2.0 * ran1(idum) - 1.0);
    // move the particle by a vector with length displ
    rx_new = rx + displ * cos(theta) * sin(phi);
    ry_new = ry + displ * sin(theta) * sin(phi);
    rz_new = rz + displ * cos(phi);
    // wrap around the periodic boundary condition
    rx_new = periodic(rx_new, box_size);
    ry_new = periodic(ry_new, box_size);
    rz_new = periodic(rz_new, box_size);
    // identify the grid that the particle is moved to, taking a grid shift into account
    xgrid_new = rx_new / grid_size + grid_shift[0];
    xgrid_new = xgrid_new % maxsite_1d;
    ygrid_new = ry_new / grid_size + grid_shift[1];
    ygrid_new = ygrid_new % maxsite_1d;
    zgrid_new = rz_new / grid_size + grid_shift[2];
    zgrid_new = zgrid_new % maxsite_1d;
    // printf("x y z: %d %d %d\n", xgrid, ygrid, zgrid);
    // printf("x_new y_new z_new: %d %d %d\n", xgrid_new, ygrid_new, zgrid_new);

    double phiA_new, phiB_new;
    double phiA_old_from, phiB_old_from, phiA_old_to, phiB_old_to;
    double phiA_new_from, phiB_new_from, phiA_new_to, phiB_new_to;
    double enrg_local_old, enrg_local_new;
    double enrg_bond_old, enrg_bond_new;
    double enrg_diff;
    int iaccept;

    // Update the new density, ONLY WHEN THE GRID CHANGES!
    if (xgrid_new != xgrid || ygrid_new != xgrid || zgrid_new != zgrid)
    {
        // record the old density of the grid that this particle is:
        // moving from (_from) and moving to (_to)
        phiA_old_from = density_grids[0][xgrid][ygrid][zgrid];
        phiB_old_from = density_grids[1][xgrid][ygrid][zgrid];
        phiA_old_to = density_grids[0][xgrid_new][ygrid_new][zgrid_new];
        phiB_old_to = density_grids[1][xgrid_new][ygrid_new][zgrid_new];

        enrg_local_old = calc_singlemesh_energy(phiA_old_from, phiB_old_from, chi, kappa, density) +
                         calc_singlemesh_energy(phiA_old_to, phiB_old_to, chi, kappa, density);
        if (atom_type == 0)
        {
            phiA_new_from = phiA_old_from - 1 / density_A_ideal;
            phiA_new_to = phiA_old_to + 1 / density_A_ideal;
            phiB_new_from = phiB_old_from;
            phiB_new_to = phiB_old_to;
        }
        else
        {
            phiB_new_from = phiB_old_from - 1 / density_B_ideal;
            phiB_new_to = phiB_old_to + 1 / density_B_ideal;
            phiA_new_from = phiA_old_from;
            phiA_new_to = phiA_old_to;
        }
        // calculate the single mesh energy
        enrg_local_new = calc_singlemesh_energy(phiA_new_from, phiB_new_from, chi, kappa, density) +
                         calc_singlemesh_energy(phiA_new_to, phiB_new_to, chi, kappa, density);
        // Calculate bond energy
        enrg_bond_old = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);

        // move the particle to the new position and calculate the new bond energy
        coords[iatom][0] = rx_new;
        coords[iatom][1] = ry_new;
        coords[iatom][2] = rz_new;
        enrg_bond_new = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);

        enrg_diff = (enrg_local_new - enrg_local_old) +
                    (enrg_bond_new - enrg_bond_old);

        // Use metropolis criteria to determine whether to accept the move
        iaccept = metro_crit(enrg_diff, idum);
        // printf("phiA_old_from = %lf, phiA_old_to = %lf\n", phiA_old_from, phiA_old_to);
        // printf("phiB_old_from = %lf, phiB_old_to = %lf\n", phiB_old_from, phiB_old_to);
        // printf("enrg_local_old = %lf, enrg_local_new = %lf\n, \t enrg_bond_old = %lf, enrg_bond_new = %lf\n, \t enrg_diff = %lf, accept = %d\n", enrg_local_old, enrg_local_new, enrg_bond_old, enrg_bond_new, enrg_diff, iaccept);
        *naccept += iaccept;
        if (iaccept == 0)
        {
            // move back the particle if not accepted
            coords[iatom][0] = rx;
            coords[iatom][1] = ry;
            coords[iatom][2] = rz;
        }
        else
        {
            // update the density profile if accepted
            if (atom_type == 0)
            {
                density_grids[0][xgrid][ygrid][zgrid] = phiA_new_from;
                density_grids[0][xgrid_new][ygrid_new][zgrid_new] = phiA_new_to;
            }
            else
            {
                density_grids[1][xgrid][ygrid][zgrid] = phiB_new_from;
                density_grids[1][xgrid_new][ygrid_new][zgrid_new] = phiB_new_to;
            }
        }
    }
    // only calculate bond energy if the grid is not changed
    else
    {
        // Calculate bond energy
        enrg_bond_old = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);

        // move the particle to the new position and calculate the new bond energy
        coords[iatom][0] = rx_new;
        coords[iatom][1] = ry_new;
        coords[iatom][2] = rz_new;
        enrg_bond_new = calc_atom_bond_energy_harmonic(coords, atom_bond_list, iatom, kspring, box_size);
        enrg_diff = enrg_bond_new - enrg_bond_old;

        // Use metropolis criteria to determine whether to accept the move
        iaccept = metro_crit(enrg_diff, idum);
        // printf("phiA_old_from = %lf, phiA_old_to = %lf\n", phiA_old_from, phiA_old_to);
        // printf("phiB_old_from = %lf, phiB_old_to = %lf\n", phiB_old_from, phiB_old_to);
        // printf("enrg_local_old = %lf, enrg_local_new = %lf\n, \t enrg_bond_old = %lf, enrg_bond_new = %lf\n, \t enrg_diff = %lf, accept = %d\n", enrg_local_old, enrg_local_new, enrg_bond_old, enrg_bond_new, enrg_diff, iaccept);
        *naccept += iaccept;
        if (iaccept == 0)
        {
            // move back the particle if not accepted
            coords[iatom][0] = rx;
            coords[iatom][1] = ry;
            coords[iatom][2] = rz;
        }
    }
}
