#include "header.h"

int main(int argc, const char *argv[])
{
    long long int natom_perchain, natom_total, natom_A, natom_B, Nbb, Nsc, f_branch;
    long long int nchain, nchain_A, nchain_B;
    int nbonds_perchain, nbonds_total;

    int pm_type;
    double kappa, kspring, chi, bondlen, bondlen_sqr;
    double f_A;
    double grid_size;
    int maxsite_1d;

    int save_every, tot_mc_cycles, equil_cycles, grid_shift_every;
    double ave_bondlen;

    time_t start_t, current_t;
    double diff_t;
    int diff_days, diff_hours, diff_minutes, diff_seconds;
    time(&start_t);

    // File Input
    FILE *inputfile;
    inputfile = fopen("Input.txt", "r");

    // fscanf(inputfile, "nchain_A = %d\n", &nchain_A);             // # of A polymers
    // fscanf(inputfile, "nchain_B = %d\n", &nchain_B);             // # of B polymers

    fscanf(inputfile, "Nbb = %lld\n", &Nbb);           // # of backbone atoms per chain
    fscanf(inputfile, "Nsc = %lld\n", &Nsc);           // # of backbone atoms per chain
    fscanf(inputfile, "f_branch = %lld\n", &f_branch); // # of backbone atoms per chain
    // fscanf(inputfile, "f_A = %lf\n", &f_A);                        // mol fraction of A
    fscanf(inputfile, "pm_type = %d\n", &pm_type);  // PM type (0- or 1-order)
    fscanf(inputfile, "bondlen = %lf\n", &bondlen); // bond length           // End-to-end distance for a polymer
    fscanf(inputfile, "chi = %lf\n", &chi);         // F-H parameter
    fscanf(inputfile, "kappa = %lf\n", &kappa);     // imcompressible parameter
    fscanf(inputfile, "kspring = %lf\n", &kspring); // spring constant

    fscanf(inputfile, "grid_size(Re) = %lf\n", &grid_size);          // mesh size
    fscanf(inputfile, "maxsite_1d = %d\n", &maxsite_1d);             // number of grids per direction
    fscanf(inputfile, "save_every = %d\n", &save_every);             // saving interval
    fscanf(inputfile, "tot_mc_cycles = %d\n", &tot_mc_cycles);       // # of steps (exclude equil. steps)
    fscanf(inputfile, "equil_cycles = %d\n", &equil_cycles);         // # of equil. steps
    fscanf(inputfile, "grid_shift_every = %d\n", &grid_shift_every); // frequency of shift grid

    // fscanf(inputfile, "Nbb = %d\n", &Nbb);
    // fscanf(inputfile, "Nsc = %d\n", &Nsc);
    fclose(inputfile);

    // nchain = nchain_A + nchain_B;          // # of B polymer
    // natom_A = nchain_A * natom_perchain;   // total # of A atoms
    // natom_B = nchain_B * natom_perchain;   // total # of B atoms
    // natom_total = nchain * natom_perchain; // total # of polymers
    // nbonds_perchain = natom_perchain - 1;  // # of bonds per chain
    // nbonds_total = natom_total - nchain;   // total # of bonds

    double density, rho_norm;
    double box_size, box_volume;

    int nsites;
    bondlen_sqr = bondlen * bondlen;
    natom_perchain = Nbb * (f_branch * Nsc + 1);
    nbonds_perchain = natom_perchain - 1; // # of bonds per chain
    double Re, ReReRe;
    Re = sqrt((Nbb - 1) * bondlen_sqr);
    ReReRe = Re * Re * Re;
    grid_size = grid_size * Re;
    // IMPORTANT NOTE:
    // "density" refers to the number density of bead, whereas
    // "rho_norm" refers to the normalized density in perfectly homogeneous melt
    // These two are DIFFERENT!!!
    // Here, assume nint = 30, density = nint / deltaL**3
    density = 30 / (grid_size * grid_size * grid_size);

    // spring constant is 1.5 / b^2
    kspring = kspring / bondlen_sqr;
    box_size = maxsite_1d * grid_size;
    box_volume = box_size * box_size * box_size;
    nsites = maxsite_1d * maxsite_1d * maxsite_1d;

    // nchain_A = f_A * density * box_volume / natom_perchain;       // total # of A atoms
    // nchain_B = (1 - f_A) * density * box_volume / natom_perchain; // total # of B atoms
    // nchain = nchain_A + nchain_B;                                 // # of B polymer
    // natom_A = nchain_A * natom_perchain;
    // natom_B = nchain_B * natom_perchain;
    nchain = density * box_volume / natom_perchain;
    natom_total = natom_perchain * nchain;
    nbonds_total = natom_total - nchain;     // total # of bonds
    rho_norm = (double)natom_total / nsites; // normalized density in perfectly homogeneous melt

    FILE *outputfile;
    outputfile = fopen("output.txt", "w");
    // fprintf(outputfile, "nchain_A = %lld\n", nchain_A);
    // fprintf(outputfile, "nchain_B = %lld\n", nchain_B);
    fprintf(outputfile, "Nbb = %lld\n", Nbb);
    fprintf(outputfile, "Nsc = %lld\n", Nsc);
    fprintf(outputfile, "f_branch = %lld\n", f_branch);
    fprintf(outputfile, "natom_perchain = %lld\n", natom_perchain);
    fprintf(outputfile, "natom_total = %lld\n", natom_total);
    fprintf(outputfile, "nbonds_total = %d\n", nbonds_total);
    fprintf(outputfile, "pm_type = %d\n", pm_type);
    fprintf(outputfile, "density = %lf\n", density);
    fprintf(outputfile, "rho_norm = %lf\n", rho_norm);
    fprintf(outputfile, "bondlen = %lf\n", bondlen);
    fprintf(outputfile, "chi = %lf\n", chi);
    fprintf(outputfile, "kappa = %lf\n", kappa);
    fprintf(outputfile, "kspring = %lf\n", kspring);
    fprintf(outputfile, "box_size = %lf\n", box_size);
    fprintf(outputfile, "grid_size = %lf\n", grid_size);
    fprintf(outputfile, "maxsite_1d = %d\n", maxsite_1d);
    fprintf(outputfile, "nsites = %d\n", nsites);

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
        // Assuming maximum 10 bonds per atom!!!
        atom_bond_list[i] = (int *)calloc(10, sizeof(int));
    }
    // Then initialize this list to -1 to make sure check can be performed
    for (int i = 0; i < natom_total; ++i)
    {
        // Assuming maximum 10 bonds per atom!!!
        for (int j = 0; j < 10; ++j)
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
    initialize_bottlebrush(nchain, 0, bondlen_sqr, Nbb, Nsc, f_branch, box_size, start_idx,
                           mol_start_idx, particle_type, molecule_id, coords, bonds, atom_bond_list, idum);

    // Write the initial configuration to a file
    FILE *init_config;
    init_config = fopen("ini.cfg", "w");

    fprintf(init_config, "#Initial configuration\n");
    fprintf(init_config, "%lld atoms\n", natom_total);
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
    get_mesh_density(natom_total, grid_size, maxsite_1d, pm_type, rho_norm,
                     grid_shift, coords, particle_type, density_grids);

    // equilibration run
    // Write the thermo info to a file
    double mesh_energy, bond_energy, total_energy;
    FILE *thermo;
    char *thermo_filename = malloc(sizeof(char) * FILENAME_MAX);
    sprintf(thermo_filename, "thermo_Nbb%lld_Nsc%lld_f%lld.txt", Nbb, Nsc, f_branch);
    thermo = fopen(thermo_filename, "w");
    fprintf(thermo, "t\tmesh_energy\tbond_energy\ttotal_energy\tacceptance_rate\tave_bondlen\n");
    fclose(thermo);

    // Write the density info to a file
    FILE *density_profile;
    char *density_filename = malloc(sizeof(char) * FILENAME_MAX);
    sprintf(density_filename, "dens_profile_Nbb%lld_Nsc%lld_f%lld.txt", Nbb, Nsc, f_branch);
    density_profile = fopen(density_filename, "w");
    fprintf(density_profile, "type\tx\ty\tz\tdensity\n");
    fclose(density_profile);

    thermo = fopen(thermo_filename, "a");
    density_profile = fopen(density_filename, "a");
    unsigned long long t;
    unsigned long long naccept = 0;
    unsigned long long nattempt = 0;
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
            grid_shift[0] = 0.5 * ran1(idum);
            grid_shift[1] = 0.5 * ran1(idum);
            grid_shift[2] = 0.5 * ran1(idum);
            get_mesh_density(natom_total, grid_size, maxsite_1d, pm_type, rho_norm,
                             grid_shift, coords, particle_type, density_grids);
        }

        // MC atom displacement move for each particle in the system
        for (int i = 0; i < natom_total; ++i)
        {
            mc_atom_displ(i, bondlen, pm_type, &naccept, density_grids, density,
                          particle_type, coords, atom_bond_list, box_size, grid_size, rho_norm, maxsite_1d, chi, kappa,
                          kspring, natom_total, grid_shift, idum);
            nattempt += 1;
        }

        // save thermo info every so often
        if (t % save_every == 0)
        {
            acc_rate = (double)naccept / nattempt;
            mesh_energy = calc_total_particlemesh_energy(density_grids, chi, kappa, density, maxsite_1d);
            bond_energy = calc_total_bond_energy_harmonic(coords, bonds, nbonds_total, kspring, box_size);
            total_energy = mesh_energy + bond_energy;
            ave_bondlen = calc_ave_bondlen(coords, bonds, nbonds_total, box_size);
            fprintf(thermo, "%llu\t%lf\t%lf\t%lf\t%lf\t%lf\n", t, mesh_energy, bond_energy, total_energy, acc_rate, ave_bondlen);
        }
    }

    time(&current_t);
    diff_t = difftime(current_t, start_t);
    diff_days = diff_t / 86400;
    diff_hours = (diff_t - diff_days * 86400) / 3600;
    diff_minutes = (diff_t - diff_days * 86400 - diff_hours * 3600) / 60;
    diff_seconds = diff_t - diff_days * 86400 - diff_hours * 3600 - diff_minutes * 60;

    printf("Time elapsed: %dd-%dh-%dm-%ds\n", diff_days, diff_hours, diff_minutes, diff_seconds);
    printf("End of Equilibration\n");
    fprintf(thermo, "\nEnd of Equilibration\n\n");

    // production run
    naccept = 0;
    nattempt = 0;

    // Write the initial configuration to xyz file
    FILE *traj;
    char *traj_name = malloc(sizeof(char) * FILENAME_MAX);
    sprintf(traj_name, "traj_N%lld_A%lld_B%lld.xyz", natom_perchain, nchain_A, nchain_B);
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
            grid_shift[0] = 0.5 * ran1(idum);
            grid_shift[1] = 0.5 * ran1(idum);
            grid_shift[2] = 0.5 * ran1(idum);
            get_mesh_density(natom_total, grid_size, maxsite_1d, pm_type, rho_norm,
                             grid_shift, coords, particle_type, density_grids);
        }

        // MC atom displacement move for each particle in the system
        for (int i = 0; i < natom_total; ++i)
        {
            mc_atom_displ(i, bondlen, pm_type, &naccept, density_grids, density,
                          particle_type, coords, atom_bond_list, box_size, grid_size, rho_norm, maxsite_1d, chi, kappa,
                          kspring, natom_total, grid_shift, idum);
            nattempt += 1;
        }

        // save trajectory, thermo and density profile every so often
        if (t % save_every == 0)
        {
            time(&current_t);
            diff_t = difftime(current_t, start_t);
            diff_days = diff_t / 86400;
            diff_hours = (diff_t - diff_days * 86400) / 3600;
            diff_minutes = (diff_t - diff_days * 86400 - diff_hours * 3600) / 60;
            diff_seconds = diff_t - diff_days * 86400 - diff_hours * 3600 - diff_minutes * 60;

            printf("Timestep %llu finished. Time elapsed: %dd-%dh-%dm-%ds.\n", t, diff_days, diff_hours, diff_minutes, diff_seconds);
            // trajectory
            fprintf(traj, "%lld\n%llu 1.500000\n", natom_total, t);
            for (int i = 0; i < natom_total; ++i)
            {
                fprintf(traj, "%d %lf %lf %lf\n", particle_type[i], coords[i][0], coords[i][1], coords[i][2]);
            }
            // thermo
            acc_rate = (double)naccept / nattempt;
            mesh_energy = calc_total_particlemesh_energy(density_grids, chi, kappa, density, maxsite_1d);
            bond_energy = calc_total_bond_energy_harmonic(coords, bonds, nbonds_total, kspring, box_size);
            total_energy = mesh_energy + bond_energy;
            ave_bondlen = calc_ave_bondlen(coords, bonds, nbonds_total, box_size);
            // density profile
            fprintf(density_profile, "t = %llu\n", t);
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

    // Write the final configuration to a file
    FILE *final_config;
    final_config = fopen("final.cfg", "w");
    fprintf(final_config, "#Final configuration\n");
    fprintf(final_config, "%lld atoms\n", natom_total);
    fprintf(final_config, "2 atom types\n");
    fprintf(final_config, "%d bonds\n", nbonds_total);
    fprintf(final_config, "1 bond types\n");
    fprintf(final_config, "0.0 %lf xlo xhi\n", box_size);
    fprintf(final_config, "0.0 %lf ylo yhi\n", box_size);
    fprintf(final_config, "0.0 %lf zlo zhi\n", box_size);
    fprintf(final_config, "0 0 0 xy xz yz\n");
    fprintf(final_config, "\n");
    fprintf(final_config, "Atoms # bond\n");
    fprintf(final_config, "\n");
    for (int i = 0; i < natom_total; ++i)
    {
        fprintf(final_config, "%d %d %d %f %f %f\n", i + 1, molecule_id[i], particle_type[i] + 1,
                coords[i][0], coords[i][1], coords[i][2]);
    }
    fprintf(final_config, "\n");
    fprintf(final_config, "Bonds\n");
    fprintf(final_config, "\n");
    for (int i = 0; i < nbonds_total; ++i)
    {
        fprintf(final_config, "%d 1 %d %d\n", i + 1, bonds[i][0] + 1, bonds[i][1] + 1);
    }
    fclose(final_config);
}
