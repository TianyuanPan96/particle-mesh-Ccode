#ifndef HEADER
#define HEADER

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

int mod(int a, int b);
float ran1(long *idum);
float gasdev(long *idum);
long initRan();
double periodic(double ri, double bl);
double dist_calc(double ra_x, double rb_x, double ra_y, double rb_y, double ra_z, double rb_z, double bl);

void initialize_polymer(long long int nchain, long long int natom_perchain, int itype, double bb, double box_size,
                        int start_idx, int mol_start_idx, int *particle_type, int *molecule_id,
                        double **coords, int **bonds, int **atom_bond_list, long *idum);
void initialize_bottlebrush(long long int nchain, int itype, double bb, int Nbb, int Nsc, double box_size,
                            int start_idx, int mol_start_idx, int *particle_type, int *molecule_id,
                            double **coords, int **bonds, int **atom_bond_list, long *idum);
int metro_crit(double enrg_diff, long *idum);
void get_mesh_density(long long int natom_total, double grid_size, int maxsite_1d, int pm_type, double rho_norm,
                      double grid_shift[3], double **coords, int *particle_type, double ****density_grids);
double calc_total_particlemesh_energy(double ****density_grids, double chi, double kappa,
                                      double density, int maxsite_1d);
double calc_singlemesh_energy(double phiA, double phiB, double chi, double kappa, double density);
double calc_total_bond_energy_harmonic(double **coords, int **bonds, int nbonds_total,
                                       double kspring, double box_size);
double calc_atom_bond_energy_harmonic(double **coords, int **atom_bond_list, int iatom,
                                      double kspring, double box_size);
void mc_atom_displ(int iatom, double displ, int pm_type, unsigned long long *naccept, double ****density_grids, double density,
                   int *particle_type, double **coords, int **atom_bond_list, double box_size, double grid_size, double rho_norm,
                   int maxsite_1d, double chi, double kappa, double kspring, long long int natom_total,
                   double grid_shift[3], long *idum);

#endif