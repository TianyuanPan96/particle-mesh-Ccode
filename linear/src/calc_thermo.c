#include "header.h"
// calculate all bonds energy
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

// calculate a single mesh's energy
double calc_singlemesh_energy(double phiA, double phiB, double chi, double kappa, double density)
{
    double res;
    res = density * (chi * phiA * phiB + kappa / 2 * (1 - phiA - phiB) * (1 - phiA - phiB));
    return res;
}

// calculate the bond energy connected to a given atom using atom_bond_list
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

// calculate the entire mesh energy using the formula
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