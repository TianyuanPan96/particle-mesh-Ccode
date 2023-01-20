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
    // assuming maximum 10 bonds per atom!!
    for (int j = 0; j < 10; ++j)
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

// calculate average bond length
double calc_ave_bondlen(double **coords, int **bonds, int nbonds_total, double box_size)
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

        res += dist;
    }
    res /= nbonds_total;
    return res;
}

// calculate the end-to-end distance of all backbones
double calc_reedsq_bb(double **coords, int Nbb, int nchain, int natom_perchain, double box_size)
{
    double res = 0;
    double rx1, ry1, rz1;
    double rx2, ry2, rz2;
    int i, bb_beg, bb_end;
    double dist;
    for (i = 0; i < nchain; ++i)
    {
        bb_beg = i * natom_perchain;
        bb_end = bb_beg + Nbb - 1;
        rx1 = coords[bb_beg][0];
        ry1 = coords[bb_beg][1];
        rz1 = coords[bb_beg][2];
        rx2 = coords[bb_end][0];
        ry2 = coords[bb_end][1];
        rz2 = coords[bb_end][2];
        dist = dist_sq_calc(rx1, rx2, ry1, ry2, rz1, rz2, box_size);

        res += dist;
    }
    res /= nchain;
    return res;
}

// calculate the end-to-end distance of all side chains
double calc_reedsq_sc(double **coords, int Nbb, int Nsc, int f_branch, int nchain, double box_size)
{
    double res = 0;
    double rx1, ry1, rz1;
    double rx2, ry2, rz2;
    int i, j, k, sc_beg, sc_end;
    int natom_perchain, nsidechain;
    double dist;
    natom_perchain = Nbb * (f_branch * Nsc + 1);
    nsidechain = nchain * Nbb * f_branch;
    for (i = 0; i < nchain; ++i)
    {
        for (j = 0; j < Nbb; ++j)
        {
            for (k = 0; k < f_branch; ++k)
            {
                sc_beg = i * natom_perchain + Nbb + (j * f_branch + k) * Nsc;
                sc_end = sc_beg + Nsc - 1;
                rx1 = coords[sc_beg][0];
                ry1 = coords[sc_beg][1];
                rz1 = coords[sc_beg][2];
                rx2 = coords[sc_end][0];
                ry2 = coords[sc_end][1];
                rz2 = coords[sc_end][2];
                dist = dist_sq_calc(rx1, rx2, ry1, ry2, rz1, rz2, box_size);

                res += dist;
            }
        }
    }
    res /= nsidechain;
    return res;
}

// calculate the radius of gyration of all backbones
double calc_rgsq_bb(double **coords, int Nbb, int nchain, int natom_perchain, double box_size)
{
    double res = 0;
    double rg, rgsq;
    double rx, ry, rz;
    double com_rx, com_ry, com_rz;
    double dx, dy, dz;
    int i, j, bb_beg;
    // a buffer variable which stores the position of beads within a molecule
    double **molbuf;
    molbuf = (double **)malloc(Nbb * sizeof(double *));
    for (int i = 0; i < Nbb; ++i)
    {
        molbuf[i] = (double *)malloc(3 * sizeof(double));
    }
    for (i = 0; i < nchain; ++i)
    {
        rg = 0;
        rgsq = 0;
        com_rx = 0;
        com_ry = 0;
        com_rz = 0;
        bb_beg = i * natom_perchain;
        // store coordinates of the backbone beads to the buffer
        for (j = 0; j < Nbb; ++j)
        {
            rx = coords[bb_beg + j][0];
            ry = coords[bb_beg + j][1];
            rz = coords[bb_beg + j][2];
            molbuf[j][0] = rx;
            molbuf[j][1] = ry;
            molbuf[j][2] = rz;
        }
        // calculate the relative distance of the beads w.r.t. the first bead
        // accounting for PBC
        for (j = 0; j < Nbb; ++j)
        {
            molbuf[j][0] = molbuf[j][0] - molbuf[0][0];
            molbuf[j][1] = molbuf[j][1] - molbuf[0][1];
            molbuf[j][2] = molbuf[j][2] - molbuf[0][2];
            molbuf[j][0] = simbox_1dgetimage(molbuf[j][0], box_size);
            molbuf[j][1] = simbox_1dgetimage(molbuf[j][1], box_size);
            molbuf[j][2] = simbox_1dgetimage(molbuf[j][2], box_size);
        }
        // calculate the center of mass of this backbone
        for (j = 0; j < Nbb; ++j)
        {
            com_rx += molbuf[j][0];
            com_ry += molbuf[j][1];
            com_rz += molbuf[j][2];
        }
        com_rx /= Nbb;
        com_ry /= Nbb;
        com_rz /= Nbb;
        // calculate Rg
        for (j = 0; j < Nbb; ++j)
        {
            dx = molbuf[j][0] - com_rx;
            dy = molbuf[j][1] - com_ry;
            dz = molbuf[j][2] - com_rz;
            rgsq += dx * dx + dy * dy + dz * dz;
        }
        rgsq /= Nbb;
        res += rgsq;
    }
    res /= nchain;
    return res;
}