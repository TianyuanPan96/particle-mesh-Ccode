#include "header.h"
// monte carlo move for a single atom displacement
void mc_atom_displ(int iatom, double displ, int pm_type, unsigned long long *naccept, double ****density_grids, double density,
                   int *particle_type, double **coords, int **atom_bond_list, double box_size, double grid_size, double rho_norm,
                   int maxsite_1d, double chi, double kappa, double kspring, long long int natom_total,
                   double grid_shift[3], long *idum)
{
    if (pm_type == 0)
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
                phiA_new_from = phiA_old_from - 1 / rho_norm;
                phiA_new_to = phiA_old_to + 1 / rho_norm;
                phiB_new_from = phiB_old_from;
                phiB_new_to = phiB_old_to;
            }
            else
            {
                phiB_new_from = phiB_old_from - 1 / rho_norm;
                phiB_new_to = phiB_old_to + 1 / rho_norm;
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
            // For debug:
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
    if (pm_type == 1)
    {
        printf("No implementation of PM1 yet! Exiting the program....\n");
        exit(0);
    }
}