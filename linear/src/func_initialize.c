#include "header.h"

// initialization of linear polymer
void initialize_polymer(long long int nchain, long long int natom_perchain, int itype, double bb, double box_size,
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

// initialization of bottlebrush polymers
void initialize_bottlebrush(long long int nchain, int itype, double bb, int Nbb, int Nsc, double box_size,
                            int start_idx, int mol_start_idx, int *particle_type, int *molecule_id,
                            double **coords, int **bonds, int **atom_bond_list, long *idum)
{
    double theta, phi;
    int index, start_thischain;
    double b = sqrt(bb);
    int ibond;
    int natom_perchain;
    ibond = start_idx - mol_start_idx;
    natom_perchain = Nbb * (Nsc + 1);
    for (int i = 0; i < nchain; ++i)
    {
        // Randomizing positions/locations of the first bead(=index 0, 1N, 2N-th...) of each polymer
        start_thischain = start_idx + i * natom_perchain;
        index = start_thischain;
        coords[index][0] = ran1(idum) * box_size;
        coords[index][1] = ran1(idum) * box_size;
        coords[index][2] = ran1(idum) * box_size;
        particle_type[index] = itype;
        molecule_id[index] = mol_start_idx + i + 1;

        // Create the backbone first
        for (int j = 1; j < Nbb; ++j)
        {
            // from index 1(=2nd bead of each chain)
            index = start_thischain + j;
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
        }

        // then for each backbone bead, attach a side-chain to it
        for (int ibb = 0; ibb < Nbb; ++ibb)
        {
            index = start_thischain + Nbb + ibb * Nsc;
            particle_type[index] = itype;
            molecule_id[index] = mol_start_idx + i + 1;
            bonds[ibond][0] = start_thischain + ibb;
            bonds[ibond][1] = index;
            ibond += 1;
            atom_bond_list[index][0] = start_thischain + ibb;
            atom_bond_list[index - 1][1] = index;
            theta = ran1(idum) * 2.0 * M_PI;
            phi = acos(2.0 * ran1(idum) - 1.0);
            coords[index][0] = coords[index - 1][0] + b * cos(theta) * sin(phi);
            coords[index][1] = coords[index - 1][1] + b * sin(theta) * sin(phi);
            coords[index][2] = coords[index - 1][2] + b * cos(phi);
            coords[index][0] = periodic(coords[index][0], box_size);
            coords[index][1] = periodic(coords[index][1], box_size);
            coords[index][2] = periodic(coords[index][2], box_size);

            for (int k = 1; k < Nsc; ++k)
            {
                index = start_thischain + Nbb + ibb * Nsc + k;
                particle_type[index] = itype;
                molecule_id[index] = mol_start_idx + i + 1;
                bonds[ibond][0] = index - 1;
                bonds[ibond][1] = index;
                ibond += 1;
                atom_bond_list[index][0] = index - 1;
                atom_bond_list[index - 1][1] = index;
                theta = ran1(idum) * 2.0 * M_PI;
                phi = acos(2.0 * ran1(idum) - 1.0);
                coords[index][0] = coords[index - 1][0] + b * cos(theta) * sin(phi);
                coords[index][1] = coords[index - 1][1] + b * sin(theta) * sin(phi);
                coords[index][2] = coords[index - 1][2] + b * cos(phi);
                coords[index][0] = periodic(coords[index][0], box_size);
                coords[index][1] = periodic(coords[index][1], box_size);
                coords[index][2] = periodic(coords[index][2], box_size);
            }
        }
    }
}