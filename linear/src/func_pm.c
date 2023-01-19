#include "header.h"
// calculation of density profile of the meshgrids
void get_mesh_density(long long int natom_total, double grid_size, int maxsite_1d, int pm_type, double rho_norm,
                      double grid_shift[3], double **coords, int *particle_type, double ****density_grids)
{
    double rx, ry, rz;
    int xgrid, ygrid, zgrid;
    int atom_type;

    // TODO: set the density grid to zero within the function
    // memset(density_grids, 0, 2 * maxsite_1d * maxsite_1d * maxsite_1d * sizeof(double));
    // Nearest Grid Point (NGC / PM0) Scheme, assign the particle to only the nearest grid
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
            density_grids[atom_type][xgrid][ygrid][zgrid] += 1 / rho_norm;
        }
    }

    // Cloud in Cell (CIC / PM1) Scheme, assign the particle to eight grids that the particle
    // is closest to
    if (pm_type == 1)
    {
        double xc, yc, zc;
        double dx, dy, dz;
        double tx, ty, tz;
        int xneigh, yneigh, zneigh;
        for (int i = 0; i < natom_total; ++i)
        {
            atom_type = particle_type[i];
            rx = coords[i][0];
            ry = coords[i][1];
            rz = coords[i][2];
            // identify the grid that the particle is in and the corresponding grid center;
            // taking a grid shift into account
            xgrid = rx / grid_size + grid_shift[0];
            xc = ((double)xgrid + 0.5 - grid_shift[0]) * grid_size;
            xgrid = mod(xgrid, maxsite_1d);
            ygrid = ry / grid_size + grid_shift[1];
            yc = ((double)ygrid + 0.5 - grid_shift[1]) * grid_size;
            ygrid = mod(ygrid, maxsite_1d);
            zgrid = rz / grid_size + grid_shift[2];
            zc = ((double)zgrid + 0.5 - grid_shift[2]) * grid_size;
            zgrid = mod(zgrid, maxsite_1d);
            // the distance of particle w.r.t the center of the grid (NORMALIZED BY GRID SIZE!!!)
            dx = (rx - xc) / grid_size;
            dy = (ry - yc) / grid_size;
            dz = (rz - zc) / grid_size;
            // For debug:
            // if (fabs(dx) > 0.5 || fabs(dy) > 0.5 || fabs(dz) > 0.5)
            // {
            //     printf("Error in distance calculation!");
            // }
            tx = 1 - fabs(dx);
            ty = 1 - fabs(dy);
            tz = 1 - fabs(dz);
            // identify the neighbor grid that the particle is contributing density to
            xneigh = xgrid + ((dx > 0) ? 1 : -1);
            xneigh = mod(xneigh, maxsite_1d);
            yneigh = ygrid + ((dy > 0) ? 1 : -1);
            yneigh = mod(yneigh, maxsite_1d);
            zneigh = zgrid + ((dz > 0) ? 1 : -1);
            zneigh = mod(zneigh, maxsite_1d);
            // taking the normalization factor into account
            density_grids[atom_type][xgrid][ygrid][zgrid] += tx * ty * tz / rho_norm;
            density_grids[atom_type][xneigh][ygrid][zgrid] += dx * ty * tz / rho_norm;
            density_grids[atom_type][xgrid][yneigh][zgrid] += tx * dy * tz / rho_norm;
            density_grids[atom_type][xgrid][ygrid][zneigh] += tx * ty * dz / rho_norm;
            density_grids[atom_type][xneigh][yneigh][zgrid] += dx * dy * tz / rho_norm;
            density_grids[atom_type][xneigh][ygrid][zneigh] += dx * ty * dz / rho_norm;
            density_grids[atom_type][xgrid][yneigh][zneigh] += tx * dy * dz / rho_norm;
            density_grids[atom_type][xneigh][yneigh][zneigh] += dx * dy * dz / rho_norm;
        }
    }
    // For debug:
    double local_des = 0;
    // printf("Grid shift: %f %f %f\n", grid_shift[0], grid_shift[1], grid_shift[2]);
    for (int itype = 0; itype < 2; ++itype)
    {
        for (int i = 0; i < maxsite_1d; ++i)
        {
            for (int j = 0; j < maxsite_1d; ++j)
            {
                for (int k = 0; k < maxsite_1d; ++k)
                {
                    local_des = density_grids[itype][i][j][k];
                    if (local_des < 0)
                        printf("Total density: %f\n", local_des);
                }
            }
        }
    }
}