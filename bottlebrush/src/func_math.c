#include "header.h"

// get the mod between two integers (guarantee only positive number result)
int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}

// calculate the position (scaler) due to periodic boundary condition
double periodic(double ri, double bl)
{
    if (ri < 0)
        ri += bl;
    else if (ri > bl)
        ri -= bl;

    return ri;
}

// distance calculation that accounts for periodic boundary condition
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

// distance (squared) that accounts for periodic boundary condition
double dist_sq_calc(double ra_x, double rb_x, double ra_y, double rb_y, double ra_z, double rb_z, double bl)
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

    return dpos[0] * dpos[0] + dpos[1] * dpos[1] + dpos[2] * dpos[2];
}

// get the periodic image of a position (1d)
double simbox_1dgetimage(double rx, double bl)
{
    double res = 0;
    res = rx - bl * round(rx / bl);
    return res;
}

// metropolis criteria determination of monte carlo moves
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