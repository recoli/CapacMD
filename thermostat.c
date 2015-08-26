/*****************************************************************************
 *  This file is part of the CapacMD program.                                *  
 *  Copyright (C) 2015 Xin Li <lixin.reco@gmail.com>                         *  
 *                                                                           *
 *  Filename:  thermostat.c                                                  *
 *  Function:  temperature coupling                                          *
 *  Version:   1.0                                                           *
 *  Updated:   2015-Jun-30                                                   *
 *  License:   GNU Public License, version 2                                 *
 *                                                                           *  
 *  This program is free software; you can redistribute it and/or modify     *  
 *  it under the terms of the GNU General Public License as published by     *  
 *  the Free Software Foundation; either version 2 of the License, or        *  
 *  (at your option) any later version.                                      *  
 *                                                                           *
 *  This program is distributed in the hope that it will be useful,          *  
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *  
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *  
 *  GNU General Public License for more details.                             *  
 *                                                                           *
 *  You should have received a copy of the GNU General Public License along  *  
 *  with this program; if not, write to the Free Software Foundation, Inc.,  *  
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.              *  
 *****************************************************************************/

#include <math.h>
#include "typedef.h"

//=================================
// compute kinetic energy and 
// update temperature
//=================================

void kinetic_energy(System *p_system, int nAtoms, Atom_Info *atom_info)
{
    // diagonal components of the p^2/m tensor
    double eKsum_xx = 0.0;
    double eKsum_yy = 0.0;
    double eKsum_zz = 0.0;
    double m;

    int i;
    for(i = 0; i < nAtoms; i ++)
    {
        m = atom_info[i].mass;
        // ignore the virtual sites
        if(m > 0.0)
        {
            eKsum_xx += m * p_system->vx[i] * p_system->vx[i];
            eKsum_yy += m * p_system->vy[i] * p_system->vy[i];
            eKsum_zz += m * p_system->vz[i] * p_system->vz[i];
        }
    }

    // update scalar p^2/m and instant temperature
    p_system->eKsum = eKsum_xx + eKsum_yy + eKsum_zz;
    p_system->inst_temp = p_system->eKsum / (K_BOLTZ * p_system->ndf);

    /*
    ! For now only NVT ensemble is supported.
    ! Pressure coupling not yet implemented.
    */
}

//===================================
// remove center of mass translation
//===================================

void remove_comm(int nAtoms, Atom_Info* atom_info, System *p_system)
{
    // find center of mass velocity
    double m_sum, vx_com, vy_com, vz_com, m;

    vx_com = 0.0;
    vy_com = 0.0;
    vz_com = 0.0;
    m_sum  = 0.0;

    int i;
    for(i = 0; i < nAtoms; ++ i) 
    {
        m = atom_info[i].mass;
        if(m > 0.0)
        {
            vx_com += p_system->vx[i] * m;
            vy_com += p_system->vy[i] * m;
            vz_com += p_system->vz[i] * m;
            m_sum  += m;
        }
    }

    vx_com /= m_sum;
    vy_com /= m_sum;
    vz_com /= m_sum;

    // zero center of mass velocity
    for(i = 0; i < nAtoms; i ++) 
    {
        m = atom_info[i].mass;
        // ignore the virtual sites
        if(m > 0.0)
        {
            p_system->vx[i] -= vx_com;
            p_system->vy[i] -= vy_com;
            p_system->vz[i] -= vz_com;
        }
    }
}


//=======================================
// Nose-Hoover chain for NVT ensemble
// see Mol. Phys, 1996, 87, 1117-1157
// DOI: 10.1080/00268979600100761
//=======================================

void nose_hoover_chain(RunSet *p_runset, System *p_system, int nAtoms)
{
    const int Nc  = 5;
    const int Nys = 5;

    //double w1 = 1.0 / (4.0 - pow(4.0, 1.0/3.0));
    //double w3 = 1.0 - 4.0 * w1;
    double w1 =  0.41449077179437571194;
    double w3 = -0.65796308717750284778;
    double w[5] = {w1, w1, w3, w1, w1};

    int    M     = p_system->num_nhc;
    double qMass = p_system->qMass;
    double eKsum = p_system->eKsum;
    double *aQ   = p_system->aQ;
    double *vQ   = p_system->vQ;

    double dt  = p_runset->dt;
    double kT  = p_runset->kT;
    double NkT = kT * p_system->ndf;

    int im, ic, iys;
    double aa;
    double scale = 1.0;

    // initialize aQ from [0] to [M-1]
    aQ[0] = (eKsum - NkT) / qMass;
    for (im = 0; im < M-1; ++ im)
    {
        aQ[im+1] = (qMass * vQ[im] * vQ[im] - kT) / qMass;
    }

    // Nose-Hoover chain
    for (ic = 0; ic < Nc; ++ ic)
    {
        for (iys = 0; iys < Nys; ++ iys)
        {
            // update vQ from [M-1] to [0]
            vQ[M-1] += aQ[M-1] * w[iys] * dt / (4 * Nc); 
            for (im = 0; im < M-1; ++ im)
            {
                aa = exp(-w[iys] * dt / (8 * Nc) * vQ[M-(im+1)]);
                vQ[M-1-(im+1)] = vQ[M-1-(im+1)] * aa * aa + 
                    w[iys] * dt / (4 * Nc) * aQ[M-1-(im+1)] * aa;
            }

            // update kinetic energy
            aa = exp(-w[iys] * dt / (2 * Nc) * vQ[0]);
            scale *= aa;
            aQ[0] = (eKsum * scale * scale - NkT) / qMass;

            // update vQ from [0] to [M-1]
            for (im = 0; im < M-1; ++ im)
            {
                aa = exp(-w[iys] * dt / (8 * Nc) * vQ[im+1]);
                vQ[im] = vQ[im] * aa * aa + 
                    w[iys] * dt / (4 * Nc) * aQ[im] * aa;
                aQ[im+1] = (qMass * vQ[im] * vQ[im] - kT) / qMass;
            }
            vQ[M-1] += aQ[M-1] * w[iys] * dt / (4 * Nc);
        }
    }

    // update velocities, kinetic energy, and temperature
    int i;
    for(i = 0; i < nAtoms; i ++)
    {
        p_system->vx[i] *= scale;
        p_system->vy[i] *= scale;
        p_system->vz[i] *= scale;
    }

    p_system->eKsum     *= scale * scale;
    p_system->inst_temp *= scale * scale;
}
