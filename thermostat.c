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


//=======================================
// Nose-Hoover coupling for NVT ensemble
// see Mol. Phys, 1996, 87, 1117-1157
// DOI: 10.1080/00268979600100761
//=======================================

void nose_hoover(RunSet *p_runset, System *p_system, int nAtoms)
{
    double aQ, scale;
    double dt_4 = 0.5 * p_runset->dt_2;
    double NkT  = p_runset->kT * p_system->ndf;

    aQ = (p_system->eKsum - NkT) / p_system->qMass;
    p_system->vQ += dt_4 * aQ;
    scale = exp(-p_runset->dt_2 * p_system->vQ);

    int i;
    for(i = 0; i < nAtoms; i ++)
    {
        p_system->vx[i] *= scale;
        p_system->vy[i] *= scale;
        p_system->vz[i] *= scale;
    }

    p_system->eKsum     *= scale * scale;
    p_system->inst_temp *= scale * scale;
    
    aQ = (p_system->eKsum - NkT) / p_system->qMass;
    p_system->vQ += dt_4 * aQ;
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


