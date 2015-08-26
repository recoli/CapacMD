/*****************************************************************************
 *  This file is part of the CapacMD program.                                *  
 *  Copyright (C) 2015 Xin Li <lixin.reco@gmail.com>                         *  
 *                                                                           *
 *  Filename:  thermostat.h                                                  *
 *  Function:  temperature and pressure coupling                             *
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

//=======================================
// compute kinetic energy and 
// update temperature and pressure
//=======================================

void kinetic_energy(System *p_system, int nAtoms, Atom_Info *atom_info);

//=======================================
// remove center of mass translation
//=======================================

void remove_comm(int nAtoms, Atom_Info* atom_info, System *p_system);

//=======================================
// Nose-Hoover chain for NVT ensemble
// see Mol. Phys, 1996, 87, 1117-1157
// DOI: 10.1080/00268979600100761
//=======================================

void nose_hoover_chain(RunSet *p_runset, System *p_system, int nAtoms);
