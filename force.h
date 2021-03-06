/*****************************************************************************
 *  This file is part of the CapacMD program.                                *  
 *  Copyright (C) 2015 Xin Li <lixin.reco@gmail.com>                         *  
 *                                                                           *
 *  Filename:  force.h                                                       *
 *  Function:  calculate bonded, nonbonded, Q-SC and CPFF forces             *
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

//=============================
// apply PBC to delta-x
//=============================

double pbc(double dx, double xbox);

//========================================
// apply PBC to the whole system
//========================================

void apply_pbc(const int nMols, const Mol_Info* mol_info, 
               double* rx, double* ry, double* rz, double* box);

//========================================
// find starting and ending atom/molecule
// on each processor
//========================================

void find_start_end(int* start_group, int* end_group, 
                    const int total_num, const int num_procs);

//========================================
// find starting and ending atomic pairs
// on each processor
//========================================

void find_start_end_long(long int* start_group, long int* end_group,
                         const long int total_num, const int num_procs);

//=============================
// sum and analyze time_used
//=============================

void sum_time_used(double** time_used, const int my_id, const int num_procs);
void analyze_time_used(double** time_used, const int num_procs);

//=============================
// sum potential energies
//=============================

void sum_potential(double* potential);

//=============================
// print potential energies
//=============================

void print_potential(double* potential);

//==================================================
// get maximal force and root mean square of force
//==================================================

void get_fmax_rms(const int nAtoms, System *p_system);

//==================================================================================
// RATTLE constraint algorithm for the 1st half step
// http://www.ccl.net/cca/software/SOURCES/FORTRAN/allen-tildesley-book/f.09.shtml
//==================================================================================

void rattle_1st(double dt, Mol_Info* mol_info, Atom_Info* atom_info,
                Topol *p_topol, System *p_system);

//==================================================================================
// RATTLE constraint algorithm for the 2nd half step
// http://www.ccl.net/cca/software/SOURCES/FORTRAN/allen-tildesley-book/f.09.shtml
//==================================================================================

void rattle_2nd(double dt, Mol_Info* mol_info, Atom_Info* atom_info,
                Topol *p_topol, System *p_system);

//=============================================
// compute all forces except those from CPIM
//=============================================

void mpi_force(Task *p_task, Topol *p_topol,
               Atom_Info* atom_info, Mol_Info* mol_info, 
               RunSet* p_runset, Metal *p_metal, System *p_system, Bicgstab *p_bicgstab,
               int my_id, int num_procs, double** time_used);
