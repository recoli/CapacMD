/*****************************************************************************
 *  This file is part of the CapacMD program.                                *  
 *  Copyright (C) 2015 Xin Li <lixin.reco@gmail.com>                         *  
 *                                                                           *
 *  Filename:  file.h                                                        *
 *  Function:  read input files, write trajectories and charge/dipoles       *
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

//==============================
// read the next line from file
//==============================

void read_next_line(FILE* file_par, char* line, char* subline);

//==================================
// read md options from input_mdset
//==================================

void read_settings(char* input_mdset, RunSet* p_runset, Metal* p_metal);

//=====================
// read input_gro file
//=====================

void read_gro(char* input_gro, System *p_system,
              int* ptr_nAtoms, Atom_Info* atom_info);

//========================
// write to traj.gro file
//========================

void write_gro(FILE* file_gro, System *p_system,
               int nAtoms, Atom_Info* atom_info, int step);

//================================================
// write metal dipole-charge vector to vec_pq.txt
//================================================

void write_vec_pq(FILE *file_pq, Metal *p_metal, int step);

//=================================
// write to binary trajectory file
//=================================

void write_binary(FILE* file_dat, System *p_system, int nAtoms, int step);

//==============================================
// read force field parameters from input_param
//==============================================

void read_param_1(char* input_param, Topol *p_topol, Metal* p_metal);

//=========================================
// read FF parameters (2) from input_param
//=========================================

void read_param_2(char* input_param, Topol *p_topol, Metal *p_metal);

//=============================================
// assign molecular and atomic indices
//=============================================

void assign_indices(Topol *p_topol, Metal *p_metal,
                    Atom_Info* atom_info, Mol_Info* mol_info);
