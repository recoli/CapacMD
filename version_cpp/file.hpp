/**
 \file file.hpp
 \author Xin Li
 \date 2015/09
 \version 1.0.2
 
 \brief read input files, write trajectories and charge/dipoles
 
 \copyright
 This file is part of the CapacMD program.
 Copyright (C) 2015 Xin Li <lixin.reco@gmail.com>
 
 \copyright
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along
 with this program; if not, see http://www.gnu.org/licenses/.
 */

#include <string>

//==============================
// read the next line from file
//==============================

void read_next_line(FILE* file_par, char* line, char* subline);

//==================================
// read md options from input_mdset
//==================================

void read_settings(std::string input_mdset, RunSet& s_runset, Metal& s_metal);

//=====================
// read input_gro file
//=====================

void read_gro(std::string input_gro, System& s_system,
              int* ptr_nAtoms, Atom_Info* atom_info);

//========================
// write to traj.gro file
//========================

void write_gro(FILE* file_gro, System& s_system,
               int nAtoms, Atom_Info* atom_info, int step);

//================================================
// write metal dipole-charge vector to vec_pq.txt
//================================================

void write_vec_pq(FILE *file_pq, Metal& s_metal, int step);

//=================================
// write to binary trajectory file
//=================================

void write_binary(FILE* file_dat, System& s_system, int nAtoms, int step);

//==============================================
// read force field parameters from input_param
//==============================================

void read_param_1(std::string input_param, Topol& s_topol, Metal& s_metal);

//=========================================
// read FF parameters (2) from input_param
//=========================================

void read_param_2(std::string input_param, Topol& s_topol, Metal& s_metal);

//=============================================
// assign molecular and atomic indices
//=============================================

void assign_indices(Topol& s_topol, Metal& s_metal,
                    Atom_Info* atom_info, Mol_Info* mol_info);
