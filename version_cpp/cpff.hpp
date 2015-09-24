/**
 \file cpff.hpp
 \author Xin Li
 \date 2015/09
 \version 1.0.2
 
 \brief capacitance-polarizability force field
 
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

//===========================================================
// constructs the CPIM external field and potential, vec_ext
//===========================================================

void mpi_cpff_vec_ext(Task& s_task, Metal& s_metal, RunSet& s_runset,
                      Atom_Info* atom_info, Topol& s_topol,
                      System& s_system, int my_id);

//==================================================
// constructs the CPIM relay matrix
//==================================================

void mpi_cpff_mat_relay_COO(Task& s_task, Metal& s_metal, System& s_system, double rCut2,
                            int my_id, int num_procs, long int *p_count_nnz);

//=======================
// calculate CPIM forces
//=======================

void mpi_cpff_force(Task& s_task, Metal& s_metal, RunSet& s_runset,
                    System& s_system, double* vec_pq,
                    int nAtoms, Atom_Info* atom_info, int my_id);
