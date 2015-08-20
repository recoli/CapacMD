/*****************************************************************************
 *  This file is part of the CapacMD program.                                *  
 *  Copyright (C) 2015 Xin Li <lixin.reco@gmail.com>                         *  
 *                                                                           *
 *  Filename:  cpff.h                                                        *
 *  Function:  capacitance-polarizability force field                        *
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

//===========================================================
// constructs the CPIM external field and potential, vec_ext
//===========================================================

void mpi_cpff_vec_ext(Task *p_task, Metal* p_metal, RunSet* p_runset, 
                      Atom_Info* atom_info, Topol *p_topol,
                      System *p_system, int my_id);

//==================================================
// constructs the CPIM relay matrix
//==================================================

void mpi_cpff_mat_relay_CRS(Task *p_task, Metal *p_metal, System *p_system, double rCut2,
                            int my_id, int num_procs, 
                            long int *p_count_size, long int incr_size, long int *p_count_nnz);

//=======================
// calculate CPIM forces
//=======================

void mpi_cpff_force(Task *p_task, Metal *p_metal, RunSet *p_runset,
                    System *p_system, double* vec_pq,
                    int nAtoms, Atom_Info* atom_info, int my_id);

//====================================================
// evaluate CPIM energies from matrices and vectors
//====================================================

void eval_cpff_pot(int n_mat, int n_NPs, double** mat_relay, double* vec_ext, double* vec_pq,
                   double* potential);
