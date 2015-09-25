/**
 \file bicg_stab.hpp
 \author Xin Li
 \date 2015/09
 \version 1.0.2
 
 \brief biconjugate gradient stablized method
 
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

//===========================================
// parallel vec-vec inner product
// the final sum is broadcasted to all procs
//===========================================

void mpi_vec_vec(int start_metal, int end_metal, int min_metal, int max_metal, int n_NPs,
                 double* a, double* b, double* ptr_aTb, int my_id, int num_procs);

//==============================================
// communicate a vector Ax among all processors
// for subsequent mat-vec multiplication
//==============================================

void mpi_comm_vec(int start_metal, int end_metal, int min_metal, int max_metal, int n_NPs,
                  double* Ax, int my_id, int num_procs);

//======================================================================
// parallel preconditioned BiCGSTAB solver for Ax=b
// for sparse matrix
// http://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
//======================================================================

void mpi_precon_bicg_stab_COO(Task& s_task, Metal& s_metal, int n_mat,
                              int my_id, int num_procs,
                              long int count_nnz, Bicgstab& s_bicgstab);
