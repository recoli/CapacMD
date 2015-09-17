/*****************************************************************************
 *  This file is part of the CapacMD program.                                *  
 *  Copyright (C) 2015 Xin Li <lixin.reco@gmail.com>                         *  
 *                                                                           *
 *  Filename:  bicg_stab.h                                                   *
 *  Function:  biconjugate gradient stablized method                         *
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

//========================================================
// Parallel matrix-vector multiplication
// Each processor has partial matrix A and full vector x
// Partial vector Ax is stored locally (no communication)
//========================================================

void mpi_mat_vec(int start_metal, int end_metal, 
                 int min_metal, int max_metal, int n_NPs,
                 double** A, double* x, double* Ax, int my_id, int num_procs);

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

void mpi_precon_bicg_stab_COO(int start_metal, int end_metal, int min_metal, int max_metal, int n_NPs,
                              int n_mat, double* diag_relay, double* vec_ext, double* vec_pq,
                              int my_id, int num_procs, Metal *p_metal, long int count_nnz,
                              Bicgstab *p_bicgstab);
