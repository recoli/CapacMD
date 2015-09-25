/**
 \file bicg_stab.cpp
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

#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <vector>

#include "typedef.hpp"
#include "file.hpp"
#include "cpff.hpp"

//========================================================
/// \brief Parallel matrix-vector multiplication
/// for sparse matrix in COO format.
/// Each processor has partial matrix A and full vector x.
/// Partial vector Ax is stored locally (no communication).
//========================================================

void mpi_mat_vec_COO(long int count_nnz, int n_mat, Metal& s_metal,
                     double* x, double* Ax)
{
    // zero the Ax vector
    for (int i = 0; i < n_mat; ++ i)
    {
        Ax[i] = 0.0;
    }

    // add contributions to the Ax vector
    for (long int idx = 0; idx < count_nnz; ++ idx)
    {
        long int irow = s_metal.vec_row_ind[idx];
        long int icol = s_metal.vec_col_ind[idx];

        Ax[irow] += s_metal.vec_val[idx] * x[icol];
    }
}

//========================================================
// Parallel matrix-vector multiplication
// Each processor has partial matrix A and full vector x
// Partial vector Ax is stored locally (no communication)
//========================================================

void mpi_mat_vec(int start_metal, int end_metal, 
                 int min_metal, int max_metal, int n_NPs,
                 double** A, double* x, double* Ax, int my_id, int num_procs)
{
    int n_metal = max_metal - min_metal + 1;
    int n_mat = n_metal*4 + n_NPs;

    int i_metal, i, j;

    for (i_metal = start_metal; i_metal <= end_metal; i_metal ++)
    {
        i = i_metal - min_metal;

        // first three for dipole, last one for charge
        Ax[i*3] = 0.0;
        Ax[i*3 + 1] = 0.0;
        Ax[i*3 + 2] = 0.0;
        Ax[n_metal*3 + i] = 0.0;

        for (j = 0; j < n_mat; j ++)
        {
            Ax[i*3] += A[i*3][j] * x[j];
            Ax[i*3+1] += A[i*3+1][j] * x[j];
            Ax[i*3+2] += A[i*3+2][j] * x[j];
            Ax[n_metal*3+i] += A[n_metal*3+i][j] * x[j];
        }
    }

    if (ROOT_PROC == my_id)
    {
        // calculate the last element on root processor
        int iNP;
        for (iNP = 0; iNP < n_NPs; iNP ++)
        {
            Ax[n_metal*4 + iNP] = 0.0;
            for (j = 0; j < n_mat; j ++) {
                Ax[n_metal*4 + iNP] += A[n_metal*4 + iNP][j] * x[j];
            }
        }
    }
}

//=================================================
/// \brief Parallel vec-vec inner product.
/// The final sum is broadcasted to all procs.
//=================================================

void mpi_vec_vec(int start_metal, int end_metal, int min_metal, int max_metal, int n_NPs,
                 double* a, double* b, double* ptr_aTb, int my_id, int num_procs)
{
    int n_metal = max_metal - min_metal + 1;

    MPI::Status status;

    int tag_51 = 51; // slave to master, partial_aTb
    int tag_52 = 52; // master to slave, aTb
    
    *ptr_aTb = 0.0;
    for (int i_metal = start_metal; i_metal <= end_metal; ++ i_metal)
    {
        int i = i_metal - min_metal;

        // first three for dipole, last one for charge
        *ptr_aTb += a[i*3    ] * b[i*3    ];
        *ptr_aTb += a[i*3 + 1] * b[i*3 + 1];
        *ptr_aTb += a[i*3 + 2] * b[i*3 + 2];
        *ptr_aTb += a[n_metal*3 + i] * b[n_metal*3 + i];
    }

    double partial_aTb;
    if (ROOT_PROC == my_id)
    {
        // get partial sum from slave processors
        for (int an_id = 1; an_id < num_procs; ++ an_id)
        {
            MPI::COMM_WORLD.Recv(&partial_aTb, 1, MPI::DOUBLE, MPI::ANY_SOURCE,
                                 tag_51, status);
            *ptr_aTb += partial_aTb;
        }

        // add the last element on root processor
        for (int iNP = 0; iNP < n_NPs; ++ iNP)
        {
            *ptr_aTb += a[n_metal*4 + iNP] * b[n_metal*4 + iNP];
        }

        // send the sum to slave processor
        for (int an_id = 1; an_id < num_procs; ++ an_id)
        {
            MPI::COMM_WORLD.Send(ptr_aTb, 1, MPI::DOUBLE, an_id, tag_52);
        }
    }
    else
    {
        // send partial sum to root processor
        MPI::COMM_WORLD.Send(ptr_aTb, 1, MPI::DOUBLE, ROOT_PROC, tag_51);

        // get the sum from root processor
        MPI::COMM_WORLD.Recv(ptr_aTb, 1, MPI::DOUBLE, ROOT_PROC,
                             tag_52, status);
    }
}


//=================================================
// communicate a vector Ax among all processors
// for subsequent mat-vec multiplication
//=================================================

void mpi_comm_vec(int start_metal, int end_metal, int min_metal, int max_metal, int n_NPs,
                  double* Ax, int my_id, int num_procs)
{
    // communicate inv_sqrt_dens among all processors
    MPI::Status status;

    int tag_91 = 91; // any to any, start
    int tag_92 = 92; // any to any, num
    int tag_93 = 93; // any to any, Ax (dipole)
    int tag_94 = 94; // any to any, Ax (charge)

    int n_metal = max_metal - min_metal + 1;

    int start = start_metal - min_metal;
    int num = end_metal - start_metal + 1;

    int s1, s2, n1, n2;

    // master: receive data
    if (ROOT_PROC == my_id)
    {
        for (int an_id = 1; an_id < num_procs; ++ an_id)
        {
            MPI::COMM_WORLD.Recv(&start, 1, MPI::INT, an_id, tag_91, status);
            MPI::COMM_WORLD.Recv(&num,   1, MPI::INT, an_id, tag_92, status);

            // Note: start and num are updated. 
            // So we need to update s1,s2,n1,n2 as well.
            s1 = start*3;
            s2 = n_metal*3 + start;
            n1 = num*3;
            n2 = num;

            MPI::COMM_WORLD.Recv(&(Ax[s1]), n1, MPI::DOUBLE, an_id, tag_93, status);
            MPI::COMM_WORLD.Recv(&(Ax[s2]), n2, MPI::DOUBLE, an_id, tag_94, status);
        }
    }
    // slave: send data
    else
    {
        s1 = start*3;
        s2 = n_metal*3 + start;
        n1 = num*3;
        n2 = num;

        MPI::COMM_WORLD.Send(&start, 1, MPI::INT, ROOT_PROC, tag_91);
        MPI::COMM_WORLD.Send(&num,   1, MPI::INT, ROOT_PROC, tag_92);

        MPI::COMM_WORLD.Send(&(Ax[s1]), n1, MPI::DOUBLE, ROOT_PROC, tag_93);
        MPI::COMM_WORLD.Send(&(Ax[s2]), n2, MPI::DOUBLE, ROOT_PROC, tag_94);
    }

    // broadcast the last element
    MPI::COMM_WORLD.Bcast(&(Ax[0]), n_metal * 4 + n_NPs, MPI::DOUBLE, ROOT_PROC);
}


//======================================================================
// parallel preconditioned BiCGSTAB solver for Ax=b
// for sparse matrix
// http://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
//======================================================================

void mpi_precon_bicg_stab_COO(Task& s_task, Metal& s_metal, int n_mat,
                              int my_id, int num_procs,
                              long int count_nnz, Bicgstab& s_bicgstab)
{
    int start_metal = s_task.start_metal[my_id];
    int end_metal   = s_task.end_metal[my_id];
    int min_metal   = s_metal.min;
    int max_metal   = s_metal.max;
    int n_NPs       = s_metal.n_NPs;

    double *diag_relay = s_metal.diag_relay;
    double *vec_ext    = s_metal.vec_ext;
    double *vec_pq     = s_metal.vec_pq;

    int n_metal = max_metal - min_metal + 1;
    
    // A == relay
    // b == vec_ext
    // x == vec_pq
    
    // At the beginning, each processor already has a copy of vec_pq

    double* Ax = s_bicgstab.Ax;
    double* r0 = s_bicgstab.r0;
    double* r  = s_bicgstab.r;
    double* p  = s_bicgstab.p;
    double* v  = s_bicgstab.v;
    double* s  = s_bicgstab.s;
    double* t  = s_bicgstab.t;
    double* y  = s_bicgstab.y;
    double* z  = s_bicgstab.z;
    double* Kt = s_bicgstab.Kt;
    double* K  = s_bicgstab.K;

    double rho_old = 1.0;
    double omega = 1.0;
    double alpha = 1.0;
    double rho, beta;

    // bTb = bT * b
    double bTb;
    mpi_vec_vec(start_metal, end_metal, min_metal, max_metal, n_NPs,
                vec_ext, vec_ext, &bTb, my_id, num_procs);

    // if vec_ext is zero vector then skip
    if (bTb < 1.0e-10) { return; }

    // convergence criteria for rTr
    double criteria = bTb * 1.0e-08;

    // Ax = A * x = relay * vec_pq
    mpi_mat_vec_COO(count_nnz, n_mat, s_metal, vec_pq, Ax);


    // set initial values
    // r0 = b - Ax,  r = r0,  v = p = 0,  K ~ A^-1
    for (int i_metal = start_metal; i_metal <= end_metal; ++ i_metal)
    {
        int i = i_metal - min_metal;

        int a1 = i*3;
        int a2 = i*3 + 1;
        int a3 = i*3 + 2;
        int a4 = n_metal*3 + i;

        r0[a1] = vec_ext[a1] - Ax[a1];
        r0[a2] = vec_ext[a2] - Ax[a2];
        r0[a3] = vec_ext[a3] - Ax[a3];
        r0[a4] = vec_ext[a4] - Ax[a4];

        r[a1] = r0[a1];
        r[a2] = r0[a2];
        r[a3] = r0[a3];
        r[a4] = r0[a4];

        v[a1] = 0.0;
        v[a2] = 0.0;
        v[a3] = 0.0;
        v[a4] = 0.0;

        p[a1] = 0.0;
        p[a2] = 0.0;
        p[a3] = 0.0;
        p[a4] = 0.0;

        K[a1] = 1.0 / diag_relay[a1];
        K[a2] = 1.0 / diag_relay[a2];
        K[a3] = 1.0 / diag_relay[a3];
        K[a4] = 1.0 / diag_relay[a4];
    }
    if (ROOT_PROC == my_id)
    {
        for (int iNP = 0; iNP < n_NPs; ++ iNP)
        {
            int a5 = n_metal*4 + iNP;
            r0[a5] = vec_ext[a5] - Ax[a5];
            r[a5] = r0[a5];
            v[a5] = 0.0;
            p[a5] = 0.0;
            K[a5] = 1.0;
        }
    }

    int iter = 0;
    while (1)
    {
        // preconditioned BiCGSTAB
        
        // rho = (r0, r)
        mpi_vec_vec(start_metal, end_metal, min_metal, max_metal, n_NPs,
                    r0, r, &rho, my_id, num_procs);

        if (rho < 1.0e-10 && iter == 0) { return; }

        // beta
        beta = (rho / rho_old) * (alpha / omega);

        // p = r + beta * (p - omega * v)
        // y = K * p
        for (int i_metal = start_metal; i_metal <= end_metal; ++ i_metal)
        {
            int i = i_metal - min_metal;

            int a1 = i*3;
            int a2 = i*3 + 1;
            int a3 = i*3 + 2;
            int a4 = n_metal*3 + i;

            p[a1] = r[a1] + beta * (p[a1] - omega * v[a1]);
            p[a2] = r[a2] + beta * (p[a2] - omega * v[a2]);
            p[a3] = r[a3] + beta * (p[a3] - omega * v[a3]);
            p[a4] = r[a4] + beta * (p[a4] - omega * v[a4]);

            y[a1] = p[a1] * K[a1];
            y[a2] = p[a2] * K[a2];
            y[a3] = p[a3] * K[a3];
            y[a4] = p[a4] * K[a4];
        }
        if (ROOT_PROC == my_id)
        {
            for (int iNP = 0; iNP < n_NPs; ++ iNP)
            {
                int a5 = n_metal*4 + iNP;
                p[a5] = r[a5] + beta * (p[a5] - omega * v[a5]);

                y[a5] = p[a5];
            }
        }

        // v = A * y
        mpi_comm_vec(start_metal, end_metal, min_metal, max_metal, n_NPs,
                     y, my_id, num_procs);

        mpi_mat_vec_COO(count_nnz, n_mat, s_metal, y, v);

        // alpha = rho / (r0, v)
        double r0_v;
        mpi_vec_vec(start_metal, end_metal, min_metal, max_metal, n_NPs,
                    r0, v, &r0_v, my_id, num_procs);
        alpha = rho / r0_v;

        // s = r - alpha * v
        // z = K * s
        for (int i_metal = start_metal; i_metal <= end_metal; ++ i_metal)
        {
            int i = i_metal - min_metal;

            int a1 = i*3;
            int a2 = i*3 + 1;
            int a3 = i*3 + 2;
            int a4 = n_metal*3 + i;

            s[a1] = r[a1] - alpha * v[a1];
            s[a2] = r[a2] - alpha * v[a2];
            s[a3] = r[a3] - alpha * v[a3];
            s[a4] = r[a4] - alpha * v[a4];

            z[a1] = s[a1] * K[a1];
            z[a2] = s[a2] * K[a2];
            z[a3] = s[a3] * K[a3];
            z[a4] = s[a4] * K[a4];
        }
        if (ROOT_PROC == my_id)
        {
            for (int iNP = 0; iNP < n_NPs; ++ iNP)
            {
                int a5 = n_metal*4 + iNP;
                s[a5] = r[a5] - alpha * v[a5];

                z[a5] = s[a5];
            }
        }

        // t = A * z
        mpi_comm_vec(start_metal, end_metal, min_metal, max_metal, n_NPs,
                     z, my_id, num_procs);

        mpi_mat_vec_COO(count_nnz, n_mat, s_metal, z, t);

        // compute Kt = K * t
        for (int i_metal = start_metal; i_metal <= end_metal; ++ i_metal)
        {
            int i = i_metal - min_metal;

            int a1 = i*3;
            int a2 = i*3 + 1;
            int a3 = i*3 + 2;
            int a4 = n_metal*3 + i;

            Kt[a1] = t[a1] * K[a1];
            Kt[a2] = t[a2] * K[a2];
            Kt[a3] = t[a3] * K[a3];
            Kt[a4] = t[a4] * K[a4];
        }
        if (ROOT_PROC == my_id)
        {
            for (int iNP = 0; iNP < n_NPs; ++ iNP)
            {
                int a5 = n_metal*4 + iNP;
                Kt[a5] = t[a5];
            }
        }

        // omega = (Kt, Ks) / (Kt, Kt),  note Ks == z
        double Kt_Ks, Kt_Kt;
        mpi_vec_vec(start_metal, end_metal, min_metal, max_metal, n_NPs,
                    Kt, z, &Kt_Ks, my_id, num_procs);
        mpi_vec_vec(start_metal, end_metal, min_metal, max_metal, n_NPs,
                    Kt, Kt, &Kt_Kt, my_id, num_procs);
        omega = Kt_Ks / Kt_Kt;

        // x = x + alpha * y + omega * z,  x == vec_pq
        // r = s - omega * t
        for (int i_metal = start_metal; i_metal <= end_metal; ++ i_metal)
        {
            int i = i_metal - min_metal;

            int a1 = i*3;
            int a2 = i*3 + 1;
            int a3 = i*3 + 2;
            int a4 = n_metal*3 + i;

            vec_pq[a1] += alpha * y[a1] + omega * z[a1];
            vec_pq[a2] += alpha * y[a2] + omega * z[a2];
            vec_pq[a3] += alpha * y[a3] + omega * z[a3];
            vec_pq[a4] += alpha * y[a4] + omega * z[a4];

            r[a1] = s[a1] - omega * t[a1];
            r[a2] = s[a2] - omega * t[a2];
            r[a3] = s[a3] - omega * t[a3];
            r[a4] = s[a4] - omega * t[a4];
        }
        if (ROOT_PROC == my_id)
        {
            for (int iNP = 0; iNP < n_NPs; ++ iNP)
            {
                int a5 = n_metal*4 + iNP;
                vec_pq[a5] += alpha * y[a5] + omega * z[a5];
                r[a5] = s[a5] - omega * t[a5];
            }
        }

        // update rho_old
        rho_old = rho;

        iter ++;

        // check convergence
        double rTr;
        mpi_vec_vec(start_metal, end_metal, min_metal, max_metal, n_NPs,
                    r, r, &rTr, my_id, num_procs);

#ifdef DEBUG
        if (ROOT_PROC == my_id)
        {
            printf("iter= %d, rTr=%e, bTb= %e\n", iter, rTr, bTb);
        }
#endif

        if (rTr < criteria)
        {
            // converged, communicate vec_pq and quit the iteration
            mpi_comm_vec(start_metal, end_metal, min_metal, max_metal, n_NPs,
                         vec_pq, my_id, num_procs);
            break; 
        }
        else if (iter > MAX_ITER)
        {
            // Loosen the convergence criteria a bit
            if (rTr < criteria * 100.0)
            {
                if (ROOT_PROC == my_id)
                {
                    printf("Note: Loosen bicg_stab convergence criteria to 1e-6 for this step.\n");
                }

                // converged, communicate vec_pq and quit the iteration
                mpi_comm_vec(start_metal, end_metal, min_metal, max_metal, n_NPs,
                             vec_pq, my_id, num_procs);
                break;
            }
            else
            {
                // not converged, terminate the program
                if (ROOT_PROC == my_id)
                {
                    printf("Error: bicg_stab not converged. rTr= %e, bTb= %e\n", rTr, bTb);
                }

                exit(-1);
            }
        }
    }
}