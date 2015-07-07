/*****************************************************************************
 *  This file is part of the CapacMD program.                                *  
 *  Copyright (C) 2015 Xin Li <lixin.reco@gmail.com>                         *  
 *                                                                           *
 *  Filename:  force.c                                                       *
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

#include <mpi.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>

#include "typedef.h"
#include "my_malloc.h"
#include "cpff.h"
#include "bicg_stab.h"

//================
// norm of vector
//================

double dist_2(const Vec_R* ptr_r)
{
    return ptr_r->x * ptr_r->x + ptr_r->y * ptr_r->y + ptr_r->z * ptr_r->z;
}

//===================
// fvec = fij * rvec
//===================

void scale_vec_1(const double ss, Vec_R* ptr_rvec)
{
	ptr_rvec->x *= ss;
	ptr_rvec->y *= ss;
	ptr_rvec->z *= ss;
}

void scale_vec_2(const double fij, const Vec_R* ptr_rvec, Vec_R* ptr_fvec)
{
	ptr_fvec->x = fij * ptr_rvec->x;
	ptr_fvec->y = fij * ptr_rvec->y;
	ptr_fvec->z = fij * ptr_rvec->z;
}

//==================================
// update virial from rvec and fvec
//==================================

void virial_rvec_fvec_half(const Vec_R* ptr_rvec, const Vec_R* ptr_fvec, 
						   double** virial)
{
    virial[0][0] -= ptr_rvec->x * ptr_fvec->x * 0.5;
    virial[0][1] -= ptr_rvec->x * ptr_fvec->y * 0.5;
    virial[0][2] -= ptr_rvec->x * ptr_fvec->z * 0.5;

    virial[1][0] -= ptr_rvec->y * ptr_fvec->x * 0.5;
    virial[1][1] -= ptr_rvec->y * ptr_fvec->y * 0.5;
    virial[1][2] -= ptr_rvec->y * ptr_fvec->z * 0.5;

    virial[2][0] -= ptr_rvec->z * ptr_fvec->x * 0.5;
    virial[2][1] -= ptr_rvec->z * ptr_fvec->y * 0.5;
    virial[2][2] -= ptr_rvec->z * ptr_fvec->z * 0.5;
}

void virial_rvec_fvec_full(const Vec_R* ptr_rvec, const Vec_R* ptr_fvec, 
						   double** virial)
{
    virial[0][0] -= ptr_rvec->x * ptr_fvec->x;
    virial[0][1] -= ptr_rvec->x * ptr_fvec->y;
    virial[0][2] -= ptr_rvec->x * ptr_fvec->z;

    virial[1][0] -= ptr_rvec->y * ptr_fvec->x;
    virial[1][1] -= ptr_rvec->y * ptr_fvec->y;
    virial[1][2] -= ptr_rvec->y * ptr_fvec->z;

    virial[2][0] -= ptr_rvec->z * ptr_fvec->x;
    virial[2][1] -= ptr_rvec->z * ptr_fvec->y;
    virial[2][2] -= ptr_rvec->z * ptr_fvec->z;
}

//======================
// apply PBC to delta-x
//======================

double pbc(double dx, double xbox)
{
    if     (dx < -0.5*xbox) {dx += xbox;}
    else if(dx >  0.5*xbox) {dx -= xbox;}
    return dx;
}

void pbc_12(Vec_R* ptr_r, const double* box)
{
    if     (ptr_r->x < -box[3]) { ptr_r->x += box[0]; }
    else if(ptr_r->x >  box[3]) { ptr_r->x -= box[0]; }
                                                      
    if     (ptr_r->y < -box[4]) { ptr_r->y += box[1]; }
    else if(ptr_r->y >  box[4]) { ptr_r->y -= box[1]; }
                                                      
    if     (ptr_r->z < -box[5]) { ptr_r->z += box[2]; }
    else if(ptr_r->z >  box[5]) { ptr_r->z -= box[2]; }
}

//===============================
// apply PBC to the whole system
//===============================

void apply_pbc(const int nMols, const Mol_Info* mol_info, 
               double* rx, double* ry, double* rz, double* box)
{
    double half_box[3];
    half_box[0] = 0.5 * box[0];
    half_box[1] = 0.5 * box[1];
    half_box[2] = 0.5 * box[2];

    // loop over molecules
    int im;
    for(im = 0; im < nMols; ++ im)
    {
        // put the center of molecule in the box
        int first_atom = mol_info[im].mini;
        int last_atom  = mol_info[im].maxi;
        const int mid = (first_atom + last_atom) / 2;

        if     (rx[mid] > box[0]) { rx[mid] -= box[0]; }
        else if(rx[mid] < 0.0   ) { rx[mid] += box[0]; }
        if     (ry[mid] > box[1]) { ry[mid] -= box[1]; }
        else if(ry[mid] < 0.0   ) { ry[mid] += box[1]; }
        if     (rz[mid] > box[2]) { rz[mid] -= box[2]; }
        else if(rz[mid] < 0.0   ) { rz[mid] += box[2]; }

        // apply PBC to other atoms in this molecule
        int i;
        for(i = first_atom; i <= last_atom; ++ i)
        {
            if(mid == i) { continue; }

            if     (rx[i] - rx[mid] >  half_box[0]) { rx[i] -= box[0]; }
            else if(rx[i] - rx[mid] < -half_box[0]) { rx[i] += box[0]; }
            if     (ry[i] - ry[mid] >  half_box[1]) { ry[i] -= box[1]; }
            else if(ry[i] - ry[mid] < -half_box[1]) { ry[i] += box[1]; }
            if     (rz[i] - rz[mid] >  half_box[2]) { rz[i] -= box[2]; }
            else if(rz[i] - rz[mid] < -half_box[2]) { rz[i] += box[2]; }
        }
    }
}


//========================================
// find starting and ending atom/molecule
// on each processor
//========================================

void find_start_end(int* start_group, int* end_group, 
                    const int total_num, const int num_procs)
{
    int avg_num_per_process = total_num / num_procs;
    int residual = total_num % num_procs;

    int an_id;
    for(an_id = 0; an_id < num_procs; ++ an_id)
    {
        int num_in_group;

        // processors with smaller id than the residual will carry one more element
        if(an_id < residual) { num_in_group = avg_num_per_process + 1;}
        else { num_in_group = avg_num_per_process; }

        if(0 == an_id) { start_group[0] = 0; }
        else { start_group[an_id] = end_group[an_id - 1] + 1; }

        end_group[an_id] = start_group[an_id] + num_in_group - 1;
    }
}

//=======================================
// find starting and ending atomic pairs
// on each processor
//=======================================

void find_start_end_long(long int* start_group, long int* end_group,
                         const long int total_num, const int num_procs)
{
    long int avg_num_per_process = total_num / num_procs;
    int residual = total_num % num_procs;

    int an_id;
    for(an_id = 0; an_id < num_procs; ++ an_id)
    {
        long int num_in_group;

        // processors with smaller id than the residual will carry one more element
        if(an_id < residual) { num_in_group = avg_num_per_process + 1; }
        else { num_in_group = avg_num_per_process; }

        if(0 == an_id) { start_group[0] = 0; }
        else { start_group[an_id] = end_group[an_id - 1] + 1; }

        end_group[an_id] = start_group[an_id] + num_in_group - 1;
    }
}


//===========================
// sum and analyze time_used
//===========================

void sum_time_used(double** time_used, const int my_id, const int num_procs)
{
    // gather time_used to root processor
    int ierr;

    const int root_process = 0;
    MPI_Status status;
    int tag_101 = 101;

    if(root_process == my_id)
    {
        int an_id;
        for(an_id = 1; an_id < num_procs; an_id++) 
        {
            ierr = MPI_Recv(&(time_used[an_id][0]), 15, MPI_DOUBLE, an_id, tag_101, 
                            MPI_COMM_WORLD, &status);
        }

        // sum up time usage
        for(an_id = 0; an_id < num_procs; ++ an_id)
        {
            time_used[an_id][0] = 0.0;

            int it;
            for(it = 1; it < 15; ++ it)
            {
                time_used[an_id][0] += time_used[an_id][it];
            }
        }
    }
    else
    {
        ierr = MPI_Send(&(time_used[my_id][0]), 15, MPI_DOUBLE, root_process, tag_101, 
                        MPI_COMM_WORLD);
    }

    if(ierr) {}
}


void analyze_time_used(double** time_used, const int num_procs)
{
    double aver_time[15];
    double max_time[15];
    double imbalance[15];

    int an_id, it;
    for(it = 0; it < 15; ++ it)
    {
        aver_time[it] = 0.0;
        max_time[it]  = 0.0;
        imbalance[it]  = 0.0;

        for(an_id = 0; an_id < num_procs; ++ an_id)
        {
            aver_time[it] += time_used[an_id][it];

            if(time_used[an_id][it] > max_time[it])
            {
                max_time[it] = time_used[an_id][it];
            }
        }

        aver_time[it] /= num_procs;

		// skip evaluation of imbalance for very short time usage
		if(aver_time[it] > 0.1)
		{
        	imbalance[it] = (max_time[it] / aver_time[it] - 1.0) * 100.0;
		}
    }

    printf ("    Total time used in force: %.3f seconds  Imb. %.1f%%\n", aver_time[0], imbalance[0]);
    printf ("    Time used in QSC density: %.3f seconds  Imb. %.1f%%\n", aver_time[1], imbalance[1]);
    printf ("    Time used in QSC force:   %.3f seconds  Imb. %.1f%%\n", aver_time[2], imbalance[2]);
    printf ("    Time used in Bonded:      %.3f seconds  Imb. %.1f%%\n", aver_time[3], imbalance[3]);
    printf ("    Time used in Nonbonded:   %.3f seconds  Imb. %.1f%%\n", aver_time[4], imbalance[4]);
    printf ("    Time used in Coulomb_LR:  %.3f seconds  Imb. %.1f%%\n", aver_time[9], imbalance[9]);
    printf ("    Time used in CPIM-matrix: %.3f seconds  Imb. %.1f%%\n", aver_time[5], imbalance[5]);
    printf ("    Time used in CPIM-vector: %.3f seconds  Imb. %.1f%%\n", aver_time[6], imbalance[6]);
    printf ("    Time used in CPIM-solve:  %.3f seconds  Imb. %.1f%%\n", aver_time[7], imbalance[7]);
    printf ("    Time used in CPIM-force:  %.3f seconds  Imb. %.1f%%\n", aver_time[8], imbalance[8]);
    printf ("\n");
    printf ("    Time used in rx communication:  %.3f seconds  Imb. %.1f%%\n", 
                            aver_time[10], imbalance[10]);
    printf ("    Time used in fx communication:  %.3f seconds  Imb. %.1f%%\n", 
                            aver_time[11], imbalance[11]);
    printf ("\n");
}


//========================
// sum potential energies
//========================

void sum_potential(double *potential)
{
    potential[0] = 0.0;

    int i;
    for(i = 1; i < 15; i++)
    {
        potential[0] += potential[i];
    }
}

//==========================
// print potential energies
//==========================

void print_potential(double* potential)
{
    printf("    Potential energies: \n");
    printf(" %15.6e   total energy                         \n", potential[0]);
    printf(" %15.6e   metal quantum Sutton-Chen energy     \n", potential[1]);
    printf(" %15.6e   non-metal bond stretching energy     \n", potential[2]);
    printf(" %15.6e   non-metal angle bending energy       \n", potential[3]);
    printf(" %15.6e   non-metal torsional energy           \n", potential[4]);
    printf(" %15.6e   Long range Coulomb                   \n", potential[6]);
    printf(" %15.6e   Coulomb energy (including 1-4)       \n", potential[7]);
    printf(" %15.6e   vdW energy (including 1-4)           \n", potential[8]);
    printf(" %15.6e   CPIM metal charge - non-metal charge \n", potential[10]);
    printf(" %15.6e   CPIM metal dipole - non-metal charge \n", potential[11]);
    printf(" %15.6e   CPIM metal charge - metal charge     \n", potential[12]);
    printf(" %15.6e   CPIM metal charge - metal dipole     \n", potential[13]);
    printf(" %15.6e   CPIM metal dipole - metal dipole     \n", potential[14]);
    printf("\n");
}


//=================================================
// get maximal force and root mean square of force
//=================================================

void get_fmax_rms(const int nAtoms, System *p_system)
{
    double fmax2 = 0.0;
    double f_rms = 0.0;
    int i;
    for (i = 0; i < nAtoms; ++ i)
    {
        double f2 = p_system->fx[i] * p_system->fx[i] + 
					p_system->fy[i] * p_system->fy[i] + 
					p_system->fz[i] * p_system->fz[i];

        f_rms += f2;

        if (f2 > fmax2) { fmax2 = f2; }
    }
    p_system->f_max = sqrt(fmax2);
    p_system->f_rms = sqrt(f_rms / nAtoms);
}


//===========================================
// Quantum Sutton-Chen density, parallelized
//===========================================

void mpi_qsc_dens(Task *p_task, const double rCut2, Metal *p_metal, 
				  System *p_system, const int my_id, const int num_procs)
{
	const int start_metal = p_task->start_metal[my_id];
	const int end_metal   = p_task->end_metal[my_id];

	Vec_R rr;
	double rij2;
    int i, j;

    // compute densities and get inv_sqrt of the densities
    for(i = start_metal; i <= end_metal; ++ i)
    {
        int i_metal = i - p_metal->min;
        double qsc_dens = 0.0;

        for(j = p_metal->min; j <= p_metal->max; ++ j)
        {
            if (j == i) { continue; }


			rr.x = p_system->rx[j] - p_system->rx[i];
			rr.y = p_system->ry[j] - p_system->ry[i];
			rr.z = p_system->rz[j] - p_system->rz[i];
			pbc_12(&rr, p_system->box);
			rij2 = dist_2(&rr);


            // using cutoff for Q-SC potential
            if(rij2 < rCut2)
            {
                qsc_dens += pow((p_metal->qsc_a / sqrt(rij2)), p_metal->qsc_m);
            }
        }

        p_metal->inv_sqrt_dens[i_metal] = 1.0 / sqrt(qsc_dens);
    }

    // communicate inv_sqrt_dens among all processors
    int an_id, ierr;
    MPI_Status status;

    int tag_21 = 21; // any to any, start
    int tag_22 = 22; // any to any, num
    int tag_23 = 23; // any to any, inv_sqrt_dens
    
    const int root_process = 0;

    int start, num;

    // master: receive data
    if (root_process == my_id)
    {
            for (an_id = 1; an_id < num_procs; ++ an_id)
            {
                    ierr = MPI_Recv(&start, 1, MPI_INT, an_id, tag_21, MPI_COMM_WORLD, &status);
                    ierr = MPI_Recv(&num,   1, MPI_INT, an_id, tag_22, MPI_COMM_WORLD, &status);
                    ierr = MPI_Recv(&(p_metal->inv_sqrt_dens[start]), num, MPI_DOUBLE, an_id, tag_23, 
                                                    MPI_COMM_WORLD, &status);
            }
    }
    // slave: send data
    else
    {
    start = start_metal - p_metal->min;
    num = end_metal - start_metal + 1;

     ierr = MPI_Send(&start, 1, MPI_INT, root_process, tag_21, MPI_COMM_WORLD);
     ierr = MPI_Send(&num,   1, MPI_INT, root_process, tag_22, MPI_COMM_WORLD);
     ierr = MPI_Send(&(p_metal->inv_sqrt_dens[start]), num, MPI_DOUBLE, root_process, tag_23, 
                                            MPI_COMM_WORLD);
    }

    MPI_Bcast(&(p_metal->inv_sqrt_dens[0]), p_metal->num, MPI_DOUBLE, root_process, MPI_COMM_WORLD);

    if(ierr) {}
}


//=========================================
// Quantum Sutton-Chen force, parallelized
//=========================================

void mpi_qsc_force(Task *p_task, const double rCut2, Metal* p_metal, 
				   System *p_system, const int my_id)
{
	const int start_metal = p_task->start_metal[my_id];
	const int end_metal   = p_task->end_metal[my_id];

    int i_atom, j_atom;
    double rij2, rij; 
    double cm_2, a_rij, a_rij_n, a_rij_m, sum_inv_sqrt, fij;

	Vec_R rvec, fvec;

    // compute forces
    cm_2 = p_metal->qsc_c * p_metal->qsc_m / 2.0;

    for(i_atom = start_metal; i_atom <= end_metal; ++ i_atom)
    {
        int i_metal = i_atom - p_metal->min;

        for(j_atom = p_metal->min; j_atom <= p_metal->max; ++ j_atom)
        {
            if (j_atom == i_atom) { continue; }

            int j_metal = j_atom - p_metal->min;


			rvec.x = p_system->rx[j_atom] - p_system->rx[i_atom];
			rvec.y = p_system->ry[j_atom] - p_system->ry[i_atom];
			rvec.z = p_system->rz[j_atom] - p_system->rz[i_atom];
			pbc_12(&rvec, p_system->box);
			rij2 = dist_2(&rvec);


            //===================================================================
            // Force calculation of the Q-SC potential
            // http://www.ens-lyon.fr/DSM/SDMsite/M2/stages_M2/Forster2012.pdf
            // see also: Surface Science 281 (1993) 191-206
            //===================================================================

            // using cutoff for Q-SC potential
            if(rij2 < rCut2)
            {
                rij = sqrt(rij2);
                a_rij = p_metal->qsc_a / rij,
                a_rij_n = pow(a_rij, p_metal->qsc_n);
                a_rij_m = pow(a_rij, p_metal->qsc_m);

                sum_inv_sqrt = p_metal->inv_sqrt_dens[i_metal] + p_metal->inv_sqrt_dens[j_metal];

                fij = p_metal->qsc_eps * 
                      (p_metal->qsc_n * a_rij_n - cm_2 * sum_inv_sqrt * a_rij_m) / rij2;


				// fvec = fij * rvec;
				scale_vec_2(fij, &rvec, &fvec);


                // Note: not using Newton's third law
                // update forces on i_atom only
                p_system->fx[i_atom] -= fvec.x;
                p_system->fy[i_atom] -= fvec.y;
                p_system->fz[i_atom] -= fvec.z;

                // virial increment scaled by 0.5
				virial_rvec_fvec_half(&rvec, &fvec, p_system->virial);


                // potential[1] = metal qsc energy
                p_system->potential[1] += p_metal->qsc_eps * a_rij_n * 0.5;
            }
        }

        p_system->potential[1] -= p_metal->qsc_eps * p_metal->qsc_c / p_metal->inv_sqrt_dens[i_metal];
    }
}


void zero_vec_nb(Vec_nb* ptr_nonbonded)
{
    ptr_nonbonded->lj.fi.x = 0.0;
    ptr_nonbonded->lj.fi.y = 0.0;
    ptr_nonbonded->lj.fi.z = 0.0;
    ptr_nonbonded->qq.fi.x = 0.0;
    ptr_nonbonded->qq.fi.y = 0.0;
    ptr_nonbonded->qq.fi.z = 0.0;

    ptr_nonbonded->lj.fj.x = 0.0;
    ptr_nonbonded->lj.fj.y = 0.0;
    ptr_nonbonded->lj.fj.z = 0.0;
    ptr_nonbonded->qq.fj.x = 0.0;
    ptr_nonbonded->qq.fj.y = 0.0;
    ptr_nonbonded->qq.fj.z = 0.0;
}

void virial_vec_nb(const double rxij, const double ryij, const double rzij,
					const Vec_nb* ptr_nonbonded, double** virial)
{
	double fxij, fyij, fzij;

	fxij = ptr_nonbonded->lj.fj.x + ptr_nonbonded->qq.fj.x;
	fyij = ptr_nonbonded->lj.fj.y + ptr_nonbonded->qq.fj.y;
	fzij = ptr_nonbonded->lj.fj.z + ptr_nonbonded->qq.fj.z;

    // update virial
    virial[0][0] -= rxij * fxij;
    virial[0][1] -= rxij * fyij;
    virial[0][2] -= rxij * fzij;

    virial[1][0] -= ryij * fxij;
    virial[1][1] -= ryij * fyij;
    virial[1][2] -= ryij * fzij;

    virial[2][0] -= rzij * fxij;
    virial[2][1] -= rzij * fyij;
    virial[2][2] -= rzij * fzij;
}


//=============================================
// Erf-vdW potential
// note: no 1-4 scaling factor, setting to 1.0
//=============================================

void compute_erf_vdw(double rxij, double ryij, double rzij, double rij2,
                     double c6, double c12, double width,
                     double* potential, double** virial, Vec_nb* ptr_nonbonded)
{
	double rij6, rij12, cr12, cr6, rR, erf_rR, rij, vij, fij;
    rij = sqrt(rij2);

    // zero nonbonded force for this i-j pair
	zero_vec_nb(ptr_nonbonded);

    // vdW contribution
    if (c6!=0.0 || c12!=0.0)
    {
        rij6  = rij2 * rij2 * rij2;
        rij12 = rij6 * rij6;

        cr12 = c12 / rij12;
        cr6  = c6  / rij6;

        rR = rij / width;
        erf_rR = erf(rR);

        vij = (cr12 - cr6) * erf_rR;
        fij = (12.0 * cr12 - 6.0 * cr6) / rij * erf_rR -
              (cr12 - cr6) * 2.0 * INV_SQRT_PI * exp(-rR * rR) / width;
        fij /= rij;

        // potential[8] = vdW energy
        potential[8] += vij;

        ptr_nonbonded->lj.fj.x = fij * rxij;
        ptr_nonbonded->lj.fj.y = fij * ryij;
        ptr_nonbonded->lj.fj.z = fij * rzij;

        ptr_nonbonded->lj.fi.x = -ptr_nonbonded->lj.fj.x;
        ptr_nonbonded->lj.fi.y = -ptr_nonbonded->lj.fj.y;
        ptr_nonbonded->lj.fi.z = -ptr_nonbonded->lj.fj.z;
    }

    // update virial
	virial_vec_nb(rxij, ryij, rzij, ptr_nonbonded, virial);
}


//=================
// Morse potential
//=================

void compute_morse(double rxij, double ryij, double rzij, double rij2,
                   double morse_D, double morse_a, double morse_R,
                   double* potential, double** virial, Vec_nb* ptr_nonbonded)
{
	double fij, rij, exp__ar, exp__2ar;
    rij  = sqrt(rij2);

	zero_vec_nb(ptr_nonbonded);

    // Morse potential
    exp__ar = exp(-morse_a * (rij - morse_R));
    exp__2ar = exp__ar * exp__ar;

    fij = 2.0 * morse_D * morse_a * (exp__2ar - exp__ar) / rij;

    // potential[8] = vdW energy
    potential[8] += morse_D * (exp__2ar - 2.0 * exp__ar);

    ptr_nonbonded->lj.fj.x = fij * rxij;
    ptr_nonbonded->lj.fj.y = fij * ryij;
    ptr_nonbonded->lj.fj.z = fij * rzij;

    ptr_nonbonded->lj.fi.x = -ptr_nonbonded->lj.fj.x;
    ptr_nonbonded->lj.fi.y = -ptr_nonbonded->lj.fj.y;
    ptr_nonbonded->lj.fi.z = -ptr_nonbonded->lj.fj.z;

    // update virial
	virial_vec_nb(rxij, ryij, rzij, ptr_nonbonded, virial);
}


//======================
// Buckingham potential
//======================

void compute_buckingham(double rxij, double ryij, double rzij, double rij2,
                        double buck_A, double buck_B, double c6,
                        double* potential, double** virial, Vec_nb* ptr_nonbonded)
{
	double rij6, rij8, rij, exp__Br, fij;
    rij  = sqrt(rij2);

	zero_vec_nb(ptr_nonbonded);

    // Buckingham potential
    rij6 = rij2 * rij2 * rij2;
    rij8 = rij6 * rij2;

    exp__Br = exp(-buck_B * rij);
    fij = buck_A * exp__Br * buck_B / rij - 6.0 * c6 / rij8;

    // potential[8] = vdW energy
    potential[8] += (buck_A * exp__Br - c6 / rij6);

    ptr_nonbonded->lj.fj.x = fij * rxij;
    ptr_nonbonded->lj.fj.y = fij * ryij;
    ptr_nonbonded->lj.fj.z = fij * rzij;

    ptr_nonbonded->lj.fi.x = -ptr_nonbonded->lj.fj.x;
    ptr_nonbonded->lj.fi.y = -ptr_nonbonded->lj.fj.y;
    ptr_nonbonded->lj.fi.z = -ptr_nonbonded->lj.fj.z;

    // update virial
	virial_vec_nb(rxij, ryij, rzij, ptr_nonbonded, virial);
}


//==========================================
// Non-bonded force between a pair of atoms
//==========================================

void compute_nonbonded(double rxij, double ryij, double rzij, double rij2,
					   RunSet* p_runset, double scaleLJ, double scaleQQ,
					   int is_14_pair,
                       double c6, double c12, double qi, double qj,
                       double* potential, double** virial, Vec_nb* ptr_nonbonded)
{
	double rij6, rij12, rij, fij;
	double f_r, f_c, erfc_arij;
    rij = sqrt(rij2);
	fij = 0.0; // initialize fij

	double rCut, inv_rc12, inv_rc6;
	rCut     = p_runset->rCut;
	inv_rc12 = p_runset->inv_rc12;
	inv_rc6  = p_runset->inv_rc6;

	double alpha, a2_sqrtPI, wolfConst, erfc_arCut;
    alpha      = p_runset->w_alpha;
    a2_sqrtPI  = p_runset->w_a2_sqrtPI;
    wolfConst  = p_runset->w_Const;
    erfc_arCut = p_runset->w_erfc_arCut;

	zero_vec_nb(ptr_nonbonded);

    // vdW contribution
    if (c6!=0.0 || c12!=0.0)
    {
        rij6 = rij2 * rij2 * rij2;
        rij12 = rij6 * rij6;

        if(1 == p_runset->use_vdw && 0 == is_14_pair)
        {
            // The shifted force method for Lennard-Jones potential
            // Ref: Toxvaerd and Dyre, J. Chem. Phys. 2011, 134, 081102
            //      dx.doi.org/10.1063/1.3558787
			
			// note: no scaleLJ for non-14 pair!

            f_r = 12.0 * c12 / rij12 - 6.0 * c6 / rij6;
            f_c = 12.0 * c12 * inv_rc12 - 6.0 * c6 * inv_rc6;

            fij = (f_r - f_c) / rij2;

            // potential[8] = vdW energy
            potential[8] += ((c12 / rij12 - c6 / rij6) +
                             (rij - rCut) * f_c / rCut -
                             (c12 * inv_rc12 - c6 * inv_rc6));
        }

        else if(0 == p_runset->use_vdw || 1 == is_14_pair)
        {
            fij = (c12 * 2.0 / rij12 - c6 / rij6) * 6.0 / rij2 * scaleLJ;

            potential[8] += (c12 / rij12 - c6 / rij6) * scaleLJ;
        }

        ptr_nonbonded->lj.fj.x = fij * rxij;
        ptr_nonbonded->lj.fj.y = fij * ryij;
        ptr_nonbonded->lj.fj.z = fij * rzij;
 
        ptr_nonbonded->lj.fi.x = -ptr_nonbonded->lj.fj.x;
        ptr_nonbonded->lj.fi.y = -ptr_nonbonded->lj.fj.y;
        ptr_nonbonded->lj.fi.z = -ptr_nonbonded->lj.fj.z;
    }

    // Coulomb contribution
    if (qi!=0.0 && qj!=0.0)
    {

        if(1 == p_runset->use_coulomb && 0 == is_14_pair)
        {
            // The damped shifted force (DSF) approach for electrostatics
            // Ref.: Fennell and Gezelter, J. Chem. Phys. 2006, 124, 234104
            //       dx.doi.org/10.1063/1.2206581

			// note: no scaleQQ for non-14 pair!
            
			erfc_arij = erfc(alpha * rij);

            fij = FQQ * qi * qj * 
                  ((erfc_arij / rij2 + 
                   a2_sqrtPI * exp(-alpha * alpha * rij2) / rij) - wolfConst) / rij;

            // potential[7] = Coulomb energy
            potential[7] += FQQ * qi * qj *
                            (erfc_arij / rij - erfc_arCut + wolfConst * (rij - rCut));
        }

        else if(0 == p_runset->use_coulomb || 1 == is_14_pair)
        {
            fij = scaleQQ * FQQ * qi * qj / (rij2 * rij);

            potential[7] += scaleQQ * FQQ * qi * qj / rij;
        }

        ptr_nonbonded->qq.fj.x = fij * rxij;
        ptr_nonbonded->qq.fj.y = fij * ryij;
        ptr_nonbonded->qq.fj.z = fij * rzij;
 
        ptr_nonbonded->qq.fi.x = -ptr_nonbonded->qq.fj.x;
        ptr_nonbonded->qq.fi.y = -ptr_nonbonded->qq.fj.y;
        ptr_nonbonded->qq.fi.z = -ptr_nonbonded->qq.fj.z;
    }

    // update virial
	virial_vec_nb(rxij, ryij, rzij, ptr_nonbonded, virial);
}

//==============================================
// bond stretching force, harmonic
//==============================================

Vec_2 compute_bond_1(double rxi, double ryi, double rzi, 
                     double rxj, double ryj, double rzj,
                     double b0, double kb, 
                     double* box, double* potential, double** virial)
{
    Vec_2 bond_force;
    double rxij, ryij, rzij, rij, fij;

    rxij = pbc(rxj - rxi, box[0]);
    ryij = pbc(ryj - ryi, box[1]);
    rzij = pbc(rzj - rzi, box[2]);
    rij  = sqrt(rxij*rxij + ryij*ryij + rzij*rzij);

    // Eij = 1/2 kb (rij-r0)^2
    // fij = -dEij/drij = -kb(rij-r0) = kb(r0-rij)
    // Vector_fij = fij * Vector_rij / rij

    fij  = kb * ( b0 - rij ) / rij;
    bond_force.fj.x = fij * rxij;
    bond_force.fj.y = fij * ryij;
    bond_force.fj.z = fij * rzij;
    bond_force.fi.x = -bond_force.fj.x;
    bond_force.fi.y = -bond_force.fj.y;
    bond_force.fi.z = -bond_force.fj.z;

    // potential[2] = non-metal bond stretching energy
    potential[2] += 0.5 * kb * (b0 - rij) * (b0 - rij);

    virial[0][0] += -rxij * bond_force.fj.x;
    virial[0][1] += -rxij * bond_force.fj.y;
    virial[0][2] += -rxij * bond_force.fj.z;
    virial[1][0] += -ryij * bond_force.fj.x;
    virial[1][1] += -ryij * bond_force.fj.y;
    virial[1][2] += -ryij * bond_force.fj.z;
    virial[2][0] += -rzij * bond_force.fj.x;
    virial[2][1] += -rzij * bond_force.fj.y;
    virial[2][2] += -rzij * bond_force.fj.z;

    return bond_force;
}

//==============================================
// bond stretching force, Morse potential
//==============================================

Vec_2 compute_bond_3(double rxi, double ryi, double rzi,
                     double rxj, double ryj, double rzj,
                     double b0, double D, double beta, 
                     double* box, double* potential, double** virial)
{
    Vec_2 bond_force;
    double rxij, ryij, rzij, rij, fij;

    rxij = pbc(rxj - rxi, box[0]);
    ryij = pbc(ryj - ryi, box[1]);
    rzij = pbc(rzj - rzi, box[2]);
    rij  = sqrt(rxij*rxij + ryij*ryij + rzij*rzij);

    // Eij = D*[1.0-exp(-beta*(rij-r0))]^2
    // fij = -dEij/drij = -D * 2.0*[1.0-exp(-beta*(rij-r0))] *
    //                    exp(-beta*(rij-r0)) * beta
    // Vector_fij = fij * Vector_rij / rij

    double exp_b_rij = exp(-beta * (rij - b0));

    fij = -2.0 * D * beta * exp_b_rij * ( 1.0 - exp_b_rij ) / rij;
    bond_force.fj.x = fij * rxij;
    bond_force.fj.y = fij * ryij;
    bond_force.fj.z = fij * rzij;
    bond_force.fi.x = -bond_force.fj.x;
    bond_force.fi.y = -bond_force.fj.y;
    bond_force.fi.z = -bond_force.fj.z;

    // potential[2] = non-metal bond stretching energy
    potential[2] += D * pow(1.0 - exp_b_rij, 2.0);

    virial[0][0] += -rxij * bond_force.fj.x;
    virial[0][1] += -rxij * bond_force.fj.y;
    virial[0][2] += -rxij * bond_force.fj.z;
    virial[1][0] += -ryij * bond_force.fj.x;
    virial[1][1] += -ryij * bond_force.fj.y;
    virial[1][2] += -ryij * bond_force.fj.z;
    virial[2][0] += -rzij * bond_force.fj.x;
    virial[2][1] += -rzij * bond_force.fj.y;
    virial[2][2] += -rzij * bond_force.fj.z;

    return bond_force;
}

//==============================================
// product between scalar and vector
//==============================================

Vec_R scalar_x_vector ( double ss, Vec_R vv )
{
    Vec_R sv;
    sv.x = ss * vv.x;
    sv.y = ss * vv.y;
    sv.z = ss * vv.z;
    return sv;
}

//==============================================
// inner product between vectors
//==============================================

double inner_product ( Vec_R ri, Vec_R rj )
{
    double inner_p;
    inner_p = ri.x*rj.x + ri.y*rj.y + ri.z*rj.z;
    return inner_p;
}

//==============================================
// outer product between vectors
// http://en.wikipedia.org/wiki/Cross_product
//==============================================

Vec_R outer_product ( Vec_R ri, Vec_R rj )
{
    Vec_R outer_p;

    outer_p.x = ri.y*rj.z - ri.z*rj.y;
    outer_p.y = ri.z*rj.x - ri.x*rj.z;
    outer_p.z = ri.x*rj.y - ri.y*rj.x;
    return outer_p;
}

//==================================================================
// angle bending force, harmonic
// J. Comput. Chem. 2000, 21, 553-561
// doi: 10.1002/(SICI)1096-987X(200005)21:7<553::AID-JCC4>3.0.CO;2-1
//==================================================================

Vec_3 compute_angle_1(double rxi, double ryi, double rzi,
                      double rxj, double ryj, double rzj,
                      double rxk, double ryk, double rzk,
                      double a0, double ka,
                      double *box, double *potential, double **virial)
{
    Vec_3 angle_force;
    Vec_R r21, r21_unit, r32, r32_unit;
    double r21_d, r32_d, cos123, theta123, sin123, fij;

    r21.x = pbc(rxi - rxj, box[0]);
    r21.y = pbc(ryi - ryj, box[1]);
    r21.z = pbc(rzi - rzj, box[2]);
    r21_d = sqrt(dist_2(&r21));
    r21_unit = scalar_x_vector(1.0 / r21_d, r21);

    r32.x = pbc(rxj - rxk, box[0]);
    r32.y = pbc(ryj - ryk, box[1]);
    r32.z = pbc(rzj - rzk, box[2]);
    r32_d = sqrt(dist_2(&r32));
    r32_unit = scalar_x_vector(1.0 / r32_d, r32);

    cos123 = -inner_product(r21_unit, r32_unit);
    theta123 = acos(cos123);
    sin123 = sin(theta123);

    double delta_a = a0 * M_PI / 180.0 - theta123; // in radian
    fij = ka * delta_a;
    angle_force.fi.x =  (cos123 * r21_unit.x + r32_unit.x) / (sin123 * r21_d) * fij;
    angle_force.fi.y =  (cos123 * r21_unit.y + r32_unit.y) / (sin123 * r21_d) * fij;
    angle_force.fi.z =  (cos123 * r21_unit.z + r32_unit.z) / (sin123 * r21_d) * fij;
    angle_force.fk.x = -(cos123 * r32_unit.x + r21_unit.x) / (sin123 * r32_d) * fij;
    angle_force.fk.y = -(cos123 * r32_unit.y + r21_unit.y) / (sin123 * r32_d) * fij;
    angle_force.fk.z = -(cos123 * r32_unit.z + r21_unit.z) / (sin123 * r32_d) * fij;
    angle_force.fj.x = -angle_force.fi.x - angle_force.fk.x;
    angle_force.fj.y = -angle_force.fi.y - angle_force.fk.y;
    angle_force.fj.z = -angle_force.fi.z - angle_force.fk.z;

    // potential[3] = non-metal angle bending energy
    potential[3] += 0.5 * ka * delta_a * delta_a;

    // scalar viral from angle-dependent potential is zero
	// but we add individual contributions anyways
	// ftp://ftp.dl.ac.uk/ccp5.newsletter/26/pdf/smith.pdf

	Vec_R rvec;

	rvec.x = r21.x;
	rvec.y = r21.y;
	rvec.z = r21.z;
	virial_rvec_fvec_full(&rvec, &(angle_force.fi), virial);

	rvec.x = -r32.x;
	rvec.y = -r32.y;
	rvec.z = -r32.z;
	virial_rvec_fvec_full(&rvec, &(angle_force.fk), virial);


    return angle_force;
}

//==============================================
// angle bending force, Urey-Bradley
//==============================================

Vec_3 compute_angle_5(double rxi, double ryi, double rzi,
                      double rxj, double ryj, double rzj,
                      double rxk, double ryk, double rzk,
                      double a0, double ka,
                      double b0_ub, double kb_ub,
                      double* box, double* potential, double** virial)
{
    Vec_3 angle_force;
    angle_force = compute_angle_1(rxi, ryi, rzi,
                                  rxj, ryj, rzj,
                                  rxk, ryk, rzk,
                                  a0, ka,
                                  box, potential, virial);

    if(kb_ub > 0.0)
    {
        Vec_2 bond_force;
        bond_force = compute_bond_1(rxi, ryi, rzi,
                                    rxk, ryk, rzk,
                                    b0_ub, kb_ub, 
                                    box, potential, virial);

        angle_force.fi.x += bond_force.fi.x;
        angle_force.fi.y += bond_force.fi.y;
        angle_force.fi.z += bond_force.fi.z;
        angle_force.fk.x += bond_force.fj.x;
        angle_force.fk.y += bond_force.fj.y;
        angle_force.fk.z += bond_force.fj.z;
    }

    return angle_force;
}


//=========================================================
// distribute forces from vitual sites to real atoms
// see GROMACS manual
//=========================================================

Vec_4 compute_vsite_4(double rxi, double ryi, double rzi,
                      double rxj, double ryj, double rzj,
                      double rxk, double ryk, double rzk,
                      double rxs, double rys, double rzs,
                      double fxs, double fys, double fzs,
                      double a, double b, double c, double* box)
{
    double xik, yik, zik, xij, yij, zij;
    Vec_4  vsite_force;

    xij = pbc(rxj - rxi, box[0]);
    yij = pbc(ryj - ryi, box[1]);
    zij = pbc(rzj - rzi, box[2]);

    xik = pbc(rxk - rxi, box[0]);
    yik = pbc(ryk - ryi, box[1]);
    zik = pbc(rzk - rzi, box[2]);

    vsite_force.fj.x =  a * fxs - c * zik * fys + c * yik * fzs;
    vsite_force.fj.y =  c * zik * fxs + a * fys - c * xik * fzs;
    vsite_force.fj.z = -c * yik * fxs + c * xik * fys + a * fzs;

    vsite_force.fk.x =  b * fxs + c * zij * fys - c * yij * fzs;
    vsite_force.fk.y = -c * zij * fxs + b * fys + c * xij * fzs;
    vsite_force.fk.z =  c * yij * fxs - c * xij * fys + b * fzs;

    vsite_force.fi.x = fxs - vsite_force.fj.x - vsite_force.fk.x;
    vsite_force.fi.y = fys - vsite_force.fj.y - vsite_force.fk.y;
    vsite_force.fi.z = fzs - vsite_force.fj.z - vsite_force.fk.z;

    return vsite_force;
}


//===================================================================
// dihedral / torsion force, RB-type
// J. Comput. Chem. 2000, 21, 553-561
// doi: 10.1002/(SICI)1096-987X(200005)21:7<553::AID-JCC4>3.0.CO;2-1
//===================================================================

Vec_4 compute_dihedral_3(double rxi, double ryi, double rzi,
                         double rxj, double ryj, double rzj,
                         double rxk, double ryk, double rzk,
                         double rxl, double ryl, double rzl,
                         double c0, double c1, double c2, 
                         double c3, double c4, double c5,
                         double *box, double *potential, double **virial)
{
    Vec_R  r21, r32, r43, r21_unit, r32_unit, r43_unit;
    double r21_d, r32_d, r43_d;
    double cos123, theta123, sin123;
    double cos234, theta234, sin234;
    Vec_R  r_21_32, r_43_32;
    double const_1, const_4;
    double cos_phi, c123, b432, fij;
    Vec_4  dihedral_force;

    r21.x = pbc(rxi - rxj, box[0]);
    r21.y = pbc(ryi - ryj, box[1]);
    r21.z = pbc(rzi - rzj, box[2]);
    r21_d = sqrt(dist_2(&r21));
    r21_unit = scalar_x_vector(1.0 / r21_d, r21);

    r32.x = pbc(rxj - rxk, box[0]);
    r32.y = pbc(ryj - ryk, box[1]);
    r32.z = pbc(rzj - rzk, box[2]);
    r32_d = sqrt(dist_2(&r32));
    r32_unit = scalar_x_vector(1.0 / r32_d, r32);

    r43.x = pbc(rxk - rxl, box[0]);
    r43.y = pbc(ryk - ryl, box[1]);
    r43.z = pbc(rzk - rzl, box[2]);
    r43_d = sqrt(dist_2(&r43));
    r43_unit = scalar_x_vector(1.0 / r43_d, r43);

    cos123 = -inner_product(r21_unit, r32_unit);
    theta123 = acos(cos123);
    sin123 = sin(theta123);

    cos234 = -inner_product(r32_unit, r43_unit);
    theta234 = acos(cos234);
    sin234 = sin(theta234);

    cos_phi = (cos123 * cos234 - inner_product(r21_unit, r43_unit))
                 / (sin123 * sin234);

    // Note: in RB function the polymer convention is used
    // psi = phi - 180(degree)
    // coefficient (-1)^n is multiplied with c_n

    // potential[4] = non-metal torsional energy
    potential[4] += ( c0 - c1 * cos_phi 
                         + c2 * pow(cos_phi, 2.0) 
                         - c3 * pow(cos_phi, 3.0) 
                         + c4 * pow(cos_phi, 4.0) 
                         - c5 * pow(cos_phi, 5.0) );

    fij = -( -c1 + c2 * 2.0 * cos_phi 
                 - c3 * 3.0 * pow(cos_phi,2.0) 
                 + c4 * 4.0 * pow(cos_phi,3.0) 
                 - c5 * 5.0 * pow(cos_phi,4.0) );

    r_21_32 = outer_product(r21_unit, r32_unit);
    r_43_32 = outer_product(r43_unit, r32_unit);

    const_1 = inner_product(r43_unit, r_21_32)
				/ (r21_d * pow(sin123,3.0) * sin234);
    const_4 = inner_product(r43_unit, r_21_32) 
				/ (r43_d * pow(sin234,3.0) * sin123);

    c123 = r21_d * cos123 / r32_d - 1.0;
    b432 = r43_d * cos234 / r32_d;

    dihedral_force.fi.x = -const_1 * r_21_32.x * fij;
    dihedral_force.fi.y = -const_1 * r_21_32.y * fij;
    dihedral_force.fi.z = -const_1 * r_21_32.z * fij;

    dihedral_force.fl.x = -const_4 * r_43_32.x * fij;
    dihedral_force.fl.y = -const_4 * r_43_32.y * fij;
    dihedral_force.fl.z = -const_4 * r_43_32.z * fij;

    dihedral_force.fj.x = c123 * dihedral_force.fi.x - b432 * dihedral_force.fl.x;
    dihedral_force.fj.y = c123 * dihedral_force.fi.y - b432 * dihedral_force.fl.y;
    dihedral_force.fj.z = c123 * dihedral_force.fi.z - b432 * dihedral_force.fl.z;

    dihedral_force.fk.x = -dihedral_force.fi.x
                         - dihedral_force.fj.x
                         - dihedral_force.fl.x;
    dihedral_force.fk.y = -dihedral_force.fi.y
                         - dihedral_force.fj.y
                         - dihedral_force.fl.y;
    dihedral_force.fk.z = -dihedral_force.fi.z
                         - dihedral_force.fj.z
                         - dihedral_force.fl.z;

    // scalar viral from angle-dependent potential is zero
	// but contribution to virial elements is not zero

	// add -ri(x)fi to virial
	// note that the ri's have been treated with PBC
	// ftp://ftp.dl.ac.uk/ccp5.newsletter/26/pdf/smith.pdf
	
	Vec_R rvec;

	// i
	rvec.x = rxi;
	rvec.y = ryi;
	rvec.z = rzi;
	virial_rvec_fvec_full(&rvec, &(dihedral_force.fi), virial);

	// j
	rvec.x -= r21.x;
	rvec.y -= r21.y;
	rvec.z -= r21.z;
	virial_rvec_fvec_full(&rvec, &(dihedral_force.fj), virial);

	// k
	rvec.x -= r32.x;
	rvec.y -= r32.y;
	rvec.z -= r32.z;
	virial_rvec_fvec_full(&rvec, &(dihedral_force.fk), virial);

	// l
	rvec.x -= r43.x;
	rvec.y -= r43.y;
	rvec.z -= r43.z;
	virial_rvec_fvec_full(&rvec, &(dihedral_force.fl), virial);


    return dihedral_force;
}

//==================================================================
// dihedral / torsion force, periodic
//==================================================================

Vec_4 compute_dihedral_9(double rxi, double ryi, double rzi,
                         double rxj, double ryj, double rzj,
                         double rxk, double ryk, double rzk,
                         double rxl, double ryl, double rzl,
                         double phi0, double kphi, int n,
                         double *box, double *potential, double **virial)
{
    double c0, c1, c2, c3, c4, c5;
    Vec_4 dihedral_force;

    c0 = 0.0; c1 = 0.0; c2 = 0.0;
    c3 = 0.0; c4 = 0.0; c5 = 0.0;

    if (n == 1)
    {
        if (phi0 == 0.0)
        {
            c0 = kphi;
            c1 = kphi;
        }
        else if (phi0 == 180.0)
        {
            c0 = kphi;
            c1 = kphi * (-1.0);
        }
        else
        {
            printf("<> Error: phi0 is not equal to 0 or pi!\n");
            printf("    phi0=%f, kphi=%f, n=%d\n", phi0, kphi, n);
            exit(1);
        }
    }
    else if (n == 2)
    {
        if (phi0 == 0.0)
        {
            c2 = kphi * 2.0;
        }
        else if (phi0 == 180.0)
        {
            c0 = kphi * 2.0;
            c2 = kphi * (-2.0);
        }
        else
        {
            printf("<> Error: phi0 is not equal to 0 or pi!\n");
            printf("    phi0=%f, kphi=%f, n=%d\n", phi0, kphi, n);
            exit(1);
        }
    }
    else if (n == 3)
    {
        if (phi0 == 0.0)
        {
            c0 = kphi;
            c1 = kphi * (-3.0);
            c3 = kphi * 4.0;
        }
        else if (phi0 == 180.0)
        {
            c0 = kphi;
            c1 = kphi * 3.0;
            c3 = kphi * (-4.0);
        }
        else
        {
            printf ("<> Error: phi0 is not equal to 0 or pi!\n");
            printf ("    phi0=%f, kphi=%f, n=%d\n", phi0, kphi, n);
            exit(1);
        }
    }
    else if (n == 4)
    {
        if (phi0 == 0.0)
        {
            c0 = kphi * 2.0;
            c2 = kphi * (-8.0);
            c4 = kphi * 8.0;
        }
        else if (phi0 == 180.0)
        {
            c2 = kphi * 8.0;
            c4 = kphi * (-8.0);
        }
        else
        {
            printf("<> Error: phi0 is not equal to 0 or pi!\n");
            printf("    phi0=%f, kphi=%f, n=%d\n", phi0, kphi, n);
            exit(1);
        }
    }
    else
    {
        printf("<> Error: n is not equal to 1, 2, 3 or 4!\n");
        printf("    phi0=%f, kphi=%f, n=%d\n", phi0, kphi, n);
        exit(1);
    }

    c1 *= -1.0;
    c3 *= -1.0;
    c5 *= -1.0;

    dihedral_force = compute_dihedral_3(rxi, ryi, rzi,
                                        rxj, ryj, rzj,
                                        rxk, ryk, rzk,
                                        rxl, ryl, rzl,
                                        c0,  c1,  c2,
                                        c3,  c4,  c5,
                                        box, potential, virial);

    return dihedral_force;
}


//===================
// print a 1-D array
//===================

void print_vector(int n, double* A)
{
    int row;
    for(row = 0; row < n; row ++)
    {
        printf("%12.6e\n", A[row]);
    }
}

//===================
// print a 2-D array
//===================

void print_matrix(int n, double** A)
{
    int row, col;
    for(row = 0; row < n; row ++)
    {
        for(col = 0; col < n; col ++)
        {
            printf("%12.3e", A[row][col]);
        }
        printf("\n");
    }
}

//==========================
// print a random 2-D array
//==========================

void print_rnd_matrix(int nrows, int ncols, double** A)
{
    int row, col;
    for(row = 0; row < nrows; row ++)
    {
        for(col = 0; col < ncols; col ++)
        {
            printf("%12.3e", A[row][col]);
        }
        printf("\n");
    }
}


//==============================================
// Non-bonded forces, parallelized w.r.t. pairs
//==============================================

void mpi_nonb_pair(Topol *p_topol, Task *p_task, Atom_Info* atom_info,
				   RunSet *p_runset, System *p_system, int my_id)
{
	long int pair_start = p_task->start_pair[my_id];
	long int pair_end   = p_task->end_pair[my_id];

    int i, j, iType, jType, mol_i, atom_i, atom_j, mol_ID_i, funct;
    double rxi, ryi, rzi, rij2, qi, qj;

	double rCut2 = p_runset->rCut2;
	double rCut  = p_runset->rCut;

    Vec_nb nonbonded;
	Vec_R rvec;

    long int count = 0;
    // intermolecular nonbonded interaction
    for(i = 0; i < p_topol->n_atoms; ++ i)
    {
        if(count > pair_end) { break; }

        if(count + i < pair_start)
        {
            count += i;
            continue;
        }

        mol_i    = atom_info[i].iMol;
        atom_i   = atom_info[i].iAtom;
        mol_ID_i = atom_info[i].molID;

        rxi = p_system->rx[i];
        ryi = p_system->ry[i];
        rzi = p_system->rz[i];

        qi    = atom_info[i].charge;
        iType = atom_info[i].atomtype;

        for(j = 0; j < i; ++ j)
        {
            if(count > pair_end)
                break;

            if(count < pair_start)
            {
                ++ count;
                continue;
            }

            ++ count;

            //mol_j  = atom_info[j].iMol;
            atom_j = atom_info[j].iAtom;

            // check exclusion; note j < i
            if (atom_info[j].molID == mol_ID_i && 
				1 == p_topol->exclude[mol_i][atom_i][atom_j])
            { 
            	continue;
            }


            // check distance_ij
			rvec.x = p_system->rx[j] - rxi;
			rvec.y = p_system->ry[j] - ryi;
			rvec.z = p_system->rz[j] - rzi;

			pbc_12(&rvec, p_system->box);
            if(rvec.x >= rCut && rvec.y >= rCut && rvec.z >= rCut ) { continue; }

			rij2 = dist_2(&rvec);
            if(rij2 >= rCut2) { continue; }


            // get nonbonded parameters
            qj     = atom_info[j].charge;
            jType = atom_info[j].atomtype;

            funct = p_topol->nonbonded_param[iType][jType].funct;


            // erf_vdw potential
            if(4 == funct)
            {
                double c6    = p_topol->nonbonded_param[iType][jType].C6;
                double c12   = p_topol->nonbonded_param[iType][jType].C12;
                double width = p_topol->nonbonded_param[iType][jType].A;

                compute_erf_vdw(rvec.x, rvec.y, rvec.z, rij2, c6, c12, width,
                                p_system->potential, p_system->virial, &nonbonded);
            }

            // Morse potential
			else if(3 == funct)
            {
                double morse_D = p_topol->nonbonded_param[iType][jType].A;
                double morse_a = p_topol->nonbonded_param[iType][jType].B;
                double morse_R = p_topol->nonbonded_param[iType][jType].C6;

                compute_morse(rvec.x, rvec.y, rvec.z, rij2, 
                              morse_D, morse_a, morse_R,
                              p_system->potential, p_system->virial, &nonbonded);
            }

            // Buckingham potential
            else if(2 == funct)
            {
                double buck_A = p_topol->nonbonded_param[iType][jType].A;
                double buck_B = p_topol->nonbonded_param[iType][jType].B;
                double c6     = p_topol->nonbonded_param[iType][jType].C6;

                compute_buckingham(rvec.x, rvec.y, rvec.z, rij2, 
                                   buck_A, buck_B, c6, 
                                   p_system->potential, p_system->virial, &nonbonded);
            }

            // Lennard-Jones potential
            else if(1 == funct)
            {
                double c6  = p_topol->nonbonded_param[iType][jType].C6;
                double c12 = p_topol->nonbonded_param[iType][jType].C12;

                compute_nonbonded(rvec.x, rvec.y, rvec.z, rij2,
								  p_runset, 1.0, 1.0, 0,
                                  c6, c12, qi, qj,
                                  p_system->potential, p_system->virial, &nonbonded);
            }

            // update forces
            p_system->fx[i] += nonbonded.lj.fi.x + nonbonded.qq.fi.x;
            p_system->fy[i] += nonbonded.lj.fi.y + nonbonded.qq.fi.y;
            p_system->fz[i] += nonbonded.lj.fi.z + nonbonded.qq.fi.z;

            p_system->fx[j] += nonbonded.lj.fj.x + nonbonded.qq.fj.x;
            p_system->fy[j] += nonbonded.lj.fj.y + nonbonded.qq.fj.y;
            p_system->fz[j] += nonbonded.lj.fj.z + nonbonded.qq.fj.z;
        }
    }
}


//=============================================================
// zero the forces, potential energies and virial
//=============================================================

void zero_force_pot_vir(int nAtoms, double* fx, double* fy, double* fz,
                        double* potential, double** virial)
{
    // zero force
    int i, j;
    for(i = 0; i < nAtoms; i++)
    {
        fx[i] = 0.0;
        fy[i] = 0.0;
        fz[i] = 0.0;
    }
    
    // zero potential
    for(i = 0; i < 15; i++)
    {
        potential[i] = 0.0;
    }
    
    // zero virial
    for(i = 0; i < 3; i++)
    {
        for(j = 0; j < 3; j++)
        {
            virial[i][j] = 0.0;
        }
    }
}



//======================================
// build virtual sites from real atoms
// for now only TIP5P type is supported
//======================================

void build_vsites(int mol_start, int mol_end, 
                  Mol_Info* mol_info, Atom_Info* atom_info,
                  int* nVSites, int** vsite_funct, VSite_4** vsite_4,
                  double* rx, double* ry, double* rz, double* box)
{
    // loop over molecules
    int im;
    for(im = mol_start; im <= mol_end; ++ im)
    {
        int first_atom = mol_info[im].mini;
        int mol = atom_info[first_atom].iMol;
    
        // build each virtual site using absolute value (a,b,c)
        int iVSite;
        for(iVSite = 0; iVSite < nVSites[mol]; ++ iVSite)
        {
            int i, j, k, s;
            double a, b, c;
            Vec_R  rij, rik, rijik;

            if ( vsite_funct[mol][iVSite] == 4 )
            {
                i = vsite_4[mol][iVSite].atom_i + first_atom;
                j = vsite_4[mol][iVSite].atom_j + first_atom;
                k = vsite_4[mol][iVSite].atom_k + first_atom;
                s = vsite_4[mol][iVSite].atom_s + first_atom;

                a = vsite_4[mol][iVSite].a;
                b = vsite_4[mol][iVSite].b;
                c = vsite_4[mol][iVSite].c;

                rij.x = pbc(rx[j] - rx[i], box[0]);
                rij.y = pbc(ry[j] - ry[i], box[1]);
                rij.z = pbc(rz[j] - rz[i], box[2]);
    
                rik.x = pbc(rx[k] - rx[i], box[0]);
                rik.y = pbc(ry[k] - ry[i], box[1]);
                rik.z = pbc(rz[k] - rz[i], box[2]);
    
                rijik = outer_product(rij, rik);
    
                rx[s] = rx[i] + a * rij.x + b * rik.x + c * rijik.x;
                ry[s] = ry[i] + a * rij.y + b * rik.y + c * rijik.y;
                rz[s] = rz[i] + a * rij.z + b * rik.z + c * rijik.z;
            }
        }
    }
}


//========================================================
// distribute all forces from virtual sites to real atoms
//========================================================

void mpi_vsites(int mol_start, int mol_end, 
                Mol_Info* mol_info, Atom_Info* atom_info,
                int* nVSites, int** vsite_funct, VSite_4** vsite_4,
                double* rx, double* ry, double* rz, double* box,
                double* fx, double* fy, double* fz)
{
    // loop over molecules
    int im;
    for(im = mol_start; im <= mol_end; ++ im)
    {
        int first_atom = mol_info[im].mini;
        int mol = atom_info[first_atom].iMol;
    
        // distribute forces for each virtual site
        int iVSite;
        for(iVSite = 0; iVSite < nVSites[mol]; ++ iVSite)
        {
            int i, j, k, s;
            Vec_4 vsite_force;
            if(vsite_funct[mol][iVSite] == 4)
            {
                i = vsite_4[mol][iVSite].atom_i + first_atom;
                j = vsite_4[mol][iVSite].atom_j + first_atom;
                k = vsite_4[mol][iVSite].atom_k + first_atom;
                s = vsite_4[mol][iVSite].atom_s + first_atom;
    
                vsite_force = compute_vsite_4 (rx[i], ry[i], rz[i],
                                               rx[j], ry[j], rz[j],
                                               rx[k], ry[k], rz[k],
                                               rx[s], ry[s], rz[s],
                                               fx[s], fy[s], fz[s],
                                               vsite_4[mol][iVSite].a,
                                               vsite_4[mol][iVSite].b,
                                               vsite_4[mol][iVSite].c,
                                               box);
    
                fx[i] += vsite_force.fi.x;
                fy[i] += vsite_force.fi.y;
                fz[i] += vsite_force.fi.z;
                fx[j] += vsite_force.fj.x;
                fy[j] += vsite_force.fj.y;
                fz[j] += vsite_force.fj.z;
                fx[k] += vsite_force.fk.x;
                fy[k] += vsite_force.fk.y;
                fz[k] += vsite_force.fk.z;
                fx[s] = 0.0;
                fy[s] = 0.0;
                fz[s] = 0.0;
            }
        }
    }
}


//=============================================================
// compute all bonded forces, parallelized w.r.t. molecules
//=============================================================

void mpi_bonded(Task *p_task, Mol_Info* mol_info, Atom_Info* atom_info,
				Topol *p_topol, RunSet* p_runset, System* p_system, int my_id)
{
	int mol_start = p_task->start_mol[my_id];
	int mol_end   = p_task->end_mol[my_id];

    // bonded interaction (including 1-4 pair)
    int im;
    for(im = mol_start; im <= mol_end; ++ im)
    {
        int first_atom = mol_info[im].mini;
        int mol = atom_info[first_atom].iMol;
    
        // bonds
        int iBond;
        for(iBond = 0; iBond < p_topol->n_bonds[mol]; ++ iBond)
        {
            int i, j;
            Vec_2 bond_force;

            i = p_topol->bond_param[mol][iBond].atom_i + first_atom;
            j = p_topol->bond_param[mol][iBond].atom_j + first_atom;

            if(p_topol->bond_param[mol][iBond].funct == 1)
            {
                bond_force = compute_bond_1(p_system->rx[i], p_system->ry[i], p_system->rz[i],
                                            p_system->rx[j], p_system->ry[j], p_system->rz[j],
                                            p_topol->bond_param[mol][iBond].b0,
                                            p_topol->bond_param[mol][iBond].kb,
                                            p_system->box, p_system->potential, p_system->virial);
            }
            else if(p_topol->bond_param[mol][iBond].funct == 3)
            {
                bond_force = compute_bond_3(p_system->rx[i], p_system->ry[i], p_system->rz[i],
                                            p_system->rx[j], p_system->ry[j], p_system->rz[j],
                                            p_topol->bond_param[mol][iBond].b0,
                                            p_topol->bond_param[mol][iBond].D,
                                            p_topol->bond_param[mol][iBond].beta,
                                            p_system->box, p_system->potential, p_system->virial);
            }
    
            p_system->fx[i] += bond_force.fi.x;
            p_system->fy[i] += bond_force.fi.y;
            p_system->fz[i] += bond_force.fi.z;
            p_system->fx[j] += bond_force.fj.x;
            p_system->fy[j] += bond_force.fj.y;
            p_system->fz[j] += bond_force.fj.z;
        }
    
        // pairs
        int iPair;
        for(iPair = 0; iPair < p_topol->n_pairs[mol]; ++ iPair)
        {
            int i, j, iType, jType;
            double qi, qj, c6, c12;
            Vec_nb nonbonded;

            i = p_topol->pair_param[mol][iPair].atom_i + first_atom;
            j = p_topol->pair_param[mol][iPair].atom_j + first_atom;
    
            double rxij = pbc(p_system->rx[j] - p_system->rx[i], p_system->box[0]);
            double ryij = pbc(p_system->ry[j] - p_system->ry[i], p_system->box[1]);
            double rzij = pbc(p_system->rz[j] - p_system->rz[i], p_system->box[2]);
            double rij2 = rxij*rxij + ryij*ryij + rzij*rzij;

            // normal 1-4 pair
            if(p_topol->pair_param[mol][iPair].funct == 1)
            {
                qi     = atom_info[i].charge;
                iType = atom_info[i].atomtype;
                qj     = atom_info[j].charge;
                jType = atom_info[j].atomtype;
                c6     = p_topol->nonbonded_param[iType][jType].C6;
                c12    = p_topol->nonbonded_param[iType][jType].C12;
    
                compute_nonbonded(rxij, ryij, rzij, rij2,
								  p_runset, p_topol->scaleLJ, p_topol->scaleQQ, 1,
                                  c6, c12, qi, qj,
                                  p_system->potential, p_system->virial, &nonbonded);
            }

            // user-defined pair
            else if(p_topol->pair_param[mol][iPair].funct == 2)
            {
                // V = sigma, W = epsilon
                qi  = p_topol->pair_param[mol][iPair].q_i;
                qj  = p_topol->pair_param[mol][iPair].q_j;

                double sig_6 = pow(p_topol->pair_param[mol][iPair].V, 6.0);
                double eps_4 = 4.0 * p_topol->pair_param[mol][iPair].W;
                c6  = eps_4 * sig_6;
                c12 = c6 * sig_6;
    
                compute_nonbonded(rxij, ryij, rzij, rij2,
								  p_runset, 1.0, p_topol->pair_param[mol][iPair].scaleQQ, 1,
                                  c6, c12, qi, qj,
                                  p_system->potential, p_system->virial, &nonbonded);
            }
    
            p_system->fx[i] += nonbonded.lj.fi.x + nonbonded.qq.fi.x;
            p_system->fy[i] += nonbonded.lj.fi.y + nonbonded.qq.fi.y;
            p_system->fz[i] += nonbonded.lj.fi.z + nonbonded.qq.fi.z;
            p_system->fx[j] += nonbonded.lj.fj.x + nonbonded.qq.fj.x;
            p_system->fy[j] += nonbonded.lj.fj.y + nonbonded.qq.fj.y;
            p_system->fz[j] += nonbonded.lj.fj.z + nonbonded.qq.fj.z;
        }
    
        // angles
        int iAngle;
        for(iAngle = 0; iAngle < p_topol->n_angles[mol]; ++ iAngle)
        {
            int i, j, k;
            Vec_3 angle_force;

            i = p_topol->angle_param[mol][iAngle].atom_i + first_atom;
            j = p_topol->angle_param[mol][iAngle].atom_j + first_atom;
            k = p_topol->angle_param[mol][iAngle].atom_k + first_atom;
    
            if(p_topol->angle_param[mol][iAngle].funct == 1)
            {
                angle_force = compute_angle_1(p_system->rx[i], p_system->ry[i], p_system->rz[i],
                                              p_system->rx[j], p_system->ry[j], p_system->rz[j],
                                              p_system->rx[k], p_system->ry[k], p_system->rz[k],
                                              p_topol->angle_param[mol][iAngle].a0,
                                              p_topol->angle_param[mol][iAngle].ka,
                                              p_system->box, p_system->potential, p_system->virial);
            }
            else if(p_topol->angle_param[mol][iAngle].funct == 5)
            {
                angle_force = compute_angle_5(p_system->rx[i], p_system->ry[i], p_system->rz[i],
                                              p_system->rx[j], p_system->ry[j], p_system->rz[j],
                                              p_system->rx[k], p_system->ry[k], p_system->rz[k],
                                              p_topol->angle_param[mol][iAngle].a0,
                                              p_topol->angle_param[mol][iAngle].ka,
                                              p_topol->angle_param[mol][iAngle].b0_ub,
                                              p_topol->angle_param[mol][iAngle].kb_ub,
                                              p_system->box, p_system->potential, p_system->virial);
            }
    
            p_system->fx[i] += angle_force.fi.x;
            p_system->fy[i] += angle_force.fi.y;
            p_system->fz[i] += angle_force.fi.z;
            p_system->fx[j] += angle_force.fj.x;
            p_system->fy[j] += angle_force.fj.y;
            p_system->fz[j] += angle_force.fj.z;
            p_system->fx[k] += angle_force.fk.x;
            p_system->fy[k] += angle_force.fk.y;
            p_system->fz[k] += angle_force.fk.z;
        }
    
        // dihedrals
        int iDihedral;
        for(iDihedral = 0; iDihedral < p_topol->n_dihedrals[mol]; ++ iDihedral)
        {
            int i, j, k, l;
            Vec_4 dihedral_force;

            i = p_topol->dihedral_param[mol][iDihedral].atom_i + first_atom;
            j = p_topol->dihedral_param[mol][iDihedral].atom_j + first_atom;
            k = p_topol->dihedral_param[mol][iDihedral].atom_k + first_atom;
            l = p_topol->dihedral_param[mol][iDihedral].atom_l + first_atom;
    
            if(p_topol->dihedral_param[mol][iDihedral].funct == 3)
            {
                dihedral_force = compute_dihedral_3(p_system->rx[i], p_system->ry[i], p_system->rz[i],
                                                    p_system->rx[j], p_system->ry[j], p_system->rz[j],
                                                    p_system->rx[k], p_system->ry[k], p_system->rz[k],
                                                    p_system->rx[l], p_system->ry[l], p_system->rz[l],
                                                    p_topol->dihedral_param[mol][iDihedral].c0,
                                                    p_topol->dihedral_param[mol][iDihedral].c1,
                                                    p_topol->dihedral_param[mol][iDihedral].c2,
                                                    p_topol->dihedral_param[mol][iDihedral].c3,
                                                    p_topol->dihedral_param[mol][iDihedral].c4,
                                                    p_topol->dihedral_param[mol][iDihedral].c5,
                                                    p_system->box, p_system->potential, p_system->virial);
            }
            else if(p_topol->dihedral_param[mol][iDihedral].funct == 1 ||
                    p_topol->dihedral_param[mol][iDihedral].funct == 4 ||
                    p_topol->dihedral_param[mol][iDihedral].funct == 9 )
            {
                dihedral_force = compute_dihedral_9(p_system->rx[i], p_system->ry[i], p_system->rz[i],
                                                    p_system->rx[j], p_system->ry[j], p_system->rz[j],
                                                    p_system->rx[k], p_system->ry[k], p_system->rz[k],
                                                    p_system->rx[l], p_system->ry[l], p_system->rz[l],
                                                    p_topol->dihedral_param[mol][iDihedral].phi0,
                                                    p_topol->dihedral_param[mol][iDihedral].kphi,
                                                    p_topol->dihedral_param[mol][iDihedral].n,
                                                    p_system->box, p_system->potential, p_system->virial);
            }
    
            p_system->fx[i] += dihedral_force.fi.x;
            p_system->fy[i] += dihedral_force.fi.y;
            p_system->fz[i] += dihedral_force.fi.z;
            p_system->fx[j] += dihedral_force.fj.x;
            p_system->fy[j] += dihedral_force.fj.y;
            p_system->fz[j] += dihedral_force.fj.z;
            p_system->fx[k] += dihedral_force.fk.x;
            p_system->fy[k] += dihedral_force.fk.y;
            p_system->fz[k] += dihedral_force.fk.z;
            p_system->fx[l] += dihedral_force.fl.x;
            p_system->fy[l] += dihedral_force.fl.y;
            p_system->fz[l] += dihedral_force.fl.z;
        }
    }
}



//==================================================================================
// RATTLE constraint algorithm for the 1st half step
// http://www.ccl.net/cca/software/SOURCES/FORTRAN/allen-tildesley-book/f.09.shtml
//==================================================================================

void rattle_1st(double dt, Mol_Info* mol_info, Atom_Info* atom_info,
				Topol *p_topol, System *p_system)
{
    //double tol = 1.0e-05;
    double tol2 = 2.0e-05;

    // loop over all molecules
    int im;
    for (im = 0; im < p_topol->n_mols; ++ im) 
    {
		int first_atom, mol, iter, done;

        first_atom = mol_info[im].mini;
        mol = atom_info[first_atom].iMol;

        iter = 0;
        done = 0;

        // start rattle_1
        while (0 == done)
        {
            done = 1;

            int iCstr;
            for (iCstr = 0; iCstr < p_topol->n_constraints[mol]; ++ iCstr)
            {
				int a, b;
				double rab, pxab, pyab, pzab, pabSQ, rabSQ, diffSQ;
				double rxab, ryab, rzab, rpab, rma, rmb, gab, dx, dy, dz;

                a = p_topol->constraint[mol][iCstr].atom_i + first_atom;
                b = p_topol->constraint[mol][iCstr].atom_j + first_atom;
                rab = p_topol->constraint[mol][iCstr].r;

                pxab = pbc(p_system->rx[a] - p_system->rx[b], p_system->box[0]);
                pyab = pbc(p_system->ry[a] - p_system->ry[b], p_system->box[1]);
                pzab = pbc(p_system->rz[a] - p_system->rz[b], p_system->box[2]);

                pabSQ = pxab * pxab + pyab * pyab + pzab * pzab;
                rabSQ = rab * rab;
                diffSQ = rabSQ - pabSQ;

                if (fabs(diffSQ) > rabSQ * tol2)
                {
                    rxab = pbc(p_system->old_rx[a] - p_system->old_rx[b], p_system->box[0]);
                    ryab = pbc(p_system->old_ry[a] - p_system->old_ry[b], p_system->box[1]);
                    rzab = pbc(p_system->old_rz[a] - p_system->old_rz[b], p_system->box[2]);
                    rpab = rxab * pxab + ryab * pyab + rzab * pzab;

                    rma = atom_info[a].inv_mass;
                    rmb = atom_info[b].inv_mass;
                    gab = diffSQ / (2.0 * (rma + rmb) * rpab);

                    dx = rxab * gab;
                    dy = ryab * gab;
                    dz = rzab * gab;

                    p_system->rx[a] += rma * dx;
                    p_system->ry[a] += rma * dy;
                    p_system->rz[a] += rma * dz;
                    p_system->rx[b] -= rmb * dx;
                    p_system->ry[b] -= rmb * dy;
                    p_system->rz[b] -= rmb * dz;

                    dx /= dt;
                    dy /= dt;
                    dz /= dt;

                    p_system->vx[a] += rma * dx;
                    p_system->vy[a] += rma * dy;
                    p_system->vz[a] += rma * dz;
                    p_system->vx[b] -= rmb * dx;
                    p_system->vy[b] -= rmb * dy;
                    p_system->vz[b] -= rmb * dz;

                    done = 0;
                }
            }

            iter ++;
            if (iter > MAX_ITER) { break; }
        }

        if (0 == done) 
        {
            printf( "RATTLE_1 not converged!\n");
            exit(1);
        }
        //printf("iter = %d in RATTLE_1\n", iter);
    }
}

//==================================================================================
// RATTLE constraint algorithm for the 2nd half step
// http://www.ccl.net/cca/software/SOURCES/FORTRAN/allen-tildesley-book/f.09.shtml
//==================================================================================

void rattle_2nd(double dt, Mol_Info* mol_info, Atom_Info* atom_info,
				Topol *p_topol, System *p_system)
{
    double tol = 1.0e-05;
    double half_dt = 0.5 * dt;

    // loop over all molecules
    int im;
    for(im = 0; im < p_topol->n_mols; ++ im) 
    {
        int first_atom = mol_info[im].mini;
        int mol = atom_info[first_atom].iMol;

        int iter = 0;
        int done = 0;

        // start rattle_2
        while(done == 0)
        {
            done = 1;

            int iCstr;
            for (iCstr = 0; iCstr < p_topol->n_constraints[mol]; ++ iCstr)
            {
                int a = p_topol->constraint[mol][iCstr].atom_i + first_atom;
                int b = p_topol->constraint[mol][iCstr].atom_j + first_atom;
                double rab = p_topol->constraint[mol][iCstr].r;

                double vxab = p_system->vx[a] - p_system->vx[b];
                double vyab = p_system->vy[a] - p_system->vy[b];
                double vzab = p_system->vz[a] - p_system->vz[b];

                double rxab = pbc(p_system->rx[a] - p_system->rx[b], p_system->box[0]);
                double ryab = pbc(p_system->ry[a] - p_system->ry[b], p_system->box[1]);
                double rzab = pbc(p_system->rz[a] - p_system->rz[b], p_system->box[2]);

                double rvab = rxab * vxab + ryab * vyab + rzab * vzab;
                double rma = atom_info[a].inv_mass;
                double rmb = atom_info[b].inv_mass;
                double rabSQ = rab * rab;
                double gab = -rvab / ((rma + rmb) * rabSQ);

                if(fabs(gab) > tol)
                {
                    double dx = rxab * gab;
                    double dy = ryab * gab;
                    double dz = rzab * gab;

                    p_system->vx[a] += rma * dx;
                    p_system->vy[a] += rma * dy;
                    p_system->vz[a] += rma * dz;
                    p_system->vx[b] -= rmb * dx;
                    p_system->vy[b] -= rmb * dy;
                    p_system->vz[b] -= rmb * dz;

                    // update virial tensor

                    dx /= half_dt;
                    dy /= half_dt;
                    dz /= half_dt;

                    p_system->virial[0][0] -= rxab * dx;
                    p_system->virial[0][1] -= rxab * dy;
                    p_system->virial[0][2] -= rxab * dz;
                    p_system->virial[1][0] -= ryab * dx;
                    p_system->virial[1][1] -= ryab * dy;
                    p_system->virial[1][2] -= ryab * dz;
                    p_system->virial[2][0] -= rzab * dx;
                    p_system->virial[2][1] -= rzab * dy;
                    p_system->virial[2][2] -= rzab * dz;

                    done = 0;
                }
            }

            iter ++;
            if(iter > MAX_ITER) { break; }
        }

        if(done == 0) 
        {
            printf( "RATTLE_2 not converged!\n");
            exit(1);
        }
        //printf("iter = %d in RATTLE_2", iter);
    }
}


//=============================================
// compute Q-SC, bonded and non-bonded forces
//=============================================

void mpi_force(Task *p_task, Topol *p_topol,
               Atom_Info* atom_info, Mol_Info* mol_info, 
			   RunSet* p_runset, Metal *p_metal, System *p_system,
               int my_id, int num_procs, double** time_used)
{
	//double rCut     = p_runset->rCut;
	double rCut2    = p_runset->rCut2;


    // monitor time usage
    double t0, t1;
    struct timeval tv;

    const int root_process = 0;
    int an_id, ierr;
    MPI_Status status;

    int tag_11 = 11; // slave to master, partial_fx, partial_fy, partial_fz 
    int tag_12 = 12; // slave to master, p_system->partial_vir
    int tag_13 = 13; // slave to master, p_system->partial_pot

    //int i_mat, j_mat;
    int n_mat = p_metal->num * 4 + p_metal->n_NPs;

    int i, j;


    //======== build virtual sites on root processor ====================

    if(my_id == root_process) 
    {
        // build all virtual sites before distributing atomic coordinates
        build_vsites(0, p_topol->n_mols-1, mol_info, atom_info,
                     p_topol->n_vsites, p_topol->vsite_funct, p_topol->vsite_4, 
					 p_system->rx, p_system->ry, p_system->rz, p_system->box);
    }


    //======== broadcast rx,ry,rz and do calculations on each processor =============

    gettimeofday(&tv, NULL);
    t0 = tv.tv_sec + tv.tv_usec * 1.0e-6;

    ierr = MPI_Bcast(&(p_system->rx[0]), p_topol->n_atoms*DIM, MPI_DOUBLE, 
					 root_process, MPI_COMM_WORLD);

    zero_force_pot_vir(p_topol->n_atoms, p_system->fx, p_system->fy, p_system->fz, 
					   p_system->potential, p_system->virial);



	// update box volume
	p_system->volume = p_system->box[0] * p_system->box[1] * p_system->box[2];
	p_system->inv_volume = 1.0 / p_system->volume;


    gettimeofday(&tv, NULL);
    t1 = tv.tv_sec + tv.tv_usec * 1.0e-6;
    time_used[my_id][10] += (t1 - t0);
    t0 = t1;


    // calculate bonded forces
    mpi_bonded(p_task, mol_info, atom_info, p_topol, p_runset, p_system, my_id);

    gettimeofday(&tv, NULL);
    t1 = tv.tv_sec + tv.tv_usec * 1.0e-6;
    time_used[my_id][3] += (t1 - t0);
    t0 = t1;

    // calculate nonbonded forces
    mpi_nonb_pair(p_topol, p_task, atom_info, p_runset, p_system, my_id);

    gettimeofday(&tv, NULL);
    t1 = tv.tv_sec + tv.tv_usec * 1.0e-6;
    time_used[my_id][4] += (t1 - t0);
    t0 = t1;


    // compute QSC densities and forces
    if(p_metal->min >=0 && p_metal->num > 0)
    {
        mpi_qsc_dens(p_task, rCut2, p_metal, p_system, my_id, num_procs);

        gettimeofday(&tv, NULL);
        t1 = tv.tv_sec + tv.tv_usec * 1.0e-6;
        time_used[my_id][1] += (t1 - t0);
        t0 = t1;

        mpi_qsc_force(p_task, rCut2, p_metal, p_system, my_id);

        gettimeofday(&tv, NULL);
        t1 = tv.tv_sec + tv.tv_usec * 1.0e-6;
        time_used[my_id][2] += (t1 - t0);
        t0 = t1;

        // CPIM model
        if(p_metal->use_cpff)
        {
            // construct p_metal->vec_ext
            mpi_cpff_vec_ext(p_task, p_metal, p_runset, atom_info, p_topol,
							 p_system, my_id);

            gettimeofday(&tv, NULL);
            t1 = tv.tv_sec + tv.tv_usec * 1.0e-6;
            time_used[my_id][6] += (t1 - t0);
            t0 = t1;

			// construct the sparse CPIM matrix
			long int count_nnz;
            mpi_cpff_mat_relay_count(p_task, p_metal, p_system, p_runset->rCut2, 
									 my_id, num_procs, &count_nnz);

			p_metal->val     = my_malloc_2(sizeof(double)   * count_nnz, "val");
			p_metal->col_ind = my_malloc_2(sizeof(long int) * count_nnz, "col_ind");
			p_metal->row_ind = my_malloc_2(sizeof(long int) * count_nnz, "row_ind");

            // construct p_metal->mat_relay
            mpi_cpff_mat_relay_CRS(p_task, p_metal, p_system, p_runset->rCut2, my_id, num_procs);


            gettimeofday(&tv, NULL);
            t1 = tv.tv_sec + tv.tv_usec * 1.0e-6;
            time_used[my_id][5] += (t1 - t0);
            t0 = t1;

            // solve Ax=b
			// using preconditioned BiCGSTAB
			{
				mpi_precon_bicg_stab_CRS(p_task->start_metal[my_id], p_task->end_metal[my_id], 
									 p_metal->min, p_metal->max, p_metal->n_NPs, n_mat, 
									 p_metal->diag_relay, p_metal->vec_ext, p_metal->vec_pq, 
									 my_id, num_procs, p_metal, count_nnz);
			}
            // at the end of mpi_bicg_stab each proc has a copy of p_metal->vec_pq


			free(p_metal->val);
			free(p_metal->col_ind);
			free(p_metal->row_ind);

			p_metal->val = NULL;
			p_metal->col_ind = NULL;
			p_metal->row_ind = NULL;


            gettimeofday(&tv, NULL);
            t1 = tv.tv_sec + tv.tv_usec * 1.0e-6;
            time_used[my_id][7] += (t1 - t0);
            t0 = t1;

            // compute CPIM forces
            mpi_cpff_force(p_task, p_metal, p_runset,
						   p_system, p_metal->vec_pq,
                           p_topol->n_atoms, atom_info, my_id);

            gettimeofday(&tv, NULL);
            t1 = tv.tv_sec + tv.tv_usec * 1.0e-6;
            time_used[my_id][8] += (t1 - t0);
            t0 = t1;
        }
    }


    //======== root processor ==========================
    if(my_id == root_process) 
    {

        // receive nonbonded & QSC forces from slave processors
        for(an_id = 1; an_id < num_procs; an_id++) 
        {
            ierr = MPI_Recv(&(p_system->partial_fx[0]), p_topol->n_atoms*DIM, MPI_DOUBLE, 
							MPI_ANY_SOURCE, tag_11, MPI_COMM_WORLD, &status);
            //sender = status.MPI_SOURCE;

            ierr = MPI_Recv(&(p_system->partial_vir[0][0]), DIM*DIM, MPI_DOUBLE, 
							MPI_ANY_SOURCE, tag_12, MPI_COMM_WORLD, &status);
            ierr = MPI_Recv(&(p_system->partial_pot[0]), 15, MPI_DOUBLE, 
							MPI_ANY_SOURCE, tag_13, MPI_COMM_WORLD, &status);

            for(i = 0; i < p_topol->n_atoms; i++)
            {
                p_system->fx[i] += p_system->partial_fx[i];
                p_system->fy[i] += p_system->partial_fy[i];
                p_system->fz[i] += p_system->partial_fz[i];
            }

            for(i = 0; i < DIM; i++)
			{
                for(j = 0; j < DIM; j++)
				{
                    p_system->virial[i][j] += p_system->partial_vir[i][j];
				}
			}

            for(i = 0; i < 15; i++)
			{
                p_system->potential[i] += p_system->partial_pot[i];
			}
        }


        // distribute forces from virtual sites to real atoms 
        // after computation of all forces
        mpi_vsites(0, p_topol->n_mols-1, mol_info, atom_info, 
				   p_topol->n_vsites, p_topol->vsite_funct, p_topol->vsite_4,
                   p_system->rx, p_system->ry, p_system->rz, p_system->box, 
				   p_system->fx, p_system->fy, p_system->fz);

    }
    //======== slave processor ==========================
    else 
    {
        // send nonbonded & QSC forces back to root processor
        ierr = MPI_Send(&(p_system->fx[0]), p_topol->n_atoms*DIM, MPI_DOUBLE, root_process, tag_11, MPI_COMM_WORLD);
        ierr = MPI_Send(&(p_system->virial[0][0]), DIM*DIM, MPI_DOUBLE, root_process, tag_12, MPI_COMM_WORLD);
        ierr = MPI_Send(&(p_system->potential[0]), 15, MPI_DOUBLE, root_process, tag_13, MPI_COMM_WORLD);
    }

    gettimeofday(&tv, NULL);
    t1 = tv.tv_sec + tv.tv_usec * 1.0e-6;
    time_used[my_id][11] += (t1 - t0);
    t0 = t1;

    if(ierr) {}
}

