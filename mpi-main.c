/*****************************************************************************
 *  This file is part of the CapacMD program.                                *  
 *  Copyright (C) 2015 Xin Li <lixin.reco@gmail.com>                         *  
 *                                                                           *
 *  Filename:  mpi-main.c                                                    *
 *  Function:  main engine of the CapacMD program                            *
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
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "typedef.h"
#include "my_malloc.h"

#include "file.h"
#include "force.h"
#include "thermostat.h"

int main(int argc, char *argv[])
{
    //============= set up MPI =================
    
    MPI_Status status;
    const int root_process = 0;
    int ierr, my_id, an_id, num_procs;

    // tag for MPI message: energy minimization convergence 
    const int tag_99 = 99;

    // initialize MPI, get my_id for this process and total number of processes
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);


    //============ get input files =============

    // default filenames
    char input_gro[MAX_STR_LEN], input_param[MAX_STR_LEN], input_mdset[MAX_STR_LEN];
    strcpy(input_gro,   "init.gro");
    strcpy(input_param, "param.txt");
    strcpy(input_mdset, "mdset.txt");

    // parse arguments (input files)
    // -gro for gro file
    // -par for parameter file
    // -mds for MD settings
    int iarg = 1;
    while (1)
    {
        if (iarg >= argc) { break; }

        if (0 == strcmp(argv[iarg], "-gro") && iarg < argc)
        {
            strcpy(input_gro, argv[iarg + 1]);
            iarg += 2;
        }
        else if (0 == strcmp(argv[iarg], "-par") && iarg < argc)
        {
            strcpy(input_param, argv[iarg + 1]);
            iarg += 2;
        }
        else if (0 == strcmp(argv[iarg], "-mds") && iarg < argc)
        {
            strcpy(input_mdset, argv[iarg + 1]);
            iarg += 2;
        }
        else
        {
            ++ iarg;
        }
    }


    //========= define variables and read MD settings ===========

    // varialbles to read from input_mdset
    RunSet *p_runset = (RunSet *)my_malloc(sizeof(RunSet));
    Metal  *p_metal  = (Metal  *)my_malloc(sizeof(Metal));

    p_runset->external_efield = (double *)my_malloc(sizeof(double) * DIM);
    p_runset->external_efield[0] = 0.0;
    p_runset->external_efield[1] = 0.0;
    p_runset->external_efield[2] = 0.0;


    // read md settings from input_mdset
    read_settings(input_mdset, p_runset, p_metal);


    // initialize timer
    time_t start_t = time(NULL);
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double start_time = (tv.tv_sec) + (tv.tv_usec) * 1.0e-6;

    // time_used[0] = total
    // time_used[1] = QSC density
    // time_used[2] = QSC force
    // time_used[3] = Bonded
    // time_used[4] = Nonbonded
    // time_used[5] = CPIM matrix
    // time_used[6] = CPIM vector
    // time_used[7] = CPIM solve
    // time_used[8] = CPIM force
    // time_used[10] = QSC density communication
    
    double **time_used = (double **)my_malloc(sizeof(double *) * num_procs);
    for(an_id = 0; an_id < num_procs; ++ an_id)
    {
        time_used[an_id] = (double *)my_malloc(sizeof(double) * 15);

        int it;
        for(it = 0; it < 15; ++ it)
        {
            time_used[an_id][it] = 0.0;
        }
    }


    // Coulomb type: cut_off, wolf_sum
    if (0 == strcmp(p_runset->coulomb_type, "cut_off" )) 
    { 
        p_runset->use_coulomb = 0; 
    }
    else if (0 == strcmp(p_runset->coulomb_type, "wolf_sum")) 
    { 
        p_runset->use_coulomb = 1; 
    }
    else 
    { 
        printf("Error: unknown coulomb_type %s!\n", p_runset->coulomb_type); 
        exit(1); 
    }


    // Damped shifted force (DSF) approach for electrostatics 
    // Ref.:      Fennell and Gezelter, J. Chem. Phys. 2006, 124, 234104
    //            dx.doi.org/10.1063/1.2206581
    // Based on:  Wolf et al., J. Chem. Phys. 1999, 110, 8254
    //            dx.doi.org/10.1063/1.478738
    //            Zahn et al., J. Phys. Chem. B 2002, 106, 10725-10732
    //            dx.doi.org/10.1021/jp025949h
    // Benchmark: McCann and Acevedo, J. Chem. Theory Comput. 2013, 9, 944-950
    //            dx.doi.org/10.1021/ct300961e

    // note: p_runset->rCut and p_runset->w_alpha were read from mdset.txt
    p_runset->rCut2 = p_runset->rCut * p_runset->rCut;
    double rCut     = p_runset->rCut;
    double rCut2    = p_runset->rCut2;

    double w_alpha         = p_runset->w_alpha;
    p_runset->w_a2_sqrtPI  = w_alpha * 2.0 * INV_SQRT_PI;
    p_runset->w_Const      = erfc(w_alpha * rCut) / rCut2 + 
                             p_runset->w_a2_sqrtPI * exp(-w_alpha * w_alpha * rCut2) / rCut;
    p_runset->w_erfc_arCut = erfc(w_alpha * rCut)/ rCut;



    // van der Waals type: cut_off or shifted
    if (0 == strcmp(p_runset->vdw_type, "cut_off")) 
    { 
        p_runset->use_vdw = 0; 
    }
    else if (0 == strcmp(p_runset->vdw_type, "shifted")) 
    { 
        p_runset->use_vdw = 1; 
    }
    else 
    { 
        printf("Error: unknown vdw_type %s!\n", p_runset->vdw_type); 
        exit(1); 
    }


    // Shifted force method for Lennard-Jones potential
    // Ref:   Toxvaerd and Dyre, J. Chem. Phys. 2011, 134, 081102
    //        dx.doi.org/10.1063/1.3558787

    p_runset->inv_rc12 = 1.0 / pow(rCut, 12.0);
    p_runset->inv_rc6  = 1.0 / pow(rCut, 6.0);


    //============== read force field parameters from param.txt =================
    
    Topol *p_topol = (Topol *)my_malloc(sizeof(Topol));

    int mol, atom;
    int number_VSites, number_Cstrs;

    // variables for bonded and nonbonded interaction


    // CPIM: capacitance-polarizability interaction model
    // Ref.: a) Jensen and Jensen, J. Phys. Chem. C, 2008, 112, 15697-15703
    //          dx.doi.org/10.1021/jp804116z
    //       b) Morton and Jensen, J. Chem. Phys., 2010, 133, 074103
    //          dx.doi.org/10.1063/1.3457365
    //
    // p_metal->cpff_polar: polarizability
    // p_metal->cpff_capac: capacitance
    // p_metal->n_NPs:      number of nanoparticles
    // p_metal->cpff_chg:   total charge of a nanoparticle
    // p_metal->start_NP:   first atom of a nanoparticle
    // p_metal->end_NP:     last atom of a nanoparticle
    

    // read parameters from input_param, step 1
    read_param_1(input_param, p_topol, p_metal);


    // allocate memory for arrays

    // CPIM charge and indices
    p_metal->cpff_chg = (double *)my_malloc(p_metal->n_NPs * sizeof(double));
    p_metal->start_NP = (int *)my_malloc(p_metal->n_NPs * sizeof(int));
    p_metal->end_NP   = (int *)my_malloc(p_metal->n_NPs * sizeof(int));


    // molecules and atom parameters
    // For a given molecule type (mol from 0 to p_topol->mol_types-1):
    // p_topol->atom_num[mol]:  its number of atoms
    // p_topol->mol_num[mol]:   the number of this type of molecule in the system
    p_topol->atom_num   = (int *)my_malloc(p_topol->mol_types * sizeof(int));
    p_topol->mol_num    = (int *)my_malloc(p_topol->mol_types * sizeof(int));
    p_topol->atom_param = (AtomParam **)my_malloc(p_topol->mol_types * sizeof(AtomParam *));


    // van der Waals interaction parameters
    NonbondedParam *data_nonbonded = 
        (NonbondedParam *)my_malloc(p_topol->n_types * p_topol->n_types * sizeof(NonbondedParam));

    p_topol->nonbonded_param =
        (NonbondedParam **)my_malloc(p_topol->n_types * sizeof(NonbondedParam *));

    int i_type;
    for(i_type = 0; i_type < p_topol->n_types; ++ i_type)
    {
        p_topol->nonbonded_param[i_type] =
            &(data_nonbonded[p_topol->n_types * i_type]);
    }


    // bonded potentials: bond, pair, angle, dihedral
    int *data_bonded = (int *)my_malloc(sizeof(int) * p_topol->mol_types * 6);
    p_topol->n_bonds       = &(data_bonded[0]);
    p_topol->n_pairs       = &(data_bonded[p_topol->mol_types]);
    p_topol->n_angles      = &(data_bonded[p_topol->mol_types * 2]);
    p_topol->n_dihedrals   = &(data_bonded[p_topol->mol_types * 3]);
    p_topol->n_vsites      = &(data_bonded[p_topol->mol_types * 4]);
    p_topol->n_constraints = &(data_bonded[p_topol->mol_types * 5]);

    p_topol->vsite_funct = (int **)my_malloc(p_topol->mol_types * sizeof(int *)) ;

    p_topol->bond_param     = (BondParam **)my_malloc(p_topol->mol_types * sizeof(BondParam *));
    p_topol->pair_param     = (PairParam **)my_malloc(p_topol->mol_types * sizeof(PairParam *)) ;
    p_topol->angle_param    = (AngleParam **)my_malloc(p_topol->mol_types * sizeof(AngleParam *));
    p_topol->dihedral_param = (DihedralParam **)my_malloc(p_topol->mol_types * sizeof(DihedralParam *)) ;
    p_topol->vsite_4        = (VSite_4 **)my_malloc(p_topol->mol_types * sizeof(VSite_4 *)) ;
    p_topol->constraint     = (Constraint **)my_malloc(p_topol->mol_types * sizeof(Constraint *)) ;

    p_topol->exclude = (int ***)my_malloc(p_topol->mol_types * sizeof(int **)) ;


    // Quantum Sutton-Chen densities for metal
    if (p_metal->min >=0 && p_metal->max >= p_metal->min)
    {
        p_metal->inv_sqrt_dens = (double *)my_malloc(sizeof(double) * p_metal->num);
    }


    // read parameters from input_param, step 2
    read_param_2(input_param, p_topol, p_metal);

    int nAtoms = p_topol->n_atoms;
    int nMols  = p_topol->n_mols;


    // count number of virtual sites and constraints
    number_VSites = 0;
    number_Cstrs  = 0;
    for(mol = 0; mol < p_topol->mol_types; ++ mol)
    {
        number_VSites += p_topol->n_vsites[mol] * p_topol->mol_num[mol];
        number_Cstrs  += p_topol->n_constraints[mol] * p_topol->mol_num[mol];
    }


    // Gaussian distribution width for capacitance-polarizability model
    // see Mayer, Phys. Rev. B 2007, 75, 045407
    // and Jensen, J. Phys. Chem. C 2008, 112, 15697

    p_metal->inv_polar = 1.0 / p_metal->cpff_polar;
    p_metal->inv_capac = 1.0 / p_metal->cpff_capac;
    
    double R_q = sqrt(2.0 / M_PI) * p_metal->cpff_capac;
    double R_p = pow(sqrt(2.0 / M_PI) * p_metal->cpff_polar / 3.0, 1.0 / 3.0);

    p_metal->inv_R_qq = 1.0 / sqrt(R_q * R_q + R_q * R_q);
    p_metal->inv_R_pq = 1.0 / sqrt(R_p * R_p + R_q * R_q);
    p_metal->inv_R_pp = 1.0 / sqrt(R_p * R_p + R_p * R_p);


    // print info
    if (root_process == my_id) 
    {
        printf("\n");
        printf("              +-----------------------------------------------------+\n");
        printf("              |             CapacMD program version 1.0             |\n");
        printf("              |         Xin Li, TheoChemBio, KTH, Stockholm         |\n");
        printf("              +-----------------------------------------------------+\n");
        printf("\n");

        printf("              .------------------ reference paper ------------------.\n");
        printf("\n");
        printf("    Molecular Dynamics Simulations using a Capacitance-Polarizability Force Field,\n");
        printf("    Xin Li and Hans Agren, J. Phys. Chem. C, 2015, 119, 19430-19437.\n");
        printf("    http://pubs.acs.org/doi/abs/10.1021/acs.jpcc.5b04347\n");
        printf("\n");
        printf("\n");

        printf("    Calculation started at %s", ctime(&start_t));
        printf("    Parallelized via MPI, number of processors = %d\n", num_procs);
        printf("\n");
        printf("\n");

        printf("              .------------------ run parameters -------------------.\n");
        printf("\n");
        printf("    run_type = %s, ensemble = %s\n", p_runset->run_type, p_runset->ensemble);
        printf("    vdw_type = %s, coulomb_type = %s\n", p_runset->vdw_type, p_runset->coulomb_type);

        printf("    rCut = %.3f nm, ", rCut);
        if (1 == p_runset->use_coulomb)
        {
            printf("alpha = %.3f nm^-1", p_runset->w_alpha);
        }
        printf("\n");

        printf("    ref_T = %.1f K\n", p_runset->ref_temp);
        printf("\n");
        printf("\n");

        printf("              .------------------- molecule info -------------------.\n");
        printf("\n");
        printf("    There are %d types of molecules.\n", p_topol->mol_types);
        printf("\n");

        for(mol = 0; mol < p_topol->mol_types; ++ mol)
        {
            printf("    Molecule[%5d], num= %5d\n", mol, p_topol->mol_num[mol]);
            for(atom = 0; atom < p_topol->atom_num[mol]; ++ atom)
            {
                printf("    Atom[%5d], charge= %8.3f, mass= %8.3f, atomtype= %5d\n",
                       atom, p_topol->atom_param[mol][atom].charge,
                       p_topol->atom_param[mol][atom].mass, 
                       p_topol->atom_param[mol][atom].atomtype);
            }
            printf("\n");
        }
        printf("\n");
    }


    //=========== distribute molecules/atoms/metals among the procs ==============

    Task *p_task = (Task *)my_malloc(sizeof(Task));

    int *data_start_end = (int *)my_malloc(sizeof(int) * num_procs * 6);
    p_task->start_mol   = &(data_start_end[0]);
    p_task->end_mol     = &(data_start_end[num_procs]);
    p_task->start_atom  = &(data_start_end[num_procs * 2]);
    p_task->end_atom    = &(data_start_end[num_procs * 3]);
    p_task->start_metal = &(data_start_end[num_procs * 4]);
    p_task->end_metal   = &(data_start_end[num_procs * 5]);

    find_start_end(p_task->start_mol,   p_task->end_mol,   nMols,        num_procs);
    find_start_end(p_task->start_atom,  p_task->end_atom,  nAtoms,       num_procs);
    find_start_end(p_task->start_metal, p_task->end_metal, p_metal->num, num_procs);

    long int *data_long_start_end = (long int *)my_malloc(sizeof(long int) * num_procs * 2);
    p_task->start_pair = &(data_long_start_end[0]);
    p_task->end_pair   = &(data_long_start_end[num_procs]);

    long int n_pairs = nAtoms * (nAtoms - 1) / 2;
    find_start_end_long(p_task->start_pair,  p_task->end_pair,  n_pairs, num_procs);


    //============= assign indices, masses and charges =======================

    // For each atom in the system (i from 0 to p_topol->n_atoms-1)
    // atom_info[i].iAtom:  the index of this atom in its molecule type
    // atom_info[i].iMol:   the index of its molecule type
    // atom_info[i].molID:  the index of its molecule in the system

    // For each molecule in the system (im from 0 to p_topol->n_mols-1)
    // mol_info[im].mini:  the index of its first atom in the system
    // mol_info[im].maxi:  the index of its last atom in the system
    
    Atom_Info *atom_info = (Atom_Info *)my_malloc(nAtoms * sizeof(Atom_Info));
    Mol_Info  *mol_info  = (Mol_Info  *)my_malloc(nMols  * sizeof(Mol_Info));

    assign_indices(p_topol, p_metal, atom_info, mol_info);


    //======= allocate memory for coordinates, velocities and forces ==========

    System *p_system = (System *)my_malloc(sizeof(System));

    // potential energy
    // p_system->potential[0] = total energy
    // p_system->potential[1] = metal quantum Sutton-Chen energy
    // p_system->potential[2] = non-metal bond stretching energy
    // p_system->potential[3] = non-metal angle bending energy
    // p_system->potential[4] = non-metal torsional energy
    // p_system->potential[5] = 
    // p_system->potential[6] = Long range Coulomb
    // p_system->potential[7] = Coulomb energy (including 1-4)
    // p_system->potential[8] = vdW energy (including 1-4)
    // p_system->potential[9] = 
    // p_system->potential[10] = CPIM metal charge - non-metal charge
    // p_system->potential[11] = CPIM metal dipole - non-metal charge
    // p_system->potential[12] = CPIM metal charge - metal charge
    // p_system->potential[13] = CPIM metal charge - metal dipole
    // p_system->potential[14] = CPIM metal dipole - metal dipole
    
    double *data_potential = (double *)my_malloc(sizeof(double) * 15 * 2);
    p_system->potential   = &(data_potential[0]);
    p_system->partial_pot = &(data_potential[15]);
    p_system->old_potential = 0.0;


    // virial tensor
    p_system->virial      = (double **)my_malloc(DIM * sizeof(double*));
    p_system->partial_vir = (double **)my_malloc(DIM * sizeof(double*));

    double *data_vir = (double *)my_malloc(DIM*2 * DIM * sizeof(double));
    int i;
    for (i = 0; i < DIM; ++ i)
    {
        p_system->virial[i]      = &(data_vir[DIM * i]);
        p_system->partial_vir[i] = &(data_vir[DIM * (DIM + i)]);
    }


    // box size
    // Note: for now we treat rectangular box only.
    // "p_system->box" has six elements
    // the first three are length in x, y, z
    // the second three are half of the length in x, y, z
    p_system->box = (double *)my_malloc(sizeof(double) * DIM*2);


    // coordinates, velocities and forces
    double *data_rvf = (double *)my_malloc(nAtoms*7 * DIM * sizeof(double));
            
    p_system->rx = &(data_rvf[0]);
    p_system->ry = &(data_rvf[nAtoms]);
    p_system->rz = &(data_rvf[nAtoms*2]);
            
    p_system->vx = &(data_rvf[nAtoms*3]);
    p_system->vy = &(data_rvf[nAtoms*4]);
    p_system->vz = &(data_rvf[nAtoms*5]);
    
    p_system->fx = &(data_rvf[nAtoms*6]);
    p_system->fy = &(data_rvf[nAtoms*7]);
    p_system->fz = &(data_rvf[nAtoms*8]);
            
    // forces from slave processors
    p_system->partial_fx = &(data_rvf[nAtoms*9]);
    p_system->partial_fy = &(data_rvf[nAtoms*10]);
    p_system->partial_fz = &(data_rvf[nAtoms*11]);
            
    // old position for RATTLE constraints
    p_system->old_rx = &(data_rvf[nAtoms*12]);
    p_system->old_ry = &(data_rvf[nAtoms*13]);
    p_system->old_rz = &(data_rvf[nAtoms*14]);
            
    // old force for CG optimization
    p_system->old_fx = &(data_rvf[nAtoms*15]);
    p_system->old_fy = &(data_rvf[nAtoms*16]);
    p_system->old_fz = &(data_rvf[nAtoms*17]);

    // direction for CG optimization
    p_system->sx = &(data_rvf[nAtoms*18]);
    p_system->sy = &(data_rvf[nAtoms*19]);
    p_system->sz = &(data_rvf[nAtoms*20]);


    //================ read input gro file ==================
    
    // vQ and vP are "velocities" of the thermostat/barostat particles
    p_system->num_nhc = 10; // number of NH-chains hard-coded as 10
    p_system->vQ = (double *)my_malloc(sizeof(double) * p_system->num_nhc);
    p_system->aQ = (double *)my_malloc(sizeof(double) * p_system->num_nhc);
    int i_nhc;
    for (i_nhc = 0; i_nhc < p_system->num_nhc; ++ i_nhc)
    {
        p_system->vQ[i_nhc] = 0.0;
        p_system->aQ[i_nhc] = 0.0;
    }
    p_system->vP = 0.0;

    int groNAtoms;
    read_gro(input_gro, p_system, &groNAtoms, atom_info);

    if (groNAtoms != nAtoms)
    {
       printf("Error: groNAtoms(%d) not equal to nAtoms(%d)!\n", 
              groNAtoms, nAtoms);
       exit(1);
    }

    // half of box length for PBC
    p_system->box[3] = p_system->box[0] * 0.5;
    p_system->box[4] = p_system->box[1] * 0.5;
    p_system->box[5] = p_system->box[2] * 0.5;

    p_system->volume = p_system->box[0] * p_system->box[1] * p_system->box[2];
    p_system->inv_volume = 1.0 / p_system->volume;



    //================== set MD variables ==========================

    // degree of freedom
    p_system->ndf = 3 * (nAtoms - number_VSites) - 3 - number_Cstrs;

    // temperature coupling; kT = kB*T, in kJ mol^-1
    p_runset->kT    = K_BOLTZ * p_runset->ref_temp;
    p_system->qMass = (double)p_system->ndf * p_runset->kT * p_runset->tau_temp * p_runset->tau_temp;
    //p_system->pMass = (double)p_system->ndf * p_runset->kT * p_runset->tau_pres * p_runset->tau_pres;


    // temperature control
    p_system->first_temp = 0.0;
    p_system->ext_temp = 0.0;


    // time step
    p_runset->dt_2 = 0.5 * p_runset->dt;


    //=================== CPIM matrix and arrays =========================
    // Relay matrix
    // external electric field and potential
    // induced dipoles and charges
    
    double *data_relay = NULL;

    // dimension of matrix: 4M + n_NPs
    int n_mat = p_metal->num * 4 + p_metal->n_NPs;

    // initialize mat_relay, vec_ext and vec_pq
    if (p_metal->min >=0 && p_metal->max >= p_metal->min)
    {
        data_relay = (double *)my_malloc(sizeof(double) * n_mat * 3);

        p_metal->vec_pq  = &(data_relay[0]);
        p_metal->vec_ext = &(data_relay[n_mat]);
        p_metal->diag_relay = &(data_relay[n_mat * 2]);

        int i_mat;
        for(i_mat = 0; i_mat < n_mat; ++ i_mat)
        {
            p_metal->vec_pq[i_mat]  = 0.0;
            p_metal->vec_ext[i_mat] = 0.0;
            p_metal->diag_relay[i_mat] = 1.0;
        }
    }


    //=================== vectors for the BiCGSTAB solver =========================

    double* data_bicgstab = (double *)my_malloc(sizeof(double) * n_mat * 11);
    Bicgstab  *p_bicgstab = (Bicgstab *)my_malloc(sizeof(Bicgstab));

    p_bicgstab->Ax = &(data_bicgstab[0]);
    p_bicgstab->r0 = &(data_bicgstab[n_mat]);
    p_bicgstab->r  = &(data_bicgstab[n_mat * 2]);
    p_bicgstab->p  = &(data_bicgstab[n_mat * 3]);
    p_bicgstab->v  = &(data_bicgstab[n_mat * 4]);
    p_bicgstab->s  = &(data_bicgstab[n_mat * 5]);
    p_bicgstab->t  = &(data_bicgstab[n_mat * 6]);
    p_bicgstab->y  = &(data_bicgstab[n_mat * 7]);
    p_bicgstab->z  = &(data_bicgstab[n_mat * 8]);
    p_bicgstab->Kt = &(data_bicgstab[n_mat * 9]);
    p_bicgstab->K  = &(data_bicgstab[n_mat * 10]);


    //================== compute forces ==========================

    mpi_force(p_task, p_topol, atom_info, mol_info,
              p_runset, p_metal, p_system, p_bicgstab,
              my_id, num_procs, time_used);


    //========== file handles: gro, vec_pq, binary dat, parameters ========
    
    FILE *file_gro, *file_pq, *file_dat;
    file_gro = NULL;
    file_pq  = NULL;
    file_dat = NULL;


    //========== adjust velocities and write trajectories =================

    if (root_process == my_id) 
    {
        // sum potential energies
        sum_potential(p_system->potential);

        // update temperature
        remove_comm(nAtoms, atom_info, p_system);
        kinetic_energy(p_system, nAtoms, atom_info);

        // save the starting temperature
        p_system->first_temp = p_system->inst_temp;
        p_system->ext_temp   = p_runset->ref_temp;

        // creat traj.gro for writing
        file_gro = fopen("traj.gro","w") ;
        if (NULL == file_gro) 
        {
            printf( "Cannot write to traj.gro!\n" ) ;
            exit(1);
        }

        // creat vec_pq.txt for writing
        file_pq = fopen("vec_pq.txt","w") ;
        if (NULL == file_pq) 
        {
            printf( "Cannot write to vec_pq.txt!\n" ) ;
            exit(1);
        }

        // creat traj.dat for writing
        file_dat = fopen("traj.dat","w") ;
        if (NULL == file_dat) 
        {
            printf( "Cannot write to traj.dat!\n" ) ;
            exit(1);
        }

        // get maximal force
        get_fmax_rms(nAtoms, p_system);

        // write the starting geometry to gro and dat file
        write_gro(file_gro, p_system, nAtoms, atom_info, 0);
        write_binary(file_dat, p_system, nAtoms, 0);

        if (p_metal->min >=0 && p_metal->max >= p_metal->min)
        {
            write_vec_pq(file_pq, p_metal, 0);
        }


        printf("              .---------------- start MD calculation ---------------.\n");
        printf("\n");

        // check p_runset->run_type
        if (0 == strcmp(p_runset->run_type, "em") || 
            0 == strcmp(p_runset->run_type, "cg"))
        {
            printf("    Step %-5d Fmax=%10.3e  E=%15.8e\n", 
                    0, p_system->f_max, p_system->potential[0]);
        }
        else if (0 == strcmp(p_runset->run_type, "md"))
        {
            printf("  %10.3f  Fmax=%.2e  E=%.6e  T=%.3e\n", 
                   0.0, p_system->f_max, p_system->potential[0], 
                   p_system->inst_temp);
        }
    }


    //========= Now start MD steps ============================

    int step = 0;

    //===================================================
    // Energy minimization using steepest descent or CG
    //===================================================

    if (0 == strcmp(p_runset->run_type, "em") || 
        0 == strcmp(p_runset->run_type, "cg"))
    {
        p_system->vQ[0] = 0.0;
        p_system->vP = 0.0;

        int converged = 0;
        double gamma = 0.0;
        double delta_pot;

        // initialize direction sx,sy,sz
        for(i = 0; i < nAtoms; ++ i)
        {
            p_system->old_fx[i] = p_system->fx[i];
            p_system->old_fy[i] = p_system->fy[i];
            p_system->old_fz[i] = p_system->fz[i];

            p_system->sx[i] = p_system->fx[i];
            p_system->sy[i] = p_system->fy[i];
            p_system->sz[i] = p_system->fz[i];
        }

        for(step = 1; step <= p_runset->em_steps && 0 == converged; ++ step) 
        {
            // update coordinates on root processor
            if (root_process == my_id) 
            {
                p_system->old_potential = p_system->potential[0];

                for(i = 0; i < nAtoms; ++ i)
                {
                    p_system->old_rx[i] = p_system->rx[i];
                    p_system->old_ry[i] = p_system->ry[i];
                    p_system->old_rz[i] = p_system->rz[i];

                    // fix metal coordinates?
                    if (1 == p_metal->fix_pos && 1 == atom_info[i].is_metal)
                    {
                        continue;
                    }

                    p_system->rx[i] += p_system->sx[i] / p_system->f_max * p_runset->em_length;
                    p_system->ry[i] += p_system->sy[i] / p_system->f_max * p_runset->em_length;
                    p_system->rz[i] += p_system->sz[i] / p_system->f_max * p_runset->em_length;
                }

                // apply constraints
                rattle_1st(p_runset->dt, mol_info, atom_info, p_topol, p_system);

                // zero velocities
                for(i = 0; i < nAtoms; ++ i)
                {
                    p_system->vx[i] = 0.0;
                    p_system->vy[i] = 0.0;
                    p_system->vz[i] = 0.0;
                }
            }

            // update forces
            mpi_force(p_task, p_topol, atom_info, mol_info,
                      p_runset, p_metal, p_system, p_bicgstab,
                      my_id, num_procs, time_used);

            if (root_process == my_id) 
            {
                // check potential and fmax on root processor
                sum_potential(p_system->potential);
                get_fmax_rms(nAtoms, p_system);
                delta_pot = p_system->potential[0] - p_system->old_potential;

                if (delta_pot <= 0.0) 
                { 
                    p_runset->em_length *= 1.2; 
                }
                else 
                { 
                    p_runset->em_length *= 0.2; 
                }
                
                // print info and write to the gro file
                printf("    Step %-5d Fmax=%10.3e  E=%15.8e\n", 
                        step, p_system->f_max, p_system->potential[0]);

                // write trajectories
                write_gro(file_gro, p_system, nAtoms, atom_info, step);
                write_binary(file_dat, p_system, nAtoms, step);

                if (p_metal->min >=0 && p_metal->max >= p_metal->min)
                {
                    write_vec_pq(file_pq, p_metal, step);
                }

                // check convergence
                if (p_system->f_max < p_runset->em_tol && 
                    p_system->f_rms < p_runset->em_tol * 0.5 && 
                    fabs(delta_pot) < p_runset->em_tol * 0.1)
                {
                    printf("\n");
                    printf("    F_max   (%13.6e) smaller than %e\n", 
                            p_system->f_max, p_runset->em_tol);
                    printf("    F_rms   (%13.6e) smaller than %e\n", 
                            p_system->f_rms, p_runset->em_tol * 0.5);
                    printf("    delta_E (%13.6e) smaller than %e\n", 
                            delta_pot, p_runset->em_tol * 0.1);
                    printf("\n");
                    printf("    ===========  Optimization converged ============\n");
                    converged = 1;
                }
                else
                {
                    // update gamma for CG optimization
                    // for steep descent, gamma = 0.0
                    if (0 == strcmp(p_runset->run_type, "cg"))
                    {
                        double g22 = 0.0;
                        double g12 = 0.0;
                        double g11 = 0.0;
                        for(i = 0; i < nAtoms; ++ i)
                        {
                            g22 += p_system->fx[i] * p_system->fx[i] + 
                                   p_system->fy[i] * p_system->fy[i] + 
                                   p_system->fz[i] * p_system->fz[i];
                            g12 += p_system->old_fx[i] * p_system->fx[i] + 
                                   p_system->old_fy[i] * p_system->fy[i] + 
                                   p_system->old_fz[i] * p_system->fz[i];
                            g11 += p_system->old_fx[i] * p_system->old_fx[i] + 
                                   p_system->old_fy[i] * p_system->old_fy[i] + 
                                   p_system->old_fz[i] * p_system->old_fz[i];
                        }
                        gamma = (g22 - g12) / g11;
                    }

                    for(i = 0; i < nAtoms; ++ i)
                    {
                        p_system->sx[i] = p_system->fx[i] + gamma * p_system->sx[i];
                        p_system->sy[i] = p_system->fy[i] + gamma * p_system->sy[i];
                        p_system->sz[i] = p_system->fz[i] + gamma * p_system->sz[i];

                        p_system->old_fx[i] = p_system->fx[i];
                        p_system->old_fy[i] = p_system->fy[i];
                        p_system->old_fz[i] = p_system->fz[i];
                    }
                }

                // communicate convergence
                for(an_id = 1; an_id < num_procs; ++ an_id) 
                {
                    ierr = MPI_Send(&converged, 1, MPI_INT, an_id, tag_99, MPI_COMM_WORLD);
                }
            }
            else
            {
                ierr = MPI_Recv(&converged, 1, MPI_INT, root_process, tag_99, 
                                MPI_COMM_WORLD, &status);
            }

            // exit loop if converged
            if (1 == converged) { break; }
        }
    }


    //===============================
    // MD with nvt ensemble
    //===============================

    else if (0 == strcmp(p_runset->run_type, "md"))
    {
        for (step = 1; step <= p_runset->nSteps; ++ step) 
        {
            if (root_process == my_id) 
            {
                // gradually increase the reference temperature
                if (step < p_runset->nHeating)
                {
                    p_runset->ref_temp = 
                        p_system->first_temp + 
                        (p_system->ext_temp - p_system->first_temp) * step / p_runset->nHeating;
                }
                else
                {
                    p_runset->ref_temp = p_system->ext_temp;
                }

                // update kT accordingly
                p_runset->kT = K_BOLTZ * p_runset->ref_temp;
            
                // thermostat for 1st half step
                if (0 == strcmp(p_runset->ensemble, "nvt"))
                {
                    nose_hoover_chain(p_runset, p_system, nAtoms);
                }

                // update velocity for 1st half step
                for(i = 0; i < nAtoms; ++ i)
                {
                    // fix metal coordinates?
                    if (1 == p_metal->fix_pos && 1 == atom_info[i].is_metal)
                    {
                        continue;
                    }

                    // for virtual sites, inv_mass == 0.0;
                    double inv_mass = atom_info[i].inv_mass;
                    p_system->vx[i] += p_runset->dt_2 * p_system->fx[i] * inv_mass;
                    p_system->vy[i] += p_runset->dt_2 * p_system->fy[i] * inv_mass;
                    p_system->vz[i] += p_runset->dt_2 * p_system->fz[i] * inv_mass;
                }

                // update position for the whole time step
                // using velocity at half time step
                for(i = 0; i < nAtoms; ++ i)
                {
                    // fix metal coordinates?
                    if (1 == p_metal->fix_pos && 1 == atom_info[i].is_metal)
                    {
                        continue;
                    }

                    p_system->old_rx[i] = p_system->rx[i];
                    p_system->old_ry[i] = p_system->ry[i];
                    p_system->old_rz[i] = p_system->rz[i];

                    p_system->rx[i] += p_runset->dt * p_system->vx[i];
                    p_system->ry[i] += p_runset->dt * p_system->vy[i];
                    p_system->rz[i] += p_runset->dt * p_system->vz[i];
                }

                // apply constraints for the 1st half
                rattle_1st(p_runset->dt, mol_info, atom_info, p_topol, p_system);
            }


            // compute forces
            mpi_force(p_task, p_topol, atom_info, mol_info,
                      p_runset, p_metal, p_system, p_bicgstab,
                      my_id, num_procs, time_used);
            

            // update velocities
            if (root_process == my_id) 
            {
                // sum potential energies
                sum_potential(p_system->potential);

                // update velocity for 2nd half step
                for(i = 0; i < nAtoms; ++ i)
                {
                    // fix metal coordinates?
                    if (1 == p_metal->fix_pos && 1 == atom_info[i].is_metal)
                    {
                        continue;
                    }

                    // for virtual sites, inv_mass == 0.0;
                    double inv_mass = atom_info[i].inv_mass;
                    p_system->vx[i] += p_runset->dt_2 * p_system->fx[i] * inv_mass;
                    p_system->vy[i] += p_runset->dt_2 * p_system->fy[i] * inv_mass;
                    p_system->vz[i] += p_runset->dt_2 * p_system->fz[i] * inv_mass;
                }

                // apply constraints for the 2nd half
                rattle_2nd(p_runset->dt, mol_info, atom_info, p_topol, p_system);

                // update temperature
                kinetic_energy(p_system, nAtoms, atom_info);


                // thermostat for 2nd half step
                if (0 == strcmp(p_runset->ensemble, "nvt"))
                {
                    nose_hoover_chain(p_runset, p_system, nAtoms);
                }


                // remove center of mass motion and update temperature
                remove_comm(nAtoms, atom_info, p_system);
                kinetic_energy(p_system, nAtoms, atom_info);


                // print information for this step
                if (0 == step % p_runset->nSave)
                {
                    // apply PBC
                    apply_pbc(nMols, mol_info, p_system->rx, p_system->ry, p_system->rz, p_system->box);
             
                    get_fmax_rms(nAtoms, p_system);
                    printf("  %10.3f  Fmax=%.2e  E=%.6e  T=%.3e (%.1f)\n", 
                           p_runset->dt * step, p_system->f_max, p_system->potential[0], 
                           p_system->inst_temp, p_runset->ref_temp);

                    // write to dat file
                    write_binary(file_dat, p_system, nAtoms, step);
                }
             
                // regularly write to gro file
                if (0 == step % (p_runset->nSave*10))
                {
                    write_gro(file_gro, p_system, nAtoms, atom_info, step);

                    if (p_metal->min >=0 && p_metal->max >= p_metal->min)
                    {
                        write_vec_pq(file_pq, p_metal, step);
                    }
                }
            }

        }
    }


    //=========== Finalize MD ====================
    
    sum_time_used(time_used, my_id, num_procs);

    if (root_process == my_id) 
    {
        fclose(file_dat);
        fclose(file_pq);
        fclose(file_gro);

#ifdef DEBUG
        int i;
        for(i = 0; i < nAtoms; ++ i)
        {
            printf("    f[%5d]= %12.5e, %12.5e, %12.5e\n", 
                    i, p_system->fx[i], p_system->fy[i], p_system->fz[i]);
        }
#endif

        printf("\n");
        printf("\n");
        printf("              .------------- final potential energies --------------.\n");
        printf("\n");

        print_potential(p_system->potential);

        time_t end_t = time(NULL);
        gettimeofday(&tv, NULL);
        double end_time = (tv.tv_sec) + (tv.tv_usec) * 1.0e-6;

        printf("    Calculation ended normally at %s", ctime(&end_t));
        printf("    %.3f seconds were used\n", end_time - start_time );
        printf("\n");
        printf("\n");

        printf("              .------------------- time usage ----------------------.\n");
        printf("\n");

        analyze_time_used(time_used, num_procs);
        printf("\n");
    }

    // free arrays
    for(an_id = 0; an_id < num_procs; ++ an_id)
    {
        free(time_used[an_id]);
    }
    free(time_used);

    free(p_runset->external_efield);

    free(p_metal->cpff_chg);
    free(p_metal->start_NP);
    free(p_metal->end_NP);

    free(data_bonded);
    data_bonded = NULL;

    for(mol = 0; mol < p_topol->mol_types; ++ mol) 
    {
        int atom_i;
        for(atom_i = 0; atom_i < p_topol->atom_num[mol]; ++ atom_i) 
        {
            free(p_topol->exclude[mol][atom_i]);
        }

        free(p_topol->exclude[mol]);

        free(p_topol->bond_param[mol]);
        free(p_topol->pair_param[mol]);
        free(p_topol->angle_param[mol]);
        free(p_topol->dihedral_param[mol]);

        free(p_topol->vsite_4[mol]);
        free(p_topol->vsite_funct[mol]);

        free(p_topol->constraint[mol]);
    }

    free(p_topol->exclude);

    free(p_topol->bond_param);
    free(p_topol->pair_param);
    free(p_topol->angle_param);
    free(p_topol->dihedral_param);
    p_topol->bond_param = NULL;
    p_topol->pair_param = NULL;
    p_topol->angle_param = NULL;
    p_topol->dihedral_param = NULL;

    free(p_topol->vsite_4);
    free(p_topol->vsite_funct);
    p_topol->vsite_4 = NULL;
    p_topol->vsite_funct = NULL;

    free(p_topol->constraint);
    p_topol->constraint = NULL;

    free(p_topol->mol_num);
    free(p_topol->atom_num);
    for(mol = 0; mol < p_topol->mol_types; ++ mol) 
    {
        free(p_topol->atom_param[mol]);
    }
    free(p_topol->atom_param);

    if (p_metal->min >=0 && p_metal->max >= p_metal->min)
    {
        free(p_metal->inv_sqrt_dens);
        p_metal->inv_sqrt_dens = NULL;

        free(data_relay);
        //free(p_metal->mat_relay);
        data_relay = NULL;
        //p_metal->mat_relay = NULL;
        p_metal->vec_pq    = NULL;
        p_metal->vec_ext   = NULL;
        p_metal->diag_relay = NULL;
    }

    free(data_bicgstab);


    free(data_start_end);
    free(data_long_start_end);
    data_start_end      = NULL;
    data_long_start_end = NULL;

    p_task->start_mol   = NULL;
    p_task->end_mol     = NULL;
    p_task->start_atom  = NULL;
    p_task->end_atom    = NULL;
    p_task->start_metal = NULL;
    p_task->end_metal   = NULL;
    p_task->start_pair  = NULL;
    p_task->end_pair    = NULL;

    free(atom_info);
    free(mol_info);
    atom_info = NULL;
    mol_info  = NULL;

    free(data_potential);
    data_potential = NULL;
    p_system->potential   = NULL;
    p_system->partial_pot = NULL;


    free(p_system->virial);
    free(p_system->partial_vir);
    free(data_vir);


    free(p_system->box);
    p_system->box = NULL;

    free(data_rvf);
    data_rvf = NULL;
    p_system->rx = NULL;
    p_system->ry = NULL;
    p_system->rz = NULL;
    p_system->vx = NULL;
    p_system->vy = NULL;
    p_system->vz = NULL;
    p_system->fx = NULL;
    p_system->fy = NULL;
    p_system->fz = NULL;
    p_system->partial_fx = NULL;
    p_system->partial_fy = NULL;
    p_system->partial_fz = NULL;
    p_system->old_rx = NULL;
    p_system->old_ry = NULL;
    p_system->old_rz = NULL;
    p_system->old_fx = NULL;
    p_system->old_fy = NULL;
    p_system->old_fz = NULL;
    p_system->sx = NULL;
    p_system->sy = NULL;
    p_system->sz = NULL;

    free(p_topol->nonbonded_param);
    free(data_nonbonded);
    p_topol->nonbonded_param = NULL;
    data_nonbonded = NULL;

    free(p_task);
    free(p_system);
    free(p_topol);
    free(p_metal);
    free(p_runset);
    p_task   = NULL;
    p_system = NULL;
    p_topol  = NULL;
    p_metal  = NULL;
    p_runset = NULL;

    ierr = MPI_Finalize();

    if (ierr) {}

    return 0;
}

