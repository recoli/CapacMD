/**
 \file main.cpp
 \author Xin Li
 \date 2015/09
 \version 1.0.2
 
 \brief main engine of the CapacMD program
 
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
#include <sys/time.h>
#include <string>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "typedef.hpp"
#include "file.hpp"
#include "force.hpp"
#include "thermostat.hpp"

int main(int argc, char *argv[])
{
    //============= set up MPI =================
    
    MPI::Status status;

    // initialize MPI, get my_id for this process and total number of processes
    MPI::Init(argc, argv);
    int my_id     = MPI::COMM_WORLD.Get_rank();
    int num_procs = MPI::COMM_WORLD.Get_size();

    // tag for MPI message: energy minimization convergence 
    const int tag_99 = 99;

    
    // initialize timer
    time_t start_t = time(NULL);
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double start_time = (tv.tv_sec) + (tv.tv_usec) * 1.0e-6;
    

    //============ get input files =============

    // default filenames
    std::string input_gro = "init.gro";
    std::string input_param = "param.txt";
    std::string input_mdset = "mdset.txt";

    // parse arguments (input files)
    // -gro for gro file
    // -par for parameter file
    // -mds for MD settings
    int iarg = 1;
    while (1)
    {
        if (iarg >= argc) { break; }

        if ((std::string("-gro") == argv[iarg]) && (iarg < argc - 1))
        {
            input_gro = argv[iarg + 1];
            iarg += 2;
        }
        else if ((std::string("-par") == argv[iarg]) && (iarg < argc - 1))
        {
            input_param = argv[iarg + 1];
            iarg += 2;
        }
        else if ((std::string("-mds") == argv[iarg]) && (iarg < argc - 1))
        {
            input_mdset = argv[iarg + 1];
            iarg += 2;
        }
        else
        {
            ++ iarg;
        }
    }

    if (ROOT_PROC == my_id) {
        printf("\n");
        printf("              +-----------------------------------------------------+\n");
        printf("              |            CapacMD program version 1.0.2            |\n");
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
        
        printf("    Job started at %s", ctime(&start_t));
        printf("    Parallelized via MPI, number of processors = %d\n", num_procs);
        printf("\n");
        printf("\n");
    }

    //========= define variables and read MD settings ===========

    // varialbles to read from input_mdset
    RunSet s_runset;
    Metal  s_metal;
    
    // read md settings from input_mdset
    read_settings(input_mdset, s_runset, s_metal, my_id);


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
    
    double **time_used = new (std::nothrow) double* [num_procs];
    assert(time_used != nullptr);
    for (int an_id = 0; an_id < num_procs; ++ an_id)
    {
        time_used[an_id] = new (std::nothrow) double [N_TIMER];
        assert(time_used[an_id] != nullptr);

        for (int it = 0; it < N_TIMER; ++ it)
        {
            time_used[an_id][it] = 0.0;
        }
    }


    // Coulomb type: cut_off, wolf_sum
    if (std::string("cut_off") == s_runset.coulomb_type)
    { 
        s_runset.use_coulomb = 0; 
    }
    else if (std::string("wolf_sum") == s_runset.coulomb_type)
    { 
        s_runset.use_coulomb = 1; 
    }
    else 
    { 
        printf("Error: unknown coulomb_type %s!\n", s_runset.coulomb_type.c_str());
        exit(-1); 
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

    // note: s_runset.rCut and s_runset.w_alpha were read from mdset.txt
    s_runset.rCut2 = s_runset.rCut * s_runset.rCut;
    double rCut    = s_runset.rCut;
    double rCut2   = s_runset.rCut2;

    double w_alpha        = s_runset.w_alpha;
    s_runset.w_a2_sqrtPI  = w_alpha * 2.0 * INV_SQRT_PI;
    s_runset.w_Const      = erfc(w_alpha * rCut) / rCut2 + 
                            s_runset.w_a2_sqrtPI *
                            exp(-w_alpha * w_alpha * rCut2) / rCut;
    s_runset.w_erfc_arCut = erfc(w_alpha * rCut)/ rCut;


    // van der Waals type: cut_off or shifted
    if (std::string("cut_off") == s_runset.vdw_type)
    {
        s_runset.use_vdw = 0; 
    }
    else if (std::string("shifted") == s_runset.vdw_type)
    { 
        s_runset.use_vdw = 1; 
    }
    else 
    { 
        printf("Error: unknown vdw_type %s!\n", s_runset.vdw_type.c_str());
        exit(-1); 
    }


    // Shifted force method for Lennard-Jones potential
    // Ref:   Toxvaerd and Dyre, J. Chem. Phys. 2011, 134, 081102
    //        dx.doi.org/10.1063/1.3558787

    s_runset.inv_rc12 = 1.0 / pow(rCut, 12.0);
    s_runset.inv_rc6  = 1.0 / pow(rCut, 6.0);


    //============== read force field parameters from param.txt =================
    
    Topol s_topol;


    // read parameters from input_param, step 1
    read_param_1(input_param, s_topol, s_metal);
    

    // CPIM: capacitance-polarizability interaction model
    // Ref.: a) Jensen and Jensen, J. Phys. Chem. C, 2008, 112, 15697-15703
    //          dx.doi.org/10.1021/jp804116z
    //       b) Morton and Jensen, J. Chem. Phys., 2010, 133, 074103
    //          dx.doi.org/10.1063/1.3457365
    //
    // s_metal.cpff_polar: polarizability
    // s_metal.cpff_capac: capacitance
    // s_metal.n_NPs:      number of nanoparticles
    // s_metal.cpff_chg:   total charge of a nanoparticle
    // s_metal.start_NP:   first atom of a nanoparticle
    // s_metal.end_NP:     last atom of a nanoparticle
    
    // CPIM charge and indices
    s_metal.cpff_chg = new (std::nothrow) double [s_metal.n_NPs];
    s_metal.start_NP = new (std::nothrow) int [s_metal.n_NPs];
    s_metal.end_NP   = new (std::nothrow) int [s_metal.n_NPs];
    assert(s_metal.cpff_chg != nullptr);
    assert(s_metal.start_NP != nullptr);
    assert(s_metal.end_NP   != nullptr);

    // molecules and atom parameters
    // For a given molecule type (mol from 0 to s_topol.mol_types-1):
    // s_topol.atom_num[mol]:  its number of atoms
    // s_topol.mol_num[mol]:   the number of this type of molecule in the system
    s_topol.atom_num   = new (std::nothrow) int [s_topol.mol_types];
    s_topol.mol_num    = new (std::nothrow) int [s_topol.mol_types];
    s_topol.atom_param = new (std::nothrow) AtomParam* [s_topol.mol_types];
    assert(s_topol.atom_num   != nullptr);
    assert(s_topol.mol_num    != nullptr);
    assert(s_topol.atom_param != nullptr);

    
    // variables for bonded and nonbonded interaction

    // van der Waals interaction parameters
    NonbondedParam *data_nonbonded = new (std::nothrow) NonbondedParam
                                     [s_topol.n_types * s_topol.n_types];
    assert(data_nonbonded != nullptr);

    s_topol.nonbonded_param = new (std::nothrow) NonbondedParam* [s_topol.n_types];
    assert(s_topol.nonbonded_param != nullptr);

    for (int i_type = 0; i_type < s_topol.n_types; ++ i_type)
    {
        s_topol.nonbonded_param[i_type] =
            &(data_nonbonded[s_topol.n_types * i_type]);
    }

    // bonded potentials: bond, pair, angle, dihedral
    int *data_bonded = new (std::nothrow) int [s_topol.mol_types * 6];
    assert(data_bonded != nullptr);
    
    s_topol.n_bonds       = &(data_bonded[0]);
    s_topol.n_pairs       = &(data_bonded[s_topol.mol_types]);
    s_topol.n_angles      = &(data_bonded[s_topol.mol_types * 2]);
    s_topol.n_dihedrals   = &(data_bonded[s_topol.mol_types * 3]);
    s_topol.n_vsites      = &(data_bonded[s_topol.mol_types * 4]);
    s_topol.n_constraints = &(data_bonded[s_topol.mol_types * 5]);
    
    s_topol.bond_param     = new (std::nothrow) BondParam*     [s_topol.mol_types];
    s_topol.pair_param     = new (std::nothrow) PairParam*     [s_topol.mol_types];
    s_topol.angle_param    = new (std::nothrow) AngleParam*    [s_topol.mol_types];
    s_topol.dihedral_param = new (std::nothrow) DihedralParam* [s_topol.mol_types];
    assert(s_topol.bond_param     != nullptr);
    assert(s_topol.pair_param     != nullptr);
    assert(s_topol.angle_param    != nullptr);
    assert(s_topol.dihedral_param != nullptr);

    // virtual sites, constraints and exclusions
    s_topol.vsite_funct    = new (std::nothrow) int*        [s_topol.mol_types];
    s_topol.vsite_4        = new (std::nothrow) VSite_4*    [s_topol.mol_types];
    s_topol.constraint     = new (std::nothrow) Constraint* [s_topol.mol_types];
    assert(s_topol.vsite_funct    != nullptr);
    assert(s_topol.vsite_4        != nullptr);
    assert(s_topol.constraint     != nullptr);

    s_topol.exclude = new (std::nothrow) int** [s_topol.mol_types];
    assert(s_topol.exclude != nullptr);

    /// \todo will improve the "have_metal" checking
    // Quantum Sutton-Chen densities for metal
    if (s_metal.min >= 0 && s_metal.max >= s_metal.min)
    {
        s_metal.inv_sqrt_dens = new (std::nothrow) double [s_metal.num];
        assert(s_metal.inv_sqrt_dens != nullptr);
    }

    // read parameters from input_param, step 2
    read_param_2(input_param, s_topol, s_metal);

    int nAtoms = s_topol.n_atoms;
    int nMols  = s_topol.n_mols;

    // count number of virtual sites and constraints
    int number_VSites = 0;
    int number_Cstrs  = 0;
    for (int mol = 0; mol < s_topol.mol_types; ++ mol)
    {
        number_VSites += s_topol.n_vsites[mol]      * s_topol.mol_num[mol];
        number_Cstrs  += s_topol.n_constraints[mol] * s_topol.mol_num[mol];
    }

    // Gaussian distribution width for capacitance-polarizability model
    // see Mayer, Phys. Rev. B 2007, 75, 045407
    // and Jensen, J. Phys. Chem. C 2008, 112, 15697
    s_metal.inv_polar = 1.0 / s_metal.cpff_polar;
    s_metal.inv_capac = 1.0 / s_metal.cpff_capac;
    
    double R_q = sqrt(2.0 / M_PI) * s_metal.cpff_capac;
    double R_p = pow(sqrt(2.0 / M_PI) * s_metal.cpff_polar / 3.0, 1.0 / 3.0);

    s_metal.inv_R_qq = 1.0 / sqrt(R_q * R_q + R_q * R_q);
    s_metal.inv_R_pq = 1.0 / sqrt(R_p * R_p + R_q * R_q);
    s_metal.inv_R_pp = 1.0 / sqrt(R_p * R_p + R_p * R_p);


    //=========== print program version, reference and running information ==============

    if (ROOT_PROC == my_id)
    {
        printf("              .------------------ run parameters -------------------.\n");
        printf("\n");
        printf("    run_type = %s, ensemble = %s, ",
               s_runset.run_type.c_str(), s_runset.ensemble.c_str());
        printf("vdw_type = %s, coulomb_type = %s,\n",
               s_runset.vdw_type.c_str(), s_runset.coulomb_type.c_str());

        printf("    rCut = %.3f nm, ", rCut);
        if (1 == s_runset.use_coulomb)
        {
            printf("alpha = %.3f nm^-1, ", s_runset.w_alpha);
        }

        printf("ref_T = %.1f K\n", s_runset.ref_temp);
        printf("\n");
        printf("\n");

        printf("              .------------------- molecule info -------------------.\n");
        printf("\n");
        printf("    There are %d types of molecules.\n", s_topol.mol_types);
        printf("\n");

        for (int mol = 0; mol < s_topol.mol_types; ++ mol)
        {
            printf("    Molecule[%5d], num= %5d\n", mol, s_topol.mol_num[mol]);
            for (int atom = 0; atom < s_topol.atom_num[mol]; ++ atom)
            {
                printf("    Atom[%5d], charge= %8.3f, mass= %8.3f, atomtype= %5d\n",
                       atom, s_topol.atom_param[mol][atom].charge,
                       s_topol.atom_param[mol][atom].mass, 
                       s_topol.atom_param[mol][atom].atomtype);
            }
            printf("\n");
        }
        printf("\n");
    }


    //=========== distribute molecules/atoms/metals among the procs ==============

    Task s_task;

    int *data_start_end = new (std::nothrow) int [num_procs * 6];
    assert(data_start_end != nullptr);
    s_task.start_mol   = &(data_start_end[0]);
    s_task.end_mol     = &(data_start_end[num_procs]);
    s_task.start_atom  = &(data_start_end[num_procs * 2]);
    s_task.end_atom    = &(data_start_end[num_procs * 3]);
    s_task.start_metal = &(data_start_end[num_procs * 4]);
    s_task.end_metal   = &(data_start_end[num_procs * 5]);

    find_start_end(s_task.start_mol,   s_task.end_mol,   nMols,       num_procs);
    find_start_end(s_task.start_atom,  s_task.end_atom,  nAtoms,      num_procs);
    find_start_end(s_task.start_metal, s_task.end_metal, s_metal.num, num_procs);

    long int *data_long_start_end = new (std::nothrow) long int [num_procs * 2];
    assert(data_long_start_end != nullptr);
    s_task.start_pair = &(data_long_start_end[0]);
    s_task.end_pair   = &(data_long_start_end[num_procs]);

    long int n_pairs = nAtoms * (nAtoms - 1) / 2;
    find_start_end(s_task.start_pair,  s_task.end_pair,  n_pairs, num_procs);


    //============= assign indices, masses and charges =======================

    // For each atom in the system (i from 0 to s_topol.n_atoms-1)
    // atom_info[i].iAtom:  the index of this atom in its molecule type
    // atom_info[i].iMol:   the index of its molecule type
    // atom_info[i].molID:  the index of its molecule in the system

    // For each molecule in the system (im from 0 to s_topol.n_mols-1)
    // mol_info[im].mini:  the index of its first atom in the system
    // mol_info[im].maxi:  the index of its last atom in the system
    
    Atom_Info *atom_info = new (std::nothrow) Atom_Info [nAtoms];
    Mol_Info  *mol_info  = new (std::nothrow) Mol_Info  [nMols];
    assert(atom_info != nullptr);
    assert(mol_info  != nullptr);

    assign_indices(s_topol, s_metal, atom_info, mol_info);


    //======= allocate memory for coordinates, velocities and forces ==========

    System s_system;

    // potential energy
    // s_system.potential[0] = total energy
    // s_system.potential[1] = metal quantum Sutton-Chen energy
    // s_system.potential[2] = non-metal bond stretching energy
    // s_system.potential[3] = non-metal angle bending energy
    // s_system.potential[4] = non-metal torsional energy
    // s_system.potential[5] = 
    // s_system.potential[6] = Long range Coulomb
    // s_system.potential[7] = Coulomb energy (including 1-4)
    // s_system.potential[8] = vdW energy (including 1-4)
    // s_system.potential[9] = 
    // s_system.potential[10] = CPIM metal charge - non-metal charge
    // s_system.potential[11] = CPIM metal dipole - non-metal charge
    // s_system.potential[12] = CPIM metal charge - metal charge
    // s_system.potential[13] = CPIM metal charge - metal dipole
    // s_system.potential[14] = CPIM metal dipole - metal dipole
    
    double *data_potential = new (std::nothrow) double [N_POT * 2];
    assert(data_potential != nullptr);
    s_system.potential   = &(data_potential[0]);
    s_system.partial_pot = &(data_potential[N_POT]);
    s_system.old_potential = 0.0;


    // coordinates, velocities and forces
    double *data_rvf = new (std::nothrow) double [nAtoms * DIM * 7];
    assert(data_rvf != nullptr);
            
    s_system.rx = &(data_rvf[0]);
    s_system.ry = &(data_rvf[nAtoms]);
    s_system.rz = &(data_rvf[nAtoms*2]);
            
    s_system.vx = &(data_rvf[nAtoms*3]);
    s_system.vy = &(data_rvf[nAtoms*4]);
    s_system.vz = &(data_rvf[nAtoms*5]);
    
    s_system.fx = &(data_rvf[nAtoms*6]);
    s_system.fy = &(data_rvf[nAtoms*7]);
    s_system.fz = &(data_rvf[nAtoms*8]);
            
    // forces from slave processors
    s_system.partial_fx = &(data_rvf[nAtoms*9]);
    s_system.partial_fy = &(data_rvf[nAtoms*10]);
    s_system.partial_fz = &(data_rvf[nAtoms*11]);
            
    // old position for RATTLE constraints
    s_system.old_rx = &(data_rvf[nAtoms*12]);
    s_system.old_ry = &(data_rvf[nAtoms*13]);
    s_system.old_rz = &(data_rvf[nAtoms*14]);
            
    // old force for CG optimization
    s_system.old_fx = &(data_rvf[nAtoms*15]);
    s_system.old_fy = &(data_rvf[nAtoms*16]);
    s_system.old_fz = &(data_rvf[nAtoms*17]);

    // direction for CG optimization
    s_system.sx = &(data_rvf[nAtoms*18]);
    s_system.sy = &(data_rvf[nAtoms*19]);
    s_system.sz = &(data_rvf[nAtoms*20]);


    //================ read input gro file ==================
    
    // vQ and vP are "velocities" of the thermostat/barostat particles
    for (int i_nhc = 0; i_nhc < s_system.num_nhc; ++ i_nhc)
    {
        s_system.vQ[i_nhc] = 0.0;
        s_system.aQ[i_nhc] = 0.0;
    }
    s_system.vP = 0.0;

    int groNAtoms;
    read_gro(input_gro, s_system, &groNAtoms, atom_info);

    if (groNAtoms != nAtoms)
    {
       printf("Error: groNAtoms(%d) not equal to nAtoms(%d)!\n", 
              groNAtoms, nAtoms);
       exit(-1);
    }

    // half of box length for PBC
    s_system.box[3] = s_system.box[0] * 0.5;
    s_system.box[4] = s_system.box[1] * 0.5;
    s_system.box[5] = s_system.box[2] * 0.5;

    s_system.volume = s_system.box[0] * s_system.box[1] * s_system.box[2];
    s_system.inv_volume = 1.0 / s_system.volume;


    //================== set MD variables ==========================

    // degree of freedom
    s_system.ndf = 3 * (nAtoms - number_VSites) - 3 - number_Cstrs;

    // temperature coupling; kT = kB*T, in kJ mol^-1
    s_runset.kT    = K_BOLTZ * s_runset.ref_temp;
    s_system.qMass = (double)s_system.ndf * s_runset.kT *
                     s_runset.tau_temp * s_runset.tau_temp;
    //s_system.pMass = (double)s_system.ndf * s_runset.kT * s_runset.tau_pres * s_runset.tau_pres;

    // temperature control
    s_system.first_temp = 0.0;
    s_system.ext_temp = 0.0;

    // time step
    s_runset.dt_2 = 0.5 * s_runset.dt;


    //=================== CPIM matrix and arrays =========================
    // Relay matrix
    // external electric field and potential
    // induced dipoles and charges
    
    double *data_relay = NULL;

    // dimension of matrix: 4M + n_NPs
    int n_mat = s_metal.num * 4 + s_metal.n_NPs;

    // initialize mat_relay, vec_ext and vec_pq
    if (s_metal.min >= 0 && s_metal.max >= s_metal.min)
    {
        data_relay = new (std::nothrow) double [n_mat * 3];
        assert(data_relay != nullptr);

        s_metal.vec_pq     = &(data_relay[0]);
        s_metal.vec_ext    = &(data_relay[n_mat]);
        s_metal.diag_relay = &(data_relay[n_mat * 2]);

        for (int i_mat = 0; i_mat < n_mat; ++ i_mat)
        {
            s_metal.vec_pq[i_mat]  = 0.0;
            s_metal.vec_ext[i_mat] = 0.0;
            s_metal.diag_relay[i_mat] = 1.0;
        }
    }


    //=================== vectors for the BiCGSTAB solver =========================

    double* data_bicgstab = new (std::nothrow) double [n_mat * 11];
    assert(data_bicgstab != nullptr);
    
    Bicgstab s_bicgstab;

    s_bicgstab.Ax = &(data_bicgstab[0]);
    s_bicgstab.r0 = &(data_bicgstab[n_mat]);
    s_bicgstab.r  = &(data_bicgstab[n_mat * 2]);
    s_bicgstab.p  = &(data_bicgstab[n_mat * 3]);
    s_bicgstab.v  = &(data_bicgstab[n_mat * 4]);
    s_bicgstab.s  = &(data_bicgstab[n_mat * 5]);
    s_bicgstab.t  = &(data_bicgstab[n_mat * 6]);
    s_bicgstab.y  = &(data_bicgstab[n_mat * 7]);
    s_bicgstab.z  = &(data_bicgstab[n_mat * 8]);
    s_bicgstab.Kt = &(data_bicgstab[n_mat * 9]);
    s_bicgstab.K  = &(data_bicgstab[n_mat * 10]);


    //================== compute forces ==========================

    mpi_force(s_task, s_topol, atom_info, mol_info,
              s_runset, s_metal, s_system, s_bicgstab,
              my_id, num_procs, time_used);


    //========== file handles: gro, vec_pq, binary dat, parameters ========
    
    FILE *file_gro, *file_pq, *file_dat;
    file_gro = NULL;
    file_pq  = NULL;
    file_dat = NULL;


    //========== adjust velocities and write trajectories =================

    if (ROOT_PROC == my_id) 
    {
        // sum potential energies
        sum_potential(s_system.potential);

        // update temperature
        remove_comm(nAtoms, atom_info, s_system);
        kinetic_energy(s_system, nAtoms, atom_info);

        // save the starting temperature
        s_system.first_temp = s_system.inst_temp;
        s_system.ext_temp   = s_runset.ref_temp;

        // creat traj.gro for writing
        file_gro = fopen("traj.gro","w") ;
        if (NULL == file_gro) 
        {
            printf( "Cannot write to traj.gro!\n" ) ;
            exit(-1);
        }

        // creat vec_pq.txt for writing
        file_pq = fopen("vec_pq.txt","w") ;
        if (NULL == file_pq) 
        {
            printf( "Cannot write to vec_pq.txt!\n" ) ;
            exit(-1);
        }

        // creat traj.dat for writing
        file_dat = fopen("traj.dat","w") ;
        if (NULL == file_dat) 
        {
            printf( "Cannot write to traj.dat!\n" ) ;
            exit(-1);
        }

        // get maximal force
        get_fmax_rms(nAtoms, s_system);

        // write the starting geometry to gro and dat file
        write_gro(file_gro, s_system, nAtoms, atom_info, 0);
        write_binary(file_dat, s_system, nAtoms, 0);

        if (s_metal.min >= 0 && s_metal.max >= s_metal.min)
        {
            write_vec_pq(file_pq, s_metal, 0);
        }


        printf("              .---------------- start MD calculation ---------------.\n");
        printf("\n");

        // check s_runset.run_type
        if (std::string("em") == s_runset.run_type ||
            std::string("cg") == s_runset.run_type)
        {
            printf("    Step %-5d Fmax=%10.3e  E=%15.8e\n", 
                    0, s_system.f_max, s_system.potential[0]);
        }
        else if (std::string("md") == s_runset.run_type)
        {
            printf("  %10.3f  Fmax=%.2e  E=%.6e  T=%.3e\n", 
                   0.0, s_system.f_max, s_system.potential[0], s_system.inst_temp);
        }
    }


    //========= Now start MD steps ============================

    int step = 0;

    //===================================================
    // Energy minimization using steepest descent or CG
    //===================================================

    if (std::string("em") == s_runset.run_type ||
        std::string("cg") == s_runset.run_type)
    {
        s_system.vQ[0] = 0.0;
        s_system.vP = 0.0;

        int converged = 0;
        double gamma = 0.0;
        double delta_pot;

        // initialize direction sx,sy,sz
        for (int i = 0; i < nAtoms; ++ i)
        {
            s_system.old_fx[i] = s_system.fx[i];
            s_system.old_fy[i] = s_system.fy[i];
            s_system.old_fz[i] = s_system.fz[i];

            s_system.sx[i] = s_system.fx[i];
            s_system.sy[i] = s_system.fy[i];
            s_system.sz[i] = s_system.fz[i];
        }

        for (step = 1; step <= s_runset.em_steps && 0 == converged; ++ step) 
        {
            // update coordinates on root processor
            if (ROOT_PROC == my_id) 
            {
                s_system.old_potential = s_system.potential[0];

                for (int i = 0; i < nAtoms; ++ i)
                {
                    s_system.old_rx[i] = s_system.rx[i];
                    s_system.old_ry[i] = s_system.ry[i];
                    s_system.old_rz[i] = s_system.rz[i];

                    // fix metal coordinates?
                    if (1 == s_metal.fix_pos && 1 == atom_info[i].is_metal)
                    {
                        continue;
                    }

                    s_system.rx[i] += s_system.sx[i] / s_system.f_max *
                                      s_runset.em_length;
                    s_system.ry[i] += s_system.sy[i] / s_system.f_max *
                                      s_runset.em_length;
                    s_system.rz[i] += s_system.sz[i] / s_system.f_max *
                                      s_runset.em_length;
                }

                // apply constraints
                rattle_1st(s_runset.dt, mol_info, atom_info, s_topol, s_system);

                // zero velocities
                for (int i = 0; i < nAtoms; ++ i)
                {
                    s_system.vx[i] = 0.0;
                    s_system.vy[i] = 0.0;
                    s_system.vz[i] = 0.0;
                }
            }

            // update forces
            mpi_force(s_task, s_topol, atom_info, mol_info,
                      s_runset, s_metal, s_system, s_bicgstab,
                      my_id, num_procs, time_used);

            if (ROOT_PROC == my_id) 
            {
                // check potential and fmax on root processor
                sum_potential(s_system.potential);
                get_fmax_rms(nAtoms, s_system);
                delta_pot = s_system.potential[0] - s_system.old_potential;

                // update EM step size
                if (delta_pot <= 0.0) 
                { 
                    s_runset.em_length *= 1.2; 
                }
                else 
                { 
                    s_runset.em_length *= 0.2; 
                }
                
                // print info and write to the gro file
                printf("    Step %-5d Fmax=%10.3e  E=%15.8e\n", 
                        step, s_system.f_max, s_system.potential[0]);

                // write trajectories
                write_gro(file_gro, s_system, nAtoms, atom_info, step);
                write_binary(file_dat, s_system, nAtoms, step);

                if (s_metal.min >=0 && s_metal.max >= s_metal.min)
                {
                    write_vec_pq(file_pq, s_metal, step);
                }

                // check convergence
                if (s_system.f_max < s_runset.em_tol && 
                    s_system.f_rms < s_runset.em_tol * 0.5 && 
                    fabs(delta_pot) < s_runset.em_tol * 0.1)
                {
                    printf("\n");
                    printf("    F_max   (%13.6e) smaller than %e\n", 
                            s_system.f_max, s_runset.em_tol);
                    printf("    F_rms   (%13.6e) smaller than %e\n", 
                            s_system.f_rms, s_runset.em_tol * 0.5);
                    printf("    delta_E (%13.6e) smaller than %e\n", 
                            delta_pot, s_runset.em_tol * 0.1);
                    printf("\n");
                    printf("    ===========  Optimization converged ============\n");
                    converged = 1;
                }
                else
                {
                    // update gamma for CG optimization
                    // for steep descent, gamma = 0.0
                    if (std::string("cg") == s_runset.run_type)
                    {
                        double g22 = 0.0;
                        double g12 = 0.0;
                        double g11 = 0.0;
                        for (int i = 0; i < nAtoms; ++ i)
                        {
                            g22 += s_system.fx[i] * s_system.fx[i] + 
                                   s_system.fy[i] * s_system.fy[i] + 
                                   s_system.fz[i] * s_system.fz[i];
                            g12 += s_system.old_fx[i] * s_system.fx[i] + 
                                   s_system.old_fy[i] * s_system.fy[i] + 
                                   s_system.old_fz[i] * s_system.fz[i];
                            g11 += s_system.old_fx[i] * s_system.old_fx[i] + 
                                   s_system.old_fy[i] * s_system.old_fy[i] + 
                                   s_system.old_fz[i] * s_system.old_fz[i];
                        }
                        gamma = (g22 - g12) / g11;
                    }

                    for (int i = 0; i < nAtoms; ++ i)
                    {
                        s_system.sx[i] = s_system.fx[i] + gamma * s_system.sx[i];
                        s_system.sy[i] = s_system.fy[i] + gamma * s_system.sy[i];
                        s_system.sz[i] = s_system.fz[i] + gamma * s_system.sz[i];

                        s_system.old_fx[i] = s_system.fx[i];
                        s_system.old_fy[i] = s_system.fy[i];
                        s_system.old_fz[i] = s_system.fz[i];
                    }
                }

                // communicate convergence
                for (int an_id = 1; an_id < num_procs; ++ an_id) 
                {
                    MPI::COMM_WORLD.Send(&converged, 1, MPI::INT, an_id, tag_99);
                }
            }
            else
            {
                MPI::COMM_WORLD.Recv(&converged, 1, MPI::INT, ROOT_PROC,
                                     tag_99, status);
            }

            // exit loop if converged
            if (1 == converged) { break; }
        }
    }


    //===============================
    // MD with nvt ensemble
    //===============================

    else if (std::string("md") == s_runset.run_type)
    {
        for (step = 1; step <= s_runset.nSteps; ++ step) 
        {
            if (ROOT_PROC == my_id) 
            {
                // gradually increase the reference temperature
                if (step < s_runset.nHeating)
                {
                    s_runset.ref_temp = 
                        s_system.first_temp + 
                        (s_system.ext_temp - s_system.first_temp) *
                        step / s_runset.nHeating;
                }
                else
                {
                    s_runset.ref_temp = s_system.ext_temp;
                }

                // update kT accordingly
                s_runset.kT = K_BOLTZ * s_runset.ref_temp;
            
                // thermostat for 1st half step
                if (std::string("nvt") == s_runset.ensemble)
                {
                    nose_hoover_chain(s_runset, s_system, nAtoms);
                }

                // update velocity for 1st half step
                for (int i = 0; i < nAtoms; ++ i)
                {
                    // fix metal coordinates?
                    if (1 == s_metal.fix_pos && 1 == atom_info[i].is_metal)
                    {
                        continue;
                    }

                    // for virtual sites, inv_mass == 0.0;
                    double inv_mass = atom_info[i].inv_mass;
                    s_system.vx[i] += s_runset.dt_2 * s_system.fx[i] * inv_mass;
                    s_system.vy[i] += s_runset.dt_2 * s_system.fy[i] * inv_mass;
                    s_system.vz[i] += s_runset.dt_2 * s_system.fz[i] * inv_mass;
                }

                // update position for the whole time step
                // using velocity at half time step
                for (int i = 0; i < nAtoms; ++ i)
                {
                    // fix metal coordinates?
                    if (1 == s_metal.fix_pos && 1 == atom_info[i].is_metal)
                    {
                        continue;
                    }

                    s_system.old_rx[i] = s_system.rx[i];
                    s_system.old_ry[i] = s_system.ry[i];
                    s_system.old_rz[i] = s_system.rz[i];

                    s_system.rx[i] += s_runset.dt * s_system.vx[i];
                    s_system.ry[i] += s_runset.dt * s_system.vy[i];
                    s_system.rz[i] += s_runset.dt * s_system.vz[i];
                }

                // apply constraints for the 1st half
                rattle_1st(s_runset.dt, mol_info, atom_info, s_topol, s_system);
            }


            // compute forces
            mpi_force(s_task, s_topol, atom_info, mol_info,
                      s_runset, s_metal, s_system, s_bicgstab,
                      my_id, num_procs, time_used);
            

            // update velocities
            if (ROOT_PROC == my_id) 
            {
                // sum potential energies
                sum_potential(s_system.potential);

                // update velocity for 2nd half step
                for (int i = 0; i < nAtoms; ++ i)
                {
                    // fix metal coordinates?
                    if (1 == s_metal.fix_pos && 1 == atom_info[i].is_metal)
                    {
                        continue;
                    }

                    // for virtual sites, inv_mass == 0.0;
                    double inv_mass = atom_info[i].inv_mass;
                    s_system.vx[i] += s_runset.dt_2 * s_system.fx[i] * inv_mass;
                    s_system.vy[i] += s_runset.dt_2 * s_system.fy[i] * inv_mass;
                    s_system.vz[i] += s_runset.dt_2 * s_system.fz[i] * inv_mass;
                }

                // apply constraints for the 2nd half
                rattle_2nd(s_runset.dt, mol_info, atom_info, s_topol, s_system);

                // update temperature
                kinetic_energy(s_system, nAtoms, atom_info);


                // thermostat for 2nd half step
                if (std::string("nvt") == s_runset.ensemble)
                {
                    nose_hoover_chain(s_runset, s_system, nAtoms);
                }


                // remove center of mass motion and update temperature
                remove_comm(nAtoms, atom_info, s_system);
                kinetic_energy(s_system, nAtoms, atom_info);


                // print information for this step
                if (0 == step % s_runset.nSave)
                {
                    // apply PBC
                    apply_pbc(nMols, mol_info, s_system.rx, s_system.ry, s_system.rz, s_system.box);
             
                    get_fmax_rms(nAtoms, s_system);
                    printf("  %10.3f  Fmax=%.2e  E=%.6e  T=%.3e (%.1f)\n", 
                           s_runset.dt * step, s_system.f_max, s_system.potential[0], 
                           s_system.inst_temp, s_runset.ref_temp);

                    // write to dat file
                    write_binary(file_dat, s_system, nAtoms, step);
                }
             
                // regularly write to gro file
                if (0 == step % (s_runset.nSave * 10))
                {
                    write_gro(file_gro, s_system, nAtoms, atom_info, step);

                    if (s_metal.min >=0 && s_metal.max >= s_metal.min)
                    {
                        write_vec_pq(file_pq, s_metal, step);
                    }
                }
            }

        }
    }


    //=========== Finalize MD ====================
    
    sum_time_used(time_used, my_id, num_procs);

    if (ROOT_PROC == my_id) 
    {
        fclose(file_dat);
        fclose(file_pq);
        fclose(file_gro);

#ifdef DEBUG
        int i;
        for (i = 0; i < nAtoms; ++ i)
        {
            printf("    f[%5d]= %12.5e, %12.5e, %12.5e\n", 
                    i, s_system.fx[i], s_system.fy[i], s_system.fz[i]);
        }
#endif

        printf("\n");
        printf("\n");
        printf("              .------------- final potential energies --------------.\n");
        printf("\n");

        print_potential(s_system.potential);

        time_t end_t = time(NULL);
        gettimeofday(&tv, NULL);
        double end_time = (tv.tv_sec) + (tv.tv_usec) * 1.0e-6;

        printf("    Job ended normally at %s", ctime(&end_t));
        printf("    %.2f seconds were used\n", end_time - start_time );
        printf("\n");
        printf("\n");

        printf("              .------------------- time usage ----------------------.\n");
        printf("\n");

        analyze_time_used(time_used, num_procs);
        printf("\n");
    }

    // free arrays
    for (int an_id = 0; an_id < num_procs; ++ an_id)
    {
        delete[] time_used[an_id];
    }
    delete[] time_used;

    for (int mol = 0; mol < s_topol.mol_types; ++ mol)
    {
        delete[] s_topol.bond_param[mol];
        delete[] s_topol.pair_param[mol];
        delete[] s_topol.angle_param[mol];
        delete[] s_topol.dihedral_param[mol];

        delete[] s_topol.vsite_funct[mol];
        delete[] s_topol.vsite_4[mol];
        delete[] s_topol.constraint[mol];
        
        for (int atom_i = 0; atom_i < s_topol.atom_num[mol]; ++ atom_i)
        {
            delete[] s_topol.exclude[mol][atom_i];
        }
        delete[] s_topol.exclude[mol];
        
        delete[] s_topol.atom_param[mol];
    }

    delete[] s_topol.bond_param;
    delete[] s_topol.pair_param;
    delete[] s_topol.angle_param;
    delete[] s_topol.dihedral_param;
    
    delete[] s_topol.vsite_funct;
    delete[] s_topol.vsite_4;
    delete[] s_topol.constraint;
    delete[] s_topol.exclude;

    delete[] s_topol.mol_num;
    delete[] s_topol.atom_num;
    delete[] s_topol.atom_param;
    delete[] s_topol.nonbonded_param;

    if (s_metal.min >=0 && s_metal.max >= s_metal.min)
    {
        delete[] s_metal.inv_sqrt_dens;
        delete[] data_relay;
    }
    delete[] s_metal.cpff_chg;
    delete[] s_metal.start_NP;
    delete[] s_metal.end_NP;
    
    delete[] data_bonded;
    delete[] data_nonbonded;
    delete[] data_bicgstab;
    delete[] data_start_end;
    delete[] data_long_start_end;
    delete[] data_potential;
    delete[] data_rvf;
    
    delete[] atom_info;
    delete[] mol_info;

    MPI::Finalize();

    return 0;
}
