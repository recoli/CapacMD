/*! \mainpage CapacMD
   
    \section Introduction

    CapacMD is a molecular dynamics simulation engine using capacitance-polarizability force field. 
    See http://www.theochem.kth.se/~lixin/capacmd/index.html

    Related paper: X. Li and H. Agren, J. Phys. Chem. C, 2015, 119, 19430-19437. 
    http://pubs.acs.org/doi/abs/10.1021/acs.jpcc.5b04347

    \section Usage
   
    \subsection subsec_1 1: Run the simulation
    mpirun -np 4 ./mpi-main -gro init.gro -par param.txt -mds mdset.txt
*/


/*! @file    typedef.h
    @brief   define constants and structure types

    This file is part of the CapacMD program.                                   
    Copyright (C) 2015 Xin Li <lixin.reco@gmail.com>                            
                                                                              
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

//================= Physical and mathematical constants =====================

/// maximal length of string
#define MAX_STR_LEN 256

/// dimension
#define DIM 3

/// PI
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/// 2.0 * M_PI
#define TWO_PI 6.28318530717959

/// 1.0 / sqrt(M_PI)
#define INV_SQRT_PI 0.564189583547756

/// 1.0 / (4.0 * M_PI * epsilon_0), in MD unit.
///
/// Physical constants from 2010 CODATA:
/// c (speed of light in vacuum)=  299792458 m s^-1;
/// mu_0 (magnetic constant)=      4*pi * 10^-7 N A^-2 (A=Ampere);
/// e (elementary charge)=         1.602176565 * 10^-19 C;
/// nA (Avogadro constant)=        6.02214129 * 10^23 mol^-1;
/// kB (Boltzmann constant)=       1.3806488 * 10^-23 J K^-1.
///
/// Coulomb constant k=1/(4*pi*epsilon_0)=c^2*mu_0/(4*pi);
/// fQQ = 2.99792458^2 * 10^16 * 10^-7 N m^2 C^-2
///     = 2.99792458^2 * 6.02214129 * 1.602176565^2 kJ mol^-1 nm e^-2
///     = 138.93545782935 [kJ mol^-1 nm e^-2].
///
/// fQQ in atomic unit is 1.0 [Eh a0 e^-2];
/// 1 Eh = 2625.49964037578 kJ/mol;
/// 1 Bohr = 0.052917721092 nm;
/// fQQ in MD unit is 138.935457696552 [kJ mol^-1 nm e^-2].
///
/// In this program we use FQQ = 138.935457696552.

#define FQQ 138.935457696552

/// 1.0 eV = 96.4853363969351 kJ/mol.
#define EV2KJMOL 96.4853363969351

/// 1.0 nm = 18.8972612456506 Bohr.
#define NM2BOHR 18.8972612456506

/// Boltzmann constant in MD unit [kJ mol^-1 K^-1].
/// kB = 1.3806488e-3 * 6.02214129 kJ mol^-1 K^-1.
#define K_BOLTZ 0.00831446214546895

/// Pressure converting factor.
/// 1 bar = 1e+05 Pa = 1e+05 N/m^2 = 1e+05 J/m^3;
/// 1 kJ mol^-1 nm^-3 = 16.6053892103219 bar.
#define P_CONST 16.6053892103219

/// maixmal number of iterations in BiCGSTAB or RATTLE.
#define MAX_ITER 1000


/// threshold for charge-charge interaction element in CPIM matrix.
#define THRSHD_Q  1.0e-15
/// threshold for dipole-charge interaction element in CPIM matrix.
#define THRSHD_M  1.0e-15
/// threshold for dipole-dipole interaction element in CPIM matrix.
#define THRSHD_P  1.0e-15


/// vector
typedef struct {
    double x, y, z;
} Vec_R;

/// bond force
typedef struct {
    Vec_R fi, fj;
} Vec_2;

/// angle force
typedef struct {
    Vec_R fi, fj, fk;
} Vec_3;

/// dihedral force
typedef struct {
    Vec_R fi, fj, fk, fl;
} Vec_4;

/// nonbonded force
typedef struct {
    Vec_2 lj, qq;
} Vec_nb;

/// nonbonded parameters
typedef struct {
    double charge, mass;
    int atomtype;
} AtomParam;

/// virtual site, type 4 (as in TIP5P water model)
typedef struct {
    int atom_i, atom_j, atom_k, atom_s;
    double a, b, c;
} VSite_4;

/// constraint
typedef struct {
    int atom_i, atom_j;
    double r;
} Constraint;


/// atomic info
typedef struct {
    double charge, mass, inv_mass;
    int atomtype, iAtom, iMol, molID, resID, is_metal;
    char atomName[5], resName[5];
} Atom_Info;


/// molecular info
typedef struct {
    int mini, maxi;
} Mol_Info;


/// bond parameters
typedef struct {
    int atom_i, atom_j, funct;
    double b0, kb, D, beta;
} BondParam;


/// angle parameters
typedef struct {
    int atom_i, atom_j, atom_k, funct;
    double a0, ka, b0_ub, kb_ub;
} AngleParam;


/// dihedral parameters
typedef struct {
    int atom_i, atom_j, atom_k, atom_l, funct;
    double c0, c1, c2, c3, c4, c5, phi0, kphi;
    int n;
} DihedralParam;


/// 1-4 pair parameters
typedef struct {
    int atom_i, atom_j, funct;
    double scaleQQ, q_i, q_j, V, W;
} PairParam;


/// Wolf summation parameters
typedef struct {
    double alpha, a2_sqrtPI, wolfConst, erfc_arCut;
} WolfParam;


/// nonbonded parameters
typedef struct {
    int funct;
    double C12, C6, A, B;
} NonbondedParam;


/// energy minimization settings
typedef struct {
    int steps;
    double length, tol;
} EM_Set;


/// MD settings
typedef struct {
    /// \name General settings
    ///@{
    char run_type[MAX_STR_LEN], ensemble[MAX_STR_LEN];
    char vdw_type[MAX_STR_LEN], coulomb_type[MAX_STR_LEN];
    int  use_vdw, use_coulomb;
    ///@}

    /// \name Energy minimization settings
    ///@{
    int    em_steps;
    double em_length, em_tol;
    ///@}

    /// \name Molecular dynamics settings
    ///@{
    int    nSteps, nSave, nHeating;
    double rCut, rCut2;
    double inv_rc12, inv_rc6, scaleLJ, scaleQQ;
    double ref_temp, tau_temp, ref_pres, tau_pres;
    ///@}

    /// kT = kB*T, in kJ mol^-1
    double kT;

    /// \name Time step
    ///@{
    double dt, dt_2;
    ///@}

    /// \name Wolf summation constants
    ///@{
    double w_alpha, w_a2_sqrtPI, w_Const, w_erfc_arCut;
    ///@}

    /// External electric field
    double *external_efield;
} RunSet;


/// metal parameters
typedef struct {
    int min, max, num;
    int fix_pos, use_cpff, n_NPs;

    /// \name quantum Sutton-Chen potential
    ///@{
    int qsc_n, qsc_m;
    double qsc_eps, qsc_c, qsc_a;
    double *inv_sqrt_dens;
    ///@}

    /// \name CPFF parameters
    ///@{
    double cpff_polar, cpff_capac; 
    double inv_polar, inv_capac;
    double inv_R_pp, inv_R_pq, inv_R_qq;
    ///@}

    /// \name nanoparticle charges and numberings
    ///@{
    double *cpff_chg;
    int *start_NP, *end_NP;
    ///@}

    /// \name CPIM arrays: external field/potential, induced dipole/charge.
    ///@{
    double *vec_ext, *vec_pq;
    ///@}

    /// \name CPIM matrix in COO format
    ///@{
    double *diag_relay;
    double *val;
    long int *col_ind, *row_ind;
    ///@}
} Metal;


/// vectors that will be used in BiCGSTAB solver
typedef struct {
    double *Ax, *r0, *r, *p, *v, *s, *t, *y, *z, *Kt, *K;
} Bicgstab;


/// force field topology
typedef struct {
    /// \name number of bonded interaction, virtual sites, constraints
    ///@{
    int ***exclude;
    int *n_bonds, *n_pairs, *n_angles, *n_dihedrals, *n_vsites, *n_constraints;
    int **vsite_funct;
    ///@}

    /// scaling factor for 1-4 VDW interaction
    double scaleLJ;
    /// scaling factor for 1-4 Coulomb interaction
    double scaleQQ;

    /// number of vdW types
    int n_types;
    /// number of molecule types
    int mol_types;

    /// number of atoms
    int n_atoms;
    /// number of molecules
    int n_mols;

    /// number of atoms in a molecule type.
    int *atom_num;
    /// number of this type of molecules in the system.
    int *mol_num;

    /// atom_param: charge, mass, atomtype.
    AtomParam      **atom_param;

    /// \name bonded parameters: bond, 1-4 pair, angle, dihedral
    ///@{
    BondParam      **bond_param;
    PairParam      **pair_param;
    AngleParam     **angle_param;
    DihedralParam  **dihedral_param;
    ///@}

    /// virtual sites (type 4, as in TIP5P)
    VSite_4        **vsite_4;
    /// constraints
    Constraint     **constraint;
    /// VDW parameters
    NonbondedParam **nonbonded_param;

} Topol;


/// distributed task on each processor
typedef struct {
    int *start_mol,   *end_mol; 
    int *start_atom,  *end_atom;
    int *start_metal, *end_metal;
    long int *start_pair, *end_pair;
} Task;


/// system information
typedef struct {
    /// \name simulation box
    ///@{
    double *box;
    double volume, inv_volume;
    ///@}

    /// \name position, velocity, force, and potential
    ///@{
    double *rx, *ry, *rz, *vx, *vy, *vz, *fx, *fy, *fz, *potential;
    ///@}

    /// \name force and potential on slave processors
    ///@{
    double *partial_fx, *partial_fy, *partial_fz, *partial_pot; 
    ///@}

    /// \name potential and forces from previous step (for energy minimization)
    ///@{
    double old_potential;
    double *old_fx, *old_fy, *old_fz;
    double *sx, *sy, *sz;
    ///@}
    /// \name position from previous step (for constraint algorithm)
    ///@{
    double *old_rx, *old_ry, *old_rz;
    ///@}

    /// \name maximal force and root mean square of force
    ///@{
    double f_max, f_rms;
    ///@}

    /// \name virial
    ///@{
    double **virial, **partial_vir;
    ///@}

    /// \name temperature & pressure coupling
    ///@{
    double eKsum, qMass, pMass, vP;
    double *vQ, *aQ;
    int    num_nhc;
    double first_temp, ext_temp, inst_temp, inst_pres;
    double pressure[3];
    int    ndf;
    ///@}
} System;

