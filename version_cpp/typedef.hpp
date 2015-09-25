/**
 \file typedef.hpp
 \author Xin Li
 \date 2015/09
 \version 1.0.2
 
 \brief define constants and structure types
 
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

#include <array>
#include <vector>
#include <string>

//================= Physical and mathematical constants =====================

/// root processor for MPI
#define ROOT_PROC 0

/// number of timers
#define N_TIMER 15

/// number of potential energy terms
#define N_POT 15

/// number of Nose-Hoover chains
#define N_NHC 10

/// maximal length of string
#define MAX_STR_LEN 256

/// dimension
#define DIM 3

/// pi
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/// 2.0 * pi
#define TWO_PI 6.28318530717959

/// 1.0 / sqrt(pi)
#define INV_SQRT_PI 0.564189583547756

/// \brief Coulomb constant
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

/// \brief Boltzmann constant in MD unit.
/// kB = 1.3806488e-3 * 6.02214129 kJ mol^-1 K^-1.
#define K_BOLTZ 0.00831446214546895

/// \brief Pressure converting factor.
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


/// \brief vector
/// \details has Cartesian x,y,z components
typedef struct {
    double x, y, z;
} Vec_R;

/// \brief bond stretching force
/// \details has two vectors fi and fj
typedef struct {
    Vec_R fi, fj;
} Vec_2;

/// \brief angle bending force
/// \details has three vectors fi, fj and fk
typedef struct {
    Vec_R fi, fj, fk;
} Vec_3;

/// \brief dihedral force
/// \details has four vectors fi, fj, fk and fl
typedef struct {
    Vec_R fi, fj, fk, fl;
} Vec_4;

/// \brief nonbonded force (VDW & Coulomb)
/// \details has two Vec_2 components lj and qq
typedef struct {
    Vec_2 lj, qq;
} Vec_nb;

/// \brief nonbonded parameters
typedef struct {
    double charge, mass;
    int atomtype;
} AtomParam;

/// \brief virtual site, type 4 (as in TIP5P water model)
/// \details four atoms, three parameters for construction of virtual sites
typedef struct {
    int atom_i, atom_j, atom_k, atom_s;
    double a, b, c;
} VSite_4;

/// \brief constraint
/// \details two atoms, one distance for a constraint
typedef struct {
    int atom_i, atom_j;
    double r;
} Constraint;


/// \brief atomic info
typedef struct {
    double charge, mass, inv_mass;
    int atomtype, iAtom, iMol, molID, resID, is_metal;
    char atomName[5], resName[5];
} Atom_Info;


/// \brief molecular info
/// \details start and end indices of atoms in a molecule
typedef struct {
    int mini, maxi;
} Mol_Info;


/// bond stretching parameters
typedef struct {
    int atom_i, atom_j, funct;
    double b0, kb, D, beta;
} BondParam;


/// angle bending parameters
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


/// MD settings
typedef struct {
    /// \name General settings
    ///@{
    std::string run_type     = "md";
    std::string ensemble     = "nvt";
    std::string vdw_type     = "shifted";
    std::string coulomb_type = "wolf_sum";
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
    std::array<double, DIM> external_efield = { {0., 0., 0.} };
} RunSet;


/// metal parameters
typedef struct {
    /// \name start and end indices
    ///@{
    int min, max;
    ///@}

    int num; ///< number of metal atoms;
    int fix_pos; ///< fix the position of metal atoms?
    int use_cpff; ///< use capacitance-polarizability force field?
    int n_NPs; ///< number of nanoparticles

    /// \name quantum Sutton-Chen potential
    ///@{
    int qsc_n, qsc_m;
    double qsc_eps, qsc_c, qsc_a;
    double *inv_sqrt_dens;
    ///@}

    /// \name capacitance-polarizability parameters
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

    /// \name CPIM vectors: external field/potential & induced dipole/charge.
    ///@{
    double *vec_ext, *vec_pq;
    ///@}

    /// \name CPIM matrix in COO format
    ///@{
    double *diag_relay;
    std::vector<double> vec_val;
    std::vector<long int> vec_row_ind, vec_col_ind;
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

    VSite_4        **vsite_4; ///< virtual sites (type 4, as in TIP5P)
    Constraint     **constraint; ///< constraint
    NonbondedParam **nonbonded_param; ///< VDW parameters
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
    // rectangular box only.
    // "s_system.box" has six elements
    // box_x, box_y, box_z, 0.5 * box_x, 0.5 * box_y, 0.5 * box_z
    double box[DIM * 2];
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
    double virial[DIM][DIM], partial_vir[DIM][DIM];
    ///@}

    /// \name temperature & pressure coupling
    ///@{
    double eKsum, qMass, pMass, vP;
    double vQ[N_NHC], aQ[N_NHC];
    int    num_nhc = N_NHC;
    double first_temp, ext_temp, inst_temp, inst_pres;
    double pressure[3];
    int    ndf;
    ///@}
} System;
