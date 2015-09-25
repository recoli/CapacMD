/**
 \file cpff.cpp
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

#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "typedef.hpp"

// include some functions from force.cpp
double dist_2(const Vec_R* ptr_r);
void scale_vec_1(const double ss, Vec_R* ptr_rvec);
void scale_vec_2(const double fij, const Vec_R* ptr_rvec, Vec_R* ptr_fvec);
double pbc(double dx, double xbox);
void pbc_12(Vec_R* ptr_r, const double* box);

//===========================================================
// constructs the CPIM external field and potential, vec_ext
//===========================================================

void mpi_cpff_vec_ext(Task& s_task, Metal& s_metal, RunSet& s_runset,
                      Atom_Info* atom_info, Topol& s_topol,
                      System& s_system, int my_id)
{
    int start_metal = s_task.start_metal[my_id];
    int end_metal   = s_task.end_metal[my_id];

    double alpha      = s_runset.w_alpha;
    double a2_sqrtPI  = s_runset.w_a2_sqrtPI;
    double wolfConst  = s_runset.w_Const;
    double erfc_arCut = s_runset.w_erfc_arCut;

    int use_coulomb = s_runset.use_coulomb;

    double rCut2 = s_runset.rCut2;
    double rCut  = s_runset.rCut;
    
    Vec_R rvec, evec;
    double rij2, rij, qi;
    double e_x_metal, e_y_metal, e_z_metal, pot_metal, e_field, v_field;

    // electric field in atomic unit: Eh/(ea0)
    // CODATA:
    //     Avogadro constant = 6.02214129e+23 mol^-1
    //     1 Eh = 4.35974434e-18 J
    //     1 Bohr = 0.52917721092e-10 m
    // 1 Eh/(e a0) = 49614.7526045429 kJ/mol/(e nm)
    double conv_fac_e = 2.01552955019358e-05;

    // 1 Eh = 2625.49964037578 kJ/mol
    double conv_fac_v = 3.80879884583368e-04;

    //=======================================================================
    // external electric field and potential from non-metallic point charges
    // !!!!!
    // NOTE: This part is evaluated in MD units (nm, kJ/mol, ps, etc.) and
    // then converted to atomic unit (Bohr, Hartree, etc.)
    // !!!!!
    //=======================================================================
    for (int i_metal = start_metal; i_metal <= end_metal; ++ i_metal)
    {
        e_x_metal = 0.0;
        e_y_metal = 0.0;
        e_z_metal = 0.0;
        pot_metal = 0.0;

        for (int i_atom = 0; i_atom < s_topol.n_atoms; ++ i_atom)
        {
            // exclude metal atoms and atoms with zero charges
            if (1 == atom_info[i_atom].is_metal) { continue; }

            // save atomic charge in qi for further use
            qi = atom_info[i_atom].charge;
            if (0.0 == qi) { continue; }

            // distance in nm
            rvec.x = s_system.rx[i_metal] - s_system.rx[i_atom];
            rvec.y = s_system.ry[i_metal] - s_system.ry[i_atom];
            rvec.z = s_system.rz[i_metal] - s_system.rz[i_atom];
            pbc_12(&rvec, s_system.box);
            rij2 = dist_2(&rvec);

            // electric field and potential in MD unit
            // E_vec = F_vec / q = -dV/dr_vec / q
            // Note: v_field, potential V = v_field
            //       e_field, electric field E_vec = e_field * r_vec
            e_field = 0.0;
            v_field = 0.0;

            // short-range contribution
            if (rij2 < rCut2)
            {
                rij = sqrt(rij2);

                // wolf dampled shifted force
                if (1 == use_coulomb)
                {
                    // Modified electric field of surrounding point charnges
                    // based on the damped shifted force (DSF) approach for electrostatics
                    // Ref.:        Fennell and Gezelter, J. Chem. Phys. 2006, 124, 234104
                    //                dx.doi.org/10.1063/1.2206581

                    // note: alpha in nm^-1, rij in nm
                    double erfc_arij = erfc(alpha * rij);
                    double exp__a2rij2 = exp(-alpha * alpha * rij2);

                    e_field += FQQ * qi * 
                               (erfc_arij / rij2 +
                                a2_sqrtPI * exp__a2rij2 / rij - wolfConst) / rij;

                    v_field += FQQ * qi * 
                               (erfc_arij / rij - erfc_arCut + wolfConst * (rij - rCut));
                }

                else if (0 == use_coulomb)
                {
                    e_field += FQQ * qi / (rij2 * rij);

                    v_field += FQQ * qi / rij;
                }
            }

            // convert electric field from MD unit to a.u.
            scale_vec_2(e_field * conv_fac_e, &rvec, &evec);
            e_x_metal += evec.x;
            e_y_metal += evec.y;
            e_z_metal += evec.z;

            // convert potential energy from MD unit to a.u.
            pot_metal += v_field * conv_fac_v;
        }

        // note: e_x_metal, e_y_metal, e_z_metal and pot_metal in a.u.
        int i = i_metal - s_metal.min;
        s_metal.vec_ext[i * 3 + 0] = e_x_metal;
        s_metal.vec_ext[i * 3 + 1] = e_y_metal;
        s_metal.vec_ext[i * 3 + 2] = e_z_metal;

        s_metal.vec_ext[s_metal.num * 3 + i] = pot_metal;

        // external electric field
        // remember to convert from MD unit to a.u.
        s_metal.vec_ext[i * 3 + 0] += s_runset.external_efield[0] * conv_fac_e;
        s_metal.vec_ext[i * 3 + 1] += s_runset.external_efield[1] * conv_fac_e;
        s_metal.vec_ext[i * 3 + 2] += s_runset.external_efield[2] * conv_fac_e;

        double pot_external = s_runset.external_efield[0] * s_system.rx[i_metal] +
                              s_runset.external_efield[1] * s_system.ry[i_metal] +
                              s_runset.external_efield[2] * s_system.rz[i_metal];
        s_metal.vec_ext[s_metal.num * 3 + i] -= pot_external * conv_fac_v;
    }

    // last element, total charge of nanoparticle
    if (ROOT_PROC == my_id)
    {
        for (int iNP = 0; iNP < s_metal.n_NPs; ++ iNP)
        {
            s_metal.vec_ext[s_metal.num * 4 + iNP] = s_metal.cpff_chg[iNP];
        }
    }

    // no need to gather vec_ext to root processor
    // vec_ext will stay distributed on each proc for subsequent BiCGSTAB
}

//==================================================
// construct sparse CPIM matrix
//==================================================

void mpi_cpff_mat_relay_COO(Task& s_task, Metal& s_metal, System& s_system, double rCut2,
                            int my_id, int num_procs, long int *p_count_nnz)
{
    // initialize nnz and size counter
    long int count_nnz  = 0;

    int start_metal = s_task.start_metal[my_id];
    int end_metal   = s_task.end_metal[my_id];

    double inv_polar = s_metal.inv_polar;
    double inv_capac = s_metal.inv_capac;
    double inv_R_pp  = s_metal.inv_R_pp;
    double inv_R_pq  = s_metal.inv_R_pq;
    double inv_R_qq  = s_metal.inv_R_qq;

    int min_metal = s_metal.min;
    int max_metal = s_metal.max;
    int n_metal   = s_metal.num;
    int n_NPs     = s_metal.n_NPs;

    // the delta function
    double delta[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

    //=======================================================================
    // Relay matrix for electrointeraction within metal nanoparticle
    // !!!!!
    // NOTE: This part is evaluated in atomic units!
    // !!!!!
    //=======================================================================
    
    // i_metal loops over local metal atoms on each processor
    for (int i_metal = start_metal; i_metal <= end_metal; ++ i_metal)
    {
        int i = i_metal - min_metal;

        // j_metal loops over all metal atoms
        for (int j_metal = min_metal; j_metal <= max_metal; ++ j_metal)
        {
            int j = j_metal - min_metal;
            
            Vec_R rvec;
            double r[3], rij2, rij;
            double exp_rijRpp2, inv_Rpp3, const_1, const_2, exp__rijRpq2, const_3;

            if (j != i)
            {
                rvec.x = s_system.rx[j_metal] - s_system.rx[i_metal];
                rvec.y = s_system.ry[j_metal] - s_system.ry[i_metal];
                rvec.z = s_system.rz[j_metal] - s_system.rz[i_metal];
                pbc_12(&rvec, s_system.box);
                scale_vec_1(NM2BOHR, &rvec); // nm to atomic unit

                rij2 = dist_2(&rvec);
                rij  = sqrt(rij2);

                r[0] = rvec.x;
                r[1] = rvec.y;
                r[2] = rvec.z;

                exp_rijRpp2 = exp(-(rij2 * inv_R_pp * inv_R_pp));
                inv_Rpp3 = pow(inv_R_pp, 3.0);
                
                const_1 = 1.0 / pow(rij, 5.0) * 
                          (erf(rij * inv_R_pp) - 
                           2.0 * INV_SQRT_PI * inv_R_pp * rij * exp_rijRpp2);
                const_2 = 4.0 * INV_SQRT_PI * inv_Rpp3 / rij2 * exp_rijRpp2;

                exp__rijRpq2 = exp(-(rij2 * inv_R_pq * inv_R_pq));
                
                const_3 = 1.0 / pow(rij, 3.0) * 
                          (erf(rij * inv_R_pq) - 
                           2.0 * INV_SQRT_PI * inv_R_pq * rij * exp__rijRpq2);
            }
    
            // submatrix [A -M 0]: row from start_metal * 3 to end_metal * 3
            for (int a = 0; a < DIM; ++ a)
            {
                int irow = i * 3 + a;

                // submatrix Aij, dimension 3x3
                for (int b = 0; b < DIM; ++ b)
                {
                    int icol = j * 3 + b;

                    if (j == i)
                    {
                        if (b == a)
                        {
                            double element = inv_polar;
                            s_metal.diag_relay[irow] = element;

                            if (fabs(element) > THRSHD_P)
                            {
                                if (count_nnz >= s_metal.vec_val.size()) {
                                    s_metal.vec_val.push_back(element);
                                    s_metal.vec_row_ind.push_back(irow);
                                    s_metal.vec_col_ind.push_back(icol);
                                } else {
                                    s_metal.vec_val[count_nnz] = element;
                                    s_metal.vec_row_ind[count_nnz] = irow;
                                    s_metal.vec_col_ind[count_nnz] = icol;
                                }
                                ++ count_nnz;
                            }
                        }
                    }
                    else
                    {
                        double element =
                            -(const_1 * (3.0 * r[a] * r[b] - delta[a][b] * rij2) -
                              const_2 * r[a] * r[b]);

                        if (fabs(element) > THRSHD_P)
                        {
                            if (count_nnz >= s_metal.vec_val.size()) {
                                s_metal.vec_val.push_back(element);
                                s_metal.vec_row_ind.push_back(irow);
                                s_metal.vec_col_ind.push_back(icol);
                            } else {
                                s_metal.vec_val[count_nnz] = element;
                                s_metal.vec_row_ind[count_nnz] = irow;
                                s_metal.vec_col_ind[count_nnz] = icol;
                            }
                            ++ count_nnz;
                        }
                    }

                }

                // submatrix -Mij, dimension 3x1
                int icol = n_metal * 3 + j;
                {
                    if (j != i)
                    {
                        // -Mij for mu_i - q_j interaction
                        double element = const_3 * r[a];

                        if (fabs(element) > THRSHD_M)
                        {
                            if (count_nnz >= s_metal.vec_val.size()) {
                                s_metal.vec_val.push_back(element);
                                s_metal.vec_row_ind.push_back(irow);
                                s_metal.vec_col_ind.push_back(icol);
                            } else {
                                s_metal.vec_val[count_nnz] = element;
                                s_metal.vec_row_ind[count_nnz] = irow;
                                s_metal.vec_col_ind[count_nnz] = icol;
                            }
                            ++ count_nnz;
                        }
                    }
                }

                // submatrix 0, do nothing
            }

            // submatrix [-M^T -C 1_or_0]: row from n_metal*3 + start_metal to n_metal*3 + end_metal
            int irow = n_metal*3 + i;
            {
                // submatrix -Mij^T, dimension 1x3
                for (int b = 0; b < 3; ++ b)
                {
                    int icol = j * 3 + b;

                    if (j != i)
                    {
                        // -Mji for mu_j - q_i interaction
                        double element = -const_3 * r[b];

                        if (fabs(element) > THRSHD_M)
                        {
                            if (count_nnz >= s_metal.vec_val.size()) {
                                s_metal.vec_val.push_back(element);
                                s_metal.vec_row_ind.push_back(irow);
                                s_metal.vec_col_ind.push_back(icol);
                            } else {
                                s_metal.vec_val[count_nnz] = element;
                                s_metal.vec_row_ind[count_nnz] = irow;
                                s_metal.vec_col_ind[count_nnz] = icol;
                            }
                            ++ count_nnz;
                        }
                    }
                }

                // submatrix -Cij, dimension 1x1
                int icol = n_metal * 3 + j;
                if (j == i)
                {
                    double element = -inv_capac;
                    s_metal.diag_relay[irow] = element;

                    if (fabs(element) > THRSHD_Q)
                    {
                        if (count_nnz >= s_metal.vec_val.size()) {
                            s_metal.vec_val.push_back(element);
                            s_metal.vec_row_ind.push_back(irow);
                            s_metal.vec_col_ind.push_back(icol);
                        } else {
                            s_metal.vec_val[count_nnz] = element;
                            s_metal.vec_row_ind[count_nnz] = irow;
                            s_metal.vec_col_ind[count_nnz] = icol;
                        }
                        ++ count_nnz;
                    }
                }
                else
                {
                    double element = -erf(rij * inv_R_qq) / rij;

                    if (fabs(element) > THRSHD_Q)
                    {
                        if (count_nnz >= s_metal.vec_val.size()) {
                            s_metal.vec_val.push_back(element);
                            s_metal.vec_row_ind.push_back(irow);
                            s_metal.vec_col_ind.push_back(icol);
                        } else {
                            s_metal.vec_val[count_nnz] = element;
                            s_metal.vec_row_ind[count_nnz] = irow;
                            s_metal.vec_col_ind[count_nnz] = icol;
                        }
                        ++ count_nnz;
                    }
                }

                // submatrix 1_or_0, Lagrangian part, dimension 1 x n_NPs
                if (min_metal == j_metal)
                {
                    for (int iNP = 0; iNP < n_NPs; ++ iNP)
                    {
                        int icol = n_metal * 4 + iNP;

                        if (i_metal >= s_metal.start_NP[iNP] && 
                            i_metal <= s_metal.end_NP[iNP])
                        {
                            double element = 1.0;

                            if (fabs(element) > THRSHD_Q)
                            {
                                if (count_nnz >= s_metal.vec_val.size()) {
                                    s_metal.vec_val.push_back(element);
                                    s_metal.vec_row_ind.push_back(irow);
                                    s_metal.vec_col_ind.push_back(icol);
                                } else {
                                    s_metal.vec_val[count_nnz] = element;
                                    s_metal.vec_row_ind[count_nnz] = irow;
                                    s_metal.vec_col_ind[count_nnz] = icol;
                                }
                                ++ count_nnz;
                            }
                        }
                    }
                }
            }
        }
    }

    // submatrix [0 1_or_0 0]
    // dimension n_NPs x 1 for each metal atom
    // saved on root processor
    if (ROOT_PROC == my_id)
    {
        for (int iNP = 0; iNP < n_NPs; ++ iNP)
        {
            int irow = n_metal * 4 + iNP;

            // note: on root processor, i_metal loops over all metal atoms
            for (int i_metal = min_metal; i_metal <= max_metal; ++ i_metal)
            {
                int i = i_metal - min_metal;
                int icol = n_metal * 3 + i;

                if (i_metal >= s_metal.start_NP[iNP] && 
                    i_metal <= s_metal.end_NP[iNP])
                {
                    double element = 1.0;

                    if (fabs(element) > THRSHD_Q)
                    {
                        if (count_nnz >= s_metal.vec_val.size()) {
                            s_metal.vec_val.push_back(element);
                            s_metal.vec_row_ind.push_back(irow);
                            s_metal.vec_col_ind.push_back(icol);
                        } else {
                            s_metal.vec_val[count_nnz] = element;
                            s_metal.vec_row_ind[count_nnz] = irow;
                            s_metal.vec_col_ind[count_nnz] = icol;
                        }
                        ++ count_nnz;
                    }
                }
            }
        }
    }

    *p_count_nnz  = count_nnz;

    // no need to gather mat_relay on root processor
    // mat_relay will stay distributed on each proc for subsequent mpi_bicg_stab
}

//=======================
// calculate CPIM forces
//=======================

void mpi_cpff_force(Task& s_task, Metal& s_metal, RunSet& s_runset,
                    System& s_system, double* vec_pq,
                    int nAtoms, Atom_Info* atom_info, int my_id)
{
    int start_metal = s_task.start_metal[my_id];
    int end_metal   = s_task.end_metal[my_id];

    double alpha      = s_runset.w_alpha;
    double a2_sqrtPI  = s_runset.w_a2_sqrtPI;
    double wolfConst  = s_runset.w_Const;
    double erfc_arCut = s_runset.w_erfc_arCut;

    int use_coulomb = s_runset.use_coulomb;

    double rCut2 = s_runset.rCut2;
    double rCut  = s_runset.rCut;

    double inv_polar = s_metal.inv_polar;
    double inv_capac = s_metal.inv_capac;
    double inv_R_pp  = s_metal.inv_R_pp;
    double inv_R_pq  = s_metal.inv_R_pq;
    double inv_R_qq  = s_metal.inv_R_qq;

    int min_metal = s_metal.min;
    int max_metal = s_metal.max;
    int n_metal   = s_metal.num;

    // the delta function
    double delta[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

    double rij2, rij, qi, qj, const_a0, const_a1, const_b1, const_c1;
    double erfc_arij, exp__a2rij2;

    //=======================================================================
    // Now that we have all the dipoles and charges for metallic atoms
    // Compute the corresponding forces
    //=======================================================================

    double alpha3   = pow(alpha, 3.0);
    double inv_Rpp3 = pow(inv_R_pp, 3.0);

    for (int i_metal = start_metal; i_metal <= end_metal; ++ i_metal)
    {
        // charge-charge and charge-dipole interaction between non-metallic
        // and metallic atoms
        // !!!NOTE: this part is evaluated in MD units

        for (int i_atom = 0; i_atom < nAtoms; ++ i_atom)
        {
            // skip metal atoms and atoms with zero charge
            if (1 == atom_info[i_atom].is_metal) { continue; }
            if (0.0 == atom_info[i_atom].charge) { continue; }

            // distances in nm
            Vec_R rvec;
            rvec.x = s_system.rx[i_metal] - s_system.rx[i_atom];
            rvec.y = s_system.ry[i_metal] - s_system.ry[i_atom];
            rvec.z = s_system.rz[i_metal] - s_system.rz[i_atom];
            pbc_12(&rvec, s_system.box);
            rij2 = dist_2(&rvec);

            if (rij2 < rCut2)
            {
                rij = sqrt(rij2);

                // qi: permanent charge of nonmetallic atom
                // qj: induced charge of metal atom
                qi = atom_info[i_atom].charge;

                const int j = i_metal - min_metal;
                qj = vec_pq[n_metal * 3 + j];

                double r[3] = {rvec.x, rvec.y, rvec.z};
                double f[3] = {0., 0., 0.};

                // damped shifted force between nonmetallic point charge
                // and metallic induced charge
                //==========================================================================
                // fij = FQQ * qi * qj * 
                //       ((erfc(alpha * rij)/ rij2 + 
                //        a2_sqrtPI * exp(-alpha * alpha * rij2) / rij) - wolfConst) / rij;
                //==========================================================================

                // const_a0 in nm^-3
                // note: const_a0 will also be used to compute the potential of
                // CPFF metal dipole - non-metal charge interaction
                const_a0 = 0.0;

                if (1 == use_coulomb)
                {
                    erfc_arij   = erfc(alpha * rij);
                    exp__a2rij2 = exp(-alpha * alpha * rij2);

                    const_a0 = (erfc_arij / rij2 + 
                                a2_sqrtPI * exp__a2rij2 / rij - wolfConst) / rij;

                    for (int a = 0; a < 3; ++ a)
                    {
                        f[a] += const_a0 * r[a] * FQQ * qi * qj;
                    }

                    // potential[10] = CPFF metal charge - non-metal charge
                    s_system.potential[10] += FQQ * qi * qj * 
                                     (erfc_arij / rij - erfc_arCut + wolfConst * (rij - rCut));
                }

                else if (0 == use_coulomb)
                {
                    const_a0 = 1.0 / (rij2 * rij);

                    for (int a = 0; a < 3; ++ a)
                    {
                        f[a] += const_a0 * r[a] * FQQ * qi * qj;
                    }

                    s_system.potential[10] += FQQ * qi * qj / rij;
                }

                // damped shifted force between nonmetallic point charge
                // and metallic induced dipole
                //=============================================================================
                // f[b] += -((3.0 * r[a] * r[b] - delta[a][b] * rij2) / pow(rij, 5.0) * 
                //           (erfc_arij + a2_sqrtPI * rij * exp__a2rij2) + 
                //           4.0 * INV_SQRT_PI * alpha3 * r[a] * r[b] / rij2 * exp__a2rij2 -
                //           wolfConst * (r[a] * r[b] - delta[a][b] * rij2) / pow(rij, 3.0)) *
                //          FQQ * qi * vec_pq[j * 3 + a] / NM2BOHR;
                //=============================================================================
                // !!!NOTE: convert dipole moment from a.u.(ea0) to (e nm)

                if (1 == use_coulomb)
                {
                    const_a1 = 1.0 / pow(rij, 5.0) *  
                               (erfc_arij + a2_sqrtPI * rij * exp__a2rij2);
                    const_b1 = 4.0 * INV_SQRT_PI * alpha3 / rij2 * exp__a2rij2;
                    const_c1 = wolfConst / pow(rij, 3.0);

                    for (int a = 0; a < 3; ++ a)
                    {
                        for (int b = 0; b < 3; ++ b)
                        {
                            f[b] += -(const_a1 * (3.0 * r[a] * r[b] - delta[a][b] * rij2) +
                                      const_b1 * r[a] * r[b] -
                                      const_c1 * (r[a] * r[b] - delta[a][b] * rij2)) *
                                     FQQ * qi * vec_pq[j * 3 + a] / NM2BOHR;
                        }
                    }
                }

                else if (0 == use_coulomb)
                {
                    for (int a = 0; a < 3; ++ a)
                    {
                        for (int b = 0; b < 3; ++ b)
                        {
                            f[b] += -((3.0 * r[a] * r[b] - delta[a][b] * rij2) / pow(rij, 5.0)) *
                                     FQQ * qi * vec_pq[j * 3 + a] / NM2BOHR;
                        }
                    }
                }

                // potential[11] = CPFF metal dipole - non-metal charge
                for (int a = 0; a < 3; ++ a)
                {
                    s_system.potential[11] += -FQQ * qi * vec_pq[j * 3 + a] / NM2BOHR *
                                      const_a0 * r[a];
                }

                s_system.fx[i_metal] += f[0];
                s_system.fy[i_metal] += f[1];
                s_system.fz[i_metal] += f[2];
                s_system.fx[i_atom]  -= f[0];
                s_system.fy[i_atom]  -= f[1];
                s_system.fz[i_atom]  -= f[2];

                // virial from metal - non-metal interaction
                for (int a = 0; a < 3; ++ a)
                {
                    for (int b = 0; b < 3; ++ b)
                    {
                        s_system.virial[a][b] += -r[a] * f[b];
                    }
                }
            }

        }

        //========= self-interaction terms ===============
        const int i = i_metal - min_metal;

        // charge-charge
        // inv_capac in Bohr^-1
        // (inv_capac * NM2BOHR) in nm^-1
        qi = vec_pq[n_metal * 3 + i];
        s_system.potential[12] += 0.5 * FQQ * qi * qi * (inv_capac * NM2BOHR);

        // dipole-dipole (only isotropic contribution)
        // inv_polar in Bohr^-3
        // (inv_polar * NM2BOHR) in Bohr^-2 nm^-1
        for (int a = 0; a < 3; ++ a)
        {
            s_system.potential[14] += 0.5 * FQQ * vec_pq[i * 3 + a] * vec_pq[i * 3 + a] * 
                             (inv_polar * NM2BOHR);
        }
        //================================================

        // charge-charge, charge-dipole and dipole-dipole interaction 
        // within metallic nanoparticle
        // !!!NOTE: this part is evaluated in atomic units
        //          and then converted to MD units
            
        // we need FQQ in atomic units...
        // CODATA:
        //     Avogadro constant = 6.02214129e+23 mol^-1
        //     1 Eh = 4.35974434e-18 J
        //     1 Bohr = 0.52917721092e-10 m
        // 1 Eh/a0 = 49614.7526045429 kJ/mol/nm
        // 1 Eh = 2625.49964037578 kJ/mol
        double fQQ_au2 = 49614.7526045429;
        double fQQ_au1 = 2625.49964037578;

        for (int j_metal = 0; j_metal <= max_metal; ++ j_metal)
        {
            const int j = j_metal - min_metal;
            if (j == i) { continue; }

            // distances 
            Vec_R rvec;
            rvec.x = s_system.rx[j_metal] - s_system.rx[i_metal];
            rvec.y = s_system.ry[j_metal] - s_system.ry[i_metal];
            rvec.z = s_system.rz[j_metal] - s_system.rz[i_metal];
            pbc_12(&rvec, s_system.box);
            scale_vec_1(NM2BOHR, &rvec); // nm to atomic unit
            
            rij2 = dist_2(&rvec);
            rij  = sqrt(rij2);

            // induced charges
            qi = vec_pq[n_metal * 3 + i];
            qj = vec_pq[n_metal * 3 + j];

            // component of i->j vector
            double r[3] = {rvec.x, rvec.y, rvec.z};

            // force
            double f[3] = {0., 0., 0.};

            // charge(i) - charge(j)
            //===========================================================================================
            // Note: T(0) contribute to the potential energy
            //
            // f[a] += r[a] / pow(rij, 3.0) *  
            //         (erf(rij * inv_R_qq) - 2.0 * INV_SQRT_PI * inv_R_qq * rij * exp___rijRqq2) * 
            //         fQQ_au2 * qi * qj;
            //===========================================================================================

            double exp___rijRqq2 = exp(-(rij2 * inv_R_qq * inv_R_qq));

            // const_A0 in Bohr^-3
            double const_A0 = 1.0 / pow(rij, 3.0) * 
                              (erf(rij * inv_R_qq) - 2.0 * INV_SQRT_PI * inv_R_qq * rij * exp___rijRqq2);

            for (int a = 0; a < 3; ++ a)
            {
                f[a] += const_A0 * r[a] * fQQ_au2 * qi * qj;
            }

            // potential[12] = CPFF metal charge - metal charge
            s_system.potential[12] += fQQ_au1 * qi * qj * erf(rij * inv_R_qq) / rij * 0.5;

            // charge - dipole
            //===========================================================================================
            // Note: -T(1) contribute to the potential energy
            //
            // f[b] += -((3.0 * r[a] * r[b] - delta[a][b] * rij2) / pow(rij, 5.0) * 
            //           (erf(rij * inv_R_pq) - 2.0 * INV_SQRT_PI * inv_R_pq * rij * exp__rijRpq2) - 
            //           4.0 * INV_SQRT_PI * inv__Rpq3 * r[a] * r[b] / rij2 * exp__rijRpq2) *
            //          fQQ_au2 * qi * vec_pq[j * 3 + a];
            //===========================================================================================

            double exp__rijRpq2 = exp(-(rij2 * inv_R_pq * inv_R_pq));
            double inv__Rpq3 = pow(inv_R_pq, 3.0);

            // const_Apq0 in Bohr^-3
            double const_Apq0 = 1.0 / pow(rij, 3.0) * 
                                (erf(rij * inv_R_pq) - 2.0 * INV_SQRT_PI * inv_R_pq * rij * exp__rijRpq2);

            // const_A1 in Bohr^-5
            // const_B1 in Bohr^-5
            double const_A1 = 1.0 / pow(rij, 5.0) *  
                              (erf(rij * inv_R_pq) - 2.0 * INV_SQRT_PI * inv_R_pq * rij * exp__rijRpq2);
            double const_B1 = 4.0 * INV_SQRT_PI * inv__Rpq3 / rij2 * exp__rijRpq2;

            for (int a = 0; a < 3; ++ a)
            {
                // !!!NOTE: convert dipole moment from a.u.(ea0) to (e nm)
                for (int b = 0; b < 3; ++ b)
                {
                    // charge(i) - dipole(j)
                    f[b] += -(const_A1 * (3.0 * r[a] * r[b] - delta[a][b] * rij2) -
                              const_B1 * r[a] * r[b]) *
                             fQQ_au2 * qi * vec_pq[j * 3 + a];

                    // charge(j) - dipole(i)
                    // note: changing sign due to i-j exchange
                    f[b] += (const_A1 * (3.0 * r[a] * r[b] - delta[a][b] * rij2) -
                             const_B1 * r[a] * r[b] ) *
                            fQQ_au2 * qj * vec_pq[i * 3 + a];
                }

                // potential[13] = CPFF metal charge - metal dipole
                s_system.potential[13] += -fQQ_au1 * qi * vec_pq[j * 3 + a] * r[a] * const_Apq0 * 0.5;

                // potential[13] = CPFF metal charge - metal dipole
                // note: changing sign due to i-j exchange
                s_system.potential[13] +=  fQQ_au1 * qj * vec_pq[i * 3 + a] * r[a] * const_Apq0 * 0.5;
            }

            // dipole(i) - dipole(j)
            //===========================================================================================
            // Note: -T(2) contribute to the potential energy
            //
            // f[c] += -((5.0 * r[a] * r[b] * r[c] - 
            //            rij2 * (r[a] * delta[b][c] + r[b] * delta[c][a] + r[c] * delta[a][b])) / 
            //           pow(rij, 7.0) *
            //           (3.0 * erf(rij * inv_R_pp) - 
            //            6.0 * INV_SQRT_PI * inv_R_pp * rij * exp_rijRpp2 -
            //            4.0 * pow(rij, 3.0) * INV_SQRT_PI * inv_Rpp3 * exp_rijRpp2) -
            //           8.0 * r[a] * r[b] * r[c] * INV_SQRT_PI * pow(inv_R_pp, 5.0) / rij2 * exp_rijRpp2) *
            //          fQQ_au2 * vec_pq[i * 3 + a] * vec_pq[j * 3 + b];
            //===========================================================================================

            double exp_rijRpp2 = exp(-(rij2 * inv_R_pp * inv_R_pp));

            // const_App1 in Bohr^-5
            // const_Bpp1 in Bohr^-5
            double const_App1 = 1.0 / pow(rij, 5.0) *  
                                (erf(rij * inv_R_pp) - 2.0 * INV_SQRT_PI * inv_R_pp * rij * exp_rijRpp2);
            double const_Bpp1 = 4.0 * INV_SQRT_PI * inv_Rpp3 / rij2 * exp_rijRpp2;

            // const_A2 in Bohr^-7
            // const_B2 in Bohr^-7
            double const_A2 = 1.0 / pow(rij, 7.0) *
                              (3.0 * erf(rij * inv_R_pp) - 
                               6.0 * INV_SQRT_PI * inv_R_pp * rij * exp_rijRpp2 -
                               4.0 * pow(rij, 3.0) * INV_SQRT_PI * inv_Rpp3 * exp_rijRpp2);
            double const_B2 = 8.0 * INV_SQRT_PI * pow(inv_R_pp, 5.0) / rij2 * exp_rijRpp2;

            for (int a = 0; a < 3; ++ a)
            {
                for (int b = 0; b < 3; ++ b)
                {
                    // !!!NOTE: convert dipole moment from a.u.(ea0) to (e nm)
                    for (int c = 0; c < 3; ++ c)
                    {
                        f[c] += -(const_A2 * (5.0 * r[a] * r[b] * r[c] - 
                                  rij2 * (r[a] * delta[b][c] + r[b] * delta[c][a] + r[c] * delta[a][b])) -
                                  const_B2 * r[a] * r[b] * r[c]) *
                                 fQQ_au2 * vec_pq[i * 3 + a] * vec_pq[j * 3 + b];
                    }

                    // potential[14] = CPFF metal dipole- metal dipole
                    s_system.potential[14] += -fQQ_au1 * vec_pq[i * 3 + a] * vec_pq[j * 3 + b] *
                                     ((3.0 * r[a] * r[b] - delta[a][b] * rij2) * const_App1 - 
                                      r[a] * r[b] * const_Bpp1) * 0.5;
                }
            }

            // Not using Newton's third law
            // also scaling virial by 0.5
            s_system.fx[i_metal] -= f[0];
            s_system.fy[i_metal] -= f[1];
            s_system.fz[i_metal] -= f[2];

            // virial from metal - metal interaction
            for (int a = 0; a < 3; ++ a)
            {
                for (int b = 0; b < 3; ++ b)
                {
                    s_system.virial[a][b] += -r[a] * f[b] * 0.5;
                }
            }
        }

        //====== external electric field ======
        
        // force on the induced charge
        double q = vec_pq[n_metal * 3 + i];
        s_system.fx[i_metal] += s_runset.external_efield[0] * q;
        s_system.fy[i_metal] += s_runset.external_efield[1] * q;
        s_system.fz[i_metal] += s_runset.external_efield[2] * q;

        // potential on the induced charge
        s_system.potential[9] -= s_runset.external_efield[0] * q * s_system.rx[i_metal] +
                                  s_runset.external_efield[1] * q * s_system.ry[i_metal] +
                                  s_runset.external_efield[2] * q * s_system.rz[i_metal];

        // force on the induced dipole is zero
        // potential on the induced dipole
        double p[3] = {vec_pq[i * 3], vec_pq[i * 3 + 1], vec_pq[i * 3 + 2]};
        s_system.potential[9] -= s_runset.external_efield[0] * p[0] +
                                  s_runset.external_efield[1] * p[1] +
                                  s_runset.external_efield[2] * p[2];
    }
}

//====================================================
// evaluate CPIM energies from matrices and vectors
//====================================================

void eval_cpff_pot(int n_mat, int n_NPs, double** mat_relay, double* vec_ext, double* vec_pq,
                   double* potential)
{
    int n_metal = (n_mat - n_NPs) / 4;
    int i, j;

    double e_qV    = 0.0;
    double e_qCq2 = 0.0;
    for (i = 0; i < n_metal; i ++)
    {
        double Vi = vec_ext[n_metal*3 + i];
        double qi = vec_pq[n_metal*3 + i];
        e_qV += qi * Vi;

        for (j = 0; j < n_metal; j ++)
        {
            double Cij = -mat_relay[n_metal*3 + i][n_metal*3 + j];
            double qj = vec_pq[n_metal*3 + j];
            e_qCq2 += 0.5 * qi * Cij * qj;
        }
    }

    double e__pE  = 0.0;
    double e_pAp2 = 0.0;
    double e__pMq = 0.0;
    for (i = 0; i < n_metal*3; i ++)
    {
        double Ei = vec_ext[i];
        double pi = vec_pq[i];
        e__pE -= pi * Ei;

        for (j = 0; j < n_metal*3; j ++)
        {
            double Aij = mat_relay[i][j];
            double pj = vec_pq[j];
            e_pAp2 += 0.5 * pi * Aij * pj;
        }

        for (j = 0; j < n_metal; j ++)
        {
            double Mij = -mat_relay[i][n_metal*3 + j];
            double qj = vec_pq[n_metal*3 + j];
            e__pMq -= pi * Mij * qj;
        }
    }

    //double Eh2kJmol = 2625.49964037578;
    //potential[10] = e_qV   * Eh2kJmol;
    //potential[11] = e__pE  * Eh2kJmol;
    //potential[12] = e_qCq2 * Eh2kJmol;
    //potential[13] = e__pMq * Eh2kJmol;
    //potential[14] = e_pAp2 * Eh2kJmol;
}