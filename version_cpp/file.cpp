/**
 \file file.cpp
 \author Xin Li
 \date 2015/09
 \version 1.0.2
 
 \brief read input files, write trajectories and charge/dipoles
 
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

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "typedef.hpp"

/// \brief substr in C
void substr(char *dest, const char *src, unsigned int start, unsigned int cnt) 
{
    strncpy(dest, src + start, cnt + 1);
    dest[cnt] = '\0';
}

/// \brief read the next line from file
void read_next_line(FILE *file_par, char *line, char *subline)
{
    while (1)
    {
        if (fgets(line, MAX_STR_LEN, file_par) != NULL)
        {
            // check if line is blank line and
            // save the first word in subline
            strcpy(subline, "");
            sscanf(line, "%s", subline);
            if (strlen(subline) > 0) { break; }
        }
        else
        {
            printf("Error while reading file_par!\n");
            exit(-1);
        }
    }
}

/// \brief read MD settings from input_mdset
void read_settings(std::string input_mdset, RunSet& s_runset, Metal& s_metal, int my_id)
{
    std::ifstream input;
    input.open(input_mdset);
    
    for (std::string line; getline(input, line); )
    {
        std::string option, value;
        std::istringstream(line) >> option >> value;
        
        if      (std::string("run_type") == option) { s_runset.run_type = value; }
        else if (std::string("ensemble") == option) { s_runset.ensemble = value; }
        
        else if (std::string("em_steps"       ) == option) { s_runset.em_steps  = std::stoi(value); }
        else if (std::string("em_length_in_nm") == option) { s_runset.em_length = std::stod(value); }
        else if (std::string("em_tolerance"   ) == option) { s_runset.em_tol    = std::stod(value); }
        
        else if (std::string("number_of_steps") == option) { s_runset.nSteps = std::stoi(value); }
        else if (std::string("time_step_in_ps") == option) { s_runset.dt     = std::stod(value); }
        else if (std::string("save_step"      ) == option) { s_runset.nSave  = std::stoi(value); }
        
        else if (std::string("vdw_type")     == option) { s_runset.vdw_type     = value; }
        else if (std::string("coulomb_type") == option) { s_runset.coulomb_type = value; }
        
        else if (std::string("cut_off_radius") == option) { s_runset.rCut    = std::stod(value); }
        else if (std::string("coulomb_alpha" ) == option) { s_runset.w_alpha = std::stod(value); }
        
        else if (std::string("heating_steps") == option) { s_runset.nHeating = std::stoi(value); }
        else if (std::string("ref_T_in_K"   ) == option) { s_runset.ref_temp = std::stod(value); }
        else if (std::string("tau_T_in_ps"  ) == option) { s_runset.tau_temp = std::stod(value); }
        
        /*
         else if (std::string("ref_P_in_bar") == option) { s_runset.ref_pres = std::stoi(value); }
         else if (std::string("tau_P_in_ps" ) == option) { s_runset.tau_pres = std::stoi(value); }
         */
        
        else if (std::string("fix_metal") == option) { s_metal.fix_pos  = std::stoi(value); }
        else if (std::string("use_cpff" ) == option) { s_metal.use_cpff = std::stoi(value); }
        
        // input external electric field in V/nm, or eV/(e nm)
        // converting to MD units kJ/mol/(e nm)
        else if (std::string("external_Ex") == option) { s_runset.external_efield[0] = std::stod(value) * EV2KJMOL; }
        else if (std::string("external_Ey") == option) { s_runset.external_efield[1] = std::stod(value) * EV2KJMOL; }
        else if (std::string("external_Ez") == option) { s_runset.external_efield[2] = std::stod(value) * EV2KJMOL; }
        
        else if (std::string("") == option) { continue; }
        else
        {
            if (ROOT_PROC == my_id)
            {
                std::cout << ">>> WARNING: unknown option " << option << " in " << input_mdset << std::endl << std::endl;
            }
        }
    }
    
    input.close();
}

/// \brief read input_gro file
void read_gro(std::string input_gro, System& s_system,
              int *ptr_nAtoms, Atom_Info *atom_info)
{
    char line[MAX_STR_LEN], subline[MAX_STR_LEN];

    FILE *gro;
    gro = fopen(input_gro.c_str(), "r") ;
    if (NULL == gro) 
    {
        printf("Cannot open %s for input!\n", input_gro.c_str()) ;
        exit (1);
    }

    // read the first line and get vQ, vP for thermo/barostats
    int step;
    char sub_2[MAX_STR_LEN], sub_3[MAX_STR_LEN];
    s_system.vQ[0] = 0.0;
    s_system.vP = 0.0;
    if (fgets(line, MAX_STR_LEN, gro) != NULL)
    {
        sscanf(line, "%s%d%s%lf%s%lf", 
               subline, &step, sub_2, &s_system.vQ[0], sub_3, &s_system.vP);
    }

    // read the second line for number of atoms
    if (fgets(line, MAX_STR_LEN, gro) != NULL)
    {
        sscanf(line, "%d", ptr_nAtoms);
    }

    // read coordinates
    int i;
    for (i = 0; i < (*ptr_nAtoms); ++ i) 
    {
        if (fgets(line, sizeof(line), gro) != NULL)
        {
            // residue ID
            substr(subline, line, 0, 5);
            sscanf(subline, "%d", &(atom_info[i].resID));

            // residue name and atom name
            substr(subline, line, 5, 5);
            sscanf(subline, "%s", atom_info[i].resName);
            substr(subline, line, 10, 5);
            sscanf(subline, "%s", atom_info[i].atomName);

            // atom coordinates and velocities
            s_system.vx[i] = 0.0;
            s_system.vy[i] = 0.0;
            s_system.vz[i] = 0.0;
            substr(subline, line, 20, MAX_STR_LEN-20);
            sscanf(subline, "%lf%lf%lf%lf%lf%lf", 
                   &s_system.rx[i], &s_system.ry[i], &s_system.rz[i], 
                   &s_system.vx[i], &s_system.vy[i], &s_system.vz[i]);
        }
    }

    // box size
    if (fgets(line, MAX_STR_LEN, gro) != NULL)
    {
        sscanf(line, "%lf%lf%lf", &s_system.box[0], &s_system.box[1], &s_system.box[2]);
    }

    fclose(gro);
}

/// \brief write to traj.gro file
void write_gro(FILE *file_gro, System& s_system,
               int nAtoms, Atom_Info *atom_info, int step)
{
    fprintf(file_gro, "step  %d  vQ  %lf  vP  %lf\n", step, s_system.vQ[0], s_system.vP);
    fprintf(file_gro, "%5d\n", nAtoms);

    for (int i = 0; i < nAtoms; ++ i)
    {
        fprintf(file_gro, "%5d%-5s%5s%5d%10.5f%10.5f%10.5f%10.6f%10.6f%10.6f\n", 
                atom_info[i].resID, atom_info[i].resName, atom_info[i].atomName, i+1,
                s_system.rx[i], s_system.ry[i], s_system.rz[i], 
                s_system.vx[i], s_system.vy[i], s_system.vz[i]);
    }

    fprintf(file_gro, "%10.5f%10.5f%10.5f\n", 
            s_system.box[0], s_system.box[1], s_system.box[2]);
}

/// \brief write metal dipole-charge vector to vec_pq.txt
void write_vec_pq(FILE *file_pq, Metal& s_metal, int step)
{
    int n_metal = s_metal.num;
    fprintf(file_pq, "step %d", step);

    // print Lagrange multipliers
    fprintf(file_pq, "  lambda");
    for (int iNP = 0; iNP < s_metal.n_NPs; ++ iNP)
    {
        fprintf(file_pq, " %.10f", s_metal.vec_pq[n_metal*4 + iNP]);
    }
    fprintf(file_pq, "\n");

    // print number of metal atoms, dipole vectors and charges
    fprintf(file_pq, "%5d\n", n_metal);

    for (int i = 0; i < n_metal; ++ i)
    {
        fprintf(file_pq, "%15.8f%15.8f%15.8f%15.8f\n", 
                s_metal.vec_pq[i*3], s_metal.vec_pq[i*3+1], s_metal.vec_pq[i*3+2], 
                s_metal.vec_pq[n_metal*3+i]);
    }
}


/// \brief write to binary trajectory file
void write_binary(FILE *file_dat, System& s_system, int nAtoms, int step)
{
    fwrite(&nAtoms, sizeof(int), 1, file_dat);
    fwrite(&step,   sizeof(int), 1, file_dat);

    fwrite(s_system.box, sizeof(double),    DIM, file_dat);
    fwrite(s_system.rx,  sizeof(double), nAtoms, file_dat);
    fwrite(s_system.ry,  sizeof(double), nAtoms, file_dat);
    fwrite(s_system.rz,  sizeof(double), nAtoms, file_dat);
}


/// \brief read force field parameters from input_param
void read_param_1(std::string input_param, Topol& s_topol, Metal& s_metal)
{
    FILE* file_par;
    char line[MAX_STR_LEN], subline[MAX_STR_LEN], my_string[MAX_STR_LEN];
    int  mol, atom;

    file_par = fopen(input_param.c_str(), "r") ;
    if (NULL == file_par)
    {
        printf("Error: Cannot open param file %s!\n", input_param.c_str()) ;
        exit(-1);
    }

    // 1-4 scaling factors
    read_next_line(file_par, line, subline);
    if (0 == strcmp(subline, "SCALE_LJ_QQ"))
    {
        sscanf(line, "%s%lf%lf", subline, &(s_topol.scaleLJ), &(s_topol.scaleQQ));
    }
    else
    {
        printf("Error in %s: expected SCALE_LJ_QQ in line\n%s", input_param.c_str(), line);
        exit(-1);
    }

    // number of molecule types
    read_next_line(file_par, line, subline);
    if (0 == strcmp(subline, "MOLECULES"))
    {
        sscanf(line, "%s%d", subline, &(s_topol.mol_types));
    }
    else
    {
        printf("Error in %s: expected MOLECULES in line\n%s", input_param.c_str(), line);
        exit(-1);
    }


    // molecule name
    char molName[s_topol.mol_types][MAX_STR_LEN];


    // molNr section: number of this type of molecule in the system
    for (mol = 0; mol < s_topol.mol_types; ++ mol) 
    {
        read_next_line(file_par, line, subline);
        sscanf(line, "%s", molName[mol]);
    }

    // atom_in_mol section: number of atoms in this type of molecule
    for (mol = 0; mol < s_topol.mol_types; ++ mol) 
    {
        read_next_line(file_par, line, subline);

        int n_atoms;
        sscanf(line, "%s%d", subline, &n_atoms);

        strcpy(my_string, molName[mol]);
        strcat(my_string, "_ATOMS");

        if (0 == strcmp(subline, my_string))
        {
            for (atom = 0; atom < n_atoms; ++ atom) 
            {
                read_next_line(file_par, line, subline);
            }
        }
        else
        {
            printf("Error in %s: expected %s in line\n%s", input_param.c_str(), my_string, line);
            exit(-1);
        }
    }


    // nonbonded parameters
    int iType, jType;

    read_next_line(file_par, line, subline);
    sscanf(line, "%s%d", subline, &(s_topol.n_types));

    if (0 == strcmp(subline, "NONBONDED"))
    {
        for (iType = 0; iType < s_topol.n_types; ++ iType)
        {
            for (jType = iType; jType < s_topol.n_types; ++ jType)
            {
                read_next_line(file_par, line, subline);
            }
        }
    }
    else
    {
        printf("Error in %s: expected NONBONDED in line\n%s", input_param.c_str(), line);
        exit(-1);
    }

    // bonded potentials: bond, pair, angle, dihedral
    for (mol = 0; mol < s_topol.mol_types; ++ mol) 
    {
        int num_bonds, num_pairs, num_angles, num_dihedrals, num_vsites, num_constraints;
        int iBond, iPair, iAngle, iDihedral, iVSite, iCstr;

        // bonds
        read_next_line(file_par, line, subline);
        sscanf(line, "%s%d", subline, &num_bonds);

        strcpy(my_string, molName[mol]);
        strcat(my_string, "_BONDS");

        if (0 == strcmp(subline, my_string))
        {
            for (iBond = 0; iBond < num_bonds; ++ iBond)
            {
                read_next_line(file_par, line, subline);
            }
        }
        else
        {
            printf("Error in %s: expected %s in line\n%s", input_param.c_str(), my_string, line);
            exit(-1);
        }

        // pairs
        read_next_line(file_par, line, subline);
        sscanf(line, "%s%d", subline, &num_pairs);

        strcpy(my_string, molName[mol]);
        strcat(my_string, "_PAIRS");

        if (0 == strcmp(subline, my_string))
        {
            for (iPair = 0; iPair < num_pairs; ++ iPair)
            {
                read_next_line(file_par, line, subline);
            }
        }
        else
        {
            printf("Error in %s: expected %s in line\n%s", input_param.c_str(), my_string, line);
            exit(-1);
        }

        // angles
        read_next_line(file_par, line, subline);
        sscanf(line, "%s%d", subline, &num_angles);

        strcpy(my_string, molName[mol]);
        strcat(my_string, "_ANGLES");

        if (0 == strcmp(subline, my_string))
        {
            for (iAngle = 0; iAngle < num_angles; ++ iAngle)
            {
                read_next_line(file_par, line, subline);
            }
        }
        else
        {
            printf("Error in %s: expected %s in line\n%s", input_param.c_str(), my_string, line);
            exit(-1);
        }

        // dihedrals
        read_next_line(file_par, line, subline);
        sscanf(line, "%s%d", subline, &num_dihedrals);

        strcpy(my_string, molName[mol]);
        strcat(my_string, "_DIHEDRALS");

        if (0 == strcmp(subline, my_string))
        {
            for (iDihedral = 0; iDihedral < num_dihedrals; ++ iDihedral)
            {
                read_next_line(file_par, line, subline);
            }
        }
        else
        {
            printf("Error in %s: expected %s in line\n%s", input_param.c_str(), my_string, line);
            exit(-1);
        }

        // virtual_sites
        read_next_line(file_par, line, subline);
        sscanf(line, "%s%d", subline, &num_vsites);

        strcpy(my_string, molName[mol]);
        strcat(my_string, "_VSITES");

        if (0 == strcmp(subline, my_string))
        {
            for (iVSite = 0; iVSite < num_vsites; ++ iVSite)
            {
                read_next_line(file_par, line, subline);
            }
        }
        else
        {
            printf("Error in %s: expected %s in line\n%s", input_param.c_str(), my_string, line);
            exit(-1);
        }

        // constraints
        read_next_line(file_par, line, subline);
        sscanf(line, "%s%d", subline, &num_constraints);

        strcpy(my_string, molName[mol]);
        strcat(my_string, "_CSTRS");

        if (0 == strcmp(subline, my_string))
        {
            for (iCstr = 0; iCstr < num_constraints; ++ iCstr)
            {
                read_next_line(file_par, line, subline);
            }
        }
        else
        {
            printf("Error in %s: expected %s in line\n%s", input_param.c_str(), my_string, line);
            exit(-1);
        }
    }

    // metal indices (min, max)
    read_next_line(file_par, line, subline);
    sscanf(line, "%s%d%d", subline, &(s_metal.min), &(s_metal.max));

    if (0 == strcmp(subline, "MM_METAL"))
    {
        if (s_metal.min >=0 && s_metal.max >= s_metal.min)
        {
            s_metal.num = s_metal.max - s_metal.min + 1;
        }
        else
        {
            s_metal.num = 0;
        }
    }
    else
    {
        printf("Error in %s: expected MM_METAL in line\n%s", input_param.c_str(), line);
        exit(-1);
    }

    // Quantum Sutton-Chen parameters
    read_next_line(file_par, line, subline);
    sscanf(line, "%s%d%d%lf%lf%lf", subline, 
           &(s_metal.qsc_n), &(s_metal.qsc_m),
           &(s_metal.qsc_eps), &(s_metal.qsc_c), &(s_metal.qsc_a));
    if (0 != strcmp(subline, "QSC_METAL"))
    {
        printf("Error in %s: expected QSC_METAL in line\n%s", input_param.c_str(), line);
        exit(-1);
    }

    // atomic polarizability and capacitance, in atomic unit
    read_next_line(file_par, line, subline);
    sscanf(line, "%s%lf%lf", subline, &(s_metal.cpff_polar), &(s_metal.cpff_capac));
    if (0 != strcmp(subline, "CPFF_METAL"))
    {
        printf("Error in %s: expected CPFF_METAL in line\n%s", input_param.c_str(), line);
        exit(-1);
    }

    // multi-NPs
    read_next_line(file_par, line, subline);
    sscanf(line, "%s%d", subline, &(s_metal.n_NPs));
    if (0 != strcmp(subline, "NANOPARTICLES"))
    {
        printf("Error in %s: expected NANOPARTICLES in line\n%s", input_param.c_str(), line);
        exit(-1);
    }

    // close param file
    fclose( file_par );
}

/// \brief read FF parameters (2) from input_param
void read_param_2(std::string input_param, Topol& s_topol, Metal& s_metal)
{
    FILE *file_par;
    char line[MAX_STR_LEN], tmp[MAX_STR_LEN];

    file_par = fopen(input_param.c_str(),"r");
    if (NULL == file_par) 
    {
        printf("Error: Cannot open param file %s!\n", input_param.c_str());
        exit(-1);
    }

    // 1-4 scaling factors
    read_next_line(file_par, line, tmp);

    // number of molecule types
    read_next_line(file_par, line, tmp);

    // s_topol.mol_num = number of this type of molecule in the system
    for (int mol = 0; mol < s_topol.mol_types; ++ mol)
    {
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d", tmp, &s_topol.mol_num[mol]);
    }

    // s_topol.atom_num = number of atoms in this type of molecule
    // s_topol.atom_param = atomic parameters (q, m, sig, eps, type)
    for (int mol = 0; mol < s_topol.mol_types; ++ mol)
    {
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d", tmp, &s_topol.atom_num[mol]);

        s_topol.atom_param[mol] = new (std::nothrow) AtomParam [s_topol.atom_num[mol]];
        assert(s_topol.atom_param[mol] != nullptr);
        
        for (int atom = 0; atom < s_topol.atom_num[mol]; ++ atom)
        {
            read_next_line(file_par, line, tmp);
            sscanf(line, "%lf%lf%d",
                   &(s_topol.atom_param[mol][atom].charge),
                   &(s_topol.atom_param[mol][atom].mass),
                   &(s_topol.atom_param[mol][atom].atomtype));
        }
    }

    // nonbonded parameters
    read_next_line(file_par, line, tmp);
    //sscanf(line, "%s%d", tmp, &s_topol.n_types);

    int iType, jType, funct;
    for (iType = 0; iType < s_topol.n_types; ++ iType)
    {
        for (jType = iType; jType < s_topol.n_types; ++ jType)
        {
            read_next_line(file_par, line, tmp);

            sscanf(line, "%d", &funct);

            // LJ
            if (1 == funct)
            {
                sscanf(line, "%d%le%le", 
                       &(s_topol.nonbonded_param[iType][jType].funct),
                       &(s_topol.nonbonded_param[iType][jType].C6),
                       &(s_topol.nonbonded_param[iType][jType].C12));
            }

            // Buckingham
            else if (2 == funct)
            {
                sscanf(line, "%d%le%le%le", 
                       &(s_topol.nonbonded_param[iType][jType].funct),
                       &(s_topol.nonbonded_param[iType][jType].A),
                       &(s_topol.nonbonded_param[iType][jType].B),
                       &(s_topol.nonbonded_param[iType][jType].C6));
            }

            // Morse
            else if (3 == funct)
            {
                sscanf(line, "%d%le%le%le", 
                       &(s_topol.nonbonded_param[iType][jType].funct),
                       &(s_topol.nonbonded_param[iType][jType].A),
                       &(s_topol.nonbonded_param[iType][jType].B),
                       &(s_topol.nonbonded_param[iType][jType].C6));
            }

            // erf_vdw
            else if (4 == funct)
            {
                sscanf(line, "%d%le%le%le", 
                       &(s_topol.nonbonded_param[iType][jType].funct),
                       &(s_topol.nonbonded_param[iType][jType].C6),
                       &(s_topol.nonbonded_param[iType][jType].C12),
                       &(s_topol.nonbonded_param[iType][jType].A));
            }

            // error
            else
            {
                printf("Error: unsupported vdW type: %d in line:\n%s\n", 
                       funct, line );
                exit(-1);
            }
        }
    }

    for (iType = 0; iType < (s_topol.n_types - 1); ++ iType)
    {
        for (jType = iType + 1; jType < s_topol.n_types; ++ jType)
        {
            s_topol.nonbonded_param[jType][iType].funct = 
                s_topol.nonbonded_param[iType][jType].funct;

            s_topol.nonbonded_param[jType][iType].A = 
                s_topol.nonbonded_param[iType][jType].A;

            s_topol.nonbonded_param[jType][iType].B = 
                s_topol.nonbonded_param[iType][jType].B;

            s_topol.nonbonded_param[jType][iType].C6 = 
                s_topol.nonbonded_param[iType][jType].C6;

            s_topol.nonbonded_param[jType][iType].C12 = 
                s_topol.nonbonded_param[iType][jType].C12;
        }
    }

    // bonded potentials: bond, pair, angle, dihedral
    for (int mol = 0; mol < s_topol.mol_types; ++ mol)
    {
        int atom_i, atom_j, atom_k, atom_l, atom_s;
        int iBond, iPair, iAngle, iDihedral, iVSite, iCstr;
        int funct;

        // exclusions
        s_topol.exclude[mol] = new (std::nothrow) int* [s_topol.atom_num[mol]];
        assert(s_topol.exclude[mol] != nullptr);
        
        for (atom_i = 0; atom_i < s_topol.atom_num[mol]; ++ atom_i) 
        {
            s_topol.exclude[mol][atom_i] = new (std::nothrow) int [s_topol.atom_num[mol]];
            assert(s_topol.exclude[mol][atom_i] != nullptr);
            
            for (atom_j = 0; atom_j < s_topol.atom_num[mol]; ++ atom_j) 
            {
                s_topol.exclude[mol][atom_i][atom_j] = 0;
            }
        }

        // bonds
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d", tmp, &s_topol.n_bonds[mol]);

        s_topol.bond_param[mol] = new (std::nothrow) BondParam [s_topol.n_bonds[mol]];
        assert(s_topol.bond_param[mol] != nullptr);

        for (iBond = 0; iBond < s_topol.n_bonds[mol]; ++ iBond)
        {
            read_next_line(file_par, line, tmp);
            sscanf(line, "%d%d%d", &atom_i, &atom_j, &funct);

            // Harmonic bond potential
            if (1 == funct)
            {
                sscanf(line, "%d%d%d%le%le",
                       &(s_topol.bond_param[mol][iBond].atom_i), 
                       &(s_topol.bond_param[mol][iBond].atom_j),
                       &(s_topol.bond_param[mol][iBond].funct),
                       &(s_topol.bond_param[mol][iBond].b0), 
                       &(s_topol.bond_param[mol][iBond].kb));
            }

            // Morse bond potential
            else if (3 == funct)
            {
                sscanf(line, "%d%d%d%le%le%le",
                       &(s_topol.bond_param[mol][iBond].atom_i), 
                       &(s_topol.bond_param[mol][iBond].atom_j),
                       &(s_topol.bond_param[mol][iBond].funct),
                       &(s_topol.bond_param[mol][iBond].b0), 
                       &(s_topol.bond_param[mol][iBond].D), 
                       &(s_topol.bond_param[mol][iBond].beta));
            }

            // bond type 5: exclusion only
            else if (5 == funct)
            {
                // do nothing
            }

            else
            {
                printf("Error: unsupported bond function type: %d in line:\n%s\n", 
                       funct, line );
                exit(-1);
            }

            // update exclusion
            s_topol.exclude[mol][atom_i][atom_j] = 1;
            s_topol.exclude[mol][atom_j][atom_i] = 1;
        }

        // pairs
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d", tmp, &s_topol.n_pairs[mol]);

        s_topol.pair_param[mol] = new (std::nothrow) PairParam [s_topol.n_pairs[mol]];
        assert(s_topol.pair_param[mol] != nullptr);

        for (iPair = 0; iPair < s_topol.n_pairs[mol]; ++ iPair)
        {
            read_next_line(file_par, line, tmp);
            sscanf(line, "%d%d%d", &atom_i, &atom_j, &funct);

            // normal 1-4 pair
            if (1 == funct)
            {
                sscanf(line, "%d%d%d",
                       &(s_topol.pair_param[mol][iPair].atom_i), 
                       &(s_topol.pair_param[mol][iPair].atom_j),
                       &(s_topol.pair_param[mol][iPair].funct));
            }

            // user-defined 1-4 pair
            else if (2 == funct)
            {
                sscanf(line, "%d%d%d%lf%lf%lf%le%le",
                       &(s_topol.pair_param[mol][iPair].atom_i), 
                       &(s_topol.pair_param[mol][iPair].atom_j),
                       &(s_topol.pair_param[mol][iPair].funct),
                       &(s_topol.pair_param[mol][iPair].scaleQQ),
                       &(s_topol.pair_param[mol][iPair].q_i), 
                       &(s_topol.pair_param[mol][iPair].q_j),
                       &(s_topol.pair_param[mol][iPair].V), 
                       &(s_topol.pair_param[mol][iPair].W));
            }

            else
            {
                printf("Error: unsupported pair function type: %d in line:\n%s\n", 
                        funct, line);
                exit(-1);
            }

            // update exclusion
            s_topol.exclude[mol][atom_i][atom_j] = 1;
            s_topol.exclude[mol][atom_j][atom_i] = 1;
        }

        // angles
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d", tmp, &s_topol.n_angles[mol]);

        s_topol.angle_param[mol] = new (std::nothrow) AngleParam [s_topol.n_angles[mol]];
        assert(s_topol.angle_param[mol] != nullptr);

        for (iAngle = 0; iAngle < s_topol.n_angles[mol]; ++ iAngle)
        {
            read_next_line(file_par, line, tmp);
            sscanf(line, "%d%d%d%d", &atom_i, &atom_j, &atom_k, &funct);

            // harmonic angle potential
            if (1 == funct)
            {
                sscanf(line, "%d%d%d%d%le%le",
                       &(s_topol.angle_param[mol][iAngle].atom_i), 
                       &(s_topol.angle_param[mol][iAngle].atom_j), 
                       &(s_topol.angle_param[mol][iAngle].atom_k),
                       &(s_topol.angle_param[mol][iAngle].funct),
                       &(s_topol.angle_param[mol][iAngle].a0), 
                       &(s_topol.angle_param[mol][iAngle].ka));
            }

            // Urey-Bradley type
            else if (5 == funct)
            {
                sscanf(line, "%d%d%d%d%le%le%le%le",
                       &(s_topol.angle_param[mol][iAngle].atom_i), 
                       &(s_topol.angle_param[mol][iAngle].atom_j), 
                       &(s_topol.angle_param[mol][iAngle].atom_k),
                       &(s_topol.angle_param[mol][iAngle].funct),
                       &(s_topol.angle_param[mol][iAngle].a0), 
                       &(s_topol.angle_param[mol][iAngle].ka),
                       &(s_topol.angle_param[mol][iAngle].b0_ub), 
                       &(s_topol.angle_param[mol][iAngle].kb_ub));
            }

            else
            {
                printf("Error: unsupported angle function type: %d in line:\n%s\n", 
                        funct, line);
                exit(-1);
            }

            // update exclusion
            s_topol.exclude[mol][atom_i][atom_j] = 1;
            s_topol.exclude[mol][atom_j][atom_i] = 1;
            s_topol.exclude[mol][atom_i][atom_k] = 1;
            s_topol.exclude[mol][atom_k][atom_i] = 1;
            s_topol.exclude[mol][atom_j][atom_k] = 1;
            s_topol.exclude[mol][atom_k][atom_j] = 1;
        }

        // dihedrals
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d", tmp, &s_topol.n_dihedrals[mol]);

        s_topol.dihedral_param[mol] = new (std::nothrow) DihedralParam
                                      [s_topol.n_dihedrals[mol]];
        assert(s_topol.dihedral_param[mol] != nullptr);
        
        for (iDihedral = 0; iDihedral < s_topol.n_dihedrals[mol]; ++ iDihedral)
        {
            read_next_line(file_par, line, tmp);
            sscanf(line, "%d%d%d%d%d", &atom_i, &atom_j, &atom_k, &atom_l, &funct);
            
            // RB type dihedral potential
            if (3 == funct)
            {
                sscanf(line, "%d%d%d%d%d%lf%lf%lf%lf%lf%lf",
                       &(s_topol.dihedral_param[mol][iDihedral].atom_i), 
                       &(s_topol.dihedral_param[mol][iDihedral].atom_j),
                       &(s_topol.dihedral_param[mol][iDihedral].atom_k), 
                       &(s_topol.dihedral_param[mol][iDihedral].atom_l),
                       &(s_topol.dihedral_param[mol][iDihedral].funct),
                       &(s_topol.dihedral_param[mol][iDihedral].c0), 
                       &(s_topol.dihedral_param[mol][iDihedral].c1),
                       &(s_topol.dihedral_param[mol][iDihedral].c2), 
                       &(s_topol.dihedral_param[mol][iDihedral].c3),
                       &(s_topol.dihedral_param[mol][iDihedral].c4), 
                       &(s_topol.dihedral_param[mol][iDihedral].c5));
            }

            // periodic dihedral potential
            else if (1 == funct || 4 == funct || 9 == funct)
            {
                sscanf(line, "%d%d%d%d%d%lf%lf%d",
                       &(s_topol.dihedral_param[mol][iDihedral].atom_i), 
                       &(s_topol.dihedral_param[mol][iDihedral].atom_j),
                       &(s_topol.dihedral_param[mol][iDihedral].atom_k), 
                       &(s_topol.dihedral_param[mol][iDihedral].atom_l),
                       &(s_topol.dihedral_param[mol][iDihedral].funct),
                       &(s_topol.dihedral_param[mol][iDihedral].phi0), 
                       &(s_topol.dihedral_param[mol][iDihedral].kphi),
                       &(s_topol.dihedral_param[mol][iDihedral].n));
            }

            else
            {
                printf("Error: unsupported dihedral function type: %d in line:\n%s\n", 
                       funct, line);
                exit(-1);
            }

            // update exclusion
            s_topol.exclude[mol][atom_i][atom_j] = 1;
            s_topol.exclude[mol][atom_j][atom_i] = 1;
            s_topol.exclude[mol][atom_i][atom_k] = 1;
            s_topol.exclude[mol][atom_k][atom_i] = 1;
            s_topol.exclude[mol][atom_i][atom_l] = 1;
            s_topol.exclude[mol][atom_l][atom_i] = 1;
            s_topol.exclude[mol][atom_j][atom_k] = 1;
            s_topol.exclude[mol][atom_k][atom_j] = 1;
            s_topol.exclude[mol][atom_j][atom_l] = 1;
            s_topol.exclude[mol][atom_l][atom_j] = 1;
            s_topol.exclude[mol][atom_k][atom_l] = 1;
            s_topol.exclude[mol][atom_l][atom_k] = 1;
        }


        // virtual_site, TIP5P type, funct=4
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d", tmp, &s_topol.n_vsites[mol]);

        s_topol.vsite_4[mol] = new (std::nothrow) VSite_4 [s_topol.n_vsites[mol]];
        s_topol.vsite_funct[mol] = new (std::nothrow) int [s_topol.n_vsites[mol]];
        assert(s_topol.vsite_4[mol] != nullptr);
        assert(s_topol.vsite_funct[mol] != nullptr);

        for (iVSite = 0; iVSite < s_topol.n_vsites[mol]; ++ iVSite)
        {
            read_next_line(file_par, line, tmp);
            sscanf(line, "%d%d%d%d%d", &atom_i, &atom_j, &atom_k, &atom_s, &funct);

            if (4 == funct)
            {
                sscanf(line, "%d%d%d%d%d%lf%lf%lf",
                       &(s_topol.vsite_4[mol][iVSite].atom_i), 
                       &(s_topol.vsite_4[mol][iVSite].atom_j),
                       &(s_topol.vsite_4[mol][iVSite].atom_k), 
                       &(s_topol.vsite_4[mol][iVSite].atom_s),
                       &(s_topol.vsite_funct[mol][iVSite]),    
                       &(s_topol.vsite_4[mol][iVSite].a), 
                       &(s_topol.vsite_4[mol][iVSite].b),      
                       &(s_topol.vsite_4[mol][iVSite].c));

                s_topol.exclude[mol][atom_s][atom_i] = 1;
                s_topol.exclude[mol][atom_s][atom_j] = 1;
                s_topol.exclude[mol][atom_s][atom_k] = 1;
                s_topol.exclude[mol][atom_i][atom_s] = 1;
                s_topol.exclude[mol][atom_j][atom_s] = 1;
                s_topol.exclude[mol][atom_k][atom_s] = 1;
            }
            else
            {
                printf("Error: unsupported virtual_sites function type: %d in line:\n%s\n", 
                       funct, line);
                exit(-1);
            }
        }


        // s_topol.constraints
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d", tmp, &s_topol.n_constraints[mol]);

        s_topol.constraint[mol] = new (std::nothrow) Constraint [s_topol.n_constraints[mol]];
        assert(s_topol.constraint[mol] != nullptr);

        for (iCstr = 0; iCstr < s_topol.n_constraints[mol]; ++ iCstr)
        {
            read_next_line(file_par, line, tmp);
            sscanf(line, "%d%d", &atom_i, &atom_j);

            sscanf(line, "%d%d%lf",
                   &(s_topol.constraint[mol][iCstr].atom_i), 
                   &(s_topol.constraint[mol][iCstr].atom_j),
                   &(s_topol.constraint[mol][iCstr].r));

            s_topol.exclude[mol][atom_i][atom_j] = 1;
            s_topol.exclude[mol][atom_j][atom_i] = 1;
        }
    }

    // metal indices and QSC parameters
    read_next_line(file_par, line, tmp);
    read_next_line(file_par, line, tmp);

    // atomic polarizability and capacitance, in atomic unit
    read_next_line(file_par, line, tmp);

    // number of nanoparticles
    read_next_line(file_par, line, tmp);
    // indices and total charge for each nanoparticle
    int i_NP;
    for (i_NP = 0; i_NP < s_metal.n_NPs; ++ i_NP)
    {
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d%d%lf", 
               tmp, &s_metal.start_NP[i_NP], &s_metal.end_NP[i_NP], &s_metal.cpff_chg[i_NP]);
    }

    // close param file
    fclose(file_par);

    // calculate total number of atoms and molecules
    s_topol.n_atoms = 0;
    s_topol.n_mols  = 0;
    for (int mol = 0; mol < s_topol.mol_types; ++ mol)
    {
        int num_atom = s_topol.atom_num[mol];
        int num_mol  = s_topol.mol_num[mol];
        s_topol.n_atoms += num_mol * num_atom;
        s_topol.n_mols  += num_mol;
    }
}

/// \brief assign molecular and atomic indices
void assign_indices(Topol& s_topol, Metal& s_metal,
                    Atom_Info* atom_info, Mol_Info* mol_info)
{
    int i, im, mol, km, atom;

    // assign atom index and mol_type index for each atom in a molecule type
    // assign starting index and ending index for each molecule in the system 
    // i = atom index, im = mol index
    i  = 0;
    im = 0;
    for (mol = 0; mol < s_topol.mol_types; ++ mol) 
    {
        for (km = 0; km< s_topol.mol_num[mol]; ++ km) 
        {
            mol_info[im].mini = i;

            for (atom = 0; atom < s_topol.atom_num[mol]; ++ atom) 
            {
                atom_info[i].iMol  = mol;
                atom_info[i].iAtom = atom;
                atom_info[i].molID = im;

                i++;
            }

            mol_info[im].maxi = i-1;

            im++;
        }
    }

    if (s_topol.n_atoms != i || s_topol.n_mols != im) 
    {
        printf("Error: Incorrect number of atoms or molecules in the system! (1st)\n");
        printf("n_atoms = %d, i = %d\n", s_topol.n_atoms, i);
        printf("n_mols = %d, im = %d\n", s_topol.n_mols, im);
        exit(-1);
    }


    // assign atom and residue, mass and charge
    for (i = 0; i < s_topol.n_atoms; ++ i)
    {
        int this_iMol  = atom_info[i].iMol;
        int this_iAtom = atom_info[i].iAtom;

        atom_info[i].mass     = s_topol.atom_param[this_iMol][this_iAtom].mass;
        atom_info[i].charge   = s_topol.atom_param[this_iMol][this_iAtom].charge;
        atom_info[i].atomtype = s_topol.atom_param[this_iMol][this_iAtom].atomtype;

        if (atom_info[i].mass > 0.0)
        {
            atom_info[i].inv_mass = 1.0 / atom_info[i].mass;
        }
        else
        {
            atom_info[i].inv_mass = 0.0;
        }

        atom_info[i].is_metal = 0;
    }


    // check metal atom
    for (i = s_metal.min; i <= s_metal.max; ++ i)
    {
        atom_info[i].is_metal = 1;
    }
}