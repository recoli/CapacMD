/*! \file file.c                                                         
    \author Xin Li
    \date 2015/09
    \version 1.0

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "typedef.h"
#include "my_malloc.h"

/// \brief substr in C
/// \fn void substr(char *dest, const char *src, unsigned int start, unsigned int cnt) 
/// \param dest destination string
/// \param src source string
/// \param start starting position in source string
/// \param cnt length of destination string
void substr(char *dest, const char *src, unsigned int start, unsigned int cnt) 
{
    strncpy(dest, src + start, cnt + 1);
    dest[cnt] = '\0';
}

/// \brief read the next line from file
/// \fn void read_next_line(FILE *file_par, char *line, char *subline)
/// \param file_par file handle of parameter file
/// \param line the line read from parameter file
/// \param subline first word in the read line
void read_next_line(FILE *file_par, char *line, char *subline)
{
    while(1)
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
            exit(1);
        }
    }
}

/// \brief read md options from input_mdset
/// \fn void read_settings(char *input_mdset, RunSet *p_runset, Metal *p_metal)
/// \param input_mdset name of the mdset input file
/// \param p_runset pointer to RunSet data structure
/// \param p_metal pointer to Metal data structure
void read_settings(char *input_mdset, RunSet *p_runset, Metal *p_metal)
{
    FILE* f_set;
    char line[MAX_STR_LEN], option[MAX_STR_LEN], value[MAX_STR_LEN];

    f_set = fopen(input_mdset, "r") ;
    if (NULL == f_set)
    {
        printf("Cannot open %s for input!\n", input_mdset);
        exit(1);
    }

    while (fgets(line, MAX_STR_LEN, f_set) != NULL)
    {
        sscanf(line, "%s%s", option, value);

        if      (0 == strcmp(option, "run_type")) { strcpy(p_runset->run_type, value); }
        else if (0 == strcmp(option, "ensemble")) { strcpy(p_runset->ensemble, value); }

        else if (0 == strcmp(option, "em_steps"       )) { p_runset->em_steps  = atoi(value); }
        else if (0 == strcmp(option, "em_length_in_nm")) { p_runset->em_length = atof(value); }
        else if (0 == strcmp(option, "em_tolerance"   )) { p_runset->em_tol    = atof(value); }

        else if (0 == strcmp(option, "number_of_steps")) { p_runset->nSteps = atoi(value); }
        else if (0 == strcmp(option, "time_step_in_ps")) { p_runset->dt     = atof(value); }
        else if (0 == strcmp(option, "save_step"      )) { p_runset->nSave  = atoi(value); }

        else if (0 == strcmp(option, "vdw_type"    )) { strcpy(p_runset->vdw_type,     value); }
        else if (0 == strcmp(option, "coulomb_type")) { strcpy(p_runset->coulomb_type, value); }

        else if (0 == strcmp(option, "cut_off_radius")) { p_runset->rCut    = atof(value); }
        else if (0 == strcmp(option, "coulomb_alpha" )) { p_runset->w_alpha = atof(value); }

        else if (0 == strcmp(option, "heating_steps")) { p_runset->nHeating = atoi(value); }
        else if (0 == strcmp(option, "ref_T_in_K"   )) { p_runset->ref_temp = atof(value); }
        else if (0 == strcmp(option, "tau_T_in_ps"  )) { p_runset->tau_temp = atof(value); }

        /*
        else if (0 == strcmp(option, "ref_P_in_bar")) { p_runset->ref_pres = atof(value); }
        else if (0 == strcmp(option, "tau_P_in_ps" )) { p_runset->tau_pres = atof(value); }
        */

        else if (0 == strcmp(option, "fix_metal")) { p_metal->fix_pos  = atoi(value); }
        else if (0 == strcmp(option, "use_cpff" )) { p_metal->use_cpff = atoi(value); }

        // input external electric field in V/nm, or eV/(e nm)
        // converting to MD units kJ/mol/(e nm)
        else if (0 == strcmp(option, "external_Ex")) { p_runset->external_efield[0] = atof(value) * EV2KJMOL; }
        else if (0 == strcmp(option, "external_Ey")) { p_runset->external_efield[1] = atof(value) * EV2KJMOL; }
        else if (0 == strcmp(option, "external_Ez")) { p_runset->external_efield[2] = atof(value) * EV2KJMOL; }

        else if (0 == strcmp(option, "")) { continue; }
        else    { printf("Warning: unkonwn option %s in %s\n", option, input_mdset); }
    }

    fclose(f_set);
}

/// \brief read input_gro file
/// \fn void read_gro(char *input_gro, System *p_system, int *ptr_nAtoms, Atom_Info *atom_info)
/// \param input_gro name of gro input file
/// \param p_system pointer to System data structure
/// \param ptr_nAtoms pointer to number of atoms
/// \param atom_info array of atomic information
void read_gro(char *input_gro, System *p_system,
              int *ptr_nAtoms, Atom_Info *atom_info)
{
    char line[MAX_STR_LEN], subline[MAX_STR_LEN];

    FILE *gro;
    gro = fopen(input_gro, "r") ;
    if (NULL == gro) 
    {
        printf("Cannot open %s for input!\n", input_gro) ;
        exit (1);
    }

    // read the first line and get vQ, vP for thermo/barostats
    int step;
    char sub_2[MAX_STR_LEN], sub_3[MAX_STR_LEN];
    p_system->vQ[0] = 0.0;
    p_system->vP = 0.0;
    if (fgets(line, MAX_STR_LEN, gro) != NULL)
    {
        sscanf(line, "%s%d%s%lf%s%lf", 
                subline, &step, sub_2, &p_system->vQ[0], sub_3, &p_system->vP);
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
            p_system->vx[i] = 0.0;
            p_system->vy[i] = 0.0;
            p_system->vz[i] = 0.0;
            substr(subline, line, 20, MAX_STR_LEN-20);
            sscanf(subline, "%lf%lf%lf%lf%lf%lf", 
                   &p_system->rx[i], &p_system->ry[i], &p_system->rz[i], 
                   &p_system->vx[i], &p_system->vy[i], &p_system->vz[i]);
        }
    }

    // box size
    if (fgets(line, MAX_STR_LEN, gro) != NULL)
    {
        sscanf(line, "%lf%lf%lf", 
                &p_system->box[0], &p_system->box[1], &p_system->box[2]);
    }

    fclose(gro);
}

/// \brief write to traj.gro file
/// \fn void write_gro(FILE *file_gro, System *p_system, int nAtoms, Atom_Info *atom_info, int step)
/// \param file_gro name of gro output file
/// \param p_system pointer to System data structure
/// \param nAtoms number of atoms
/// \param atom_info array of atomic information
/// \param step number of MD step
void write_gro(FILE *file_gro, System *p_system,
               int nAtoms, Atom_Info *atom_info, int step)
{
    fprintf(file_gro, "step  %d  vQ  %lf  vP  %lf\n", step, p_system->vQ[0], p_system->vP);
    fprintf(file_gro, "%5d\n", nAtoms);

    int i;
    for (i = 0; i < nAtoms; ++ i) 
    {
        fprintf(file_gro, "%5d%-5s%5s%5d%10.5f%10.5f%10.5f%10.6f%10.6f%10.6f\n", 
                atom_info[i].resID, atom_info[i].resName, atom_info[i].atomName, i+1,
                p_system->rx[i], p_system->ry[i], p_system->rz[i], 
                p_system->vx[i], p_system->vy[i], p_system->vz[i]);
    }

    fprintf(file_gro, "%10.5f%10.5f%10.5f\n", 
            p_system->box[0], p_system->box[1], p_system->box[2]);
}

/// \brief write metal dipole-charge vector to vec_pq.txt
/// \fn void write_vec_pq(FILE *file_pq, Metal *p_metal, int step)
/// \param file_pq name of dipole/charge output file
/// \param p_metal pointer to Metal data structure
/// \param step number of MD step
void write_vec_pq(FILE *file_pq, Metal *p_metal, int step)
{
    int n_metal = p_metal->num;
    fprintf(file_pq, "step %d", step);

    // print Lagrange multipliers
    int iNP;
    fprintf(file_pq, "  lambda");
    for (iNP = 0; iNP < p_metal->n_NPs; ++ iNP)
    {
        fprintf(file_pq, " %.10f", p_metal->vec_pq[n_metal*4 + iNP]);
    }
    fprintf(file_pq, "\n");

    // print number of metal atoms, dipole vectors and charges
    fprintf(file_pq, "%5d\n", n_metal);

    int i;
    for (i = 0; i < n_metal; ++ i)
    {
        fprintf(file_pq, "%15.8f%15.8f%15.8f%15.8f\n", 
                p_metal->vec_pq[i*3], p_metal->vec_pq[i*3+1], p_metal->vec_pq[i*3+2], 
                p_metal->vec_pq[n_metal*3+i]);
    }
}


/// \brief write to binary trajectory file
/// \fn void write_binary(FILE *file_dat, System *p_system, int nAtoms, int step)
/// \param file_dat name of binary output trajectory file
/// \param p_system pointer to System data structure
/// \param nAtoms number of atoms
/// \param step number of MD step
void write_binary(FILE *file_dat, System *p_system, int nAtoms, int step)
{
    fwrite(&nAtoms, sizeof(int), 1, file_dat);
    fwrite(&step,   sizeof(int), 1, file_dat);

    fwrite(p_system->box, sizeof(double),      3, file_dat);
    fwrite(p_system->rx,  sizeof(double), nAtoms, file_dat);
    fwrite(p_system->ry,  sizeof(double), nAtoms, file_dat);
    fwrite(p_system->rz,  sizeof(double), nAtoms, file_dat);
}


/// \brief read force field parameters from input_param
/// \fn void read_param_1(char *input_param, Topol *p_topol, Metal* p_metal)
/// \param input_param name of input parameter file
/// \param p_topol pointer to Topol data structure
/// \param p_metal pointer to Metal data structure
void read_param_1(char *input_param, Topol *p_topol, Metal* p_metal)
{
    FILE* file_par;
    char line[MAX_STR_LEN], subline[MAX_STR_LEN], my_string[MAX_STR_LEN];
    int  mol, atom;

    file_par = fopen(input_param, "r") ;
    if (NULL == file_par)
    {
        printf("Error: Cannot open param file %s!\n", input_param) ;
        exit(1);
    }

    // 1-4 scaling factors
    read_next_line(file_par, line, subline);
    if (0 == strcmp(subline, "SCALE_LJ_QQ"))
    {
        sscanf(line, "%s%lf%lf", subline, &(p_topol->scaleLJ), &(p_topol->scaleQQ));
    }
    else
    {
        printf("Error in %s: expected SCALE_LJ_QQ in line\n%s", input_param, line);
        exit(1);
    }

    // number of molecule types
    read_next_line(file_par, line, subline);
    if (0 == strcmp(subline, "MOLECULES"))
    {
        sscanf(line, "%s%d", subline, &(p_topol->mol_types));
    }
    else
    {
        printf("Error in %s: expected MOLECULES in line\n%s", input_param, line);
        exit(1);
    }


    // molecule name
    char molName[p_topol->mol_types][MAX_STR_LEN];


    // molNr section: number of this type of molecule in the system
    for (mol = 0; mol < p_topol->mol_types; ++ mol) 
    {
        read_next_line(file_par, line, subline);
        sscanf(line, "%s", molName[mol]);
    }

    // atom_in_mol section: number of atoms in this type of molecule
    for (mol = 0; mol < p_topol->mol_types; ++ mol) 
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
            printf("Error in %s: expected %s in line\n%s", input_param, my_string, line);
            exit(1);
        }
    }


    // nonbonded parameters
    int iType, jType;

    read_next_line(file_par, line, subline);
    sscanf(line, "%s%d", subline, &(p_topol->n_types));

    if (0 == strcmp(subline, "NONBONDED"))
    {
        for (iType = 0; iType < p_topol->n_types; ++ iType)
        {
            for (jType = iType; jType < p_topol->n_types; ++ jType)
            {
                read_next_line(file_par, line, subline);
            }
        }
    }
    else
    {
        printf("Error in %s: expected NONBONDED in line\n%s", input_param, line);
        exit(1);
    }

    // bonded potentials: bond, pair, angle, dihedral
    for (mol = 0; mol < p_topol->mol_types; ++ mol) 
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
            printf("Error in %s: expected %s in line\n%s", input_param, my_string, line);
            exit(1);
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
            printf("Error in %s: expected %s in line\n%s", input_param, my_string, line);
            exit(1);
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
            printf("Error in %s: expected %s in line\n%s", input_param, my_string, line);
            exit(1);
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
            printf("Error in %s: expected %s in line\n%s", input_param, my_string, line);
            exit(1);
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
            printf("Error in %s: expected %s in line\n%s", input_param, my_string, line);
            exit(1);
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
            printf("Error in %s: expected %s in line\n%s", input_param, my_string, line);
            exit(1);
        }
    }

    // metal indices (min, max)
    read_next_line(file_par, line, subline);
    sscanf(line, "%s%d%d", subline, &(p_metal->min), &(p_metal->max));

    if (0 == strcmp(subline, "MM_METAL"))
    {
        if (p_metal->min >=0 && p_metal->max >= p_metal->min)
        {
            p_metal->num = p_metal->max - p_metal->min + 1;
        }
        else
        {
            p_metal->num = 0;
        }
    }
    else
    {
        printf("Error in %s: expected MM_METAL in line\n%s", input_param, line);
        exit(1);
    }

    // Quantum Sutton-Chen parameters
    read_next_line(file_par, line, subline);
    sscanf(line, "%s%d%d%lf%lf%lf", subline, 
           &(p_metal->qsc_n), &(p_metal->qsc_m),
           &(p_metal->qsc_eps), &(p_metal->qsc_c), &(p_metal->qsc_a));
    if (0 != strcmp(subline, "QSC_METAL"))
    {
        printf("Error in %s: expected QSC_METAL in line\n%s", input_param, line);
        exit(1);
    }

    // atomic polarizability and capacitance, in atomic unit
    read_next_line(file_par, line, subline);
    sscanf(line, "%s%lf%lf", subline, &(p_metal->cpff_polar), &(p_metal->cpff_capac));
    if (0 != strcmp(subline, "CPFF_METAL"))
    {
        printf("Error in %s: expected CPFF_METAL in line\n%s", input_param, line);
        exit(1);
    }

    // multi-NPs
    read_next_line(file_par, line, subline);
    sscanf(line, "%s%d", subline, &(p_metal->n_NPs));
    if (0 != strcmp(subline, "NANOPARTICLES"))
    {
        printf("Error in %s: expected NANOPARTICLES in line\n%s", input_param, line);
        exit(1);
    }

    // close param file
    fclose( file_par );
}

/// \brief read FF parameters (2) from input_param
/// \fn void read_param_2(char *input_param, Topol *p_topol, Metal *p_metal)
/// \param input_param name of input parameter file
/// \param p_topol pointer to Topol data structure
/// \param p_metal pointer to Metal data structure
void read_param_2(char *input_param, Topol *p_topol, Metal *p_metal)
{
    FILE *file_par;
    char line[MAX_STR_LEN], tmp[MAX_STR_LEN];

    file_par = fopen(input_param,"r");
    if (NULL == file_par) 
    {
        printf("Error: Cannot open param file %s!\n", input_param);
        exit(1);
    }

    // 1-4 scaling factors
    read_next_line(file_par, line, tmp);

    // number of molecule types
    read_next_line(file_par, line, tmp);

    // p_topol->mol_num = number of this type of molecule in the system
    int mol;
    for (mol = 0; mol < p_topol->mol_types; ++ mol) 
    {
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d", tmp, &p_topol->mol_num[mol]);
    }

    // p_topol->atom_num = number of atoms in this type of molecule
    // p_topol->atom_param = atomic parameters (q, m, sig, eps, type)
    for (mol = 0; mol < p_topol->mol_types; ++ mol) 
    {
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d", tmp, &p_topol->atom_num[mol]);

        p_topol->atom_param[mol] = my_malloc(p_topol->atom_num[mol] * sizeof(AtomParam));
        int atom;
        for (atom = 0; atom < p_topol->atom_num[mol]; ++ atom) 
        {
            read_next_line(file_par, line, tmp);
            sscanf(line, "%lf%lf%d",
                   &(p_topol->atom_param[mol][atom].charge),
                   &(p_topol->atom_param[mol][atom].mass),
                   &(p_topol->atom_param[mol][atom].atomtype));
        }
    }

    // nonbonded parameters
    read_next_line(file_par, line, tmp);
    //sscanf(line, "%s%d", tmp, &p_topol->n_types);

    int iType, jType, funct;
    for (iType = 0; iType < p_topol->n_types; ++ iType)
    {
        for (jType = iType; jType < p_topol->n_types; ++ jType)
        {
            read_next_line(file_par, line, tmp);

            sscanf(line, "%d", &funct);

            // LJ
            if (1 == funct)
            {
                sscanf(line, "%d%le%le", 
                       &(p_topol->nonbonded_param[iType][jType].funct),
                       &(p_topol->nonbonded_param[iType][jType].C6),
                       &(p_topol->nonbonded_param[iType][jType].C12));
            }

            // Buckingham
            else if (2 == funct)
            {
                sscanf(line, "%d%le%le%le", 
                       &(p_topol->nonbonded_param[iType][jType].funct),
                       &(p_topol->nonbonded_param[iType][jType].A),
                       &(p_topol->nonbonded_param[iType][jType].B),
                       &(p_topol->nonbonded_param[iType][jType].C6));
            }

            // Morse
            else if (3 == funct)
            {
                sscanf(line, "%d%le%le%le", 
                       &(p_topol->nonbonded_param[iType][jType].funct),
                       &(p_topol->nonbonded_param[iType][jType].A),
                       &(p_topol->nonbonded_param[iType][jType].B),
                       &(p_topol->nonbonded_param[iType][jType].C6));
            }

            // erf_vdw
            else if (4 == funct)
            {
                sscanf(line, "%d%le%le%le", 
                       &(p_topol->nonbonded_param[iType][jType].funct),
                       &(p_topol->nonbonded_param[iType][jType].C6),
                       &(p_topol->nonbonded_param[iType][jType].C12),
                       &(p_topol->nonbonded_param[iType][jType].A));
            }

            // error
            else
            {
                printf("Error: unsupported vdW type: %d in line:\n%s\n", 
                       funct, line );
                exit(1);
            }
        }
    }

    for (iType = 0; iType < (p_topol->n_types - 1); ++ iType)
    {
        for (jType = iType + 1; jType < p_topol->n_types; ++ jType)
        {
            p_topol->nonbonded_param[jType][iType].funct = 
                p_topol->nonbonded_param[iType][jType].funct;

            p_topol->nonbonded_param[jType][iType].A = 
                p_topol->nonbonded_param[iType][jType].A;

            p_topol->nonbonded_param[jType][iType].B = 
                p_topol->nonbonded_param[iType][jType].B;

            p_topol->nonbonded_param[jType][iType].C6 = 
                p_topol->nonbonded_param[iType][jType].C6;

            p_topol->nonbonded_param[jType][iType].C12 = 
                p_topol->nonbonded_param[iType][jType].C12;
        }
    }

    // bonded potentials: bond, pair, angle, dihedral
    for (mol = 0; mol < p_topol->mol_types; ++ mol)
    {
        int atom_i, atom_j, atom_k, atom_l, atom_s;
        int iBond, iPair, iAngle, iDihedral, iVSite, iCstr;
        int funct;

        // exclusions
        p_topol->exclude[mol] = my_malloc(p_topol->atom_num[mol] * sizeof(int *));
        for (atom_i = 0; atom_i < p_topol->atom_num[mol]; ++ atom_i) 
        {
            p_topol->exclude[mol][atom_i] = my_malloc(p_topol->atom_num[mol] * sizeof(int));
            for (atom_j = 0; atom_j < p_topol->atom_num[mol]; ++ atom_j) 
            {
                p_topol->exclude[mol][atom_i][atom_j] = 0;
            }
        }

        // bonds
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d", tmp, &p_topol->n_bonds[mol]);

        p_topol->bond_param[mol] = my_malloc(p_topol->n_bonds[mol] * sizeof(BondParam));

        for (iBond = 0; iBond < p_topol->n_bonds[mol]; ++ iBond)
        {
            read_next_line(file_par, line, tmp);
            sscanf(line, "%d%d%d", &atom_i, &atom_j, &funct);

            // Harmonic bond potential
            if (1 == funct)
            {
                sscanf(line, "%d%d%d%le%le",
                       &(p_topol->bond_param[mol][iBond].atom_i), 
                       &(p_topol->bond_param[mol][iBond].atom_j),
                       &(p_topol->bond_param[mol][iBond].funct),
                       &(p_topol->bond_param[mol][iBond].b0), 
                       &(p_topol->bond_param[mol][iBond].kb));
            }

            // Morse bond potential
            else if (3 == funct)
            {
                sscanf(line, "%d%d%d%le%le%le",
                       &(p_topol->bond_param[mol][iBond].atom_i), 
                       &(p_topol->bond_param[mol][iBond].atom_j),
                       &(p_topol->bond_param[mol][iBond].funct),
                       &(p_topol->bond_param[mol][iBond].b0), 
                       &(p_topol->bond_param[mol][iBond].D), 
                       &(p_topol->bond_param[mol][iBond].beta));
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
                exit(1);
            }

            // update exclusion
            p_topol->exclude[mol][atom_i][atom_j] = 1;
            p_topol->exclude[mol][atom_j][atom_i] = 1;
        }

        // pairs
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d", tmp, &p_topol->n_pairs[mol]);

        p_topol->pair_param[mol] = my_malloc(p_topol->n_pairs[mol] * sizeof(PairParam));

        for (iPair = 0; iPair < p_topol->n_pairs[mol]; ++ iPair)
        {
            read_next_line(file_par, line, tmp);
            sscanf(line, "%d%d%d", &atom_i, &atom_j, &funct);

            // normal 1-4 pair
            if (1 == funct)
            {
                sscanf(line, "%d%d%d",
                       &(p_topol->pair_param[mol][iPair].atom_i), 
                       &(p_topol->pair_param[mol][iPair].atom_j),
                       &(p_topol->pair_param[mol][iPair].funct));
            }

            // user-defined 1-4 pair
            else if (2 == funct)
            {
                sscanf(line, "%d%d%d%lf%lf%lf%le%le",
                       &(p_topol->pair_param[mol][iPair].atom_i), 
                       &(p_topol->pair_param[mol][iPair].atom_j),
                       &(p_topol->pair_param[mol][iPair].funct),
                       &(p_topol->pair_param[mol][iPair].scaleQQ),
                       &(p_topol->pair_param[mol][iPair].q_i), 
                       &(p_topol->pair_param[mol][iPair].q_j),
                       &(p_topol->pair_param[mol][iPair].V), 
                       &(p_topol->pair_param[mol][iPair].W));
            }

            else
            {
                printf("Error: unsupported pair function type: %d in line:\n%s\n", 
                        funct, line);
                exit(1);
            }

            // update exclusion
            p_topol->exclude[mol][atom_i][atom_j] = 1;
            p_topol->exclude[mol][atom_j][atom_i] = 1;
        }

        // angles
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d", tmp, &p_topol->n_angles[mol]);

        p_topol->angle_param[mol] = my_malloc(p_topol->n_angles[mol] * sizeof(AngleParam));

        for (iAngle = 0; iAngle < p_topol->n_angles[mol]; ++ iAngle)
        {
            read_next_line(file_par, line, tmp);
            sscanf(line, "%d%d%d%d", &atom_i, &atom_j, &atom_k, &funct);

            // harmonic angle potential
            if (1 == funct)
            {
                sscanf(line, "%d%d%d%d%le%le",
                       &(p_topol->angle_param[mol][iAngle].atom_i), 
                       &(p_topol->angle_param[mol][iAngle].atom_j), 
                       &(p_topol->angle_param[mol][iAngle].atom_k),
                       &(p_topol->angle_param[mol][iAngle].funct),
                       &(p_topol->angle_param[mol][iAngle].a0), 
                       &(p_topol->angle_param[mol][iAngle].ka));
            }

            // Urey-Bradley type
            else if (5 == funct)
            {
                sscanf(line, "%d%d%d%d%le%le%le%le",
                       &(p_topol->angle_param[mol][iAngle].atom_i), 
                       &(p_topol->angle_param[mol][iAngle].atom_j), 
                       &(p_topol->angle_param[mol][iAngle].atom_k),
                       &(p_topol->angle_param[mol][iAngle].funct),
                       &(p_topol->angle_param[mol][iAngle].a0), 
                       &(p_topol->angle_param[mol][iAngle].ka),
                       &(p_topol->angle_param[mol][iAngle].b0_ub), 
                       &(p_topol->angle_param[mol][iAngle].kb_ub));
            }

            else
            {
                printf("Error: unsupported angle function type: %d in line:\n%s\n", 
                        funct, line);
                exit(1);
            }

            // update exclusion
            p_topol->exclude[mol][atom_i][atom_j] = 1;
            p_topol->exclude[mol][atom_j][atom_i] = 1;
            p_topol->exclude[mol][atom_i][atom_k] = 1;
            p_topol->exclude[mol][atom_k][atom_i] = 1;
            p_topol->exclude[mol][atom_j][atom_k] = 1;
            p_topol->exclude[mol][atom_k][atom_j] = 1;
        }

        // dihedrals
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d", tmp, &p_topol->n_dihedrals[mol]);

        p_topol->dihedral_param[mol] = my_malloc(p_topol->n_dihedrals[mol] * sizeof(DihedralParam));

        for (iDihedral = 0; iDihedral < p_topol->n_dihedrals[mol]; ++ iDihedral)
        {
            read_next_line(file_par, line, tmp);
            sscanf(line, "%d%d%d%d%d", &atom_i, &atom_j, &atom_k, &atom_l, &funct);
            
            // RB type dihedral potential
            if (3 == funct)
            {
                sscanf(line, "%d%d%d%d%d%lf%lf%lf%lf%lf%lf",
                       &(p_topol->dihedral_param[mol][iDihedral].atom_i), 
                       &(p_topol->dihedral_param[mol][iDihedral].atom_j),
                       &(p_topol->dihedral_param[mol][iDihedral].atom_k), 
                       &(p_topol->dihedral_param[mol][iDihedral].atom_l),
                       &(p_topol->dihedral_param[mol][iDihedral].funct),
                       &(p_topol->dihedral_param[mol][iDihedral].c0), 
                       &(p_topol->dihedral_param[mol][iDihedral].c1),
                       &(p_topol->dihedral_param[mol][iDihedral].c2), 
                       &(p_topol->dihedral_param[mol][iDihedral].c3),
                       &(p_topol->dihedral_param[mol][iDihedral].c4), 
                       &(p_topol->dihedral_param[mol][iDihedral].c5));
            }

            // periodic dihedral potential
            else if (1 == funct || 4 == funct || 9 == funct)
            {
                sscanf(line, "%d%d%d%d%d%lf%lf%d",
                       &(p_topol->dihedral_param[mol][iDihedral].atom_i), 
                       &(p_topol->dihedral_param[mol][iDihedral].atom_j),
                       &(p_topol->dihedral_param[mol][iDihedral].atom_k), 
                       &(p_topol->dihedral_param[mol][iDihedral].atom_l),
                       &(p_topol->dihedral_param[mol][iDihedral].funct),
                       &(p_topol->dihedral_param[mol][iDihedral].phi0), 
                       &(p_topol->dihedral_param[mol][iDihedral].kphi),
                       &(p_topol->dihedral_param[mol][iDihedral].n));
            }

            else
            {
                printf("Error: unsupported dihedral function type: %d in line:\n%s\n", 
                       funct, line);
                exit(1);
            }

            // update exclusion
            p_topol->exclude[mol][atom_i][atom_j] = 1;
            p_topol->exclude[mol][atom_j][atom_i] = 1;
            p_topol->exclude[mol][atom_i][atom_k] = 1;
            p_topol->exclude[mol][atom_k][atom_i] = 1;
            p_topol->exclude[mol][atom_i][atom_l] = 1;
            p_topol->exclude[mol][atom_l][atom_i] = 1;
            p_topol->exclude[mol][atom_j][atom_k] = 1;
            p_topol->exclude[mol][atom_k][atom_j] = 1;
            p_topol->exclude[mol][atom_j][atom_l] = 1;
            p_topol->exclude[mol][atom_l][atom_j] = 1;
            p_topol->exclude[mol][atom_k][atom_l] = 1;
            p_topol->exclude[mol][atom_l][atom_k] = 1;
        }


        // virtual_site, TIP5P type, funct=4
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d", tmp, &p_topol->n_vsites[mol]);

        p_topol->vsite_4[mol] = my_malloc(p_topol->n_vsites[mol] * sizeof(VSite_4));
        p_topol->vsite_funct[mol] = my_malloc(p_topol->n_vsites[mol] * sizeof(int));

        for (iVSite = 0; iVSite < p_topol->n_vsites[mol]; ++ iVSite)
        {
            read_next_line(file_par, line, tmp);
            sscanf(line, "%d%d%d%d%d", &atom_i, &atom_j, &atom_k, &atom_s, &funct);

            if (4 == funct)
            {
                sscanf(line, "%d%d%d%d%d%lf%lf%lf",
                       &(p_topol->vsite_4[mol][iVSite].atom_i), 
                       &(p_topol->vsite_4[mol][iVSite].atom_j),
                       &(p_topol->vsite_4[mol][iVSite].atom_k), 
                       &(p_topol->vsite_4[mol][iVSite].atom_s),
                       &(p_topol->vsite_funct[mol][iVSite]),    
                       &(p_topol->vsite_4[mol][iVSite].a), 
                       &(p_topol->vsite_4[mol][iVSite].b),      
                       &(p_topol->vsite_4[mol][iVSite].c));

                p_topol->exclude[mol][atom_s][atom_i] = 1;
                p_topol->exclude[mol][atom_s][atom_j] = 1;
                p_topol->exclude[mol][atom_s][atom_k] = 1;
                p_topol->exclude[mol][atom_i][atom_s] = 1;
                p_topol->exclude[mol][atom_j][atom_s] = 1;
                p_topol->exclude[mol][atom_k][atom_s] = 1;
            }
            else
            {
                printf("Error: unsupported virtual_sites function type: %d in line:\n%s\n", 
                       funct, line);
                exit(1);
            }
        }


        // p_topol->constraints
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d", tmp, &p_topol->n_constraints[mol]);

        p_topol->constraint[mol] = my_malloc(p_topol->n_constraints[mol] * sizeof(Constraint));

        for (iCstr = 0; iCstr < p_topol->n_constraints[mol]; ++ iCstr)
        {
            read_next_line(file_par, line, tmp);
            sscanf(line, "%d%d", &atom_i, &atom_j);

            sscanf(line, "%d%d%lf",
                   &(p_topol->constraint[mol][iCstr].atom_i), 
                   &(p_topol->constraint[mol][iCstr].atom_j),
                   &(p_topol->constraint[mol][iCstr].r));

            p_topol->exclude[mol][atom_i][atom_j] = 1;
            p_topol->exclude[mol][atom_j][atom_i] = 1;
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
    for (i_NP = 0; i_NP < p_metal->n_NPs; ++ i_NP)
    {
        read_next_line(file_par, line, tmp);
        sscanf(line, "%s%d%d%lf", 
               tmp, &p_metal->start_NP[i_NP], &p_metal->end_NP[i_NP], &p_metal->cpff_chg[i_NP]);
    }

    // close param file
    fclose(file_par);

    // calculate total number of atoms and molecules
    p_topol->n_atoms = 0;
    p_topol->n_mols  = 0;
    for (mol = 0; mol < p_topol->mol_types; ++ mol) 
    {
        int num_atom = p_topol->atom_num[mol];
        int num_mol  = p_topol->mol_num[mol];
        p_topol->n_atoms += num_mol * num_atom;
        p_topol->n_mols  += num_mol;
    }
}

/// \brief assign molecular and atomic indices
/// \fn void assign_indices(Topol *p_topol, Metal *p_metal, Atom_Info* atom_info, Mol_Info* mol_info)
/// \param p_topol pointer to Topol data structure
/// \param p_metal pointer to Metal data structure
/// \param atom_info array of atomic information
/// \param mol_info array of molecular information
void assign_indices(Topol *p_topol, Metal *p_metal,
                    Atom_Info* atom_info, Mol_Info* mol_info)
{
    int i, im, mol, km, atom;

    // assign atom index and mol_type index for each atom in a molecule type
    // assign starting index and ending index for each molecule in the system 
    // i = atom index, im = mol index
    i  = 0;
    im = 0;
    for (mol = 0; mol < p_topol->mol_types; ++ mol) 
    {
        for (km = 0; km< p_topol->mol_num[mol]; ++ km) 
        {
            mol_info[im].mini = i;

            for (atom = 0; atom < p_topol->atom_num[mol]; ++ atom) 
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

    if (p_topol->n_atoms != i || p_topol->n_mols != im) 
    {
        printf("Error: Incorrect number of atoms or molecules in the system! (1st)\n");
        printf("n_atoms = %d, i = %d\n", p_topol->n_atoms, i);
        printf("n_mols = %d, im = %d\n", p_topol->n_mols, im);
        exit(1);
    }


    // assign atom and residue, mass and charge
    for (i = 0; i < p_topol->n_atoms; ++ i)
    {
        int this_iMol  = atom_info[i].iMol;
        int this_iAtom = atom_info[i].iAtom;

        atom_info[i].mass     = p_topol->atom_param[this_iMol][this_iAtom].mass;
        atom_info[i].charge   = p_topol->atom_param[this_iMol][this_iAtom].charge;
        atom_info[i].atomtype = p_topol->atom_param[this_iMol][this_iAtom].atomtype;

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
    for (i = p_metal->min; i <= p_metal->max; ++ i)
    {
        atom_info[i].is_metal = 1;
    }
}
