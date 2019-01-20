#include "atom.h"
#include "molecule.h"
#include "emptyBin.h"
#include "voids.h"
#include "settings.h"
#include<iostream>
#include<sstream>
#include<fstream>

using namespace std;

//const int max_neigh =8 ;

#ifndef FIND_NO_LINES_H
#define FIND_NO_LINES_H

void find_no_lines(int &nlines1, int &natoms1, int &nframes1);

#endif

#ifndef FIND_NO_FFLINES_H
#define FIND_NO_FFLINES_H

int find_no_fflines();

#endif

#ifndef FIND_NO_CONNT_LINES_H
#define FIND_NO_CONNT_LINES_H

int find_no_connt_lines(istream &myfile);

#endif

#ifndef READ_INPUT_XMOLOUT_H
#define READ_INPUT_XMOLOUT_H

int read_input_xmolout(istream &myfile,int nframes1, int natoms1,int periodic, atom **at_list1, string *cell_parameter_lines1,
                        double **lx1, double **ly1, double **lz1, double *alpha1, double *beta1, double *gama1);                        

#endif

#ifndef FIND_CELL_PARAMETERS_H
#define FIND_CELL_PARAMETERS_H

void find_cell_parameters(string s1, int frame_no, double **lx1, double **ly1, double **lz1, double *alpha1, double *beta1, double *gama1);

#endif

#ifndef GET_SETTINGS_H
#define GET_SETTINGS_H

void get_settings(settings& settings_list1);

#endif

#ifndef FIND_DISTANCE_H
#define FIND_DISTANCE_H

double find_distance(int nframes1,int target_atoms1,int object_atoms1,
                   int *t_atoms1, int *o_atoms1,atom **at_list1,
                   double **lx1,double **ly1,double **lz1,double *alpha1,
                   double *beta1,double *gama1, int frame1);

#endif

#ifndef STRING_TO_INTEGER_H
#define STRING_TO_INTEGER_H

int string_to_integer(char *s1);

#endif

#ifndef STRING_TO_DOUBLE_H
#define STRING_TO_DOUBLE_H

double string_to_double(char *s1);

#endif

#ifndef FLOAT_TO_STRING_H
#define FLOAT_TO_STRING_H

string float_to_string(int a);

#endif

#ifndef IDENTIFY_ATOM_TYPES_H
#define IDENTIFY_ATOM_TYPES_H

void identify_atom_types(string *line1, int start_line_per_section[], int num_atom_types, string *atom_types1);

#endif

#ifndef IDENTIFY_BOND_PARAMETERS_H
#define IDENTIFY_BOND_PARAMETERS_H

void identify_bond_parameters(double **bond_parameters1, string *lines1, int start_line_per_section1[], int no_lines);

#endif


#ifndef FIND_NO_TYPES_INPUT_ATOMS_H
#define FIND_NO_TYPES_INPUT_ATOMS_H

int find_no_types_input_atoms(atom **at_list, int natoms1, string input_atom_types1[]);

#endif

#ifndef WRITE_LAMMPS_INPUT_H
#define WRITE_LAMMPS_INPUT_H

void write_lammps_input(int natoms, int no_atom_types, int lines_per_section, atom **at_list,double **lx, double **ly, double **lz,
                    double *alpha, double *beta, double *gama, string *atom_types, string *ip_atom_types, double *mass);

#endif

#ifndef LAMMPS_OUTPUT_H
#define LAMMPS_OUTPUT_H

void lammps_output(string *ip_atom_types);

#endif

#ifndef FIND_VALENCE_PARAMS_H
#define FIND_VALENCE_PARAMS_H

void find_valence_params(string *line1, int start_line_per_section[], int num_atom_types, double *mass,
                         double *valency, double *cov_r1, double *cov_r2, double *cov_r3);

#endif

#ifndef MAKE_CONNECTION_TABLE_H
#define MAKE_CONNECTION_TABLE_H

void make_connection_table(ofstream &file1, int nframes, int natoms, atom **at_list, string *atom_types, int size, int size1, double **lx, double **ly, double **lz, 
                           double *alpha, double *beta, double *gama, double **bond_parameters,
                           double *valency, double *cov_r1, double *cov_r2,double *cov_r3, int ****grid_list1, int ***grid_atom_count1, int nx1, int ny1, int nz1);
                          
#endif

#ifndef IDENTIFY_MOLECULES_H
#define IDENTIFY_MOLECULES_H

void identify_molecules(int nframes, int act_frame,int natoms, int no_atom_types, atom **at_list, string *atom_types,
                            double **lx, double **ly, double **lz, double *cov_r2,double *cov_r3, 
                            int no_bond_parameters, double **bond_params, int size1, double *mass);

#endif

#ifndef IDENTIFY_MOLECULE_NAME_H
#define IDENTIFY_MOLECULE_NAME_H

void identify_molecule_name(int nframes, int natoms,int *number_molecules, int no_atom_types2, atom **atom_list1,
                             molecule **molecule_list, string *atom_type1, double *m1);

#endif

#ifndef IDENTIFY_CH3_H
#define IDENTIFY_CH3_H

void identify_ch3(int nframes, int natoms,int number_molecules, atom **at_list, molecule **molecule_list, 
                  int temp_neighbours[], int neighbours_index[][20]);

#endif

#ifndef IDENTIFY_C3_MOLECULES_H
#define IDENTIFY_C3_MOLECULES_H

void identify_c3_molecules(int nframes, int act_frame,int natoms,int number_molecules, atom **at_list, molecule **molecule_list);

#endif

#ifndef IDENTIFY_C3_H
#define IDENTIFY_C3_H

void identify_c3(int nframes, int act_frame,int natoms,atom **at_list);

#endif

#ifndef IDENTIFY_C2_H
#define IDENTIFY_C2_H

void identify_c2(int act_frame,int natoms,atom **at_list, int *c2_flag);

#endif

#ifndef IDENTIFY_CARBON_3_NEIGHBORS_H
#define IDENTIFY_CARBON_3_NEIGHBORS_H

void identify_carbon_3_neighbors(int nframes, int natoms,int number_molecules, atom **at_list, molecule **molecule_list, 
                  int temp_neighbours[], int neighbours_index[][20]);

#endif

#ifndef IDENTIFY_DOUBLE_BONDS_H
#define IDENTIFY_DOUBLE_BONDS_H

void identify_double_bonds(int nframes, int natoms,int number_molecules, atom **at_list, molecule **molecule_list, 
                  int temp_neighbours[], int neighbours_index[][20], double **lx, double **ly, double **lz, 
                  double *cov_r2, double *cov_r3, int no_bond_params1, double **double_params1, string *atom_types2, int size1);

#endif

#ifndef IDENTIFY_HEADER_LINES_H
#define IDENTIFY_HEADER_LINES_H

int identify_header_lines(istream &myfile);

#endif

#ifndef INPUT_CONNECTION_TABLE_H
#define INPUT_CONNECTION_TABLE_H

void input_connection_table(istream &myfile, int frmae_no, int natoms, int hlines, atom **at_list);

#endif

#ifndef READ_TABLE_H
#define READ_TABLE_H

void read_table(istream& myfile, int frame_no, int natoms1, int hlines, double cut_off,atom **at_list1);

#endif

#ifndef READ_REAX_TABLE_H
#define READ_REAX_TABLE_H

void read_reax_table(istream& myfile, int frame_no, int natoms1, double cut_off, atom **at_list1);

#endif


#ifndef IDENTIFY_ATOMS_IN_MOLECULE_H
#define IDENTIFY_ATOMS_IN_MOLECULE_H

void identify_atoms_in_molecule(int natoms, atom **at_list, molecule **mol_list1, int nmols);

#endif

#ifndef WRITE_REACTION_INPUT_FILE_H
#define WRITE_REACTION_INPUT_FILE_H

void write_reaction_input_file(int nframes1, int natoms1,int number_molecules, atom **at_list1, molecule **molecule_list1);

#endif

#ifndef BIN_ATOMS_H
#define BIN_ATOMS_H

void bin_atoms(int nframes1, int natoms1, double *grid_dr1, int grid_max_atoms1, atom **at_list1, double **lx1, double **ly1, double **lz1,int ****grid_list, int ***grid_atom_count, int nx1, int ny1, int nz1, settings settings_list1);

#endif

#ifndef RING_ANALYSIS_H
#define RING_ANALYSIS_H
void ring_analysis(int act_frame1, int natoms1, atom **at_list1, int nx1, int ny1, int nz1, double **lx1, double **ly1, double **lz1,settings settings_list1);
#endif

#ifndef IDENTIFY_VOIDS_H
#define IDENTIFY_VOIDS_H
void identify_voids(int act_frame1, int natoms1,double *grid_dr1,int no_empty_bins, emptyBin *bin_list1, int ***grid_atom_count, int ****grid_list1, atom **at_list1, double **lx1, double **ly1, double **lz1, int nx1, int ny1, int nz1);
#endif

#ifndef EQUATION_OF_PLANE_H
#define EQUATION_OF_PLANE_H
void equation_of_plane(double (&points)[3][3], double (&params)[4]);
#endif
