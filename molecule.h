#ifndef MOLECULE_H
#define MOLECULE_H
#include <string>
#include "atom.h"

using namespace std;
const int max_mol_types = 1000;

class molecule: public atom{
      protected:
              int no_atoms;
              int mol_number;
              double mol_weight;
              string mol_name;
      public:
             molecule(){
                        no_atoms = 0;
                        mol_number = 0;
                        mol_weight = 0.0;
                        mol_name = "None";
                        }
             void set_molecule_number( int n1){ mol_number = n1;}
             
             void increase_atom_nos(){
                  no_atoms = no_atoms + 1;
                  }
                  
             int return_no_atoms_in_molecule(){ return no_atoms;}
             
             void set_molecule_name(string s1){mol_name = s1;}
             void set_molecule_mass(double m1){mol_weight = m1;}
             
             string return_molecule_name(){ return mol_name;}
             double return_molecule_weight(){return mol_weight;}
             
             int return_molecule_number(){return mol_number;}
      
      };

#endif
