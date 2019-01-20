#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "atom.h"
#include "molecule.h"
#include "functions.h"

using namespace std;

void identify_atoms_in_molecule(int natoms1, atom **at_list1, molecule **mol_list1, int nmols){
     int i,j,k;
     int frame;
     int num1,num2,num3;
     string name;
     string s1,s2;
     
     cout << "Enter the name of the molecule (in all caps)"<< endl;
     cin >> name;
     frame=0;
     for(i=0;i<natoms1;i++){
                            num1 = at_list1[frame][i].return_mol_no() - 1 ;
                            if(num1 > nmols-1){
                                    cout << "something is wrong " << endl;
                                    //exit(1);
                                    //system("Pause");
                                    }
                            s1.clear();
                            s1 = mol_list1[frame][num1].return_molecule_name();
                            if(s1 == name){
                                  cout << i+1 << endl; 
                                  }
                            }
     }
