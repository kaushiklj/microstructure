#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "atom.h"
#include "molecule.h"
#include "functions.h"

using namespace std;

void identify_molecule_name(int nframes, int natoms,int *number_molecules, int no_atom_types2, atom **atom_list1,
                             molecule **molecule_list, string *atom_type1, double *mass1){
    int i1,i2,i3,i4; 
    int temp1,temp2;
    int *count_every_atom_type;
    double m1;
    string line1,line2,line3;
    ofstream myfile1,myfile2;
    
    count_every_atom_type = new int[no_atom_types2];
    
    /*myfile2.open("atom_molecule.txt");
    for (i1=0;i1<natoms;i1++){
       myfile2 << i1+1 << " " << atom_list1[0][i1].return_mol_no()<< endl;
    }
    myfile2.close();*/
    for (i1=0;i1<nframes;i1++){//start iterating over every frame
        temp1 = number_molecules[i1];
        
        for (i2=0;i2<temp1;i2++){//start iterating over number of molecules in that frame
        for (i3=0;i3<no_atom_types2;i3++){
                                         count_every_atom_type[i3] = 0;
                                         //cout << " printing from count loop " << count_every_atom_type[i3]<< endl;
                                         }
            for (i3=0;i3<natoms;i3++){//loop for checking which atoms belong to i2th molecule
                temp2 = atom_list1[i1][i3].return_mol_no();
                if(temp2 == (i2+1)){
                         //cout << " molecule number matched" << endl;
                         //cout << count_every_atom_type[0]<< endl;
                         line1 = atom_list1[i1][i3].return_atomname();
                         for (i4=0;i4<no_atom_types2;i4++){
                             if (line1 == atom_type1[i4]){
                                       count_every_atom_type[i4] = count_every_atom_type[i4] + 1; //keep track of number of atoms of every type 
                                       }
                             }   
                   }
             /*if(temp2 == 903){
                   cout << i3 << endl;
              }*/   
            }
            /*if (i2 == 902){
              for (i4=0;i4<no_atom_types2;i4++){
                 cout << count_every_atom_type[i4]<< endl;
              }
             }*/
            //cout << "in molecule  " << i2 << " number of  " << atom_type_name[1] << " atoms are   " << count_every_atom_type[1] << endl;
            line1.clear();
            m1 =0;
            for(i4 =0;i4<no_atom_types2;i4++){
                   if(count_every_atom_type[i4] != 0){//throws away those atom types that are not present in the molecule
                                                if (line1.empty()){//check if the line is empty
                                                                   line1 =atom_type1[i4];
                                                                   line2 = float_to_string(count_every_atom_type[i4]);
                                                                   line1.append(line2);
                                                                   }else{
                                                                         line2 = atom_type1[i4];
                                                                         line3 = float_to_string(count_every_atom_type[i4]);
                                                                         line1.append(line2);
                                                                         line1.append(line3);
                                                                         }
                                                m1 = m1 + mass1[i4]* count_every_atom_type[i4];
                                                }
                   }
            molecule_list[i1][i2].set_molecule_name(line1); //finally assign that molecule the name
            molecule_list[i1][i2].set_molecule_mass(m1);
             
            }
        }
    /*myfile1.open("molecules.txt");    
    for(i1=0;i1<nframes;i1++){
                              for(i2=0;i2<number_molecules[i1];i2++){
                                                                     //cout <<i2 <<" " <<  molecule_list[i1][i2].return_molecule_name() << endl;
                                                                       myfile1 << i2 <<" " <<  molecule_list[i1][i2].return_molecule_name() << endl;
                                                                     }
                              }
     myfile1.close();*/
    //cout << molecule_list[0][0].return_molecule_name() << endl;
    delete[] count_every_atom_type;
}
