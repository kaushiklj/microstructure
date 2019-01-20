#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include "atom.h"
#include "molecule.h"
#include "functions.h"

using namespace std;
const int max_neigh_size = 7000;

void identify_molecules(int nframes1, int act_frame, int natoms1, int no_atom_types1, atom **at_list1, string *atom_types1,
                            double **lx1, double **ly1, double **lz1, double *cov_r2, double *cov_r3, 
                            int no_bond_params1, double **bond_params1, int size1, double *mass1){
     cout << "starting molecule identification " << endl;
     int i,j,k;
     int l,m,n;
     int count1,count2,tmp_no,no_cycles;
     int no1,no2,no3,no4;
     int max_neighs;
     double t_mass;
     string no_types_mols[max_mol_types];         
     double mass_mols[max_mol_types];
     int pop_mols[max_mol_types];
     int no_mols[nframes1];
     //int temp_neighbours[natoms1];  //20 is assumed as maximum no of neighbours
     int temp_neighbors;
     cout << "starting molecule identification " << endl;
     //int neighbours_index[natoms1][max_neigh];
     int first_neighbours[max_neigh_size];
     int second_neighbours[max_neigh_size];
     string *mol_names;
     string sub,des;
     string flag;
     string s1,s2;
     string lines;
     string temp_lines[50];
     int flags[natoms1]; //0 means yes, 1 means NO
     cout << "starting molecule identification " << endl;
     molecule **mol_list;
     ifstream myfile;
     ofstream myfile1,myfile2;
     istringstream iss;
     
     //initialize the molecule no status of every atom to No
     for(i=0;i<nframes1;i++){
                             for(j=0;j<natoms1;j++){
                                                    at_list1[i][j].set_atom_state(1);
                                                    }
                             //cout << i << " frame is done " << endl;
                             }
     
     //cout << "finished setting atom state" << endl;
     //myfile.open("temp_connection_table.txt");
     
     for(i=0;i<nframes1;i++){
                             no_mols[i] = 0;
                             
                             //store all the neighbour information of ith frame
                             /*for(j=0;j<natoms1+1;j++){
                                                      lines.clear();
                                                      getline(myfile,lines);
                                                      if (j > 0){
                                                        iss.clear();
                                                        iss.str(lines);
                                                        count1 = 0;
                                                        while(iss >> sub){
                                                                  temp_lines[count1] = sub;
                                                                  count1 = count1 + 1;
                                                                  }
                                                      
                                                        n = atoi(temp_lines[1].c_str()); //no of neighbours
                                                        if (n > max_neigh){
                                                           cout << "number of neighbors " << n << "exceeds maximum allowed "<< max_neigh << "\n";
                                                           exit(EXIT_FAILURE);
                                                        }
                                                        temp_neighbours[j-1] = n;
                                                      
                                                        if(n > 0){
                                                             count1 = 0;
                                                             count2 = 3;
                                                             while(count1 < n){
                                                                          neighbours_index[j-1][count1] = atoi(temp_lines[count2].c_str());
                                                                          count1 = count1 + 1;
                                                                          count2 = count2 + 2;
                                                                          }
                                                             }
                                                        }
                                                      
                                                      }*/
                             //cout << neighbours_index[0][0] << endl;
                             //eliminate all those atoms which have zero neighbors
                             for(j=0;j<natoms1;j++){
                                                    if(at_list1[i][j].num_neighbors == 0){
                                                                       at_list1[i][j].set_atom_state(0);
                                                                       no_mols[i] = no_mols[i] + 1;
                                                                       at_list1[i][j].set_mol_no(no_mols[i]);
                                                                       }
                                                    }
                             //start looking for atom of every molecule
                             cout << "Number of unimolecular pieces are " << no_mols[i] << endl;
                             for(j=0;j<natoms1;j++){
                                 temp_neighbors = at_list1[i][j].num_neighbors;
                                 if(at_list1[i][j].return_atom_state() == 1 && temp_neighbors != 0){
                                    no_mols[i] = no_mols[i] + 1;
                                    //identify first immediate neighbours
                                    no1= temp_neighbors;
                                    //cout << "No of neighbours are " << no1 << endl;
                                    //count1 = 0;
                                    for(k=0;k<no1;k++){
                                        first_neighbours[k] = at_list1[i][j].neigh_indexes[k]-1;
                                        //cout << "First neighbors " << first_neighbours[k] << endl;
                                    }
                                                                                         
                                    flag = "Yes";
                                    no_cycles = 0;
                                    while(flag == "Yes"){
                                          for(k=0;k<natoms1;k++){
                                              flags[k] = 1;
                                          }
                                          count2 = 0;
                                          count1 = 0;
                                          //cout << "for cycle " << no_cycles + 1 << " count1 no is " << count1 << endl;
                                          //cout << no_cycles+1 << " started " << endl;
                                          for(k=0;k<no1;k++){
                                              no2= first_neighbours[k];
                                              if(at_list1[i][no2].return_atom_state() == 1){
                                                 at_list1[i][no2].set_mol_no(no_mols[i]);
                                                 at_list1[i][no2].set_atom_state(0);
                                              }else{
                                                 count2 = count2 + 1;
                                              }
                                                                                                            
                                              no3 = at_list1[i][no2].num_neighbors; //temp_neighbours[no2];
                                              //cout << "k with no2 " << k << " with no2 equal to " << no2 << endl;
                                                                                                            
                                              for(l=0;l<no3;l++){
                                                  no4 = at_list1[i][no2].neigh_indexes[l]-1; //neighbours_index[no2][l] - 1;
                                                  if(at_list1[i][no4].return_atom_state() == 1 && flags[no4] == 1){
                                                     second_neighbours[count1] = at_list1[i][no2].neigh_indexes[l]-1; //neighbours_index[no2][l] - 1;
                                                     flags[no4] = 0;
                                                     count1 = count1 + 1;
                                                  }
                                              }
                                              //cout << "k is " << k << "  no2 is " << no2 << "  no4 is " << no4 << endl;
                                              //cout << at_list1[i][no2].return_atomname() << " " << no2 << " done" << endl;
                                              //cout << count1 << endl;
                                          }
                                          //cout << no_cycles + 1 << " reached until here " << endl;
                                          //cout << count1 << endl;
                                          if (count1 > max_neigh_size ){
                                              cout << "Increase size of second_neighbour array" << endl;
                                              //system("Pause");
                                              exit(1);
                                          }
                                          if(count2 == no1){
                                             flag = "No";
                                             //break;
                                          }else{//transfer second neighbour into the first neighbour
                                             flag = "Yes";
                                             no1 = count1;
                                             for(k=0;k<no1;k++){
                                                 first_neighbours[k] = second_neighbours[k];
                                             }
                                          }
                                          no_cycles = no_cycles + 1;
                                          //cout << "Cycle number is " << no_cycles << " with count1 number "<< count1 << endl;
                                    }
                                                                  
                                    //cout << "atom " << j << " is done " << endl;
                                    }//end of if 
                             }//end of natoms1 for loop
                             
                             cout << "No molecules in frame " << i+1 << " is " << no_mols[i] << endl;
                 } //for loop of nframes1 ends here
     
     //myfile.close();
     
     //cout << "atom with names " << at_list1[0][24].return_atomname() << " has mol no " << at_list1[0][24].return_mol_no() << endl;
     for (i=0;i<25;i++){
         //cout << i << " atom has mol number " << at_list1[nframes1-1][i].return_mol_no() << endl;
     }
     
     //allocate array for molecules
     mol_list = new molecule *[nframes1];
     if (mol_list == NULL){
                   cout <<"could not allocate memory for molecule list \n";
                   //system("Pause");
                   exit(1);
                   }
     for (i=0;i<nframes1;i++){
        mol_list[i] = new molecule[no_mols[i]];
        }
     
     cout << "Begining identifying molecule names" << endl;
     //identify name of every molecule
     identify_molecule_name(nframes1, natoms1, no_mols, no_atom_types1, at_list1, mol_list, atom_types1,mass1);
     cout << "Done identifying molecule names" << endl;
     
     //molecular population analysis
     for(i=0;i<max_mol_types;i++){
                       pop_mols[i]=0;
                       }
     count1=0;
     for(i=0;i<no_mols[0];i++){
                               s2= "N";
                               if(i==0){
                                       no_types_mols[count1]= mol_list[0][i].return_molecule_name();
                                       mass_mols[count1] = mol_list[0][i].return_molecule_weight();
                                       pop_mols[count1] = pop_mols[count1] + 1;
                                       count1 = count1 + 1;
                                       s2 = "Y";
                                       }else{
                                             s1 = mol_list[0][i].return_molecule_name();
                                             for(j=0;j<count1;j++){
                                                                   if(s1 == no_types_mols[j]){
                                                                         pop_mols[j] = pop_mols[j] + 1;
                                                                         s2 = "Y";
                                                                         break;
                                                                         }
                                                                   }
                                             }
                               if(s2 == "N"){
                                     no_types_mols[count1] = mol_list[0][i].return_molecule_name();
                                     mass_mols[count1] = mol_list[0][i].return_molecule_weight();
                                     pop_mols[count1] = pop_mols[count1] + 1;
                                     count1 = count1 + 1;
                                     }
                               }
    
    myfile1.open("molfra.out",ios::app);
    for(i=0;i<count1;i++){
                          //cout << pop_mols[i] << "  " << no_types_mols[i] << endl;
                          myfile1 << act_frame << "   " << pop_mols[i] << "  x  " << no_types_mols[i] << "          " << mass_mols[i] << endl;
                          }
    myfile1 << "Total number of molecules:           " << no_mols[0] << endl;
    myfile1 << "Total number of atoms:         " << natoms1 << endl;
    t_mass = 0.0;
    for(i=0;i<count1;i++){
        t_mass = mass_mols[i]*pop_mols[i] + t_mass;
       }
    myfile1 << "Total system mass:   " << t_mass << endl;
    myfile1.close();

    myfile1.open("modified_molfra.out",ios::app);
    myfile1 << "#" << endl;
    myfile1 << act_frame << "  " << count1 << endl;
    for (i=0;i<count1;i++){
        myfile1 << no_types_mols[i] << "             " << pop_mols[i] << "  " <<  "      " << mass_mols[i] << endl;
    }
    myfile1.close();
    
    //identify CH3 groups
    
     //identify_ch3(nframes1,natoms1,no_mols[0],at_list1,mol_list,temp_neighbours,neighbours_index);
     
     identify_c3_molecules(nframes1,act_frame,natoms1,no_mols[0],at_list1,mol_list);
     
     //identify_carbon_3_neighbors(nframes1,natoms1,no_mols[0],at_list1,mol_list,temp_neighbours,neighbours_index);
     
     //identify_double_bonds(nframes1, natoms1, no_mols[0], at_list1, mol_list, temp_neighbours,neighbours_index,lx1,ly1, lz1, 
     //                      cov_r2,cov_r3, no_bond_params1,bond_params1,atom_types1,size1);
     
     //identify_atoms_in_molecule(natoms1,at_list1,mol_list,no_mols[0]);
     

     
     write_reaction_input_file(nframes1,natoms1,no_mols[0],at_list1,mol_list);

     //delete dynamic memory
     for(i=0;i<nframes1;i++){
                           delete[] mol_list[i];
                           mol_list[i] = NULL;
                           }
     delete[] mol_list;
     }
