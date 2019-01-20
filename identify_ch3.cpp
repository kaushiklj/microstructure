#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "atom.h"
#include "molecule.h"
#include "functions.h"

using namespace std;

void identify_ch3(int nframes1, int natoms1,int number_molecules, atom **at_list1, molecule **molecule_list1, int temp_neighbours1[], int neighbours_index1[][20]){
     int i,j,k;
     int frames;
     int n1,n2,n3;
     int count1, count2, count3;
     int no_groups_per_molecule[number_molecules];
     string s1,s2,s3;
     string sub;
     istringstream iss;
     ifstream myfile;
     ofstream file1;
     
     cout << "starting ch3 identification " << endl;
     
     file1.open("ch3_fragments.out",ios::app);
     for(i=0;i<number_molecules;i++){
                                 no_groups_per_molecule[i] = 0;
                                 }
     
     frames = 0;
     s1 = "C";
     s2 = "H";
     for(i=0;i<natoms1;i++){
                            count2 = 0;
                            if(at_list1[frames][i].return_atomname() == s1){
                                                                     count1 = 0;
                                                                     n1 = temp_neighbours1[i];
                                                                     for(j=0;j<n1;j++){
                                                                                       n2 = neighbours_index1[i][j] - 1;
                                                                                       s3 = at_list1[frames][n2].return_atomname();
                                                                                       if (s3 == s2){
                                                                                              count1 = count1 + 1;
                                                                                          }       
                                                                                       }
                                                                     if(count1 == 3){ //found ch3 group
                                                                               count2 = count2 + 1;
                                                                               n3 = at_list1[frames][i].return_mol_no() - 1;
                                                                               no_groups_per_molecule[n3] = no_groups_per_molecule[n3] + 1;                                        
                                                                               }
                                                                     }
                            
                            }//for loop of natoms1 ends
     
     for(i=0;i<number_molecules;i++){
                                     //cout << "No of ch3 groups in molecule " << i+1 << " is " << no_groups_per_molecule[i] << endl;
                                     file1 << molecule_list1[frames][i].return_molecule_name()<< " " << no_groups_per_molecule[i] << endl;
                                     }
     file1.close();
     }
