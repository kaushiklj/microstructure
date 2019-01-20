#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "atom.h"
#include "molecule.h"
#include "functions.h"

using namespace std;

void identify_carbon_3_neighbors(int nframes1, int natoms1,int number_molecules, atom **at_list1, molecule **molecule_list1, 
                     int temp_neighbours1[], int neighbours_index1[][20]){

     int i,j,k;
     int frames;
     int n1,n2,n3;
     int count1, count2, count3;
     int *temp_store;
     int no_c3_per_molecule[number_molecules];
     string s1,s2,s3;
     string sub;
     ofstream file1;
     
     cout << "starting carbon atoms with 3 neighbors identification " << endl;
     
     frames = 0;
     s1 = "C";
     s2 = "C8H17";
     n1 = 3;
     
     file1.open("3_neighbor_carbons.out",ios::app);
     
     count1 = 0;
     for(i=0;i<natoms1;i++){
                            n2 = at_list1[frames][i].return_mol_no() - 1;
                            s3 = molecule_list1[frames][n2].return_molecule_name();
                            if(at_list1[frames][i].return_atomname() == s1 && temp_neighbours1[i] == n1 && s3 == s2){
                                                                     count1 = count1 + 1;
                                                                     }
                            }
     cout << count1 << endl;
     file1 << count1 << endl;
     count2 = 0;
     for(i=0;i<natoms1;i++){
                            n2 = at_list1[frames][i].return_mol_no() - 1;
                            s3 = molecule_list1[frames][n2].return_molecule_name();
                            if(at_list1[frames][i].return_atomname() == s1 && temp_neighbours1[i] == n1 && s3 == s2){
                                                                     //temp_store[count2] = i+1;
                                                                     file1 << i+1 << endl;
                                                                     //count2 = count2 + 1;
                                                                     }
                            }
     
     file1.close();
     
     
     }
