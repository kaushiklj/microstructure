#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "atom.h"
#include "molecule.h"
#include "functions.h"

using namespace std;

void identify_c3_molecules(int nframes1, int act_frame1,int natoms1,int number_molecules, atom **at_list1, molecule **molecule_list1){

     int i,j,k;
     int frames;
     int n1,n2,n3;
     int count1,count3,global_count;
     int *c3_per_molecule;
     int c_neigh_pop[12];
     string s1,s2,s3;
     string sub;
     ofstream file1,file2;
     
     cout << "starting c3 identification " << endl;

     c3_per_molecule = new int[number_molecules];
     
     for(i=0;i<number_molecules;i++){
         c3_per_molecule[i] = 0;
     }

     for (i=0;i<6;i++){
         c_neigh_pop[i] = 0;
     }
     
     file1.open("c3_fragments.out",ios::app);
     if (act_frame1 == 0){
         file2.open("c3_pop.dat");
     }else{
         file2.open("c3_pop.dat",ios::app);
     }

     frames = 0;
     s1 = "C";
     
     //cout << "reached here " << endl;
     //cout << c3_per_molecule[0] << endl;
     for(i=0;i<natoms1;i++){
        //cout << i << endl; 
        if(at_list1[frames][i].return_atomname() == s1){
           count1 = 0;
           n1 = at_list1[frames][i].num_neighbors; //temp_neighbours1[i];
           //cout << "number of neighbors is " << n1 << endl;
           for(j=0;j<n1;j++){
               n2 = at_list1[frames][i].neigh_indexes[j]-1; //neighbours_index1[i][j] - 1;
               s2 = at_list1[frames][n2].return_atomname();
               if(s2 == s1){
                  count1 = count1 + 1;
               }
           }
           //cout <<  c3_per_molecule[0] << endl;
           if (count1 > 12){
               cout << "More than six carbon neighbors for atom " << i+1 << endl;
               system ("Pause");
               exit(1);
           }else{
               c_neigh_pop[count1-1] = c_neigh_pop[count1-1] + 1;
           }
           if(count1 == 3){
              //cout << "reached here" << endl;
              n3 = at_list1[frames][i].return_mol_no() - 1;
              c3_per_molecule[n3] = c3_per_molecule[n3] + 1;                                        
           }
         }
         //cout << i+1 << endl;
     }
     
     for(i=0;i<number_molecules;i++){
         //cout << i << endl;
         //cout << no_c3_per_molecule[i] << endl;
         //cout << "No of c3 groups in molecule " << i+1 << " is " << no_c3_per_molecule[i] << endl;
         //file1 << molecule_list1[frames][i].return_molecule_name()<< " " << no_c3_per_molecule[i] << endl;
     }
     file2 << act_frame1+1 << "  " ;
     for (i=0;i<12;i++){
          file2 << c_neigh_pop[i] << " ";
     }
     file2 << endl;
     file1.close();
     file2.close();
     delete[] c3_per_molecule; 
     cout << "Finished identify_c3 function" << endl;
}

void identify_c3(int nframes1, int act_frame1,int natoms1, atom **at_list1){

     int i,j,k;
     int frames;
     int n1,n2,n3;
     int count1,count3,global_count;
     int c_neigh_pop[12];
     string s1,s2,s3;
     string sub;
     ofstream file1,file2;
     
     cout << "starting c3 identification " << endl;

     for (i=0;i<6;i++){
         c_neigh_pop[i] = 0;
     }
     
     if (act_frame1 == 0){
         file2.open("c3_pop.dat");
     }else{
         file2.open("c3_pop.dat",ios::app);
     }

     frames = 0;
     s1 = "C";
     
     cout << "reached here " << endl;
     //global_count = 0;
     for(i=0;i<natoms1;i++){
        if(at_list1[frames][i].return_atomname() == s1){
           count1 = 0;
           n1 = at_list1[frames][i].num_neighbors; //temp_neighbours1[i];
           c_neigh_pop[n1-1] = c_neigh_pop[n1-1] + 1;
           //cout << "number of neighbors is " << n1 << endl;
           /*for(j=0;j<n1;j++){
               n2 = at_list1[frames][i].neigh_indexes[j]-1; //neighbours_index1[i][j] - 1;
               s2 = at_list1[frames][n2].return_atomname();
               if(s2 == s1){
                  count1 = count1 + 1;
               }
           }*/
           /*if (count1 > 12){
               cout << "More than six carbon neighbors for atom " << i+1 << endl;
               system ("Pause");
               exit(1);
           }else{
               c_neigh_pop[count1-1] = c_neigh_pop[count1-1] + 1;
           }*/
         }
         //cout << i+1 << endl;
     }
     
     file2 << act_frame1+1 << "  " ;
     for (i=0;i<12;i++){
          file2 << c_neigh_pop[i] << " ";
     }
     file2 << endl;
     file2.close();
     cout << "Finished identify_c3 function" << endl;
}

//Following function will be used only for saturating c2 carbon atoms with H atoms in a ring
void identify_c2(int act_frame1,int natoms1, atom **at_list1, int *c2_flag){

     int i,j,k;
     int frames;
     int n1,n2,n3;
     int count1,count3,global_count;
     int c_neigh_pop[12];
     string s1,s2,s3;
     string sub;
     
     cout << "starting c3 identification " << endl;

     frames = 0;
     s1 = "C";
     
     //cout << "reached here " << endl;
     //global_count = 0;
     for(i=0;i<natoms1;i++){
        if(at_list1[frames][i].return_atomname() == s1){
           count1 = 0;
           n1 = at_list1[frames][i].num_neighbors; //temp_neighbours1[i];
           //cout << "number of neighbors is " << n1 << endl;
           for(j=0;j<n1;j++){
               n2 = at_list1[frames][i].neigh_indexes[j]-1; //neighbours_index1[i][j] - 1;
               s2 = at_list1[frames][n2].return_atomname();
               if(s2 == s1){
                  count1 = count1 + 1;
               }
           }
           if (count1 == 2){
               c2_flag[i] = 1;
           }
         }
         //cout << i+1 << endl;
     }
     
}
