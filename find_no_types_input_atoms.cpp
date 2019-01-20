#include <fstream>
#include <iostream>
#include "functions.h"
#include "atom.h"

using namespace std;

int find_no_types_input_atoms(atom **at_list, int natoms1, string input_atom_types1[]){
    int i,j,k;
    int count1,count2;
    string s1,s2,decision;
    string temp_name[100];
    
    //cout << "finding number of input atom types " << endl;
    count1 = 0;
    count2 =0;
    for(i=0;i<natoms1;i++){
                          s1 = at_list[0][i].return_atomname();
                          if(i==0){
                                   temp_name[count1] = s1;
                                   }else{
                                         count2 = count1+1;
                                         decision = "no";
                                         for(j=0;j<count2;j++){
                                                               if (s1 == temp_name[j]){
                                                                            decision = "yes";
                                                                            break;
                                                                            }
                                                               }
                                         if (decision == "no"){
                                                      count1 = count1 + 1;
                                                      temp_name[count1] = s1;
                                                      }
                                         }
                          }
    
    
    cout << "Total number of atom types are " << count1+1<< endl;
    
    //input_atom_types1 = new string[count1 + 1];
    
    for(i=0;i<count1+1;i++){
                                 input_atom_types1[i] = temp_name[i]; 
                                 cout << "atom type " << i << " is " << input_atom_types1[i] << endl;
                                 }
    
    return count1+1;
}
