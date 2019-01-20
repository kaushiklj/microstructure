#include <iostream>
#include <fstream>
#include <cstring>
#include "functions.h"

using namespace std;

void find_valence_params(string *line1, int start_line_per_section[], int num_atom_types, double *mass1,
                 double *valency1, double *cov_r1, double *cov_r2, double *cov_r3){
     int i,j,k;
     int num1;
     double dans;
     string temp1,sub;
     char *temp2;
     istringstream iss;
     
     k = start_line_per_section[1];
     //cout << "starting valence parameters " << endl;
     for (i=0;i<num_atom_types;i++){
                                   temp1 = line1[k];
                                   iss.str(temp1);
                                   j=0;
                                   //cout << temp1 << endl;
                                   while (iss >> sub){
                                         //cout << sub << endl;
                                         
                                         //find valency
                                         if(j==2){
                                                  temp2 = new char[sub.size()+1];
                                                  strcpy(temp2, sub.c_str());
                                                  dans = string_to_double(temp2);
                                                  //cout << dans << endl;
                                                  valency1[i] = dans;
                                                  //j= j +1;
                                                  delete[] temp2;
                                                  }
                                         //find first covalent radius
                                         if(j==1){
                                                  temp2 = new char[sub.size()+1];
                                                  strcpy(temp2, sub.c_str());
                                                  dans = string_to_double(temp2);
                                                  //cout << dans << endl;
                                                  cov_r1[i] = dans;
                                                  //j= j +1;
                                                  delete[] temp2;
                                                  }
                                         //find second covalent radius
                                         if(j==7){
                                                  temp2 = new char[sub.size()+1];
                                                  strcpy(temp2, sub.c_str());
                                                  dans = string_to_double(temp2);
                                                  //cout << dans << endl;
                                                  cov_r2[i] = dans;
                                                  //j= j +1;
                                                  delete[] temp2;
                                                  }
                                         //find atomic mass 
                                         if(j==3){
                                                  temp2 = new char[sub.size()+1];
                                                  strcpy(temp2, sub.c_str());
                                                  dans = string_to_double(temp2);
                                                  //cout << dans << endl;
                                                  mass1[i] = dans;
                                                  //j= j +1;
                                                  delete[] temp2;
                                                  }
                                         j =j +1;
                                         //cout << j << endl;
                                         }
                                     
                                   //cout << j << endl;
                                   temp1.clear();
                                   iss.clear();
                                   
                                   temp1 = line1[k+2];
                                   iss.str(temp1);
                                   j=0;
                                   while(iss >> sub){
                                             if(j==0){
                                                  temp2 = new char[sub.size()+1];
                                                  strcpy(temp2, sub.c_str());
                                                  dans = string_to_double(temp2);
                                                  //cout << dans << endl;
                                                  cov_r3[i] = dans;
                                                  //j= j +1;
                                                  delete[] temp2;
                                                  }
                                             j = j + 1;
                                             }
                                   iss.clear();
                                   k = k + 4;                                  
                                   }
     
     }
