#include <iostream>
#include <fstream>
#include "functions.h"

using namespace std;

void identify_atom_types(string *line1, int start_line_per_section[], int num_atom_types, string *atom_types1){
     int i,j,k;
     string temp1,sub;
     istringstream iss;
     
     k = start_line_per_section[1];
     for (i=0;i<num_atom_types;i++){
                                   temp1 = line1[k];
                                   iss.str(temp1);
                                   j=0;
                                   while (iss >> sub && j <=0){
                                         //cout << sub << endl;
                                         atom_types1[i] = sub;
                                         j= j +2;
                                         } 
                                   k = k + 4;                                  
                                   iss.clear();
                                   }
     
     }
