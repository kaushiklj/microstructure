#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include "functions.h"

using namespace std;

void identify_bond_parameters(double **bond_parameters1, string *lines1, int start_line_per_section1[], int no_lines1){
     int i,j,k;
     double dans;
     string sub,sub1;
     string line1,line2;
     istringstream iss,iss1;
     char *temp1;
     
     i= start_line_per_section1[2];
     k=0;
     //cout << start_line_per_section1[3]-1 << endl;
     while (k < no_lines1){
           //cout << k << endl;
           j=0;
           line1 = lines1[i];
           line2 = lines1[i+1];
           iss.clear();
           iss1.clear();
           iss.str(line1);
           iss1.str(line2);
           while (iss >> sub){
                 temp1 = new char[sub.size()+1];
                 strcpy(temp1, sub.c_str());
                 dans = string_to_double(temp1);
                 bond_parameters1[k][j] = dans;
                 j= j +1;
                 delete[] temp1;
                 } 
           while (iss1 >> sub1){
                 temp1 = new char[sub1.size()+1];
                 strcpy(temp1, sub1.c_str());
                 dans = string_to_double(temp1);
                 bond_parameters1[k][j] = dans;
                 j= j +1;
                 delete[] temp1;
                 }
           //cout << "parameter completed is " << k << endl;
           i = i +2;
           k = k +1;            
          }
     
     //cout << bond_parameters1[16][0] << " " << bond_parameters1[16][1] << endl;
     }
