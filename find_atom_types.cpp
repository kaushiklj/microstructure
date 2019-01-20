#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include "atom.h"
#include "functions.h"

using namespace std;

void find_atom_types(istream& myfile, string atom_types1, string input_atom_types1){
     
     int i,j,k;
     int n1,n2,n3;
     int frames;
     int count1,count2;
     string s1,s2,s3,sub;
     string temp[50];
     istringstream iss;
     char *c1;
     
     
     while(!feof(myfile)){
                          s2.clear();
                          getline (myfile,s2);
                          iss.clear();
                          iss.str(s2);
                          count2 = 0;
                          while(iss >> sub){
                                           temp[count2] = sub;
                                           count2 = count2 + 1;
                                           }
                          if (temp[0] == "pair_coeff"){
                                                       break;
                                                      }
                          }//main for loop ends here
     for (i=4;i<count2;i++){
         c1 = new char[temp[i].length() +1]; 
         strcpy(c1,str.c_str());i
         if (isdigit(c1[0])){
             //Its a number
            }else{
             // Its a characher
            }
         delete[] c1;
         }            
}
