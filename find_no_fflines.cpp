#include <iostream>
#include <fstream>
#include "functions.h"

using namespace std;

int find_no_fflines(){
     using namespace std;
     
     int i,j,k;
     int a,nlines1;
     string line;
     string filename;
     ifstream myfile;
     istringstream instream;
     
     //cout << "Counting number of lines in file" << endl;
     
     nlines1 =0;
     myfile.open("ffield");
     if(myfile){
                while (getline(myfile, line))
                      {
                       nlines1 = nlines1 + 1;
                       }
                 }else{
                       cout << "file not found" << endl;
                       //system("Pause");
                       //exit(1);
                       }
     myfile.clear();
     myfile.seekg(0, ios::beg);
     myfile.close();
     
     //cout <<"Number of lines in force-field file are " << nlines1 << endl;
     
     return nlines1;
}
