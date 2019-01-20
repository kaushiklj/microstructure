#include <iostream>
#include <list>
#include <sstream>
#include "function.h"
#include "graph.h"

using namespace std;

const max_cycle_vertices = 20;
const max_neighs = 16;

void cycle_analysis(int nframes1, int act_frame, int natoms1){
   int i,j,k;
   int count1,count2,neigh_index;
   bool *visited;
   string lines,sub;
   string temp_lines[50];
   istringstream iss;
   ifstream ifile1;
   ofstream ofile1;

   //open output file for writing cycle results
   if (act_frame == 0){
      ofile1.open("cycles.dat")
   }else{
      ofile1.open("c3_pop.dat",ios::app);
   }
   
   graph g(natoms1);
   //read temp_connection table and populate neigh_list
   //this part can be very inefficient. Should be removed in the future
   ifile1.open("temp_connection_table.txt");
   for(i=0;i<natoms1+1;i++){
       lines.clear();
       getline(ifile1,lines);
       if (i > 0){
           iss.clear();
           iss.str(lines);
           count1 = 0;
           while(iss >> sub){
               temp_lines[count1] = sub;
               count1 = count1 + 1;
           }
           n = atoi(temp_lines[1].c_str());
           if (n > max_neighs){
               cout << " number of neighbors " << n << "exceeded maximum allowed " << max_neighs << endl;
               exit(EXIT_FAILURE);
           }
           if (n > 0){
               count1 = 0;
               count2 = 3;
               while(count1 < n){
                   neigh_index = atoi(temp_lines[count2].c_str());
                   g1.addEdge(i,neigh_index);
                   count1 = count1 + 1;
                   count2 = count2 + 2;
               }
           }
      }
   }
   ifile1.close();

   //mark all vertices as not visited at the before starting cycle search
   visited = new bool[natoms1];
   for (i=0;i<natoms;i++) visited[i] = false;

  //start cycle search. call recursive helper function defined in graph class
  for (i=0;i<natoms;i++){
       if (!visited[i]){
          if (isCyclicUtil(i,visited,-1)){
             //cycle is found
          }
       }
  }
   ofile1.close();
}
