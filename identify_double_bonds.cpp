#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include "atom.h"
#include "molecule.h"
#include "functions.h"

using namespace std;

void identify_double_bonds(int nframes1, int natoms1,int number_molecules, atom **at_list1, molecule **molecule_list1, 
                     int temp_neighbours1[], int neighbours_index1[][20], double **lx2, double **ly2, double **lz2, 
                     double *cov_r2, double *cov_r3, int no_bond_params2, double **bond_params2, string *atom_types2, int size2){

     int i,j,k;
     int frames,valency1,valency2;
     int n1,n2,n3;
     int n4,n5,n6;
     int num2,m;
     int count1, count2, count3;
     int at_no1,at_no2,i1,i2;
     int no_db_per_molecule[number_molecules];
     double dx,dy,dz,dr;
     double pi1,d1,d2,d3,d4;
     string s1,s2,s3;
     string s4,s5,s6;
     string sub;
     string status[natoms1];
     ofstream file;
     
     cout << "starting double bond identification " << endl;
     
     for(i=0;i<number_molecules;i++){
                                 no_db_per_molecule[i] = 0;
                                 }
     
     for(i=0;i<natoms1;i++){
                            status[i] = "N";
                            }
     frames = 0;
     s1 = "C";
     s2 = "C";
     
     for(i=0;i<natoms1;i++){
                            n1 = temp_neighbours1[i];
                            s3 = at_list1[frames][i].return_atomname();
                            
                            if(s1 == s3 && status[i] == "N"){
                                  for(k=0;k<size2;k++){
                                       if(s3 == atom_types2[k]){
                                               at_no1 = k+1;
                                               break;
                                               }
                                       }
                                       
                                  for(j=0;j<n1;j++){
                                              n2 = neighbours_index1[i][j] - 1;
                                              s4 = at_list1[frames][n2].return_atomname();
                                              
                                              if(s4 == s2){
                                                    for(k=0;k<size2;k++){
                                                                   if(s4 == atom_types2[k]){
                                                                         at_no2 = k+1;
                                                                         break;
                                                                          }
                                                                   }
                                                                   
                                                    dx = at_list1[frames][i].return_xcord() - at_list1[frames][n2].return_xcord();
                                                    dy = at_list1[frames][i].return_ycord() - at_list1[frames][n2].return_ycord();
                                                    dz = at_list1[frames][i].return_zcord() - at_list1[frames][n2].return_zcord();
                                                    
                                                    //check periodicity
                                                    //check for periodic images in x direction
                                                    if (dx > lx2[frames][2]/2.0){
                                                       dx = dx - lx2[frames][2];
                                                       }
                                                    if(dx <= (-1.0)*(lx2[frames][2]/2.0)){
                                                          dx = dx + lx2[frames][2];
                                                          }
                                                    //check for periodic images in y direction
                                                    if (dy > ly2[frames][2]/2.0){
                                                           dy = dy - ly2[frames][2];
                                                           }
                                                    if(dy <= (-1.0)*(ly2[frames][2]/2.0)){
                                                          dy = dy + ly2[frames][2];
                                                          }
                                                    //check for periodic images in z direction
                                                    if (dz > lz2[frames][2]/2.0){
                                                           dz = dz - lz2[frames][2];
                                                           }
                                                    if(dz <= (-1.0)*(lz2[frames][2]/2.0)){
                                                          dz = dz + lz2[frames][2];
                                                          }
                                                          
                                                    dr = sqrt(dx*dx + dy*dy + dz*dz);
                                                    
                                                    for (m=0;m<no_bond_params2;m++){
                                                                                   if(bond_params2[m][0] == at_no1 && bond_params2[m][1] == at_no2){
                                                                                                        num2 = m;
                                                                                                        //cout << num2 << endl;
                                                                                                        break;
                                                                                                        }
                                                                                    if(bond_params2[m][0] == at_no2 && bond_params2[m][1] == at_no1){
                                                                                                         num2 = m;
                                                                                                         break;
                                                                                                         }
                                                                                     }
                                                    if(num2 > no_bond_params2){
                                                            cout << "something is wrong." << endl;
                                                            //system("Pause");
                                                            //exit(1);
                                                            }
                                                    
                                                    //first pi bond
                                                    pi1 = 0.0;
                                                    d1 = bond_params2[num2][11];
                                                    d2 = bond_params2[num2][12];
                                                    
                                                    i1 = int(at_no1-1);
                                                    i2 = int(at_no2-1);
                            
                                                    if(cov_r3[i1] <=0.0 || cov_r3[i2] <=0.0){
                                                                   pi1 = 0.0;
                                                                   }else{
                                                                         d3 = (cov_r2[i1]+cov_r2[i2])/(2.0);
                                                                         d4 = dr/d3;
                                                                         pi1 = exp(d1*(pow(d4,d2)));
                                                                         }
                                                    if(pi1 > 0.5){
                                                           n6 = at_list1[frames][i].return_mol_no() - 1;
                                                           no_db_per_molecule[n6] = no_db_per_molecule[n6] + 1;
                                                           status[i] = "Y";
                                                           status[n2] = "Y";
                                                           }
                                                    }
                                              }
                                  }
                            }
     
     
     for(i=0;i<number_molecules;i++){
                                 //no_db_per_molecule[i] = no_db_per_molecule[i]/2;    
                                 cout << "No of double bonds in molecule " << i+1 << " is " << no_db_per_molecule[i] << endl;;
                                 }
     file.open("double_bonds");
     for(i=0;i<natoms1;i++){
                            if(status[i] == "Y"){
                                         file << i+1 << endl;
                                         }
                            }
     file.close();
}
