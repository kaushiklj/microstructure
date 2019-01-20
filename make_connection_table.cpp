#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <iomanip>
#include <cstdlib>
#include "atom.h"
#include "functions.h"

//using namespace std;

void make_connection_table(ofstream &file1,int nframes1, int natoms1, atom **at_list1, string *atom_types1,
                           int size1, int no_bond_params,double **lx1, double **ly1, double **lz1,double *alpha1, double *beta1, 
                           double *gama1, double **bond_params, double *valency,double *cov_r1, double *cov_r2, double *cov_r3, 
                           int ****grid_list1, int ***grid_atom_count1,int nx1,int ny1,int nz1){
                                  

          int i,j,k,m,n;
          int p,q,r;
          int plow,phigh,qlow,qhigh,rlow,rhigh;
          int nid,bin_pop;
          int num2,count;
          int i1, i2;
          int at_no1,at_no2;
          //int neighbours[15];
          int no_neighbours;
          int check_flag,flag1;
          double bond_orders[15];  //maximum 15 neighbours are assumed
          double dx,dy,dz,dr;
          double x1,y1,z1,x2,y2,z2;
          double total_bond_order;
          int xbin1,ybin1,zbin1,xbin2,ybin2,zbin2;
          int bin_indexes[6][3];
          double sigma, pi1,pi2,bo;
          double d1,d2,d3,d4;
          int file_write_flag = 1; //1 means write file1 and file2
          string temp,cond;
          string filename;
          ostringstream convert;
          ofstream file2,file3;
                  
          //cout << "starting connection table " << endl;
          convert << (nframes1+1);
          filename = convert.str();
          //file3.open((filename+"_connection_table.txt").c_str());
          if (file_write_flag == 1) file2.open("temp_connection_table.txt");
          //identify atom type number
          for (i=0;i<1;i++){
              //write connection table
              if (file_write_flag == 1){
                  file1 << "Frame no " << nframes1+1 << endl;
                  file2 << "Frame no " << nframes1+1 << endl;
              }
              //initialize neighbours arrays and bond_order array
              //for(j=0;j<natoms1;j++){
              //                       no_neighbours[j] = 0;
              //                       for(k=0;k<15;k++){
              //                                         neighbours[j][k]  = 0;
              //                                         bond_orders[j][k] = 0;
              //                                         }
              //                       }
              
              //cout << "for fram " << i << " last atom has atom no " << at_no[natoms1-1] << endl;
              for (j=0;j<natoms1;j++){
                  check_flag = 0;
                  
                  //initialize all neighbour and bond order related information
                  //total_bond_order1[i][j] = 0.0;
                  no_neighbours = 0;
                  total_bond_order = 0.0;
                  /*for(k=0;k<15;k++){
                                    neighbours[k] = 0;
                                    neighbours[k] = 0;
                                    }*/
                  
                  temp.clear();
                  temp = at_list1[i][j].return_atomname();
                  
                  for(k=0;k<size1;k++){
                                       if(temp == atom_types1[k]){
                                               at_no1 = k+1;
                                               break;
                                               }
                                       }
                  x1 = at_list1[i][j].return_xcord();
                  y1 = at_list1[i][j].return_ycord();
                  z1 = at_list1[i][j].return_zcord();

                  xbin1 = at_list1[i][j].return_xbin();
                  ybin1 = at_list1[i][j].return_ybin();
                  zbin1 = at_list1[i][j].return_zbin();

                  plow = qlow = rlow = -1;
                  phigh = qhigh = rhigh = 2;
                  
                  //following part is needed becuase some cells near pbs can have only 1-2 atoms
                  if (xbin1 == 0 || xbin1 == 1)         plow = -2;
                  if (xbin1 == nx1-1 || xbin1 == nx1-2) phigh = 3;
                  if (ybin1 == 0 || ybin1 == 1)         qlow = -2;
                  if (ybin1 == ny1-1 || ybin1 == ny1-2) qhigh = 3;
                  if (zbin1 == 0 || zbin1 == 1)         rlow = -2;
                  if (zbin1 == nz1-1 || zbin1 == nz1-2) rhigh = 3;
                  count = 0;
                  if (xbin1 == 0 || xbin1 == nx1-1 || ybin1 == 0 || ybin1 == ny1-1 || zbin1 == 0 || zbin1 == nz1-1) check_flag = 1;
                  
                  for (p=plow;p<phigh;p++){
                       xbin2 = xbin1 + p;
                       if (xbin2 < 0) xbin2 = nx1 + xbin2;
                       if (xbin2 > nx1-1) xbin2 = xbin2 - nx1;
                       for (q=qlow;q<qhigh;q++){
                            ybin2 = ybin1 + q;
                            if (ybin2 < 0) ybin2 = ny1 + ybin2;
                            if (ybin2 > ny1-1) ybin2 = ybin2 - ny1;
                           for (r=rlow;r<rhigh;r++){
                                zbin2 = zbin1 + r;
                                if (zbin2 < 0) zbin2 = nz1 + zbin2 ;
                                if (zbin2 > nz1-1) zbin2 = zbin2 - nz1;
                                bin_pop = grid_atom_count1[xbin2][ybin2][zbin2];
                                for (n=0;n<bin_pop;n++){
                                     //if (j == 1841) cout << xbin2 << " " << ybin2 << " " << zbin2 << endl;
                                     nid = grid_list1[xbin2][ybin2][zbin2][n];
                                     temp.clear();
                                     temp = at_list1[i][nid].return_atomname();
                                     for(m=0;m<size1;m++){
                                           if(temp == atom_types1[m]){
                                               at_no2 = m+1;
                                               break;
                                           }
                                     }
                                     x2 = at_list1[i][nid].return_xcord();
                                     y2 = at_list1[i][nid].return_ycord();
                                     z2 = at_list1[i][nid].return_zcord();
                                     //if (j == 3965 && nid == 3964) cout << "found the pair" << endl;
                                     dx = x2-x1;
                                     dy = y2-y1;
                                     dz = z2-z1;
                      
                                     //check periodicity
                                     //check for periodic images in x direction
                                     if (dx > lx1[i][2]/2.0){
                                        dx = dx - lx1[i][2];
                                     }
                                     if(dx <= (-1.0)*(lx1[i][2]/2.0)){
                                        dx = dx + lx1[i][2];
                                     }
                                     //check for periodic images in y direction
                                     if (dy > ly1[i][2]/2.0){
                                         dy = dy - ly1[i][2];
                                     }
                                     if(dy <= (-1.0)*(ly1[i][2]/2.0)){
                                         dy = dy + ly1[i][2];
                                     }
                                     //check for periodic images in z direction
                                     if (dz > lz1[i][2]/2.0){
                                         dz = dz - lz1[i][2];
                                     } 
                                     if(dz <= (-1.0)*(lz1[i][2]/2.0)){
                                        dz = dz + lz1[i][2];
                                     }
                      
                                     dr = sqrt(dx*dx + dy*dy + dz*dz);
                                     num2 = 200;
                                     if(dr < 4.0 && j != nid){
                            
                                     //cout << "selected atoms are " << at_no[j] << " " << at_no[k] << endl;
                                     for (m=0;m<no_bond_params;m++){
                                          if(bond_params[m][0] == at_no1 && bond_params[m][1] == at_no2){
                                                     num2 = m;
                                                     //cout << num2 << endl;
                                                     break;
                                          }
                                          if(bond_params[m][0] == at_no2 && bond_params[m][1] == at_no1){
                                                     num2 = m;
                                                     break;
                                          }
                                     }   
                            
                                     if(num2 > no_bond_params){
                                        cout << "something is wrong." << endl;
                                        //system("Pause");
                                        //exit(1);
                                     }
                                     //sigma bond
                                     sigma = 0.0;
                                     d1 = bond_params[num2][14];
                                     d2 = bond_params[num2][15];
                                     i1 = int(at_no1-1);
                                     i2 = int(at_no2-1);
                                     d3 = (cov_r1[i1]+cov_r1[i2])/(2.0);
                                     d4 = dr/d3;

                                     sigma = exp(d1*(pow(d4,d2)));
                            
                                     //first pi bond
                                     pi1 = 0.0;
                                     d1 = bond_params[num2][11];
                                     d2 = bond_params[num2][12];
                            
                                     if(cov_r2[i1] <=0.0 || cov_r2[i2] <=0.0){
                                           pi1 = 0.0;
                                     }else{
                                           d3 = (cov_r2[i1]+cov_r2[i2])/(2.0);
                                           d4 = dr/d3;
                                           pi1 = exp(d1*(pow(d4,d2)));
                                     }
                                     //cout << "for atoms " << j << " and " << k << " line is " << pi1 << endl;                            
                                     //double/2nd pi bond
                                     pi2 = 0.0;
                                     d1 = bond_params[num2][6];
                                     d2 = bond_params[num2][8];
                            
                                     if(cov_r3[i1] <=0.0 || cov_r3[i2]<=0.0){
                                           pi2 = 0.0;
                                     }else{
                                           d3 = (cov_r3[i1]+cov_r3[i2])/(2.0);
                                           d4 = dr/d3;
                                           pi2 = exp(d1*(pow(d4,d2)));
                                     }
                            
                                     //total bond order
                                     bo = sigma + pi1 + pi2;
                            
                                     //cout << " BO between " << j+1 << "  "<< k+1 << " " << sigma << " " << pi1 << " " << pi2 << " " << bo << endl;
                            
                                     if( bo > 0.70){
                                     //if( dr < 1.6){
                                         //if edge cell, then check to make sure this atom is not already in the neighbor list
                                         // this is required for very small grid like only 2 cells in any direction 
                                         flag1 = 0;
                                         if (check_flag == 1){
                                             for (m=0;m < no_neighbours;m++){
                                                  if (at_list1[i][j].neigh_indexes[m] == nid+1) flag1 = 1;
                                             }
                                         }
                                         if (flag1 == 0){
                                             count = no_neighbours;
                                             at_list1[i][j].set_neigh_index(count,nid+1);
                                             //neighbours[count] = nid+1;
                                             //bond_orders[count] = bo;
                                             no_neighbours = no_neighbours + 1;
                                             total_bond_order = total_bond_order + bo;
                                             if (no_neighbours > max_neigh){
                                                 cout << "atom " << j+1 << " has more than 16 neighbors. Change neighbors array size " << endl;
                                                 exit(EXIT_FAILURE);
                                             }
                                             //total_bond_order1[i][j] = total_bond_order1[i][j] + bo;
                            
                                             //count = no_neighbours[k];
                                             //neighbours[k][count] = j+1;
                                             //bond_orders[k][count] = bo;
                                             //no_neighbours[k]= no_neighbours[k] + 1;
                                         }
                                     } //bond order if ends here
                            
                                  }//dr if ends here
                                }
                           }
                       }
                  }
               at_list1[i][j].set_num_neighbors(no_neighbours);   
               at_list1[i][j].set_total_bo(total_bond_order);
               if (file_write_flag == 1){
                  file1 << setw(8) << j+1 << " " << setw(3) <<at_list1[i][j].num_neighbors << " " ;
                  file2 << setw(8) << j+1 << " " << setw(3) <<at_list1[i][j].num_neighbors << " " ;
                  //file3 << setw(8) << j+1 << " " << setw(3) <<no_neighbours << " " << endl ;
                  for(k=0;k<no_neighbours;k++){
                                              n = at_list1[i][j].neigh_indexes[k];
                                              //file1 << setw(2) << at_list1[i][n-1].return_atomname() << setw(5) 
                                              //<< neighbours[k] << " "<< setw(10)<< bond_orders[k] << " " ;
                                              file1 << setw(2) << at_list1[i][n-1].return_atomname() << " " << setw(8) 
                                              << n << " " ;
                                              file2 << setw(2) << at_list1[i][n-1].return_atomname() << " " << setw(8) 
                                              << n << " " ;
                                              }
                  file1 << endl;   
                  file2 << endl;
               }
              } // j for loop ends here
          
          
          }//no of frames loop
          
          if (file_write_flag == 1){
              file1.close();
              file2.close();
          }
          //file3.close();
          cout << "finished preparing connection table " << endl;
          
          //cout << at_no[nframes1-1][natoms1-1] << endl;

}
