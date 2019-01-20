#include <iostream>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include "functions.h"
#include "graph.h"

using namespace std;

//Find the equation of plane consisting three point in 3x3 array points1
//equation of plane is assumed as ax+by+cz+d=0; ,<a,b,c> is the vector 
//normal to the plane. result is stored in 1x4 params array.
//currently function does not re-maps coordinates at periodic boundary. This should
//not matter for equation of plane.
void equation_of_plane(double (&points1)[3][3], double (&params)[4]){
     int i,j,k;
     double vec1[3],vec2[3];

     for (i=0;i<3;i++){
          vec1[i] = points1[1][i] - points1[0][i];
          vec2[i] = points1[2][i] - points1[0][i];
     }
 
     params[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1];
     params[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2];
     params[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0];

     params[3] = -1*(params[0]*points1[0][0] + params[1]*points1[0][1] + params[2]*points1[0][2]);
}
