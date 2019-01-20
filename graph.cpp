#include "graph.h"
#include <list>

using namespace std;

graph::graph(int V){
    this->V= V;
    neigh = new list<int>[V];
}

void graph::addEdge(int v, int w){
     //since its undirected graph, we have update neigh_list of both vertices
     adj[v].push_back(w);
     adj[w].push_back(v);
}

bool graph:isCyclicUtil(int v, bool visited[], int parent){
   visited[v] = true;

   list<int>::iterator i;
   for(i=neigh[v].begin(); i!=neigh[v].end(); ++i){
       if (!visited[*i]){
           if (isCyclicUtil(*i,visited,v)){
               return true;
           }
       }else if(*i != parent) return true;
   }
   return false;
}

/*bool graph:: isCyclic(){
   bool *visited = new bool[V];
   for (int i= 0; i< V;i++) visited[i] = false;

   for (int i=0;i<V;i++){
       if (!visited[i])
           if (isCyclicUtil(i, visited,-1)) return true;
   }
}*/
