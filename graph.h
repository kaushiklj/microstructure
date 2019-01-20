#ifndef GRAPH_H
#define GRAPH_H
#include <list>

using namespace std;

class graph{
    int V; //number of vertices
    int cycles; //number of cycles
    list<int> *neigh; //neighbor list of vertices
    bool isCyclicUtil(int v,bool visited[],int parent);
 public:
    graph(int V);
    void addEdge(int v, int w);
    bool isCyclic(); //returns true if there is a cycle
};
#endif
