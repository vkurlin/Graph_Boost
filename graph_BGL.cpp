#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include<boost/graph/graph_concepts.hpp>
using namespace std;
using namespace boost;
typedef adjacency_list< listS, vecS, directedS > digraph;
int main()
{
    int vno,eno,first,second;
    cout<<"Enter no. of vertex";
    cin>>vno;
   digraph g(vno);
   cout<<"Enter no. of edges";
   cin>>eno;
   for(int counter=1;counter<=eno;counter++)
   {
        cout<<"Enter pair of edges";
        cin>>first>>second;
        if(first<vno && second<vno)
        {
            add_edge(first,second,g);
        }
   }
   write_graphviz(cout,g);
   return 0;
}
