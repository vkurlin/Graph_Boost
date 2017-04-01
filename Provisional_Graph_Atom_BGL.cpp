#include <iostream>
#include<vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include<boost/graph/graph_concepts.hpp>
class Atom
{
public:
    char element; // C for carbon, O for oxygen
    int index; // C1, O2
    double x,y,z;
    Atom() {}
    Atom (char e, int i, double x_, double y_, double z_) { element =
e; index = i; x = x_; y=y_; z=z_; }
};
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS,Atom> molecule;

void print_adjacency_vertex (molecule::vertex_descriptor v_, molecule m_)
{
    typename boost::graph_traits<molecule>::out_edge_iterator ei, ei_end;
    for(boost::tie(ei, ei_end) = boost::out_edges(v_, m_); ei != ei_end; ++ei)
    {
        std::cout << m_[boost::source(*ei, m_)].element
                  << " --> "
                  << m_[boost::target(*ei, m_)].element
                  << std::endl;
    }
}
int main()
{

    molecule m;
    molecule::vertex_descriptor v[3];
    char * molecule_element="COO";
    for(int i=0;i<3;i++)
    {
    v[i] = boost::add_vertex(m);
    m[v[i]].element=molecule_element[i];
    m[v[i]].index=1;
    m[v[i]].x=m[v[i]].y=m[v[i]].z=1;
    }
    boost::add_edge (v[0],v[1], m);
    boost::add_edge (v[0],v[2], m);
    int index;
    std::cout<<"Enter index of vertex to retrieve adjacency verties:";
    std::cin>>index;
    print_adjacency_vertex (v[index], m);
    return 0;
}

