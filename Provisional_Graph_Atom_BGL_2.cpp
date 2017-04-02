#include <iostream>
#include<vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include<boost/graph/graph_concepts.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/graph_traits.hpp>
class Atom
{
public:

    char element; // C for carbon, O for oxygen
    int index; // C1, O2
    double x,y,z;
    Atom () {}
    Atom (char e, int i, double x_, double y_, double z_) { element =
e; index = i; x = x_; y=y_; z=z_; }

};
class Molecule
{
public:
  double x = 0, y = 0, z = 0;
  std::vector<Atom> atoms;
  Molecule () {}
  Molecule (int x_, int y_, int z_, std::vector<Atom> atoms_)
{
    x = x_;
    y = y_;
    z = z_;
    atoms = atoms_;
}
};
void operator<<(std::ostream & out,Atom & atom)
{
    out<<atom.element<<atom.index;
}
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS,Atom,boost::no_property,Molecule> molecule_graph;
void print_adjacency_vertex (molecule_graph::vertex_descriptor v_ , molecule_graph m_)
{
    typename boost::graph_traits<molecule_graph>::out_edge_iterator ei, ei_end;
    for(boost::tie(ei, ei_end) = boost::out_edges(v_, m_); ei != ei_end; ++ei)
    {
        std::cout << m_[boost::source(*ei, m_)].element << m_[boost::source(*ei, m_)].index
                  << " --> "
                  << m_[boost::target(*ei, m_)].element << m_[boost::target(*ei, m_)].index
                  << std::endl;

    }
}
int main()
{
    std::vector<Atom> atoms_vec;
    char * molecule_element = "COO";
    Atom *atoms;

    for(int i = 0;i < strlen(molecule_element);i++)
    {
        atoms = new Atom(molecule_element[i], i,0.5+i, 0.5+i, 0.5+i);
        // I have set dummy value for x y and z plane
        atoms_vec.push_back (*atoms);
    }
    Molecule obj_molecule(5, 5, 5, atoms_vec);
    obj_molecule.print ();
    molecule_graph m;
    molecule_graph::vertex_descriptor v[3];
    for(int i = 0;i < 3;i++)
    {
    v[i] = boost::add_vertex(m);
    m[v[i]] = obj_molecule.atoms[i];
    //std::cout<<m[v[i]];
    }
    boost::add_edge (v[0], v[1], m);
    boost::add_edge (v[0], v[2], m);
//    m[graph_bundle].x=obj_molecule.x; //Getting error during setting of property for graph

    int index;
    std::cout<<"Enter index of vertex to retrieve adjacency verties:";
    std::cin>>index;
    print_adjacency_vertex (v[index], m);
    return 0;
}

