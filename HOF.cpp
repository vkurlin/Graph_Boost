#include <fstream>
#include <iostream>
#include <deque>
#include <set>
#include <string>
#include <iterator>
#include <utility>
#include <algorithm>
#include <limits>

//Boost
#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/graph/graphviz.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graph_traits.hpp>
//#include "boost/filesystem/operations.hpp"
//#include "boost/filesystem/path.hpp"
#include <boost/progress.hpp>
#include <boost/config.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/directed_graph.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/filtered_graph.hpp>

//Eigen
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/SparseCholesky>
typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;

// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace cv;

// Constants
const Point3d origin(0, 0, 0);
const double M_PI=3.1472;
//const bool debug = false;
const double distance_error = 1e-2;
const double length_NH = 1.00577;
const double length_HO = 2.62; //1.97436;
const double length_ON = 3.07;
//const double angle_NHO = 100; // degrees
//const double distance_ON = 2.82822; // min distance between N/O atoms in adjacent molecules 
const double big_constant = 1e+8;

Point3d V (Vector v) { return Point3d( v(0), v(1), v(2) ); }
Vector V (Point3d v) { Vector vector(3); vector << v.x, v.y, v.z; return vector; }

class Index_Value
{
public:
    int index;
    double value;
    Index_Value () {}
    Index_Value (int i, double v) { index=i; value=v; }
};
bool Decreasing_Values (Index_Value const& p1, Index_Value const& p2) { return p1.value > p2.value; }
bool Increasing_Values (Index_Value const& p1, Index_Value const& p2) { return p1.value < p2.value; }

double Angle_Positive (Point3d v1, Point3d v2) { return acos( v1.dot( v2 ) / ( norm( v1 ) * norm( v2 ) ) ) * 180 / M_PI; }
double Det (Point2d v1, Point2d v2) { return v1.x * v2.y - v1.y * v2.x; }
double Angle_Signed (Point2d v1, Point2d v2)
{
    double dot_product = v1.dot( v2 ) / ( norm( v1 ) * norm( v2 ) );
    double det = Det( v1, v2 );
    double angle = acos( dot_product ) * 180 / M_PI;
    if ( det >= 0 ) return angle; // range[0,pi)
    else return -angle; // range [-pi,0)
}

Matrix Cross_Product (Point3d v)
{
    Matrix m(3,3);
    m <<
    0,    -v.z, v.y,
    v.z,  0,    -v.x,
    -v.y, v.x,  0;
    return m;
}

Matrix Tensor_Product (Point3d v)
{
    Matrix m(3,3);
    m <<
    v.x * v.x, v.x * v.y, v.x * v.z,
    v.y * v.x, v.y * v.y, v.y * v.z,
    v.z * v.x, v.z * v.y, v.z * v.z;
    return m;
}

class Box
{
public:
    double a = 0, b = 0, c = 0; // sides
    double alpha = 0, beta = 0, gamma = 0; //angles
    double delta = 0, nu = 0; // delta is the nu is the vertical angle
    Matrix matrix = Eigen::ArrayXXd::Zero( 3, 3 );
    int n_molecules = 0;
    Box () { matrix << 0, 0, 0, 0, 0, 0, 0, 0, 0; }
    void Print()
    {
        std::cout<<" a="<<a<<" b="<<b<<" c="<<c<<" al="<<alpha<<" be="<<beta<<" ga="<<gamma<<" de="<<delta<<" nu="<<nu<<" m="<<n_molecules;
    }
    void Find_Matrix ()
    {
        double al = alpha * M_PI / 180;
        double be = beta * M_PI / 180;
        double ga = gamma * M_PI / 180;
        double cos_delta = ( cos(al) - cos(be) * cos(ga) ) / ( sin(be) * sin(ga) );
        delta = acos( cos_delta );
        matrix(0,0) = a;
        matrix(0,1) = b * cos(ga);
        matrix(1,1) = b * sin(ga);
        matrix(0,2) = c * cos(be);
        matrix(1,2) = c * sin(be) * cos_delta;
        matrix(2,2) = c * sin(be) * sin( delta );
        nu = acos( sin(be) * sin( delta ) ) * 180 / M_PI;
        delta *= 180 / M_PI;
    }
    Point3d Abs_Position (Point3d point) { return V( matrix * V( point ) ); }
};

class Atom
{
public:
    char element; // C for carbon, O for oxygen
    int index; // C1, O2
    cv::Point3d point_a; // coordinates in the orthogonal system
    cv::Point3d point_b; // fractional coordinates in the box
    Atom() {}
    Atom (char e, int i, double x, double y, double z) { element = e; index = i; point_b = Point3d( x, y, z ); }
    void Print() { std::cout<<"\n"<<element<<index<<" p_a="<<point_a<<" p_b="<<point_b; }
};

class Molecule
{
public:
    Point3d centre = Point3d(0,0,0), c = Point3d(0,0,0), carbon_axis = Point3d(0,0,0);
    std::map< char, std::vector<Atom> > atoms;
    std::vector<int> indices_NH;
    bool internal = false; // flag of a molecule connected to all neighbours
    double carbon_angle_ver = 0, carbon_angle_hor = 0, oxygen_angle = 0;
    std::vector<Point3d> oxygen_rays;

    Molecule () {}
    Molecule (Molecule& molecule_old, Point3d shift)
    {
        // Copy angles
        carbon_angle_ver = molecule_old.carbon_angle_ver;
        carbon_angle_hor = molecule_old.carbon_angle_hor;
        oxygen_angle = molecule_old.oxygen_angle;
        carbon_axis = molecule_old.carbon_axis;
        // Shift atoms
        centre = molecule_old.centre + shift;
        atoms = molecule_old.atoms;
        indices_NH = molecule_old.indices_NH;
        for ( auto it = atoms.begin(); it != atoms.end(); it++ )
            for ( int i = 0; i < it->second.size(); i++ )
                it->second[i].point_a += shift;
    }

    void Shift (Point3d shift, Molecule& molecule_new)
    {
        // Copy angles
        molecule_new.carbon_angle_ver = carbon_angle_ver;
        molecule_new.carbon_angle_hor = carbon_angle_hor;
        molecule_new.oxygen_angle = oxygen_angle;
        // Shift atoms
        molecule_new.centre = centre + shift;
        molecule_new.atoms = atoms;
        for ( auto it = molecule_new.atoms.begin(); it != molecule_new.atoms.end(); it++ )
            for ( int i = 0; i < it->second.size(); i++ )
                it->second[i].point_a += shift;
    }
    void Print()
    {
        std::cout<<"\ncentre="<<centre<<" C_v="<<carbon_angle_ver<<" C_h="<<carbon_angle_hor<<" O_a="<<oxygen_angle;
        for ( int i = 0; i < indices_NH.size(); i++ ) std::cout<<" N"<<i+1<<"H"<<indices_NH[i]+1;
    }

    void Find_Angles ()
    {
        carbon_axis = atoms['C'][0].point_a - atoms['C'][13].point_a;
        if ( carbon_axis.z < 0 or ( carbon_axis.z == 0 and carbon_axis.y < 0 )
            or ( carbon_axis.z == 0 and carbon_axis.y == 0 and carbon_axis.x < 0 )) carbon_axis *= -1; // goes upwards
        //std::cout<<"\ncarbon_axis="<<carbon_axis;
        carbon_angle_ver = Angle_Positive( cv::Point3d(0,0,1), carbon_axis );
        if ( carbon_angle_ver > 0 ) // only if carbon_axis can be projected to xy-plane
        {
            Point2d carbon_axis_xy( carbon_axis.x, carbon_axis.y );
            carbon_angle_hor = Angle_Signed( Point2d(1,0), carbon_axis_xy );
        }
        // Find oxygen_rays; 
        for ( int i = 0; i < atoms['O'].size(); i++ ) oxygen_rays.push_back( atoms['O'][i].point_a - centre );
        Matrix rotation = Eigen::MatrixXd::Identity(3,3);
        double c = cos( carbon_angle_ver * M_PI / 180 );
        rotation *= c;
        // Find the rotation matrix to make the carbon axis vertical 
        if ( carbon_angle_ver > 0 )
        {
            cv::Point3d a( 0, 0, 0 ); // rotation_axis
            a.x = carbon_axis.y;
            a.y = -carbon_axis.x; // a = rotation_axis is orthogonal to carbon_axis
            a *= 1.0 / norm( a );
            double s = sin( carbon_angle_ver * M_PI / 180 );
            rotation += s * Cross_Product( a ) + (1-c) * Tensor_Product( a ); //std::cout<<"\nr="<<rotation;
        }
        std::vector<Vector> oxygen_vectors( 3 );
        for ( int i = 0; i < 3; i++ ) oxygen_vectors[ i ] = rotation * V( oxygen_rays[ i ] );
        // Find 3 angles with oxygen rays
        int min_index = -1;
        double min_abs_angle = 180;
        std::vector<double> oxygen_angles( 3 );
        for ( int i = 0; i < 3; i++ )
        {
            if ( fabs( oxygen_vectors[ i ][2] ) > distance_error )
                std::cout<<"\nError in Find_Angles: oxygen_vectors[ i ][2]="<<oxygen_vectors[ i ][2];
            oxygen_angles[ i ] = Angle_Signed( Point2d(1,0), Point2d( oxygen_vectors[ i ][0], oxygen_vectors[ i ][1] ) );
            //std::cout<<"\nray"<<i<<"="<<Point2d( oxygen_rays[ i ][0], oxygen_rays[ i ][1] )<<" angle="<<oxygen_angles[ i ];
            if ( min_abs_angle > fabs( oxygen_angles[ i ] ) ) { min_abs_angle = fabs( oxygen_angles[ i ] ); min_index = i; }
        }
        if ( min_index < 0 )
        {
            std::cout<<"\nError in Find_Angles: min_index="<<min_index<<" C_axis="<<carbon_axis; //<<" rot_axis="<<a<<" mat="<<rotation;
            for ( int i = 0; i < 3; i++ ) std::cout<<" ray="<<atoms['O'][i].point_a - centre<<"->"<<oxygen_rays[ i ];
        }
        else oxygen_angle = oxygen_angles[ min_index ];
    }

    bool Find_indices_NH()
    {
        indices_NH.assign( atoms['N'].size(), -2 );
        for ( int i = 0; i < atoms['N'].size(); i++ )
        {
            for ( int j = 0; j < atoms['H'].size(); j++ )
                if ( indices_NH[i] < 0 ) // H neighbour not found yet
                {
                    double distance_NH = norm( atoms['N'][i].point_a - atoms['H'][j].point_a );
                    if ( fabs( distance_NH - length_NH ) < distance_error ) indices_NH[i] = j; //atoms['H'][j].index - 1;
                }
            if ( indices_NH[i] < 0 ) { std::cout<<"\nError in Find_indices_NH: no H neighbour of N"<<atoms['N'][i].index; return false; }
        }
        return true;
    }
};

Point3d Cross_Product (Point3d v0, Point3d v1)
{
    Point3d v;
    v.x = v0.y * v1.z - v1.y * v0.z;
    v.y = v0.z * v1.y - v1.z * v0.y;
    v.z = v0.x * v1.y - v1.x * v0.y;
    return v;
}

bool Find_Rotation (Point3d v0, Point3d v1, Matrix& rotation)  //What is the purpose of calculation of rotation matrix w.r.t verties 
{
    double angle = Angle_Positive( v0, v1 );
    rotation = Eigen::MatrixXd::Identity(3,3);
    if ( angle == 0 ) return true; // identity rotation
    double c = cos( angle * M_PI / 180 );
    double s = sin( angle * M_PI / 180 );
    rotation *= c;
    if ( angle == M_PI ) return true; 
    Point3d axis = Cross_Product( v0, v1 );
    axis *= 1.0 / norm( axis );
    rotation += s * Cross_Product( axis ) + (1-c) * Tensor_Product( axis ); //std::cout<<"\nr="<<rotation;
    return true;
}

class Edge
{
public:
    int first, second; // endpoint indices
    double length = 0;
    Point3d centre;
    Edge (double l) { length = l; }
    Edge (Point3d c) { centre = c; }
};

typedef boost::property<boost::edge_weight_t, Edge> Edge_Property;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Molecule> Graph;
typedef boost::graph_traits<Graph>::vertex_iterator Vertex_it;
typedef boost::property_map<Graph, boost::vertex_index_t>::type Graph_Map;
typedef boost::graph_traits < Graph >::adjacency_iterator Adjacency_it;

class Vertex_Point
{
public:
    Adjacency_it iter;
    Point3d point;
    Vertex_Point () {}
    Vertex_Point (Adjacency_it it, Point3d p) { iter = it; point = p; }
};
bool Increasing_Points (Vertex_Point const& v1, Vertex_Point const& v2)
{
    if ( v1.point.x < v2.point.x ) return true;
    if ( v1.point.x > v2.point.x ) return false;
    if ( v1.point.y < v2.point.y ) return true;
    if ( v1.point.y > v2.point.y ) return false;
    if ( v1.point.z < v2.point.z ) return true;
    return false;
}

bool Read_Until_Section (std::fstream& file, int section)
{
    std::string line, col1 = "", col2 = "";
    while ( col1 != "#" or col2 != std::to_string( section ) + "." )
    {
        if ( ! std::getline( file, line ) ) return false;
        std::istringstream stream( line );
        stream >> col1 >> col2;

    }
    return true;
}

bool Read_Until_String (std::fstream& file, std::string s)
{
    std::string line, column = "";
    while ( column != s )
    {
        if ( ! std::getline( file, line ) ) return false;
        std::istringstream stream( line );
        stream >> column;
    }
    return true;
}

bool Read_String_Value (std::fstream& file, std::string s, double& v)
{
    std::string line, col1 = "", col2;
    while ( col1 != s )
    {
        if ( ! std::getline( file, line ) ) return false;
        std::istringstream stream( line );
        stream >> col1 >> col2;
    }
    v = atof( col2.c_str() );
    return true;
}

bool Read_String_Value (std::fstream& file, std::string s, int& v)
{
    double value;
    Read_String_Value( file, s, value );
    v = int( value );
    return true;
}

bool Read_2columns (std::fstream& file, std::string& col1, std::string& col2)
{
    std::string line;
    col1 = ""; col2 = "";
    std::getline( file, line );
    std::istringstream stream( line );
    stream >> col1 >> col2;
    return true;
}

bool Read_idxyz (std::fstream& file, std::string& id, double& x, double& y, double& z)
{
    std::string line, col2, col3, col4, col5, col6;
    id = "";
    if ( !getline( file, line ) ) return false;
    std::istringstream stream( line );
    stream >> id >> col2 >> col3 >> col4;
    if ( id == "" ) return false;
    x = atof( col4.c_str() );
    stream >> col5;
    y = atof( col5.c_str() );
    stream >> col6;
    z = atof( col6.c_str() );
    return true;
}


bool Read_cif (std::string name, Box& box, std::vector<Molecule>& molecules)
{
    std::fstream file;
    file.open( name.c_str(), std::ios::in );
    if ( !file.is_open() )
    {
        std::cout<<"\nFile "<<name<<" not found";
        return false;
    }
    std::string line, id, col1, col2, col3, col4;
    Atom atom;
    char element;
    int ind_molecule = 0, ind_atom = 0;
    double x,y,z;
    //bool box_found = false, box_saved = false, atoms_found = false, atoms_saved = false;
    // initialise molecule parameters
    std::map<char,int> indices, max_indices;
    max_indices.insert( std::make_pair( 'O', 3 ) );
    max_indices.insert( std::make_pair( 'N', 6 ) );
    max_indices.insert( std::make_pair( 'C', 23 ) );
    max_indices.insert( std::make_pair( 'H', 14 ) );
    for ( auto i : max_indices ) indices.insert( std::make_pair( i.first, 0 ) );

    // Read box parameters
    Read_Until_Section( file, 6 );
    Read_String_Value( file, "_cell_length_a", box.a );
    Read_String_Value( file, "_cell_length_b", box.b );
    Read_String_Value( file, "_cell_length_c", box.c );
    Read_String_Value( file, "_cell_angle_alpha", box.alpha );
    Read_String_Value( file, "_cell_angle_beta", box.beta );
    Read_String_Value( file, "_cell_angle_gamma", box.gamma );
    box.Find_Matrix();
    Read_Until_String( file, "_atom_site_refinement_flags" );
    // Read x,y,z coordinates of all atoms
    while ( Read_idxyz( file, id, x, y, z ) )
    {
        element = id[0];
        ind_molecule = int( indices[ element ] / max_indices[ element ] );
        if ( ind_molecule > molecules.size() )
        { std::cout<<"\nError in Read_cif: i_mol="<<ind_molecule<<" mol_size="<<molecules.size(); return false; }
        if ( ind_molecule == molecules.size() ) molecules.push_back( Molecule() );
        ind_atom = ++indices[ element ];
        Atom atom( element, ind_atom, x, y ,z );
        atom.point_a = box.Abs_Position( atom.point_b );
        molecules[ ind_molecule ].atoms[ element ].push_back( atom );
        if ( element == 'C' )
        {
            molecules[ ind_molecule ].centre += atom.point_b;
            int ind = ind_atom - ind_molecule * max_indices[ element ];
            if ( ind == 1 or ind == 14 ) molecules[ ind_molecule ].c += atom.point_b;
        }
    }
    box.n_molecules = (int)molecules.size();
    box.Print();
    return true;
}

/*
bool Linked_ON (Molecule& m0, Molecule& m1, double distance_ON)
{
    bool linked = false;
    for ( auto O : m0.atoms['O'] )
        for ( auto N : m1.atoms['N'] )
        {
            //if ( fabs( norm( O.point_a - N.point_a ) - distance_ON ) < 1e-2 ) { return true; }
            if ( norm( O.point_a - N.point_a ) < distance_ON ) { std::cout<<" l="<<norm( O.point_a - N.point_a ); linked = true; }
            //else std::cout<<" ON="<<norm( O.point_a - N.point_a );
        }
    for ( auto O : m1.atoms['O'] )
        for ( auto N : m0.atoms['N'] )
        {
            //std::cout<<" ON="<<norm( O.point_a - N.point_a );
            if ( norm( O.point_a - N.point_a ) < distance_ON ) { std::cout<<" l="<<norm( O.point_a - N.point_a ); linked = true;  }
            //else std::cout<<" ON="<<norm( O.point_a - N.point_a );
        }
    return linked;
}*/

bool Linked_NHO (Molecule& m0, Molecule& m1, double& gap_HO)  // I could not understand the operations within this function.
{
    gap_HO = big_constant;
    bool linked = false;
    for ( int i = 0; i < m0.atoms['N'].size(); i++ )
    {
        int index_H = m0.indices_NH[ i ];
        for ( auto j = 0; j < m1.atoms['O'].size(); j++ )
        {
            Point3d vector_HO = m1.atoms['O'][j].point_a - m0.atoms['H'][ index_H ].point_a;
            double distance_HO = norm( vector_HO );
            if ( distance_HO > length_HO )
            {
                if ( gap_HO > distance_HO ) gap_HO = distance_HO;
                continue;
            }
            //std::cout<<"\nN"<<m0.atoms['N'][i].index<<"H"<<index_H+1<<"O"<<m1.atoms['O'][j].index<<"="<<distance_HO;
            double distance_ON = norm( m0.atoms['N'][i].point_a - m1.atoms['O'][j].point_a );
            if ( distance_ON > length_ON ) continue; //std::cout<<"\nPotential Error: distance_ON="<<distance_ON;
            double angle_NHO = Angle_Positive( m0.atoms['N'][i].point_a - m0.atoms['H'][ index_H ].point_a, vector_HO );
            if ( angle_NHO < 100 ) continue; //std::cout<<"\nPotential Error: angle_NHO="<<angle_NHO;
            linked = true;
            break;
        }
    }
    return linked;
}

void Order_Neighbours (Graph& s, Vertex_it v, std::vector<Vertex_Point>& c)
{
    c.clear();
    for( auto neighbours = boost::adjacent_vertices( *v, s ); neighbours.first != neighbours.second; ++neighbours.first )
        c.push_back( Vertex_Point( neighbours.first, s[ *neighbours.first ].centre ) );
    sort( c.begin(), c.end(), Increasing_Points );
    std::cout<<"\nSorted:"; for ( auto n : c ) std::cout<<" "<<n.point;
}

bool Neighbours_Equal (std::vector<Vertex_Point>const& neighbours0, std::vector<Vertex_Point>const& neighbours1) // Is it returns equal neighbours within same molecule structure 
{
    if ( neighbours0.size() != neighbours1.size() ) { std::cout<<"\nDifferent numbers of neighbours"; return false; }
    for ( int i = 0; i < neighbours0.size(); i++ )
        if ( norm( neighbours0[i].point - neighbours1[i].point ) > distance_error )
        {
            std::cout<<"\nDifferent: "<<neighbours0[i].point<<"!="<<neighbours1[i].point;
            return false;
        }
    std::cout<<"\nEqual neighbours";
    return true;
}

bool Structures_Equal (Graph& s0, Vertex_it v0, Graph& s1, Vertex_it v1) // How it checks equalities of structure through this condition
{
    if ( norm( s0[ *v0 ].centre - s1[ *v1 ].centre ) > distance_error )   
    {
        std::cout<<"\nDifferent: "<<s0[ *v0 ].centre<<s1[ *v1 ].centre; // initial vertices
        return false;
    }
    std::cout<<"\nEqual initial vertices: "<<s0[ *v0 ].centre<<s1[ *v1 ].centre;
    std::vector<Vertex_Point> neighbours0, neighbours1;
    Order_Neighbours( s0, v0, neighbours0 );
    Order_Neighbours( s1, v1, neighbours1 );
    if ( ! Neighbours_Equal( neighbours0, neighbours1 ) ) return false;

    for ( int i = 0; i < neighbours0.size(); i++ )
    {

    }
    return true;
}

int main()
{
    std::string data_folder = "E:/Research/molGeom/";
    Box box;
    std::vector<Molecule> molecules;
    int num_structures = 2; //5688
    int max_neighbours = 6;
    std::vector<Graph> structures( num_structures );
    for ( int ind_structure = 0; ind_structure < num_structures; ind_structure++ )
    {
        molecules.clear();
        std::cout<<"\n"<<ind_structure + 1;
        Read_cif( data_folder + "T2_" + std::to_string( ind_structure + 1 ) + "_num_molGeom.cif", box, molecules ); // T2_1_num_molGeom
        // Find centres and angles in the orthogonal system
        for ( int i = 0; i < molecules.size(); i++ )
        {
            molecules[i].internal = true;
            molecules[i].centre *= 1.0 / molecules[i].atoms['C'].size();
            molecules[i].centre = box.Abs_Position( molecules[i].centre );
            molecules[i].Find_Angles();
            molecules[i].Find_indices_NH();
        }

        // Add molecules by 26 shifts
        std::vector<Point3d> shifts;  // why this vector have their value in range of -1,0,1 for shift?
        for ( int i = -1; i <=1; i++ )
            for ( int j = -1; j <=1; j++ )
                for ( int k = -1; k <=1; k++ )
                    if ( i!=0 or j!=0 or k!=0)
                        shifts.push_back( Point3d( i, j, k) );
        std::vector<Molecule> molecules_new;
        for ( auto s : shifts)
        {
            for ( auto m : molecules )
                molecules_new.push_back( Molecule( m, box.Abs_Position( s ) ) );
        }
        molecules.insert( molecules.end(), molecules_new.begin(), molecules_new.end() );

        // Build a structure of molecules from 26 boxes
        std::vector<Graph::vertex_descriptor> vertices( molecules.size() );
        for ( int i = 0; i < molecules.size(); i++ )
        {
            //molecules[i].Print(); std::cout<<" m"<<i;
            vertices[i] = boost::add_vertex( structures[ ind_structure ] ); // vertices indexed as a vector
            structures[ ind_structure ][ vertices[i] ] = molecules[i];
        }
        //
        std::cout<<" n=";
        Point3d c0, c1, c;
        std::vector<Index_Value> gaps_HO;   //What is the use of gaps_HO vector?
        for ( int i = 0; i < box.n_molecules; i++ )
        {
            gaps_HO.clear();
            int neighbours = 0;
            for ( int j = 0; j < molecules.size(); j++ )
            {
                if ( i == j ) continue;
                c0 = structures[ ind_structure ][ vertices[i] ].centre;   
                c1 = structures[ ind_structure ][ vertices[j] ].centre;
                c = c1 - c0;
                double d = norm( c ), gap_HO;
                if ( ! Linked_NHO( molecules[i], molecules[j], gap_HO ) ) { gaps_HO.push_back( Index_Value( j, gap_HO ) ); continue; }
                neighbours++;
                if ( ! boost::edge( vertices[i], vertices[j], structures[ ind_structure ] ).second )
                    add_edge( vertices[i], vertices[j], structures[ ind_structure ] ); // Edge( c ),
            }
            std::cout<<neighbours<<",";
            if ( neighbours < max_neighbours ) //std::cout<<"\nError in molecule"<<i<<" neighbours="<<neighbours;
            {
                sort( gaps_HO.begin(), gaps_HO.end(), Increasing_Values );
                //std::cout<<"\nextra:";
                for ( int k = 0; k < max_neighbours - neighbours and k < gaps_HO.size(); k++ )
                {
                    c = structures[ ind_structure ][ vertices[ gaps_HO[ k ].index ] ].centre - c0;
                    //std::cout<<" gap="<<gaps_HO[ k ].value;
                    add_edge( vertices[i], vertices[ gaps_HO[ k ].index ], structures[ ind_structure ] ); //  Edge( c ),
                }
            }

        }
    }

    /* Double lopop over structures
    Matrix rotation_carbon, rotation_oxygen, rotation;
    //for ( int i = 0; i < num_structures; i++ )
    for ( int i = 0; i < 1; i++ )
        for ( int j = i+1; j < 2; j++ ) // Two structures fixed
        {
            // Double loop over internal molecules
            bool structures_equal = false;
            for ( auto v0 = boost::vertices( structures[ i ] ).first; v0 != boost::vertices( structures[ i ] ).second; v0++ )
                if ( structures[ i ][*v0].internal and ! structures_equal )
                    for ( auto v1 = boost::vertices( structures[ j ] ).first; v1 != boost::vertices( structures[ j ] ).second; v1++ )
                        if ( structures[ j ][*v1].internal and ! structures_equal ) // Two internal molecules fixed
                        {
                            // Transform the structure j molecule v1 to the structure i molecule v0
                            std::cout<<"\nCompare ";
                            structures[ i ][*v0].Print();
                            structures[ j ][*v1].Print();
                            //Point3d shift = structures[ i ][*v0].centre - structures[ j ][*v1].centre;
                            Find_Rotation( structures[ i ][*v0].carbon_axis, structures[ j ][*v1].carbon_axis, rotation_carbon );
                            Find_Rotation( structures[ i ][*v0].oxygen_rays[0], structures[ j ][*v1].oxygen_rays[0], rotation_oxygen );

                            // Build the new structure applying the found transformation to all structure j molecules
                            Graph structure;
                            boost::copy_graph( structures[ j ], structure );
                            Vertex_it v = boost::vertices( structure ).first, vc;
                            for ( ; v != boost::vertices( structure ).second; v++ )
                            {
                                //std::cout<<"\nold="<<structure[*v].centre;
                                structure[*v].centre -= structures[ j ][*v1].centre; // v1.centre at the origin
                                if ( structure[*v].centre == origin ) vc = v;
                                structure[*v].centre = V ( rotation_oxygen * rotation_carbon * V ( structure[*v].centre ) );
                                structure[*v].centre += structures[ i ][*v0].centre;
                                //std::cout<<" new="<<structure[*v].centre;

                            }
                            // Try 3 rotations and 1 symmetry: 6 transformations of the new structure with structure i

                            if ( Structures_Equal( structures[ i ], v0, structure, vc ) ) structures_equal = true;



                        }


        }
    */

    std::cout<<"\n";
    return 0;
}
