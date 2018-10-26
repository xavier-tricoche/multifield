/*#define CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE 6
void Surface_simplification_external_trace( std::string s )
{
   static std::ofstream out("log.txt");
   out << s << std::endl ;
} */

/*#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>*/

//typedef CGAL::Simple_cartesian<double> Kernel;
//typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Surface; 

//namespace SMS = CGAL::Surface_mesh_simplification ;

void MyCGALSimplifyMesh(char* input_file, char* output_file, double perc_reduction);

// EOF //