/*************************************************************************
multifield:   Surface-based Structural Analysis and Visualization 
              of Multifield Datasets

Author: Samer S. Barakat

Copyright (c) 2010-2012, Purdue University

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
**************************************************************************/
#include "MyCGALSimplifyMesh.h"

void MyCGALSimplifyMesh(char* input_file, char* output_file, double perc_reduction) 
{
	
	
  /*Surface surface; 
  
  printf("Reading file %s\n", infile);
  
  std::ifstream is(infile) ; is >> surface ;

  // This is a stop predicate (defines when the algorithm terminates).
  // In this example, the simplification stops when the number of undirected edges
  // left in the surface drops below the specified number (1000)
  double rat = perc;
  printf("Input ratio %f\n", rat);
  //SMS::Count_ratio_stop_predicate<Surface> stop(rat);
  SMS::Count_stop_predicate<Surface> stop(1000);
     
  //internal::cgal_enable_ecms_trace = true ; 
  
  // This the actual call to the simplification algorithm.
  // The surface and stop conditions are mandatory arguments.
  // The index maps are needed because the vertices and edges
  // of this surface lack an "id()" field.
  int r = SMS::edge_collapse
            (surface
            ,stop
            ,CGAL::vertex_index_map(boost::get(CGAL::vertex_external_index,surface)) 
                  .edge_index_map  (boost::get(CGAL::edge_external_index  ,surface)) 
            );
  
  std::cout << "\nFinished...\n" << r << " edges removed.\n" 
            << (surface.size_of_halfedges()/2) << " final edges.\n" ;
        
  std::ofstream os( outfile ) ; os << surface ;*/
  
  return ;      
}