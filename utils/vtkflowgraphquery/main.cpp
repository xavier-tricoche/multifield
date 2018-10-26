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
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <bitset>
#include <iostream>
#include <fstream>
#include <sstream>
#include <teem/nrrd.h>
#include <vector_types.h>
#include <vector_functions.h>
#include <cutil_inline.h>
#include <cutil_math.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

using namespace std;


int ts_numb = 0;
int fs_numb = 0;
char** meshes_files;

class Node
{ 
public:
	int id;
	int timestep;
	vector<vtkIdType> vlist;
	float* tag;
	
	Node() 
	{
		// initialize the tag array
		tag = (float*) malloc(fs_numb * sizeof(float));
		for (int i = 0; i < fs_numb; i++)
			tag[i] = 0.0;
	}
	
	bool operator<(const Node& node) const
	{
		if (this->timestep < node.timestep)
			return true;
			
		if ((this->timestep == node.timestep) && (this->id < node.id))
			return true;
			
		return false;
	};

	bool operator==(const Node& node) const
	{
		if ((this->timestep == node.timestep) && (this->id == node.id))
			return true;
		
		return false;
	}
	
	void UpdateTag(vtkIntArray* arr, vtkIdType id)
	{
		// loop on points summing the tag values
		int t = arr->GetValue(id);
		
		for (int i = 0; i < fs_numb; i++)
		{
			if ((t & 1) == 1)
				tag[i] += 1.0;
				
			t = t >> 1;
		}
	}
	
	void AverageTag()
	{
		// average on the number of points
		for (int i = 0; i < fs_numb; i++)
			tag[i] /= vlist.size();
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////
////// Definitions of graphs types
////////////////////////////////////////////////////////////////////////////////////////////////////
typedef boost::adjacency_list
  <
    //Store all edges as a std::vector
    boost::setS,
    //Store all vertices in a std::vector
    boost::setS,
    //Relations are both ways (in this example)
    //(note: but you can freely change it to boost::directedS)
    boost::undirectedS,
    //All vertices are person names of type std::string
    Node,
    //All edges are weights equal to the encounter frequencies
    boost::property<boost::edge_weight_t,double>,
    //Graph itself has a std::string name
    boost::property<boost::graph_name_t,std::string>
  > graph_t;

typedef boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
typedef boost::graph_traits<graph_t>::edge_descriptor edge_t;

////////////////////////////////////////////////////////////////////////////////////////////////////
////// Simple graph of nodes for meshes vertices
////////////////////////////////////////////////////////////////////////////////////////////////////

map<int, map<int, vertex_t> > tsid2vertex;

vertex_t GetNode(int ts, int id, graph_t& g)
{
	map<int, vertex_t>* lids = &(tsid2vertex[ts]);
	map<int, vertex_t>::iterator it = lids->find(id);
	vertex_t v;
	if (it == lids->end())
	{
		v = boost::add_vertex(g);
		g[v].id = id;
		g[v].timestep = ts;
		tsid2vertex[ts][id] = v;
	}
	else
	{
		v = it->second;
	}
	
	return v;
}

void FindVertices(vector<set<int>> steps, char** files)
{
	vtkPolyData* m1;
	vtkPolyData* m2;
	
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader;
	vtkPointData* pointsData1;
	vtkPointData* pointsData2;
	vtkIntArray* arr1;
	vtkIntArray* arr2;
	vtkIntArray* arr;
	vtkIntArray* m1t2;
	vtkIntArray* m2t1;
	boost::property_map<graph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, g);
	
	set<int> S1;
	set<int> S2;
	
	// create the graph
	for (int t = 0; t < steps->size(); t++)
	{
		// read mesh 1
		if (t == 0) 
		{
			vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
			vtkMeshReader->SetFileName(meshes_files[t]);
			vtkMeshReader->Update();
			m1 = vtkMeshReader->GetOutput();
		}
		else
			m1 = m2;
		
		// read mesh 2
		vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
		vtkMeshReader->SetFileName(meshes_files[t + 1]);
		vtkMeshReader->Update();
		m2 = vtkMeshReader->GetOutput();
		
		// get the segments arrays
		pointsData1 = m1->GetPointData();
		arr1 = (vtkIntArray*) pointsData1->GetArray("SEGTBL");
		if (arr1 == NULL) printf("Table not found!\n");
		pointsData2 = m2->GetPointData();
		arr2 = (vtkIntArray*) pointsData2->GetArray("SEGTBL");
		if (arr2 == NULL) printf("Table not found!\n");
		
		// matching tables
		m1t2 = (vtkIntArray*) pointsData1->GetArray("MATCHTBL1");
		if (m1t2 == NULL) printf("Table not found!\n");
		m2t1 = (vtkIntArray*) pointsData2->GetArray("MATCHTBL0");		
		if (m2t1 == NULL) printf("Table not found!\n");
		
		// loop on S1 finding matches in m2
		for (set<int>::iterator it = S1.begin(); it != S1.end(); it++)
		{
			vtkIdType vid1 = *it;
			vtkIdType vid2 = m1t2->GetValue(vid1);
			
			if (vid2 == -1) continue;
			
			int sid1 = arr1->GetValue(vid1);
			int sid2 = arr2->GetValue(vid2);
			
			if (!steps[t+1].contains(sid2)) continue;
			
			vertex_t v1 = GetNode(t, vid1, g);
			vertex_t v2 = GetNode(t+1, vid2, g);
			
			boost::add_edge(v1, v1, 1.0, g);
			S2.insert(vid2);
		}
		
		// loop on vertices of m2 finding the match in S1
		for (vtkIdType i = 0; i < m2->GetNumberOfPoints(); i++)
		{
			vtkIdType vid2 = i;
			vtkIdType vid1 = m2t1->GetValue(vid2);
			
			if (!S1.contains(vid1)) continue;
			
			int sid1 = arr1->GetValue(vid1);
			int sid2 = arr2->GetValue(vid2);
			
			if (!steps[t+1].contains(sid2)) continue;
			
			vertex_t v1 = GetNode(t, vid1, g);
			vertex_t v2 = GetNode(t+1, vid2, g);
			
			boost::add_edge(v1, v1, 1.0, g);
			S2.insert(vid2);			
		}
		
		// swicth sets
		S1.clear();
		S1 = S2;
		S2.clear();
	}
	
	
	// delete the bad vertices
	set<vertex_t> del;
	bool cont = true;
	while (cont)
	{
		cont = false;
		
		// loop on all vertices		
		boost::graph_traits<graph_t>::vertex_iterator vi, vend;
		for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi) 
		{
			// last step do not need to be matched
			if (g[*vi].timestep == steps.size() - 1)
				continue;
			
			bool found = false;
			// loop on all edges of each vertex
			boost::graph_traits<graph_t>::in_edge_iterator ei, ei_end;
			for (boost::tie(ei, ei_end) = in_edges(*vi, g); ei != ei_end; ++ei) 
			{
				vertex_t v1 = source(*ei, g);
				vertex_t v2 = target(*ei, g);
				
				if ((g[v1].timestep == g[*v].timestep + 1) || (g[v2].timestep == g[*v].timestep + 1))
				{
					found = true;
					break;
				}
			}
			
			if (!found)
				del.insert(*vi);
		}
		
		// remove the nodes
		for (set<vertex_t>::iterator it = del.begin(); it != del.end(); it++)
			remove_vertex(*del, g);
		
		if (del.size() > 0)
		{
			del.clear();
			cont = true;
		}
	}
	
}

int main (int argc, char *argv[])
{
	if (argc != 3)
	{
		printf("Example: vtkflowgraphquery meshes.txt network.txt\n");
		return 0;
	}

	
	
	return 0;
}