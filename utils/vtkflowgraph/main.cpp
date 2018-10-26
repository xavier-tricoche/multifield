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

int ts_count = 0;
int ts_str = 0;
int ts_end = 0;
int fs_count = 0;
char** meshes_files;

class Node
{ 
public:
	vtkIdType id;
	int timestep;
	vector<vtkIdType> vlist;
	float* tag;
	
	Node() 
	{
		// initialize the tag array
		tag = (float*) malloc(fs_count * sizeof(float));
		for (int i = 0; i < fs_count; i++)
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
		
		for (int i = 0; i < fs_count; i++)
		{
			if ((t & 1) == 1)
				tag[i] += 1.0;
				
			t = t >> 1;
		}
	}
	
	void AverageTag()
	{
		// average on the number of points
		for (int i = 0; i < fs_count; i++)
			tag[i] /= vlist.size();
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////
////// Definitions of graphs types
////////////////////////////////////////////////////////////////////////////////////////////////////

// (1) directed graph
typedef boost::adjacency_list
  <
    //Store all edges as a std::vector
    boost::setS,
    //Store all vertices in a std::vector
    boost::setS,
    //Relations are both ways (in this example)
    //(note: but you can freely change it to boost::directedS)
    boost::directedS,
    //All vertices are person names of type std::string
    Node,
    //All edges are weights equal to the encounter frequencies
    boost::property<boost::edge_weight_t,double>,
    //Graph itself has a std::string name
    boost::property<boost::graph_name_t,std::string>
  > graph_t;

typedef boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
typedef boost::graph_traits<graph_t>::edge_descriptor edge_t;

// (2) undirected graph
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
  > ugraph_t;

typedef boost::graph_traits<ugraph_t>::vertex_descriptor uvertex_t;
typedef boost::graph_traits<ugraph_t>::edge_descriptor uedge_t;

// (3) graph to be written in file
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, 
	  Node,
	  boost::property< boost::edge_weight_t, double >,
	  boost::property<boost::graph_name_t,std::string>
> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
////////////////////////////////////////////////////////////////////////////////////////////////////
////// Simple graph of nodes for meshes vertices
////////////////////////////////////////////////////////////////////////////////////////////////////
vector<vector<uvertex_t> > smap;

void CreateDenseGraph(ugraph_t& g)
{
	vtkPolyData* m;
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader;
	
	// add nodes
	cout << "Create nodes for the vertices.\n";
	for (int t = ts_str; t <= ts_end; t++)
	{
		vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
		vtkMeshReader->SetFileName(meshes_files[t]);
		vtkMeshReader->Update();
		m = vtkMeshReader->GetOutput();
		
		// create a node for each vertex in the mesh
		vector<uvertex_t> tmp;
		for (vtkIdType i = 0; i < m->GetNumberOfPoints(); i++)
		{
			uvertex_t v = boost::add_vertex(g);
			g[v].id = i;
			g[v].timestep = t;
			tmp.push_back(v);
		}
		smap.push_back(tmp);
	}
	
	// add links
	cout << "Create the links in the vertices graph.\n";
	vtkPointData* pointsData = m->GetPointData();
	vtkIntArray* barr = (vtkIntArray*) pointsData->GetArray("MATCHTBL0");
	vtkIntArray* farr = (vtkIntArray*) pointsData->GetArray("MATCHTBL1");
	if (barr == NULL) printf("Table not found!\n");
	if (farr == NULL) printf("Table not found!\n");
	for (int t = ts_str; t <= ts_end; t++)
	{
		vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
		vtkMeshReader->SetFileName(meshes_files[t]);
		vtkMeshReader->Update();
		m = vtkMeshReader->GetOutput();
		
		// create a node for each vertex in the mesh
		for (vtkIdType i = 0; i < m->GetNumberOfPoints(); i++)
		{
			vtkIdType b = barr->GetValue(i);
			vtkIdType f = farr->GetValue(i);
		
			if (b != -1)
				boost::add_edge(smap[t-1][b], smap[t][i], 1.0, g);
			if (f != -1)
				boost::add_edge(smap[t][i], smap[t+1][f], 1.0, g);
		}
	}
	
	// now iterate fixing the weights
	cout << "Iterate fixing the weights in the vertices graph.\n";
	boost::property_map<ugraph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, g);
	bool cont = true;
	while (cont)
	{
		cont = false;
		for (int t = ts_str; t <= ts_end; t++)
		{
			vector<uvertex_t>* l = &(smap[t]);
			for (vector<uvertex_t>::iterator it = l->begin(); it != l->end(); it++)
			{
				uvertex_t v = *it;
				
				double sumi = 0.0;
				double sumo = 0.0;
				
				// loop on all edges
				boost::graph_traits<ugraph_t>::in_edge_iterator ei, ei_end;
				for (boost::tie(ei, ei_end) = in_edges(v, g); ei != ei_end; ++ei) 
				{
					uvertex_t v1 = source(*ei, g);
					uvertex_t v2 = target(*ei, g);
					
					if (g[v1].timestep < t)
					{
						sumi += get(weightmap, *ei); 
					}
					else if (g[v1].timestep > t)
					{
						sumo += get(weightmap, *ei); 
					}
					
					if (g[v2].timestep < t)
					{
						sumi += get(weightmap, *ei); 
					}
					else if (g[v2].timestep > t)
					{
						sumo += get(weightmap, *ei); 
					}
				}
				
				double si = max(sumi, sumo) / sumi;
				double so = max(sumi, sumo) / sumo;
				if (sumi == sumo)
					continue;
					
				cont = true;
				for (boost::tie(ei, ei_end) = in_edges(v, g); ei != ei_end; ++ei) 
				{
					uvertex_t v1 = source(*ei, g);
					uvertex_t v2 = target(*ei, g);
					
					int mint = min(g[v1].timestep, g[v2].timestep);
					int maxt = max(g[v1].timestep, g[v2].timestep);
					
					if (mint < t)
					{
						put(weightmap, *ei, get(weightmap, *ei) * si);
					}
					else if (maxt > t)
					{
						put(weightmap, *ei, get(weightmap, *ei) * so);
					}					
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////// Now create the big graph
////////////////////////////////////////////////////////////////////////////////////////////////////
map<int, map<int, vertex_t> > tsid2vertex;

void AddNodes(int ts, graph_t& g, vtkPolyData* m)
{
	vtkPointData* pointsData = m->GetPointData();
	vtkIntArray* tagarr = (vtkIntArray*) pointsData->GetArray("TAGTBL");
	vtkIntArray* segarr = (vtkIntArray*) pointsData->GetArray("SEGTBL");
	if (tagarr == NULL) printf("Table not found!\n");
	if (segarr == NULL) printf("Table not found!\n");
	
	for (vtkIdType i = 0; i < m->GetNumberOfPoints(); i++)
	{
		int seg = segarr->GetValue(i);
		map<int, vertex_t>* lids = &(tsid2vertex[ts]);
		map<int, vertex_t>::iterator it = lids->find(seg);
		vertex_t v;
		if (it == lids->end())
		{
			v = boost::add_vertex(g);
			g[v].id = seg;
			g[v].timestep = ts;
			tsid2vertex[ts][seg] = v;
		}
		else
		{
			v = it->second;
		}
		
		g[v].vlist.push_back(i);
		g[v].UpdateTag(tagarr, i);
	}
}

void CreateGraph(graph_t& g, ugraph_t& ug)
{
	vtkPolyData* m1;
	vtkPolyData* m2;
	
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader;
	vtkPointData* pointsData1;
	vtkPointData* pointsData2;
	vtkIntArray* arr1;
	vtkIntArray* arr2;
	vtkIntArray* arr;
	boost::property_map<graph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, g);
	
	cout << "Create the segments graph.\n";
	for (int t = ts_str; t < ts_end; t++)
	{
		cout << "time step " << t << "\n";
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
		
		// add nodes for meshes
		if (t == 0) AddNodes(t, g, m1);
		AddNodes(t + 1, g, m2);
	
		// get the segments arrays
		pointsData1 = m1->GetPointData();
		arr1 = (vtkIntArray*) pointsData1->GetArray("SEGTBL");
		if (arr1 == NULL) printf("Table not found!\n");
		pointsData2 = m2->GetPointData();
		arr2 = (vtkIntArray*) pointsData2->GetArray("SEGTBL");
		if (arr2 == NULL) printf("Table not found!\n");
		
		// loop on edges from t to t+1
		boost::property_map<ugraph_t, boost::edge_weight_t>::type uweightmap = get(boost::edge_weight, ug);
		boost::property_map<graph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, g);
		boost::graph_traits<ugraph_t>::edge_iterator ei, ei_end;
		for (boost::tie(ei, ei_end) = edges(ug); ei != ei_end; ++ei) 
		{
			uvertex_t uv1 = source(*ei, ug);
			uvertex_t uv2 = target(*ei, ug);
			
			int mint = min(ug[uv1].timestep, ug[uv2].timestep);
			int maxt = max(ug[uv1].timestep, ug[uv2].timestep);
			
			if ((mint != t) || (maxt != t+1))
				continue;
				
			if (ug[uv1].timestep != t)
			{
				uvertex_t tmp;
				tmp = uv1;
				uv1 = uv2;
				uv2 = tmp;
			}
				
			// points ids
			int id1 = ug[uv1].id;
			int id2 = ug[uv2].id;
			
			// segment ids
			id1 = arr1->GetValue(id1);
			id2 = arr2->GetValue(id2);
			
			// vertices
			vertex_t v1 = tsid2vertex[t][id1];
			vertex_t v2 = tsid2vertex[t+1][id2];
			
			edge_t e12;
			bool found;
			tie(e12, found) = edge(v1, v2, g);
			if (found)
			{
				put(weightmap, e12, get(weightmap, e12) + get(uweightmap, *ei));
			}
			else
			{
				boost::add_edge(v1, v2, get(uweightmap, *ei), g);
			}
		}
	}
	
	boost::graph_traits<graph_t>::vertex_iterator vi, vend;
	for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi) 
	{
		g[*vi].AverageTag();
	}
}

struct my_node_writer {
	// my_node_writer() {}
	my_node_writer(Graph* g_) : g (g_) {};
	void operator()(std::ostream& out, const Vertex& v) {
		out << "[label=\"" << (*g)[v].timestep << "/" << (*g)[v].id << "\"]";
	};
	// bleah.  I can't get references right...
	// according to http://www.knowledgesearch.org/doc/examples.html
	// it should be a reference here, not the type itself.
	// but g++ either barfs, or the program segfaults.
	Graph* g;
};
  
int main (int argc, char *argv[])
{
	if (argc != 5)
	{
		printf("Example: vtkflowgraph meshes.txt 50 87 network.txt\n");
		return 0;
	}

	cout << "============\n";
	cout << "vtkflowgraph\n";
	cout << "============\n";
	
	ts_str = atoi(argv[2]);
	ts_end = atoi(argv[3]);
	
	// read the file for information
	fstream input(argv[1], std::ios::in);
	input >> ts_count;
	input >> fs_count;
	meshes_files = (char**) malloc(ts_count * sizeof(char*));
	for (int i = 0; i < ts_count; i++)
	{
		string buffer;
		getline(input, buffer);
		meshes_files[i] = new char[buffer.size()+1];
		strcpy(meshes_files[i], buffer.c_str());
	}
	input.close();
	
	// graphs used
	graph_t g;
	ugraph_t ug;
	
	// create vertex graph
	CreateDenseGraph(ug);
	
	// create segments graph
	CreateGraph(g, ug);

	// initialize Graph
	Graph ng;

	// copy the vertices
	map<vtkIdType, Vertex> g2ng;
	boost::graph_traits<graph_t>::vertex_iterator vi, vend;
	for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi) 
	{
		Vertex v;
		v = boost::add_vertex(ng);
		ng[v].id = g[*vi].id;
		ng[v].timestep = g[*vi].timestep;
		g2ng[ng[v].id] = v;
	}
	
	// copy the edges
	boost::property_map<graph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, g);
	boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
	for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) 
	{
		vertex_t gv1 = source(*ei, g);
		vertex_t gv2 = target(*ei, g);
		Vertex v1 = g2ng[g[gv1].id];
		Vertex v2 = g2ng[g[gv2].id];
		boost::add_edge(v1, v2, get(weightmap, *ei), ng);
	}
		
	// write the graph to a file
	std::ofstream output(argv[4], std::ios::out);
	boost::write_graphviz(output, ng, my_node_writer(&ng));
	output.close();

	
	return 0;
}