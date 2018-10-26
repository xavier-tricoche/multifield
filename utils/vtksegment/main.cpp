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
#include <queue>
#include <map>
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
#include <vtkPolyDataReader.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

using namespace std;

struct Node
{ 
	int id;
	vector<vtkIdType> vlist;
	float* tag;
};

bool operator<(const Node& i, const Node& j) { return (i.vlist.size() < j.vlist.size());}

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


bool HasSameTag(vtkIntArray* arr, const vtkIdType& p1, const vtkIdType& p2)
{
	if (arr->GetValue(p1) != arr->GetValue(p2))
		return false;
	
	return true;
}

double FindTagDist(float** tags, int numv, vtkIdType p1, vtkIdType p2)
{
	double sum = 0.0;
	for (int i = 0; i < numv; i++)
	{
		double dif = tags[p1][i] - tags[p2][i];
		sum += dif * dif;
	}
	
	return sqrt(sum);
}

double FindTagDist(Node& n1, Node& n2, const int& numv)
{
	double sum = 0.0;
	for (int i = 0; i < numv; i++)
	{
		double dif = n1.tag[i] - n2.tag[i];
		sum += dif * dif;
	}
	
	return sqrt(sum);
}

void ComputeTag(float** tags, const int& numv, Node& n)
{
	n.tag = (float*) malloc(numv * sizeof(float));
	for (int i = 0; i < numv; i++)
		n.tag[i] = 0.0;
		
	for (vector<vtkIdType>::iterator it = n.vlist.begin(); it != n.vlist.end(); it++)
	{
		for (int i = 0; i < numv; i++)
		{
			n.tag[i] += tags[*it][i];
		}
	}
	
	for (int i = 0; i < numv; i++)
		n.tag[i] /= n.vlist.size();
}

float* TagFromInt(int numv, int tag)
{
	float* ptr = (float*) malloc(numv * sizeof(float));
	for (int i = 0; i < numv; i++)
	{
		if (tag&1 == 1)
			ptr[i] = 1;
		else
			ptr[i] = 0;
			
		tag = tag >> 1;
	}
}

int main (int argc, char *argv[])
{
	if (argc != 4)
	{
		printf("Example: vtksegment 10 mesh.vtk 0.2\n");
		return 0;
	}
	
	cout << "==========\n";
	cout << "vtksegment\n";
	cout << "==========\n";

	printf("Mesh file is %s\n", argv[2]);
	int numv = atoi(argv[1]);
	
	// read the mesh file
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader->SetFileName(argv[2]);
    vtkMeshReader->Update();
	vtkPolyData* vtkMesh = vtkMeshReader->GetOutput();
	
	// start component analysis
	vtkPoints* points = vtkMesh->GetPoints();
	char* visited = (char*) malloc(points->GetNumberOfPoints() * sizeof(char));
	memset(visited, 0, points->GetNumberOfPoints() * sizeof(char));
	char* msk = (char*) malloc(points->GetNumberOfPoints() * sizeof(char));
	memset(msk, 1, points->GetNumberOfPoints() * sizeof(char));
	
	// get data array 
	vtkPointData* pointsData = vtkMesh->GetPointData();
	vtkIntArray* arr = (vtkIntArray*) pointsData->GetArray("TAGTBL");
	if (arr == NULL) printf("Table not found!\n");
	
	// iterate on components
	char buffer[100];
	graph_t g;
	vector<vertex_t> nodes;
	int totalcomp = 0;
	int* assi = (int*) malloc(points->GetNumberOfPoints() * sizeof(int));
	printf("The number of points is %d\n", points->GetNumberOfPoints());
	printf("The number of polys is %d\n", vtkMesh->GetNumberOfPolys());
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	{	
		// if any point was visited then the whole component was
		if (visited[i]) continue;
		
		// read the complete component
		float maxv = numeric_limits<float>::min();
		vector<vtkIdType> comp;
		queue<vtkIdType> cp;
		cp.push(i);
		comp.push_back(i);
		visited[i] = 1;
		assi[i] = totalcomp;
		while (!cp.empty())
		{
			// get point
			vtkIdType p = cp.front();
			cp.pop();
			
			// add non visited neihbors
			vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
			vtkMesh->GetPointCells(p, cells);
			for (int k = 0; k < cells->GetNumberOfIds(); k++)
			{
				vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
				vtkMesh->GetCellPoints(cells->GetId(k), ptIds);
				for (int j = 0; j < ptIds->GetNumberOfIds(); j++)
				{
					int id = ptIds->GetId(j);
					if ((!visited[id]) && (HasSameTag(arr, id, p)))
					{
						cp.push(id);
						comp.push_back(id);
						visited[id] = 1;
						assi[id] = totalcomp;
					}
				}
			}
		}
		
		// create the vertex in the graph
		vertex_t v = boost::add_vertex(g);
		g[v].id = totalcomp++;
		g[v].vlist = comp; 
		g[v].tag = TagFromInt(numv, arr->GetValue(comp.front()));
		nodes.push_back(v);
	}
	printf("Components are %d\n", totalcomp);
	
	// loop over all points and create links between the different components
	for (vtkIdType p = 0; p < points->GetNumberOfPoints(); p++)
	{
		vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
		vtkMesh->GetPointCells(p, cells);
		for (int k = 0; k < cells->GetNumberOfIds(); k++)
		{
			vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
			vtkMesh->GetCellPoints(cells->GetId(k), ptIds);
			for (int j = 0; j < ptIds->GetNumberOfIds(); j++)
			{
				int id = ptIds->GetId(j);
				if (assi[p] != assi[id])
				{
					// parallel edges are not allowed in the graph
					boost::add_edge(nodes[assi[p]], nodes[assi[id]], FindTagDist(tags, numv, p, id), g);
				}
			}
		}
	}
	
	// for evry small component find the closest neighbor that it can be merged with given radius
/*	double limit = atof(argv[3]);
	while (1)
	{	
		// get edge with minimum weight
		edge_t mwe;
		double minw = numeric_limits<double>::max();
		boost::property_map<graph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, g);
		boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
		for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) 
		{
			double w = get(weightmap, *ei);
			if (w < minw)
			{
				minw = w;
				mwe = *ei;
			}
		}
		if (minw > limit) break;
		vertex_t v1 = source(mwe, g);
		vertex_t v2 = target(mwe, g);

		// find all neighbors of v1 and v2
		vector<vertex_t> nbrs;
		boost::graph_traits<graph_t>::adjacency_iterator neighbourIt, neighbourEnd;
		for (boost::tie(neighbourIt, neighbourEnd) = adjacent_vertices(v1, g); neighbourIt != neighbourEnd; ++neighbourIt) 
			nbrs.push_back(*neighbourIt);
		for (boost::tie(neighbourIt, neighbourEnd) = adjacent_vertices(v2, g); neighbourIt != neighbourEnd; ++neighbourIt) 
			nbrs.push_back(*neighbourIt);
		
		// add new vertex
		vertex_t v = boost::add_vertex(g);
		g[v].id = totalcomp++;
		g[v].vlist.insert(g[v].vlist.end(), g[v1].vlist.begin(), g[v1].vlist.end()); 
		g[v].vlist.insert(g[v].vlist.end(), g[v2].vlist.begin(), g[v2].vlist.end()); 
		ComputeTag(tags, numv, g[v]);
		
		// add new edges
		for (vector<vertex_t>::iterator it = nbrs.begin(); it != nbrs.end(); it++)
		{
			if ((*it != v1) && (*it != v2))
				boost::add_edge(v, *it, FindTagDist(g[v], g[*it], numv), g);
		}
		
		// delete v1 and v2
		remove_vertex(v1, g);
		remove_vertex(v2, g);
	}
	
	// remove reprtitions
	//hmm ?
	
	// create the new table and the segment file
	int ids = 0;
	vtkSmartPointer<vtkIntArray> newScalars = vtkSmartPointer<vtkIntArray>::New();
	newScalars->SetName("SEGTBL");
	boost::graph_traits<graph_t>::vertex_iterator vi, vend;
	for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi) 
	{
		Node n = g[*vi];
		for (vector<vtkIdType>::iterator it = n.vlist.begin(); it != n.vlist.end(); it++)
		{
			newScalars->InsertValue(*it, ids);
		}
		ids++;
	}
	pointsData->AddArray(newScalars);
	printf("Number of segments is %d\n", ids);
	
	// write file
	vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkMeshWriter->SetFileName(argv[2]);
	vtkMeshWriter->SetInput(vtkMesh);
	vtkMeshWriter->Write();
	
	// free memory
	free(visited);
	free(msk);

	*/
	return 0;
}