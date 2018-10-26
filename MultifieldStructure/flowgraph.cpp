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
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <teem/nrrd.h>
#include <vector_types.h>
#include <vector_functions.h>
#include <cutil_inline.h>
#include <cutil_math.h>
#include <math/fixed_vector.hpp>
#include <math/bezier_spline.hpp>
#include <boost/config.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/read_dimacs.hpp>
#include <boost/graph/graph_utility.hpp>
#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

#include "FlowGrid.h"

using namespace std;

typedef unsigned int uint;
typedef unsigned char uchar;

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_d;
typedef K::Vector_3 Vector_d;
typedef CGAL::Search_traits_3<K> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;

double INTSTEP;

class Trajectory
{
public:
	static int ids;
	
	int id;	
	vector<float4> points;
	int status;
	
	Trajectory()
	{
		id = ids++;
	}
	
	void AddPoint(int t, float x, float y, float z, int s)
	{
		points.insert(points.end(), make_float4(x,y,z, t));
		status = s;
	}
	
	bool HasPointAtTime(int t)
	{
		int diff = t - points[0].w;
		if (diff < points.size())
			return true;
		return false;
	}
	
	float3 GetPointAtTime(int t)
	{
		int diff = t - points[0].w;
		float4 pt = points[diff];
		return make_float3(pt.x, pt.y, pt.z);
	}
	
	friend bool operator== (Trajectory &cP1, Trajectory &cP2)
	{
		if (cP1.id == cP2.id)
			return true;
		return false;
	}
};

int Trajectory::ids = 0;

class Node
{
public:
	static int ids;
	
	int id;
	int ts;
	uchar* type;
	int vn;
	set<Trajectory*> particles;
	
	Node(int _ts, int _vn)
	{
		id = ids++;
		ts = _ts;
		vn = _vn;
		type = (uchar*) malloc(vn);
		memset(type, 0, vn);
	}
	
	~Node()
	{
		free(type);
	}
	
	void AddParticle(Trajectory* pid)
	{
		particles.insert(particles.end(), pid);
	}
	
	friend bool operator== (Node &cP1, Node &cP2)
	{
		for (int i = 0; i < cP1.vn; i++)
		{
			if (cP1.type[i] != cP2.type[i])
				return false;
		}
		if (cP1.ts != cP2.ts)
			return false;
		return true;
	}
};

int Node::ids = 0;

Nrrd *readNrrd(const char* filename)
{
    Nrrd *nin = nrrdNew();
    if (nrrdLoad(nin, filename, NULL)) {
        printf("Error reading %s\n",filename);
                return NULL;
    }
    return nin;
}

vtkPolyData* loadMesh(string filename)
{
	vtkPolyDataReader* pdr = vtkPolyDataReader::New();
	pdr->SetFileName(filename.c_str());
	pdr->Update();
	vtkPolyData* ptr = pdr->GetOutput();
	pdr->Delete();
	return ptr;
}

void loadMeshes(string dp, int vn, int t, vtkPolyData** meshes)
{
	char ts[10];
	sprintf(ts, "%d", t);
	char buf[256];
	for (int i = 0; i < vn; i++)
	{
		sprintf(buf, "%d", i);
		meshes[i * 2] = loadMesh(dp + "/mesh_" + string(ts) + "_" + string(buf) + "_0.vtk");
		meshes[i * 2 + 1] = loadMesh(dp + "/mesh_" + string(ts) + "_" + string(buf) + "_1.vtk");
	}
}

void assignParticleToNode(int vn, int t, Trajectory* particle, vtkPolyData** meshes, set<Node*>* nodes, set<int>* marks)
{
	// if does not have a location at t exit
	if (!particle->HasPointAtTime(t))
		return ;
		
	// find node in graph
	float3 pt = particle->GetPointAtTime(t);
	double ptd[3]; ptd[0] = pt.x; ptd[1] = pt.y; ptd[2] = pt.z;
	Node* node = new Node(t, vn);
	for (int i = 0; i < 2 * vn; i++)
	{
		vtkIdType id = meshes[i]->FindPoint(ptd);
		double ptm[3];
		meshes[i]->GetPoint(id, ptm);
		
		ptm[0] -= ptd[0];
		ptm[1] -= ptd[1];
		ptm[2] -= ptd[2];
		
		if (sqrt(ptm[0]*ptm[0] + ptm[1]*ptm[1] + ptm[2]*ptm[2]) < 5)
		{
			node->type[i / 2] |= i % 2 + 1;				
			marks[i].insert(id);
		}
	}
	
	// check node exist
	if (nodes->find(node) != nodes->end())
	{
		delete node;
		node = *(nodes->find(node));
	}
	else
	{
		nodes->insert(node);
	}
	
	node->particles.insert(particle);
}

void findNewParticles(int vn, int t, vtkPolyData** meshes, set<Trajectory*>* particles, set<Node*>* nodes, set<int>* marks)
{
	Tree tree;
	
	// build the search tree
	for (set<Trajectory*>::iterator it=particles->begin(); it!=particles->end(); it++)
	{
		Trajectory* particle = *it;
		if (!particle->HasPointAtTime(t))
			continue;
		
		float3 pt = particle->GetPointAtTime(t);
		Point_d np(pt.x, pt.y, pt.z);
		tree.insert(np);
	}
	
	// loop on all meshes
	for (int i = 0; i < 2 * vn; i++)
	{
		vtkPolyData* mesh = meshes[i];
		
		// loop on vertices
		for (int j = 0; j < mesh->GetNumberOfPoints(); j++)
		{
			// ignore marked points
			if (marks[i].find(j) != marks[i].end()) continue;
			
			// get the point coordinates
			double pt[3];
			mesh->GetPoint(j, pt);
			
			// check point distance from tree
			double mind = numeric_limits<double>::max();
			Point_d query_point(pt[0], pt[1], pt[2]);
			Neighbor_search result(tree, query_point, 1);
			for(Neighbor_search::iterator it = result.begin(); it != result.end(); ++it)
			{
				mind = min(mind, it->second);
			}
			
			// if point is far enough add it
			if (sqrt(mind) < 5)
			{
				// add new trajectory for point
				Trajectory* particle = new Trajectory();
				particle->AddPoint(t, pt[0], pt[1], pt[2], 0);
				particles->insert(particle);
				
				// make node assignments
				assignParticleToNode(vn, t, particle, meshes, nodes, marks);
				
				// insert in the distance tree
				Point_d np(pt[0], pt[1], pt[2]);
				tree.insert(np);
			}
		}		
	}
}

float3 getFlow(double tcur, float3 y, FlowGrid* flow_0, FlowGrid* flow_1, bool& failed)
{
	float3 V0 = flow_0->GetFlowAt(y, failed);
	float3 V1 = flow_1->GetFlowAt(y, failed);
	return TimeInterpolate(flow_0->time, flow_1->time, V0, V1, tcur);
}

inline float CopySign(float x, float y) {
	if (y >= 0)
		return fabsf(x);
	else
		return fabsf(x) * -1;
}

int Integrate(float3 &y, float tstrt, float tout, FlowGrid* flow_0, FlowGrid* flow_1)
{
	volatile float h = INTSTEP;
	float tcur = tstrt;
	float tmp = tout - tcur;

	// Set stepsize for integration in the direction from t to tout
	h = CopySign(h, tmp);

	// Step by step integration - as an endless loop over steps
	while(1) {
		// Adjust stepsize if necessary to hit the output point.
		// Look ahead two steps to avoid drastic changes in the stepsize and
		// thus lessen the impact of output points on the code.
		tmp = tout - tcur;
		if (fabsf(tmp) <= fabsf(h))
			h = tmp;

		// Advance an approximate solution over one step of length h
		/***********************************************************************/
		float3 v1, v2, v3, v4;
		bool failed = false;

		v1 = h * getFlow(tcur, y, flow_0, flow_1, failed);
		v2 = h * getFlow(tcur + 0.5f * h, y + 0.5f * v1, flow_0, flow_1, failed);
		y = y + v2;
		/***********************************************************************/

		if (fabsf(tmp) <= fabsf(h)) {
			// integration complete
			return 0;
		}
		
		if (failed) {
			// left volume boundaries
			return 1;
		}

		tcur += h;
	}
}

void integrateParticles(int vn, int t, FlowGrid* flow_0, FlowGrid* flow_1, set<Trajectory*>* particles, set<Trajectory*>* particles_done, set<Node*>* nodes, vtkPolyData** meshes, set<int>* marks)
{	
	// integrate the particles
	for (set<Trajectory*>::iterator it=particles->begin(); it!=particles->end(); ++it)
	{
		Trajectory* particle = *it;
		
		// find point
		if (!particle->HasPointAtTime(t)) 
		{
			printf("This should never happen.\n");
			continue;
		}
		float3 y = particle->GetPointAtTime(t);
		
		// integrate the point
		int status = Integrate(y, flow_0->time, flow_1->time, flow_0, flow_1);
		particle->AddPoint(t+1, y.x, y.y, y.z, status);
	}
	
	// assign particles to nodes
	for (set<Trajectory*>::iterator it=particles->begin(); it!=particles->end(); ++it)
	{
		Trajectory* particle = *it;
		
		// find point
		if (!particle->HasPointAtTime(t)) 
		{
			printf("This should never happen.\n");
			continue;
		}
		
		// check the status
		if (particle->status == 1)
		{
			particles->erase(particle);
			particles_done->insert(particle);
		}
		
		// make node assignments
		assignParticleToNode(vn, t+1, particle, meshes, nodes, marks);
	}
}

int main(int argc, char *argv[])
{
	string dp(argv[1]);
	int vn = atoi(argv[2]);
	int sn = atoi(argv[3]);
	double td = atof(argv[4]);
	INTSTEP = td / 50;
	
	printf("Directory: %s\n", dp.c_str());
	printf("Number of variables: %d\n", vn);
	printf("Number of time steps: %d\n", sn);
	printf("Time step: %lf\n", td);
	printf("Integration step: %lf\n", INTSTEP);
	
	// all particles 
	set<Trajectory*> particles_done;
	set<Trajectory*> particles;
	
	// graph nodes
	set<Node*> nodes;
	
	// nrrd files for the flow
	Nrrd* nrrdflow_0 = readNrrd((dp + "/flow_0.nrrd").c_str());
	Nrrd* nrrdflow_1 = readNrrd((dp + "/flow_0.nrrd").c_str());
	
	// vtk mesh objects
	vtkPolyData** meshes_0 = (vtkPolyData**) malloc(2 * vn * sizeof(vtkPolyData*));
	vtkPolyData** meshes_1 = (vtkPolyData**) malloc(2 * vn * sizeof(vtkPolyData*));
	loadMeshes(dp, vn, 0, meshes_0);
	loadMeshes(dp, vn, 0, meshes_1);
	
	// marks for meshes' triangles
	set<int>* marks = (set<int>*) malloc(2 * vn * sizeof(set<int>)); 
	
	// loop on the time steps
	char tstr[1024]; 
	for (int t = 0; t < sn; t++)
	{
		// load flow t and t+1
		sprintf(tstr, "%d", t+1);
		nrrdNuke(nrrdflow_0);
		nrrdflow_0 = nrrdflow_1;
		nrrdflow_1 = readNrrd((dp + "/flow_" + string(tstr) + ".nrrd").c_str());
		FlowGrid* flow_0 = new FlowGrid(nrrdflow_0, t * td);
		FlowGrid* flow_1 = new FlowGrid(nrrdflow_1, (t+1) * td);
		
		// meshes t <- meshes t+1
		for (int i = 0; i < 2 * vn; i++)
		{
			meshes_0[i] = meshes_1[i];
		}
		
		// load meshes for t+1
		loadMeshes(dp, vn, t+1, meshes_1);
		
		// integrate all particles and mark vertices
		integrateParticles(vn, t, flow_0, flow_1, &particles, &particles_done, &nodes, meshes_1, marks);
		
		// find particles from mesh t
		findNewParticles(vn, t+1, meshes_1, &particles, &nodes, marks);
		
		// clear marks
		for (int i = 0; i < 2 * vn; i++)
			marks[i].clear();
	}
	
	// free memory
	free(meshes_0);
	free(meshes_1);
	
	return 0;
}