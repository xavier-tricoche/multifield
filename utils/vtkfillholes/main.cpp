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
#include <set>
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
#include <vtkFillHolesFilter.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include "triangulate.h"

using namespace std;

bool HasValidCell(vtkPolyData* vtkMesh, char* msk, vtkIdType id1, vtkIdType id2)
{
	vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
	vtkMesh->GetPointCells(id1, cells);
	for (int k = 0; k < cells->GetNumberOfIds(); k++)
	{
		// get cell information
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
		vtkMesh->GetCellPoints(cells->GetId(k), ptIds);
		
		vtkIdType i1 = ptIds->GetId(0);
		vtkIdType i2 = ptIds->GetId(1);
		vtkIdType i3 = ptIds->GetId(2);
		
		// not shared triangle
		if ((id2 != i1) && (id2 != i2) && (id2 != i3)) continue;
		
		// not valid
		if ((msk[i1] == 0) || (msk[i2] == 0) || (msk[i3] == 0)) continue;
		
		return true;
	}
	
	return false;
}

void FillCycle(vtkPolyData* vtkMesh, char* msk, vector<vtkIdType>* vert)
{
	int readded = 0;
	
	queue<vtkIdType> cp;
	set<vtkIdType> cycle;
	for (vector<vtkIdType>::iterator it = vert->begin() ; it < vert->end(); it++)
	{
		cycle.insert(cycle.end(), *it);
		cp.push(*it);
	}
	
	// read the complete component
	while (!cp.empty())
	{
		// get point
		vtkIdType p = cp.front();
		cp.pop();
		
		// find all cells of the point
		vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
		vtkMesh->GetPointCells(p, cells);
		for (int k = 0; k < cells->GetNumberOfIds(); k++)
		{
			// get cell information
			vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
			vtkMesh->GetCellPoints(cells->GetId(k), ptIds);
			
			// 2 vertices should be in the cylce list
			int incycle = 0;
			for (int j = 0; j < ptIds->GetNumberOfIds(); j++)
			{
				vtkIdType id = ptIds->GetId(j);
				if (cycle.find(id) != cycle.end())
					incycle++;
			}
			if (incycle != 2) continue;
			
			// element removed should be readded
			for (int j = 0; j < ptIds->GetNumberOfIds(); j++)
			{
				vtkIdType id = ptIds->GetId(j);
				if (msk[id] == 0)
				{
					cp.push(id);
					msk[id] = 1;
					cycle.insert(cycle.end(), id);
					readded++;
					
					if (readded > 300) return;
				}
			}
		}
	}
	
	printf("Readded is %d points\n", readded);
}

bool operator<(const int2 e1, const int2 e2)
{
	if (e1.x < e2.x)
	{
		return true;
	}
	else if (e1.x > e2.x)
	{
		return false;
	}
	else
	{
		if (e1.y < e2.y)
		{
			return true;
		}
		else if (e1.y > e2.y)
		{
			return false;
		}
		else
		{
			return false;
		}
	}
}

bool operator==(const int2 e1, const int2 e2)
{
	if ((e1.x == e2.x) && (e1.y == e2.y))
		return true;
	return false;
}

bool TesselateCell(vtkPoints* points, vector<vtkIdType>* cycle, vector<int3>* tris)
{
	// first find the best fit plane
	double A[9];
	for (vector<vtkIdType>::iterator it = cycle->begin(); it != cycle->end(); it++)
	{
		double* point = points->GetPoint(*it);
		
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				A[j + 3 * i] += point[i] * point[j];	
	}
	
	double evalr[3];
	double evecr[3][3];

	gsl_matrix_view m = gsl_matrix_view_array (A, 3, 3);
	gsl_vector *eval = gsl_vector_alloc (3);
	gsl_matrix *evec = gsl_matrix_alloc (3, 3);
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
	gsl_eigen_symmv (&m.matrix, eval, evec, w);
	gsl_eigen_symmv_free (w);
	gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);

	evalr[0] = gsl_vector_get (eval, 0);
	evalr[1] = gsl_vector_get (eval, 1);
	evalr[2] = gsl_vector_get (eval, 2);
	
	gsl_vector_view evec_0 = gsl_matrix_column (evec, 0);
	gsl_vector_view evec_1 = gsl_matrix_column (evec, 1);
	gsl_vector_view evec_2 = gsl_matrix_column (evec, 2);

	evecr[0][0] = gsl_vector_get (&evec_0.vector, 0);
	evecr[0][1] = gsl_vector_get (&evec_0.vector, 1);
	evecr[0][2] = gsl_vector_get (&evec_0.vector, 2);

	evecr[1][0] = gsl_vector_get (&evec_1.vector, 0);
	evecr[1][1] = gsl_vector_get (&evec_1.vector, 1);
	evecr[1][2] = gsl_vector_get (&evec_1.vector, 2);
	
	evecr[2][0] = gsl_vector_get (&evec_2.vector, 0);
	evecr[2][1] = gsl_vector_get (&evec_2.vector, 1);
	evecr[2][2] = gsl_vector_get (&evec_2.vector, 2);
     
	gsl_vector_free (eval);
	gsl_matrix_free (evec);
	
	float3 PLn = normalize(make_float3(evecr[0][0], evecr[0][1], evecr[0][2]));
	float3 PLV0 = make_float3(0.0,0.0,0.0);
	
	// PLn
	float3 vert[3];
	int nvert = 0;
	int jump = floor(cycle->size()/3.0);
	for (vector<vtkIdType>::iterator it = cycle->begin(); it != cycle->end(); it+=jump)
	{
		double* point = points->GetPoint(*it);
		
		vert[nvert++] = make_float3(point[0], point[1], point[2]);
		if (nvert == 3) break;
	}
	PLn = normalize(cross(vert[1] - vert[0], vert[2] - vert[0]));
	
	// axis of the 2D plane
	float sb, sn, sd;
	float3 P = make_float3(1, 1, 1);		
	sn = -dot( PLn, (P - PLV0));
	sd = dot(PLn, PLn);
	sb = sn / sd;
	P = P + sb * PLn; // arbitrary point on the plane
	float3 axis_0 = normalize(P - PLV0);
	float3 axis_1 = normalize(cross(axis_0, PLn));
		
	// project points to a plane
	//cout << "Query points " << cycle->size() << " :\n";
	Vector2dVector a;
	for (vector<vtkIdType>::iterator it = cycle->begin(); it != cycle->end(); it++)
	{
		double* point = points->GetPoint(*it);
			
		P = make_float3(point[0], point[1], point[2]);
		sn = -dot( PLn, (P - PLV0));
		sd = dot(PLn, PLn);
		sb = sn / sd;
		P = P + sb * PLn;
		
		a.push_back( Vector2d(dot(axis_0, P), dot(axis_1, P)));
		//printf("%f %f %f -> %f %f ,", (float) point[0], (float) point[1], (float) point[2], dot(axis_0, P), dot(axis_1, P));
	}
	//printf("\n");
	
	// Triangulate in 2D
	vector<int3> result;
	bool ret = Triangulate::Process(a, result);
	//printf("Process succeded %s\n", ret? "true":"false");
	int tcount = result.size();
	//printf("Solution %d :\n", tcount);
	for (int i=0; i<tcount; i++)
	{
		//printf(" %d %d %d, ", (*cycle)[result[i].x], (*cycle)[result[i].y], (*cycle)[result[i].z]);
		tris->push_back(make_int3((*cycle)[result[i].x], (*cycle)[result[i].y], (*cycle)[result[i].z]));
	}
	//printf("\n");
	
	return ret;
}

int main (int argc, char *argv[])
{
	printf("============\n");
	printf("vtkfillholes\n");
	printf("============\n");
	
	if (argc != 4)
	{
		printf("Example: vtkfillholes mesh.vtk out.vtk 500\n");
		return 0;
	}

	printf("Mesh file is %s\n", argv[1]);	
	
	// read the mesh file
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader->SetFileName(argv[1]);
    vtkMeshReader->Update();
	vtkPolyData* vtkMesh = vtkMeshReader->GetOutput();
	int limit = atoi(argv[3]);
	
	// get data array 
	vtkPointData* pointsData = vtkMesh->GetPointData();	
	vtkSmartPointer<vtkFloatArray> holesScalars = vtkSmartPointer<vtkFloatArray>::New();
	holesScalars->SetName("HOLESTBL");
	
	// set values
	vtkPoints* points = vtkMesh->GetPoints();
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	{
		holesScalars->InsertValue(i, 0.0);
	}
	
	// find all cells boundary edges
	vtkCellArray* vtkCells = vtkMesh->GetPolys();
	vtkCells->InitTraversal();
	vtkIdType npts;
	vtkIdType* pts;
	multimap<vtkIdType,vtkIdType> mymultimap;
	int count = 0;
	map<int2, int> edges;
	while (vtkCells->GetNextCell(npts, pts) != 0)
	{
		int2 edge1 = make_int2(min(pts[0], pts[1]), max(pts[0], pts[1]));
		int2 edge2 = make_int2(min(pts[0], pts[2]), max(pts[0], pts[2]));
		int2 edge3 = make_int2(min(pts[1], pts[2]), max(pts[1], pts[2]));
		
		map<int2, int>::iterator it;
		
		if ((it = edges.find(edge1)) != edges.end())
		{
			int scr = it->second + 1;
			edges.erase(it);
			edges.insert(pair<int2, int>(edge1, scr));
		}
		else
			edges.insert(pair<int2, int>(edge1, 1));
			
		if ((it = edges.find(edge2)) != edges.end())
		{	
			int scr = it->second + 1;
			edges.erase(it);
			edges.insert(pair<int2, int>(edge2, scr));
		}
		else
			edges.insert(pair<int2, int>(edge2, 1));
			
		if ((it = edges.find(edge3)) != edges.end())
		{	
			int scr = it->second + 1;
			edges.erase(it);
			edges.insert(pair<int2, int>(edge3, scr));
		}
		else
			edges.insert(pair<int2, int>(edge3, 1));
			
	}
	
	for (map<int2, int>::iterator it = edges.begin(); it != edges.end(); it++ )
	{
		if (it->second == 1)
		{
			count++;
			mymultimap.insert(pair<vtkIdType,vtkIdType>((it->first).x, (it->first).y));
			mymultimap.insert(pair<vtkIdType,vtkIdType>((it->first).y, (it->first).x));
		}
		else if (it->second != 2)
		{
			cout << "edge " << (it->first).x << "," << (it->first).y << " has " << it->second << " triangles!\n";
		}
	}
	
	cout << "Number of edges is " << count << "\n";
	
	// find all closed rings
	int ncycles = 0;
	multimap<vtkIdType,vtkIdType>::iterator it;
	vector<int3> tris;
	int tess_s = 0;
	int tess_f = 0;
	while (!mymultimap.empty())
	{
		it = mymultimap.begin();
		vector<vtkIdType> cycle;
		
		// find chain
		vtkIdType f = (*it).first;	
		vtkIdType s = (*it).second;
		cycle.push_back(f);	
		cycle.push_back(s);	
		mymultimap.erase(it);		
		while ((it = mymultimap.find(s)) != mymultimap.end())
		{
			if (((*it).second) == f)
			{
				mymultimap.erase(it);
				continue;
			}
				
			f = (*it).first;			
			s = (*it).second;
			cycle.push_back(s);
						
			mymultimap.erase(it);
		}
		
		// check cycle
		if ((cycle.size() > 0) && (cycle.size() < limit) && (cycle.front() == cycle.back()))
		{
			ncycles++;
			//printf("Cycle number %d :\n", ncycles);
			for (vector<vtkIdType>::iterator it=cycle.begin() ; it < cycle.end(); it++ )
			{
				//cout << " " << *it;
				holesScalars->SetValue(*it, -500);
			}
			//printf("\n");
			
			if (TesselateCell(points, &cycle, &tris))
				tess_s++;
			else
				tess_f++;
		}
	}
	cout << "Number of cycles is " << ncycles << "\n";
	cout << "Success rate is " << tess_s << "/" << (tess_s+tess_f) << "\n";
	
	// insert new cells
	pts = (vtkIdType*) malloc(3 * sizeof(vtkIdType));
	for(vector<int3>::iterator it = tris.begin(); it != tris.end(); it++)
	{
		int3 tri = *it;
		
		pts[0] = tri.x;
		pts[1] = tri.y;
		pts[2] = tri.z;
		
		vtkMesh->GetPolys()->InsertNextCell(3, pts);
	}	
	free(pts);

	// tmp for values
	pointsData->AddArray(holesScalars);
	
	// write file
	vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkMeshWriter->SetFileName(argv[2]);
	vtkMeshWriter->SetInput(vtkMesh);
	vtkMeshWriter->Write();
	
	return 0;
}