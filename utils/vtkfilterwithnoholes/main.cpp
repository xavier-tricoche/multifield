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

int main (int argc, char *argv[])
{
	printf("====================\n");
	printf("vtkfilterwithnoholes\n");
	printf("====================\n");
	
	if (argc != 6)
	{
		printf("Example: vtkfilterwithnoholes mesh.vtk out.vtk TABLE -1 130\n");
		return 0;
	}

	printf("Mesh file is %s\n", argv[1]);
	
	// read the mesh file
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader->SetFileName(argv[1]);
    vtkMeshReader->Update();
	vtkPolyData* vtkMesh = vtkMeshReader->GetOutput();
	
	// start component analysis
	vtkPoints* points = vtkMesh->GetPoints();
	char* msk = (char*) malloc(points->GetNumberOfPoints() * sizeof(char));
	memset(msk, 1, points->GetNumberOfPoints() * sizeof(char));
	
	// get data array 
	vtkPointData* pointsData = vtkMesh->GetPointData();
	vtkFloatArray* arr = (vtkFloatArray*) pointsData->GetArray(argv[3]);
	if (arr == NULL) printf("Table not found!\n");
	float limit = atof(argv[4]);
	
	// mark points above ridges
	int removed = 0;
	printf("The number of points is %d\n", points->GetNumberOfPoints());
	printf("The number of polys is %d\n", vtkMesh->GetNumberOfPolys());
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	{	
		if (arr->GetValue(i) > limit)
		{
			removed++;
			msk[i] = 0;
		}		
	}
	printf("Removed %d points\n", removed);
	
	// find all cells boundary edges
	vtkCellArray* vtkCells = vtkMesh->GetPolys();
	vtkCells->InitTraversal();
	vtkIdType npts;
	vtkIdType* pts;
	multimap<vtkIdType,vtkIdType> mymultimap;
	int idx = 0;
	vtkIdType ids[3];
	int count = 0;
	vector<int2> edges;
	while (vtkCells->GetNextCell(npts, pts) != 0)
	{
		idx = 0;
		if (msk[pts[0]] != 0) ids[idx++] = pts[0];
		if (msk[pts[1]] != 0) ids[idx++] = pts[1];
		if (msk[pts[2]] != 0) ids[idx++] = pts[2];
		
		if (idx == 2)
		{
			edges.insert(edges.begin(), make_int2(ids[0],ids[1]));
		}
	}
	
	for (vector<int2>::iterator it = edges.begin(); it < edges.end(); it++ )
	{
		if (HasValidCell(vtkMesh, msk, (*it).x, (*it).y))
		{
			count++;
			mymultimap.insert(pair<vtkIdType,vtkIdType>((*it).x, (*it).y));
			mymultimap.insert(pair<vtkIdType,vtkIdType>((*it).y, (*it).x));
		}
	}
	
	cout << "Number of edges is " << count << "\n";

	// set values
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	{
		arr->SetValue(i, 0.0);
	}
	
	// find all closed rings
	int ncycles = 0;
	multimap<vtkIdType,vtkIdType>::iterator it;
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
		if ((cycle.size() > 1) && (cycle.front() == cycle.back()))
		{
			ncycles++;
			for (vector<vtkIdType>::iterator it=cycle.begin() ; it < cycle.end(); it++ )
			{
				cout << " " << *it;
				arr->SetValue(*it, -500);
			}
			printf("\n");
			
			FillCycle(vtkMesh, msk, &cycle);
		}
	}
	cout << "Number of cycles is " << ncycles << "\n";
	
	// new mesh points
	vtkSmartPointer<vtkPolyData> vtkMeshNew = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> pointsNew = vtkSmartPointer<vtkPoints>::New();
	map<vtkIdType, vtkIdType> old2new;
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	{
		if (msk[i] == 1)
		{
			double* point = points->GetPoint(i);
			vtkIdType newid = pointsNew->InsertNextPoint(point[0], point[1], point[2]);
			old2new.insert(pair<vtkIdType, vtkIdType>(i, newid));
		}
	}
	vtkMeshNew->SetPoints(pointsNew);
	
	// new mesh polys
	vtkSmartPointer<vtkCellArray> cellsNew = vtkSmartPointer<vtkCellArray>::New();
	vtkCells->InitTraversal();
	while (vtkCells->GetNextCell(npts, pts) != 0)
	{	
		bool skip = false;
		for (int j = 0; j < npts; j++)
		{
			if (old2new.find(pts[j]) == old2new.end())
			{
				skip = true;
				break;
			}
			
			pts[j] = old2new.find(pts[j])->second;
		}
		
		if (npts != 3) 
		{
			printf(".");
			continue;
		}
		
		if (skip == false)
			cellsNew->InsertNextCell(npts, pts);
	}
	vtkMeshNew->SetPolys(cellsNew);
	
	// tmp for values
	vtkPointData* pD = vtkMeshNew->GetPointData();
	vtkSmartPointer<vtkFloatArray> newScalars = vtkSmartPointer<vtkFloatArray>::New();
	newScalars->SetName(argv[3]);
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	{
		if (old2new.find(i) == old2new.end()) continue;
		int j = old2new.find(i)->second;
		newScalars->InsertValue(j, arr->GetValue(i));
	}
	pD->AddArray(newScalars);
	
	// write file
	vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkMeshWriter->SetFileName(argv[2]);
	vtkMeshWriter->SetInput(vtkMeshNew);
	vtkMeshWriter->Write();
	
	// free memory
	free(msk);

	
	return 0;
}