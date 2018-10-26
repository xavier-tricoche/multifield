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

using namespace std;

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

int main (int argc, char *argv[])
{
	printf("==================\n");
	printf("vtkremovebranching\n");
	printf("==================\n");
	
	if (argc != 3)
	{
		printf("Example: vtkremovebranching mesh.vtk out.vtk\n");
		return 0;
	}

	printf("Mesh file is %s\n", argv[1]);	
	
	// read the mesh file
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader->SetFileName(argv[1]);
    vtkMeshReader->Update();
	vtkPolyData* vtkMesh = vtkMeshReader->GetOutput();
	
	// find all cells boundary edges
	vtkCellArray* vtkCells = vtkMesh->GetPolys();
	vtkCells->InitTraversal();
	vtkIdType npts;
	vtkIdType* pts;
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
	
	// find branching edges
	int count = 0;
	vtkPoints* points = vtkMesh->GetPoints();
	char* msk = (char*) malloc(points->GetNumberOfPoints() * sizeof(char));
	memset(msk, 1, points->GetNumberOfPoints() * sizeof(char));
	for (map<int2, int>::iterator it = edges.begin(); it != edges.end(); it++ )
	{
		if (it->second > 2)
		{
			msk[(it->first).x] = 0;
			msk[(it->first).y] = 0;
			count++;
		}
	}
	
	cout << "Number of branching edges is " << count << "\n";
	
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
	
	// write file
	vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkMeshWriter->SetFileName(argv[2]);
	vtkMeshWriter->SetInput(vtkMeshNew);
	vtkMeshWriter->Write();
	
	// free memory
	free(msk);	
	
	return 0;
}