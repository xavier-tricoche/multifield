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

using namespace std;

float3 getNormal(vtkPoints* points, vtkSmartPointer<vtkIdList> ptIds, vtkIdType center)
{
	if (ptIds->GetNumberOfIds() != 3)
		printf("Not a triangle!\n");
		
	vtkIdType ids[3];
	float3 pts[3];
	
	for (int j = 0; j < ptIds->GetNumberOfIds(); j++)
	{
		ids[j] = ptIds->GetId(j);
		double* point = points->GetPoint(ids[j]);
		pts[j] = make_float3((float) point[0], (float) point[1], (float) point[2]);
	}
	
	/*if (ids[1] == center)
	{
		vtkIdType tmpid = ids[0];
		float3 tmppt = pts[0];
		ids[0] = ids[1];
		pts[0] = pts[1];
		ids[1] = tmpid;
		pts[1] = tmppt;		
	}
	else if if (ids[2] == center)
	{
		vtkIdType tmpid = ids[0];
		float3 tmppt = pts[0];
		ids[0] = ids[2];
		pts[0] = pts[2];
		ids[2] = tmpid;
		pts[2] = tmppt;		
	}*/
	
	float3 v1 = pts[2] - pts[0];
	float3 v2 = pts[1] - pts[0];
	
	return normalize(cross(v1,v2));
}

int main (int argc, char *argv[])
{
	printf("=======================\n");
	printf("vtkbreakdiscontinuities\n");
	printf("=======================\n");
	
	if (argc != 4)
	{
		printf("Example: vtkbreakdiscontinuities mesh.vtk out.vtk 20\n");
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
	
	// iterate on components
	printf("The number of points is %d\n", points->GetNumberOfPoints());
	printf("The number of polys is %d\n", vtkMesh->GetNumberOfPolys());
	int removed = 0;
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	{	
		// get point cells
		vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
		vtkMesh->GetPointCells(i, cells);
		
		// loop on all combinations of 2
		for (int j = 0; j < cells->GetNumberOfIds(); j++)
		{
			for (int k = 0; k < cells->GetNumberOfIds(); k++)
			{
				// get frist triangle
				vtkSmartPointer<vtkIdList> ptIds1 = vtkSmartPointer<vtkIdList>::New();
				vtkMesh->GetCellPoints(cells->GetId(j), ptIds1);
				
				// get second triangle
				vtkSmartPointer<vtkIdList> ptIds2 = vtkSmartPointer<vtkIdList>::New();
				vtkMesh->GetCellPoints(cells->GetId(k), ptIds2);
				
				// check normals
				float3 n1 = getNormal(points, ptIds1, i);
				float3 n2 = getNormal(points, ptIds2, i);
				
				if ((abs(dot(n1, n2)) < cos(atof(argv[3]))) && (msk[i] != 0))
				{
					msk[i] = 0;
					removed++;
				}
				
				if (msk[i] == 0)
					break;
			}
			
			if (msk[i] == 0)
				break;
		}
	}
	
	printf("The number of points removed is %d\n", removed);
	
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
	vtkCellArray* vtkCells = vtkMesh->GetPolys();
	vtkCells->InitTraversal();
	vtkIdType npts;
	vtkIdType* pts;
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