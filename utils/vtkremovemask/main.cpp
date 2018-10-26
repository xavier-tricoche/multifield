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

int main (int argc, char *argv[])
{
	if (argc != 4)
	{
		printf("Example: vtkremovemask mesh.vtk out.vtk MSKTBL\n");
		return 0;
	}

	cout << "=============\n";
	cout << "vtkremovemask\n";
	cout << "=============\n";
	
	printf("Mesh file is %s\n", argv[1]);
	
	// read the mesh file
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader->SetFileName(argv[1]);
    vtkMeshReader->Update();
	vtkPolyData* vtkMesh = vtkMeshReader->GetOutput();
	vtkPoints* points = vtkMesh->GetPoints();
	
	// new mesh points
	vtkPointData* pointsData = vtkMesh->GetPointData();
	vtkFloatArray* array = (vtkFloatArray*) pointsData->GetArray(argv[3]);
	vtkSmartPointer<vtkPolyData> vtkMeshNew = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> pointsNew = vtkSmartPointer<vtkPoints>::New();
	map<vtkIdType, vtkIdType> old2new;
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	{
		if (array->GetValue(i) < 0.1)
		{
			double* point = points->GetPoint(i);
			vtkIdType newid = pointsNew->InsertNextPoint(point[0], point[1], point[2]);
			old2new.insert(pair<vtkIdType, vtkIdType>(i, newid));
		}
	}
	vtkMeshNew->SetPoints(pointsNew);
	printf("Original number of vertices is %d\n", points->GetNumberOfPoints());
	printf("Modified number of vertices is %d\n", vtkMeshNew->GetNumberOfPoints());
	
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
		
	return 0;
}