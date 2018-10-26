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
		printf("Example: vtkcomponents mesh.vtk out.vtk CCTBL\n");
		return 0;
	}

	cout << "=============\n";
	cout << "vtkcomponents\n";
	cout << "=============\n";
	
	printf("Mesh file is %s\n", argv[1]);
	
	// read the mesh file
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader->SetFileName(argv[1]);
    vtkMeshReader->Update();
	vtkPolyData* vtkMesh = vtkMeshReader->GetOutput();
	
	// start component analysis
	vtkPoints* points = vtkMesh->GetPoints();
	char* visited = (char*) malloc(points->GetNumberOfPoints() * sizeof(char));
	memset(visited, 0, points->GetNumberOfPoints() * sizeof(char));
	char* msk = (char*) malloc(points->GetNumberOfPoints() * sizeof(char));
	memset(msk, 1, points->GetNumberOfPoints() * sizeof(char));
	
	// points with no cells
	int ss = 0;
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	{	
		vtkSmartPointer<vtkIdList> pcells = vtkSmartPointer<vtkIdList>::New();
		vtkMesh->GetPointCells(i, pcells);
		if (pcells->GetNumberOfIds() == 0)
		{
			ss++;
		}
	}
	printf("Points with no cells are %d\n", ss);
	
	// create the array to fill
	vtkPointData* pointsData = vtkMesh->GetPointData();
	vtkSmartPointer<vtkIntArray> newScalars = vtkSmartPointer<vtkIntArray>::New();
	newScalars->SetName(argv[3]);
	
	// iterate on components
	int totalcomp = 0;
	int sum = 0;	
	cout << "The number of points is " << points->GetNumberOfPoints() << "\n";
	cout << "The number of polys is " << vtkMesh->GetNumberOfPolys() << "\n";
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	{	
		// if any point was visited then the whole component was
		if (visited[i]) continue;
		
		// read the complete component
		queue<vtkIdType> toprocess;
		queue<vtkIdType> comp;
		toprocess.push(i);
		comp.push(i);
		visited[i] = 1;
		sum++;
		newScalars->InsertValue(i, totalcomp);
		
		while (!toprocess.empty())
		{
			// get point
			vtkIdType p = toprocess.front();
			toprocess.pop();
			
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
					if (!visited[id])
					{
						toprocess.push(id);
						comp.push(id);
						visited[id] = 1;
						sum++;
						newScalars->InsertValue(id, totalcomp);
					}
				}
			}
		}
		
		printf("Component %d contains %d points\n", totalcomp, comp.size());
		
		totalcomp++;
	}
	pointsData->AddArray(newScalars);
	
	printf("\nPoints visited is %d\n" ,sum);
	printf("Number of components is %d\n", totalcomp);
	
	// write the mesh to file
	vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkMeshWriter->SetFileName(argv[2]);
	vtkMeshWriter->SetInput(vtkMesh);
	vtkMeshWriter->Write();
	
	// free memory
	free(visited);
	free(msk);

	
	return 0;
}