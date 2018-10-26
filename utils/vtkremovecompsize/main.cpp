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
		printf("Example: vtkremovecompsize mesh.vtk out.vtk 0.02\n");
		return 0;
	}

	cout << "=================\n";
	cout << "vtkremovecompsize\n";
	cout << "=================\n";
	
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
	int sizelimit = atof(argv[3]) * points->GetNumberOfPoints();
	printf("Size limit is %d\n", sizelimit);
	
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
	
	// iterate on components
	int totalcomp = 0;
	int removedcomp = 0;
	int sum = 0;	
	cout << "The number of points is " << points->GetNumberOfPoints() << "\n";
	cout << "The number of polys is " << vtkMesh->GetNumberOfPolys() << "\n";
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	{	
		// if any point was visited then the whole component was
		if (visited[i]) continue;
		
		// read the complete component
		queue<vtkIdType> comp;
		queue<vtkIdType> cp;
		cp.push(i);
		comp.push(i);
		visited[i] = 1;
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
					if (!visited[id])
					{
						cp.push(id);
						comp.push(id);
						visited[id] = 1;
					}
				}
			}
		}
		
		totalcomp++;
		//printf("%d ", comp.size());
		
		// evaluate conditions to delete component
		sum += comp.size();
		if (comp.size() < sizelimit)
		{
			removedcomp ++;
			while (!comp.empty())
			{
				// get point
				vtkIdType p = comp.front();
				comp.pop();
				
				msk[p] = 0;
			}
		}
	}
	
	printf("\nPoints visited is %d\n" ,sum);
	printf("Removed %d components of %d\n", removedcomp, totalcomp);
	
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
	free(visited);
	free(msk);

	
	return 0;
}