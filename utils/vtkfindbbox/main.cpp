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

using namespace std;

int main (int argc, char *argv[])
{
	if (argc < 2)
	{
		printf("Example: vtkfindbbox mesh.vtk\n");
		return 0;
	}

	cout << "===========\n";
	cout << "vtkfindbbox\n";
	cout << "===========\n";
	
	printf("Mesh file is %s\n", argv[1]);
	
	// read the mesh file
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader->SetFileName(argv[1]);
    vtkMeshReader->Update();
	vtkPolyData* vtkMesh = vtkMeshReader->GetOutput();
	
	// get the points
	float3 minb = make_float3(numeric_limits<float>::max());
	float3 maxb = make_float3(numeric_limits<float>::min());
	vtkPoints* points = vtkMesh->GetPoints();
	cout << "Number of points is " << points->GetNumberOfPoints() << "\n";
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	{
		double* point = points->GetPoint(i);
		
		minb.x = min(minb.x, (float)point[0]);
		minb.y = min(minb.y, (float)point[1]);
		minb.z = min(minb.z, (float)point[2]);
		
		maxb.x = max(maxb.x, (float)point[0]);
		maxb.y = max(maxb.y, (float)point[1]);
		maxb.z = max(maxb.z, (float)point[2]);
	}
 
	printf("Min: %f %f %f\n", minb.x, minb.y, minb.z);
	printf("Max: %f %f %f\n", maxb.x, maxb.y, maxb.z);
	
	return 0;
}