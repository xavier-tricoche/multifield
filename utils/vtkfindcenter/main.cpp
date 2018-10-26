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
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkPolyDataWriter.h>
//#include <vtkCenterOfMass.h>

using namespace std;

int main (int argc, char *argv[])
{
	if (argc != 2)
	{
		printf("Example: vtkfindcenter mesh.vtk \n");
		return 0;
	}

	cout << "=============\n";
	cout << "vtkfindcenter\n";
	cout << "=============\n";
	
	printf("Mesh file is %s\n", argv[1]);
	
	// read the mesh file
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader->SetFileName(argv[1]);
    vtkMeshReader->Update();
	vtkPolyData* vtkMesh = vtkMeshReader->GetOutput();

	// get the points
	vtkPoints* points = vtkMesh->GetPoints();
	cout << "Number of points is " << points->GetNumberOfPoints() << "\n";
	
	// Compute the center of mass
	//vtkSmartPointer<vtkCenterOfMass> centerOfMassFilter = vtkSmartPointer<vtkCenterOfMass>::New();
	//centerOfMassFilter->SetInputData(vtkMesh);
	//centerOfMassFilter->Update();
	
	double bounds[6];
	vtkMesh->GetBounds(bounds);
	
	double center[3];
	//centerOfMassFilter->GetCenter(center);
	center[0] = 0.5 * (bounds[0]+bounds[1]);
	center[1] = 0.5 * (bounds[2]+bounds[3]);
	center[2] = 0.5 * (bounds[4]+bounds[5]);
	
	printf("Mesh Bounds %f %f %f %f %f %f\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);
	printf("Center at %lf %lf %lf\n", center[0], center[1], center[2]);
	
	return 0;
}