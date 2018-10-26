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


using namespace std;


int main (int argc, char *argv[])
{
	if (argc != 3)
	{
		printf("Example: vtktagdocument mesh.vtk document.txt\n");
		return 0;
	}

	cout << "==============\n";
	cout << "vtktagdocument\n";
	cout << "==============\n";
	
	printf("Mesh file is %s\n", argv[1]);
	
	// read the mesh file
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader->SetFileName(argv[1]);
    vtkMeshReader->Update();
	vtkPolyData* vtkMesh = vtkMeshReader->GetOutput();
	
	// get the points
	vtkPoints* points = vtkMesh->GetPoints();
	cout << "Number of points is " << points->GetNumberOfPoints() << "\n";
	vtkPointData* pointsData = vtkMesh->GetPointData();
	
	// create file
	ofstream myfile;
	myfile.open(argv[2]);
	
	// find and copy existing table
	vtkIntArray* arr = (vtkIntArray*) pointsData->GetArray("TAGTBL");
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	{
		int tag = arr->GetValue(i);
		myfile << tag << " ";
	}
	
	// close file
	myfile.close();
	
	return 0;
}