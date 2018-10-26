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
#include <vtkStructuredPointsReader.h>
#include <vtkFloatArray.h>
#include <vtkContourFilter.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>

using namespace std;

int main (int argc, char *argv[])
{
	if (argc != 4)
	{
		printf("Example: vtkisocontour grid.vtk mesh.vtk 4.5\n");
		return 0;
	}

	printf("Grid file is %s\n", argv[1]);
	printf("Mesh file is %s\n", argv[2]);
	printf("Iso value is %f\n", atof(argv[3]));
	
	// read the grid file
	vtkSmartPointer<vtkStructuredPointsReader> vtkMeshReader = vtkSmartPointer<vtkStructuredPointsReader>::New();
	vtkMeshReader->SetFileName(argv[1]);
    vtkMeshReader->Update();
	
	// extract the iso contour
	vtkSmartPointer<vtkContourFilter> isocontour = vtkSmartPointer<vtkContourFilter>::New();
	isocontour->SetInputConnection(vtkMeshReader->GetOutputPort());
	isocontour->SetValue(0,atof(argv[3]));
	
	// write file
	vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkMeshWriter->SetFileName(argv[2]);
	vtkMeshWriter->SetInput(isocontour->GetOutput());
	vtkMeshWriter->Write();

	
	return 0;
}