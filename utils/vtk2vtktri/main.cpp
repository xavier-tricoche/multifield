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
#include <bitset>
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
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkTriangleFilter.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

using namespace std;

int main( int argc, char** argv ) 
{
	if (argc != 3)
	{
		printf("Example: vtk2vtktri int.vtk out.vtk\n");
		return 0;
	}

	cout << "==========\n";
	cout << "vtk2vtktri\n";
	cout << "==========\n";
	
	//read file
	vtkPolyData* m;	
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader;
	vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader->SetFileName(argv[1]);
	vtkMeshReader->Update();
	m = vtkMeshReader->GetOutput();	
	cout << "Original: " << m->GetNumberOfPoints() << " points, " << m->GetNumberOfPolys() << " triangles\n";
	
	// decimation
	vtkSmartPointer<vtkTriangleFilter> deci;
	deci = vtkSmartPointer<vtkTriangleFilter>::New();
	deci->SetInput(m);
	deci->Update();
	cout << "Modified: " << deci->GetOutput()->GetNumberOfPoints() << " points, " << deci->GetOutput()->GetNumberOfPolys() << " triangles\n";
	
	// write file
	vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkMeshWriter->SetFileName(argv[2]);
	vtkMeshWriter->SetInput(deci->GetOutput());
	vtkMeshWriter->Write();
	

	return 0 ;      
}
