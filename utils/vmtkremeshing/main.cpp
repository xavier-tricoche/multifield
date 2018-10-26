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
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>

#include "vtkvmtkPolyDataSurfaceRemeshing.h"

using namespace std;

int main( int argc, char** argv ) 
{
	if (argc != 5)
	{
		printf("Example: vmtkremeshing int.vtk out.vtk 5 0.005\n");
		return 0;
	}

	cout << "=============\n";
	cout << "vmtkremeshing\n";
	cout << "=============\n";
	
	int iteration = atoi(argv[3]);
	double edgelen = atof(argv[4]);
	
	//read file
	vtkPolyData* m;	
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader;
	vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader->SetFileName(argv[1]);
	vtkMeshReader->Update();
	m = vtkMeshReader->GetOutput();	
	cout << "Original: " << m->GetNumberOfPoints() << " points, " << m->GetNumberOfPolys() << " triangles\n";
	double bounds[6];
	m->GetBounds(bounds);
	double TargetEdgeLength = edgelen * max((bounds[1] - bounds[0]), max((bounds[3] - bounds[2]),(bounds[5] - bounds[4])));
	printf("Target length is %lf\n", TargetEdgeLength);
	
	// make sure all are triangles
	vtkTriangleFilter* triangleFilter = vtkTriangleFilter::New();
    triangleFilter->SetInput(m);
    triangleFilter->Update();

	// clean the mesh 
	vtkCleanPolyData* cleaner = vtkCleanPolyData::New();
	cleaner->SetInput(triangleFilter->GetOutput());
	cleaner->PointMergingOff();
    cleaner->Update();
  
	// target general area
	double TargetArea = 0.25 * pow(3.0, 0.5) * pow(TargetEdgeLength,2);
			
	// remesh when necessary
	vtkvmtkPolyDataSurfaceRemeshing* remesh = vtkvmtkPolyDataSurfaceRemeshing::New();
	remesh->SetInput(triangleFilter->GetOutput());
	remesh->SetNumberOfIterations(iteration);
	remesh->SetNumberOfConnectivityOptimizationIterations(2 * iteration);
	//remesh->SetCollapseAngleThreshold(0.2);
	remesh->SetElementSizeModeToTargetArea();
	remesh->SetTargetArea(TargetArea);
	remesh->Update();

	// clean the mesh 
	vtkCleanPolyData* cleanerAfter = vtkCleanPolyData::New();
    cleanerAfter->SetInput(remesh->GetOutput());
	cleanerAfter->PointMergingOff();
    cleanerAfter->Update();
	cout << "Modified: " << cleanerAfter->GetOutput()->GetNumberOfPoints() << " points, " << cleanerAfter->GetOutput()->GetNumberOfPolys() << " triangles\n";

	// write file
	vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkMeshWriter->SetFileName(argv[2]);
	vtkMeshWriter->SetInput(cleanerAfter->GetOutput());
	vtkMeshWriter->Write();

	return 0;      
}
