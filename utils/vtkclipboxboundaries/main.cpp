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
#include <vtkDecimatePro.h>
#include <vtkPolyDataWriter.h>
#include <vtkPlane.h>
#include <vtkClipPolyData.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

using namespace std;

int main( int argc, char** argv ) 
{
	if (argc != 4)
	{
		printf("Example: vtkclipboxboundaries int.vtk out.vtk 0.01\n");
		return 0;
	}

	cout << "====================\n";
	cout << "vtkclipboxboundaries\n";
	cout << "====================\n";
	
	//read file
	vtkPolyData* m;	
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader->SetFileName(argv[1]);
	vtkMeshReader->Update();
	m = vtkMeshReader->GetOutput();	
	double bbox[6];
	m->GetBounds(bbox);
	cout << "Bounding box: [" << bbox[0] << ", " << bbox[1] << "] [" << bbox[2] << ", " << bbox[3] << "] [" << bbox[4] << ", " << bbox[5] << "]\n";
	cout << "Original: " << m->GetNumberOfPoints() << " points, " << m->GetNumberOfPolys() << " triangles\n";
	
	// create all the clipping planes
	double r = atof(argv[3]);
	vtkSmartPointer<vtkPlane> clipplane[6];
	vtkSmartPointer<vtkClipPolyData> clipper[6];
	for (int i = 0; i < 6; i++)
	{
		double x = (i < 2)? 1:0;
		double y = (i >= 2) && (i < 4)? 1:0; 
		double z = (i >= 4)? 1:0;
		float3 normal = make_float3(x, y, z);
		if (i%2 > 0)
		{
			normal *= -1;
		}
		
		x *= (normal.x > 0)? bbox[1] - r * (bbox[1] - bbox[0]) : bbox[0] + r * (bbox[1] - bbox[0]);
		y *= (normal.y > 0)? bbox[3] - r * (bbox[3] - bbox[2]) : bbox[2] + r * (bbox[3] - bbox[2]);
		z *= (normal.z > 0)? bbox[5] - r * (bbox[5] - bbox[4]) : bbox[4] + r * (bbox[5] - bbox[4]);
		
		clipplane[i] = vtkSmartPointer<vtkPlane>::New();
		clipplane[i]->SetOrigin(x, y, z);
		clipplane[i]->SetNormal(-normal.x, -normal.y, -normal.z);

		clipper[i] = vtkSmartPointer<vtkClipPolyData>::New();
		clipper[i]->SetInput(m);
		clipper[i]->SetClipFunction(clipplane[i]);
		clipper[i]->Update();
		m = clipper[i]->GetOutput();
	}
	cout << "Final: " << m->GetNumberOfPoints() << " points, " << m->GetNumberOfPolys() << " triangles\n";
	
	// write file
	vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkMeshWriter->SetFileName(argv[2]);
	vtkMeshWriter->SetInput(m);
	vtkMeshWriter->Write();


	return 0 ;      
}
