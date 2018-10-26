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
#include <vtkPolyDataWriter.h>

using namespace std;

Nrrd *readNrrd(const char* filename)
{
    Nrrd *nin = nrrdNew();
    if (nrrdLoad(nin, filename, NULL)) {
        printf("Error reading %s\n",filename);
                return NULL;
    }
    return nin;
}

int main (int argc, char *argv[])
{
	if (argc != 3)
	{
		printf("Example: vtkadjustunitbbox from.vtk file.nrrd\n");
		return 0;
	}
	
	cout << "=================\n";
	cout << "vtkadjustunitbbox\n";
	cout << "=================\n";
	
	Nrrd* nrrdfile = readNrrd(argv[2]);
	
	float3 mint = make_float3(0.0);
	float3 maxt = make_float3((nrrdfile->axis[0].size - 1) * nrrdfile->axis[0].spacing,
							  (nrrdfile->axis[1].size - 1) * nrrdfile->axis[1].spacing,
							  (nrrdfile->axis[2].size - 1) * nrrdfile->axis[2].spacing);
							  
	float maxs = max(fabs(maxt.x), max(fabs(maxt.y), fabs(maxt.z)));
	float3 minf = make_float3(-fabs(maxt.x) / maxs, -fabs(maxt.y) / maxs, -fabs(maxt.z) / maxs);
	float3 maxf = make_float3(fabs(maxt.x) / maxs, fabs(maxt.y) / maxs, fabs(maxt.z) / maxs);
	
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
	
	// adjust the values
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	{
		double* point = points->GetPoint(i);
		
		point[0] = mint.x + (maxt.x - mint.x) * (point[0] - minf.x) / (maxf.x - minf.x);
		point[1] = mint.y + (maxt.y - mint.y) * (point[1] - minf.y) / (maxf.y - minf.y);
		point[2] = mint.z + (maxt.z - mint.z) * (point[2] - minf.z) / (maxf.z - minf.z);
		
		points->SetPoint(i, point[0], point[1], point[2]);
	}
	
	// write file
	vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkMeshWriter->SetFileName(argv[1]);
	vtkMeshWriter->SetInput(vtkMesh);
	vtkMeshWriter->Write();
	
	nrrdNuke(nrrdfile);
	return 0;
}