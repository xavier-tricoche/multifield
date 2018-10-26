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
#include <vtkCellLocator.h>

#include "RegularGrid.h"

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

void writeNrrd(void* data, const string& filename, int data_type, const vector< size_t >& dims, const vector< double >& spacing)
{
	Nrrd *nout = nrrdNew();

	if (nrrdWrap_nva(nout, data, data_type, dims.size(), &dims[0])) {
		cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
		exit(-1);
	}
	if (spacing.size() == dims.size()) {
		nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, (const void*)&spacing[0]);
	}
	if (nrrdSave(filename.c_str(), nout, NULL)) {
		cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
		exit(-1);
	}
}

// Load raw data from disk
void writeRawFile(float* data,const char *filename, int width, int height, int depth)
{
	vector< double > spacing;
	spacing.push_back(1.0);
	spacing.push_back(1.0);
	spacing.push_back(1.0);

	vector< size_t > dims;
	dims.push_back(width);
	dims.push_back(height);
	dims.push_back(depth);

	string file_string(filename);
	writeNrrd(data, file_string, nrrdTypeFloat, dims, spacing);

	printf("Write '%s'\n", filename);
}

double3 normalize(double3 v)
{
	double l = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
	v.x /= l;
	v.y /= l;
	v.z /= l;
	return v;
}

double dot(double3 v1, double3 v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

double3 cross(double3 a, double3 b)
{
	double3 r;
	r.x = (a.y * b.z - a.z * b.y);
	r.y = (a.z * b.x - a.x * b.z);
	r.z = (a.x * b.y - a.y * b.x);
	return r;
}

int main (int argc, char *argv[])
{
	if (argc != 5)
	{
		printf("Example: vtkridgesmooth mesh.vtk in.nrrd out.nrrd 0.2\n");
		return 0;
	}

	printf("Mesh file is %s\n", argv[1]);
	printf("Minimum distance factor is %f\n", atof(argv[4]));

	// read nrrd file and create grid
	Nrrd* nrrdfile = readNrrd(argv[2]);
	RegularGrid* dataset = new RegularGrid(nrrdfile, 
										   nrrdfile->axis[0].size, nrrdfile->axis[1].size, nrrdfile->axis[2].size, 
										   nrrdfile->axis[0].spacing, nrrdfile->axis[1].spacing, nrrdfile->axis[2].spacing);
										   
	// read the mesh file
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader->SetFileName(argv[1]);
    vtkMeshReader->Update();
	vtkPolyData* vtkMesh = vtkMeshReader->GetOutput();
	vtkPoints* points = vtkMesh->GetPoints();
	
	// create tree for the distance
	vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
	cellLocator->SetDataSet(vtkMesh);
	cellLocator->BuildLocator();
	
	// set values for the grid
	int size = nrrdfile->axis[0].size * nrrdfile->axis[1].size * nrrdfile->axis[2].size;
	float dec = atof(argv[4]);
	float* values = (float*) malloc(size * sizeof(float));
	//#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < size; i++)
	{
		if ((i%100000) == 0) 
			printf(".");
		
		int3 gpt = dataset->Addr2Coord(i);
		float3 spt = dataset->Grid2Space(make_float3(gpt.x, gpt.y, gpt.z));
		
		//Find the closest points to TestPoint
		double testPoint[3] = {spt.x, spt.y, spt.z};
		double closestPoint[3];//the coordinates of the closest point will be returned here
		double closestPointDist2; //the squared distance to the closest point will be returned here
		vtkIdType cellId; //the cell id of the cell containing the closest point will be returned here
		int subId; //this is rarely used (in triangle strips only, I believe)
		//cellLocator->FindClosestPoint(testPoint, closestPoint, cellId, subId, closestPointDist2);
		vtkIdType ret = cellLocator->FindClosestPointWithinRadius(testPoint, 4.0 * dataset->msp / dec, closestPoint, cellId, subId, closestPointDist2);
		if (ret == 0)
		{
			values[i] = 0.0;
			continue;
		}
		
		// get the cell normal vector
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
		vtkMesh->GetCellPoints(cellId, ptIds);
		if (ptIds->GetNumberOfIds() != 3)
			printf("Error: Cell is not a triangle!\n");
		double p1[3];
		double p2[3];
		double p3[3];
		double* p;
		p = points->GetPoint(ptIds->GetId(0));
		p1[0] = p[0]; p1[1] = p[1]; p1[2] = p[2];
		p = points->GetPoint(ptIds->GetId(1));
		p2[0] = p[0]; p2[1] = p[1]; p2[2] = p[2];
		p = points->GetPoint(ptIds->GetId(2));
		p3[0] = p[0]; p3[1] = p[1]; p3[2] = p[2];
		double3 v1 = make_double3((p2[0] - p1[0]), (p2[1] - p1[1]), (p2[2] - p1[2]));
		double3 v2 = make_double3((p3[0] - p1[0]), (p3[1] - p1[1]), (p3[2] - p1[2]));
		double3 n = normalize(cross(v1, v2));
		double3 v3 = normalize(make_double3((testPoint[0] - p1[0]), (testPoint[1] - p1[1]), (testPoint[2] - p1[2])));
		double fac = abs(dot(n, v3));
		
		if (isnan(fac) || isinf(fac))
		{
			printf("Error: %d\n ", i);
			fac = 0.5;
		}
		
		// compute the final value
		values[i] = (float) ((dec + (1.0 - dec) * fac) * sqrt(closestPointDist2));
		if (isnan(values[i]) || isinf(values[i]))
		{
			values[i] = 0.0;
		}
		
		if (values[i] > 3.0 * dataset->msp)
			values[i] = 3.0 * dataset->msp;
	}
	
	vector< double > spacing;
	spacing.push_back(nrrdfile->axis[0].spacing);
	spacing.push_back(nrrdfile->axis[1].spacing);
	spacing.push_back(nrrdfile->axis[2].spacing);

	vector< size_t > dims;
	dims.push_back(nrrdfile->axis[0].size);
	dims.push_back(nrrdfile->axis[1].size);
	dims.push_back(nrrdfile->axis[2].size);

	string file_string(argv[3]);
	writeNrrd(values, file_string, nrrdTypeFloat, dims, spacing);

	printf("Write '%s'\n", argv[3]);

	return 0;
}