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
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
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

double DistPoint3Triangle3Squared(float3 mTriangle[3], float3 mPoint)
{
    float3 diff = mTriangle[0] - mPoint;
    float3 edge0 = mTriangle[1] - mTriangle[0];
    float3 edge1 = mTriangle[2] - mTriangle[0];
    double a00 = dot(edge0, edge0);
    double a01 = dot(edge0, edge1);
    double a11 = dot(edge1, edge1);
    double b0 = dot(diff, edge0);
    double b1 = dot(diff, edge1);
    double c = dot(diff, diff);
    double det = abs(a00*a11 - a01*a01);
    double s = a01*b1 - a11*b0;
    double t = a01*b0 - a00*b1;
    double sqrDistance;

    if (s + t <= det)
    {
        if (s < (double)0)
        {
            if (t < (double)0)  // region 4
            {
                if (b0 < (double)0)
                {
                    t = (double)0;
                    if (-b0 >= a00)
                    {
                        s = (double)1;
                        sqrDistance = a00 + ((double)2)*b0 + c;
                    }
                    else
                    {
                        s = -b0/a00;
                        sqrDistance = b0*s + c;
                    }
                }
                else
                {
                    s = (double)0;
                    if (b1 >= (double)0)
                    {
                        t = (double)0;
                        sqrDistance = c;
                    }
                    else if (-b1 >= a11)
                    {
                        t = (double)1;
                        sqrDistance = a11 + ((double)2)*b1 + c;
                    }
                    else
                    {
                        t = -b1/a11;
                        sqrDistance = b1*t + c;
                    }
                }
            }
            else  // region 3
            {
                s = (double)0;
                if (b1 >= (double)0)
                {
                    t = (double)0;
                    sqrDistance = c;
                }
                else if (-b1 >= a11)
                {
                    t = (double)1;
                    sqrDistance = a11 + ((double)2)*b1 + c;
                }
                else
                {
                    t = -b1/a11;
                    sqrDistance = b1*t + c;
                }
            }
        }
        else if (t < (double)0)  // region 5
        {
            t = (double)0;
            if (b0 >= (double)0)
            {
                s = (double)0;
                sqrDistance = c;
            }
            else if (-b0 >= a00)
            {
                s = (double)1;
                sqrDistance = a00 + ((double)2)*b0 + c;
            }
            else
            {
                s = -b0/a00;
                sqrDistance = b0*s + c;
            }
        }
        else  // region 0
        {
            // minimum at interior point
            double invDet = ((double)1)/det;
            s *= invDet;
            t *= invDet;
            sqrDistance = s*(a00*s + a01*t + ((double)2)*b0) +
                t*(a01*s + a11*t + ((double)2)*b1) + c;
        }
    }
    else
    {
        double tmp0, tmp1, numer, denom;

        if (s < (double)0)  // region 2
        {
            tmp0 = a01 + b0;
            tmp1 = a11 + b1;
            if (tmp1 > tmp0)
            {
                numer = tmp1 - tmp0;
                denom = a00 - ((double)2)*a01 + a11;
                if (numer >= denom)
                {
                    s = (double)1;
                    t = (double)0;
                    sqrDistance = a00 + ((double)2)*b0 + c;
                }
                else
                {
                    s = numer/denom;
                    t = (double)1 - s;
                    sqrDistance = s*(a00*s + a01*t + ((double)2)*b0) +
                        t*(a01*s + a11*t + ((double)2)*b1) + c;
                }
            }
            else
            {
                s = (double)0;
                if (tmp1 <= (double)0)
                {
                    t = (double)1;
                    sqrDistance = a11 + ((double)2)*b1 + c;
                }
                else if (b1 >= (double)0)
                {
                    t = (double)0;
                    sqrDistance = c;
                }
                else
                {
                    t = -b1/a11;
                    sqrDistance = b1*t + c;
                }
            }
        }
        else if (t < (double)0)  // region 6
        {
            tmp0 = a01 + b1;
            tmp1 = a00 + b0;
            if (tmp1 > tmp0)
            {
                numer = tmp1 - tmp0;
                denom = a00 - ((double)2)*a01 + a11;
                if (numer >= denom)
                {
                    t = (double)1;
                    s = (double)0;
                    sqrDistance = a11 + ((double)2)*b1 + c;
                }
                else
                {
                    t = numer/denom;
                    s = (double)1 - t;
                    sqrDistance = s*(a00*s + a01*t + ((double)2)*b0) +
                        t*(a01*s + a11*t + ((double)2)*b1) + c;
                }
            }
            else
            {
                t = (double)0;
                if (tmp1 <= (double)0)
                {
                    s = (double)1;
                    sqrDistance = a00 + ((double)2)*b0 + c;
                }
                else if (b0 >= (double)0)
                {
                    s = (double)0;
                    sqrDistance = c;
                }
                else
                {
                    s = -b0/a00;
                    sqrDistance = b0*s + c;
                }
            }
        }
        else  // region 1
        {
            numer = a11 + b1 - a01 - b0;
            if (numer <= (double)0)
            {
                s = (double)0;
                t = (double)1;
                sqrDistance = a11 + ((double)2)*b1 + c;
            }
            else
            {
                denom = a00 - ((double)2)*a01 + a11;
                if (numer >= denom)
                {
                    s = (double)1;
                    t = (double)0;
                    sqrDistance = a00 + ((double)2)*b0 + c;
                }
                else
                {
                    s = numer/denom;
                    t = (double)1 - s;
                    sqrDistance = s*(a00*s + a01*t + ((double)2)*b0) +
                        t*(a01*s + a11*t + ((double)2)*b1) + c;
                }
            }
        }
    }

    // Account for numerical round-off error.
    if (sqrDistance < (double)0)
    {
        sqrDistance = (double)0;
    }

	//mClosestPoint0 = mPoint;
    //mClosestPoint1 = mTriangle[0] + s*edge0 + t*edge1;
    //mTriangleBary[1] = s;
    //mTriangleBary[2] = t;
    //mTriangleBary[0] = (double)1 - s - t;
    return sqrDistance;
}

void computeDistanceField_kernel(RegularGrid* dataset, float* field, int c, bool& changed)
{
	int3 gpt = dataset->Addr2Coord(c);
	float3 cpt = dataset->Grid2Space(make_float3(gpt.x, gpt.y, gpt.z));

	int x = gpt.x;
	int y = gpt.y;
	int z = gpt.z;
	float nval;
	for (int i = -1; i <= 1; i++)
	{
		if ((x + i < 0) || (x + i > dataset->width - 1))
			continue;
		for (int j = -1; j <= 1; j++)
		{
			if ((y + j < 0) || (y + j > dataset->height - 1))
				continue;
			for (int k = -1; k <= 1; k++)
			{
				if ((z + k < 0) || (z + k > dataset->depth - 1))
					continue;
					
				int idx = dataset->Coord2Addr(x + i, y + j, z + k);
				
				if (field[idx] == numeric_limits<float>::max())
					continue;
				
				nval = field[idx] + length(make_float3(i * dataset->xspc, j * dataset->yspc, k * dataset->zspc));
				if (field[c] > nval)
				{
					field[c] = nval;
					changed = true;
				}
			}
		}
	}
}

void ComputeDistanceField1(RegularGrid* dataset, float* field)
{
	bool changed = true;
	while (changed)
	{
		printf(".\n");
		changed = false;
		
		#pragma omp parallel for schedule(dynamic, 10000)
		for (int c = 0; c < dataset->width * dataset->height * dataset->depth; c++)
		{
			// check change
			computeDistanceField_kernel(dataset, field, c, changed);
		}
	}
}

void ComputeDistanceField2(RegularGrid* dataset, float* field)
{
	// array to mark the voxels to update
	char* mark = (char*) malloc(dataset->width * dataset->height * dataset->depth * sizeof(char));
	memset(mark, 0, dataset->width * dataset->height * dataset->depth * sizeof(char));
	
	// array to mark the voxels to update
	char* marked = (char*) malloc(dataset->width * dataset->height * dataset->depth * sizeof(char));
	memset(marked, 1, dataset->width * dataset->height * dataset->depth * sizeof(char));
	
	bool changed = true;
	while (changed)
	{
		printf(".\n");
		changed = false;
		
		#pragma omp parallel for schedule(dynamic, 10000)
		for (int c = 0; c < dataset->width * dataset->height * dataset->depth; c++)
		{
			if (!marked[c]) continue;
			marked[c] = 0;
			
			// check change
			computeDistanceField_kernel(dataset, field, c, changed);
			
			// now mark the neighbors if a change occured
			if (changed)
			{
				int3 cpt = dataset->Addr2Coord(c);
				
				for (int i = -1; i <= 1; i++)
				{
					for (int j = -1; j <= 1; j++)
					{
						for (int k = -1; k <= 1; k++)
						{
							if (!dataset->IsValid(cpt.x+i, cpt.y+j, cpt.z+k)) continue;

							int idx = dataset->Coord2Addr(cpt.x+i, cpt.y+j, cpt.z+k);
							mark[idx] = 1;
						}
					}
				}
			}
		}
		
		// switch 
		char* tmp = marked;
		marked = mark;
		mark = tmp;
	}
	
	free(mark);
	free(marked);
}

void ComputeDistanceField3(RegularGrid* dataset, float* field)
{
	// array to mark the voxels to update
	char* mrk = (char*) malloc(dataset->width * dataset->height * dataset->depth * sizeof(char));
	memset(mrk, 0, dataset->width * dataset->height * dataset->depth * sizeof(char));
	
	// list to check for update
	vector<int> mq;
	
	// initialize the list
	for (int c = 0; c < dataset->width * dataset->height * dataset->depth; c++)
	{
		int3 cpt = dataset->Addr2Coord(c);
		if (field[c] < numeric_limits<float>::max())
		{
			for (int i = -1; i <= 1; i++)
			{
				for (int j = -1; j <= 1; j++)
				{
					for (int k = -1; k <= 1; k++)
					{
						if (!dataset->IsValid(cpt.x+i, cpt.y+j, cpt.z+k)) continue;

						int idx = dataset->Coord2Addr(cpt.x+i, cpt.y+j, cpt.z+k);
						mq.insert(mq.end(), idx);
						mrk[idx] = 1;
					}
				}
			}
		}
	}
	printf("Initialization done\n");
	
	// now loop on the list to checck for updates
	int count = 0;
	while(!mq.empty())
	{
		count ++;		
		if (count % 100000 == 0) cout << mq.size() << "\n";
		
		// get one item from the list
		int c = mq.back();
		mq.pop_back();
		if (mrk[c] == 0) continue;
		mrk[c] = 0;

		// check the item
		bool changed = false;
		computeDistanceField_kernel(dataset, field, c, changed);
				
		// now mark the neighbors if a change occured
		if (changed)
		{
			int3 cpt = dataset->Addr2Coord(c);
			
			for (int i = -1; i <= 1; i++)
			{
				for (int j = -1; j <= 1; j++)
				{
					for (int k = -1; k <= 1; k++)
					{
						if (!dataset->IsValid(cpt.x+i, cpt.y+j, cpt.z+k)) continue;

						int idx = dataset->Coord2Addr(cpt.x+i, cpt.y+j, cpt.z+k);
						if (mrk[idx] == 0)
						{
							mq.insert(mq.end(), idx);
							mrk[idx] = 1;
						}
					}
				}
			}
		}
	}
	
	free(mrk);
}

int main (int argc, char *argv[])
{
	if (argc != 5)
	{
		printf("Example: vtkstrengthfield mesh.vtk in.nrrd out.nrrd 3.0\n");
		return 0;
	}

	printf("Mesh file is %s\n", argv[1]);
	
	cout << "================\n";
	cout << "vtkstrengthfield\n";
	cout << "================\n";

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
	
	// initialize
	int size = nrrdfile->axis[0].size * nrrdfile->axis[1].size * nrrdfile->axis[2].size;
	float* values = (float*) malloc(size * sizeof(float));
	//#pragma omp parallel for
	for (int c = 0; c < dataset->width * dataset->height * dataset->depth; c++)
	{
		values[c] = dataset->msp * atof(argv[4]);//numeric_limits<float>::max();
	}
	
	// loop on cells and update distance at neighbors
	cout << "Number of polys is " << vtkMesh->GetNumberOfPolys() << "\n";
	vtkCellArray* vtkCells = vtkMesh->GetPolys();
	vtkCells->InitTraversal();
	vtkIdType npts;
	vtkIdType* pts;
	int count = 0;
	while (vtkCells->GetNextCell(npts, pts) != 0)
	{
		if (count % 10000 == 0) printf("%d\n",count);
		count++;
		
		double* point;
		float3 minb = make_float3(numeric_limits<float>::max(),numeric_limits<float>::max(),numeric_limits<float>::max());
		float3 maxb = make_float3(numeric_limits<float>::min(),numeric_limits<float>::min(),numeric_limits<float>::min());
		float3 pt[3];
		
		// find the cell bounding box
		for (int i = 0; i < npts; i++)
		{
			point = points->GetPoint(pts[i]);
			pt[i] = make_float3((float)point[0], (float)point[1], (float)point[2]);
			
			//printf("%f %f %f\n", pt[i].x, pt[i].y, pt[i].z);
			
			minb.x = min(minb.x, pt[i].x);
			minb.y = min(minb.y, pt[i].y);
			minb.z = min(minb.z, pt[i].z);
			
			maxb.x = max(maxb.x, pt[i].x);
			maxb.y = max(maxb.y, pt[i].y);
			maxb.z = max(maxb.z, pt[i].z);
		}
		
		if (npts != 3) 
			printf("Cell is not a triangle and has been skipped!\n");
		
		// find grid bounding box
		minb = dataset->Space2Grid(minb);
		minb.x = floor(minb.x);
		minb.y = floor(minb.y);
		minb.z = floor(minb.z);		
		maxb = dataset->Space2Grid(maxb);
		maxb.x = ceil(maxb.x);
		maxb.y = ceil(maxb.y);
		maxb.z = ceil(maxb.z);
		
		//printf("%f %f %f -> %f %f %f\n", minb.x, minb.y, minb.z, maxb.x, maxb.y, maxb.z);
		
		//#pragma omp parallel for schedule(dynamic)
		for (int x = cint(minb.x); x <= cint(maxb.x); x++)
		{
			for (int y = cint(minb.y); y <= cint(maxb.y); y++)
			{
				for (int z = cint(minb.z); z <= cint(maxb.z); z++)
				{
					if (!dataset->IsValid(x, y, z)) continue;
					
					float3 spt = dataset->Grid2Space(make_float3(x,y,z));
					int idx = dataset->Coord2Addr(x,y,z);
					float dist = sqrt(DistPoint3Triangle3Squared(pt, spt));
					//#pragma omp critical
					{
						values[idx] = min(values[idx], dist);
					}
				}
			}
		}
	}
	
	// grow distance field
	printf("Propagate\n");
	ComputeDistanceField1(dataset, values);
	
	// write the output file
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