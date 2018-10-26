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
#include <limits.h>
#include <string.h>
#include <cmath>
#include <time.h>
#include <vector>
#include <queue>
#include <map>
#include <bitset>
#include <algorithm>
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

#include "RegularGrid.h"

using namespace std;

int mesh_count = 0;
char** meshes_files;
vtkPolyData** meshes;
vtkSmartPointer<vtkPolyDataReader>* vtkMeshReader;
vector< vector<vtkIdType> > sortedvertices;
char** visited;
RegularGrid* distfield;
int** rank;
Nrrd* nrrdfile;

class myvertexsort {
public:
	vtkFloatArray* arr;
	
	myvertexsort(vtkPolyData* vtkMesh)
	{
		vtkPointData* pointsData = vtkMesh->GetPointData();
		arr = (vtkFloatArray*) pointsData->GetArray("STRTBL");
		if (arr == NULL) printf("Table not found!\n");
	}
	
	bool operator() (int i,int j) 
	{ 
		return (-arr->GetValue(i) < -arr->GetValue(j));
	}	
};

bool operator<(const int2 e1, const int2 e2)
{
	if (e1.x < e2.x)
	{
		return true;
	}
	else if (e1.x > e2.x)
	{
		return false;
	}
	else
	{
		if (e1.y < e2.y)
		{
			return true;
		}
		else if (e1.y > e2.y)
		{
			return false;
		}
		else
		{
			return false;
		}
	}
}

bool operator==(const int2 e1, const int2 e2)
{
	if ((e1.x == e2.x) && (e1.y == e2.y))
		return true;
	return false;
}

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

void ComputeDistanceField(RegularGrid* dataset, float* field)
{
	bool changed = true;
	while (changed)
	{
		printf(".");
		fflush(stdout);
		changed = false;
		
		#pragma omp parallel for schedule(dynamic, 10000)
		for (int c = 0; c < dataset->width * dataset->height * dataset->depth; c++)
		{
			// check change
			computeDistanceField_kernel(dataset, field, c, changed);
		}
	}
}

void LoadMeshesAndLists()
{	
	meshes = (vtkPolyData**) malloc(mesh_count * sizeof(vtkPolyData*));
	sortedvertices.resize(mesh_count);
	visited = (char**) malloc(mesh_count * sizeof(char*));
	rank = (int**) malloc(mesh_count * sizeof(int*));
	vtkMeshReader = (vtkSmartPointer<vtkPolyDataReader>*) malloc(mesh_count * sizeof(vtkSmartPointer<vtkPolyDataReader>));
	for (int i = 0; i < mesh_count; i++)
	{
		vtkMeshReader[i] = vtkSmartPointer<vtkPolyDataReader>::New();
		vtkMeshReader[i]->SetFileName(meshes_files[i]);
		vtkMeshReader[i]->Update();
		meshes[i] = vtkMeshReader[i]->GetOutput();
		printf("Mesh %s loaded\n", meshes_files[i]);
	
		for (vtkIdType j = 0; j < meshes[i]->GetNumberOfPoints(); j++)
		{
			sortedvertices[i].push_back(j);
		}
	
		sort(sortedvertices[i].begin(), sortedvertices[i].end(), myvertexsort(meshes[i]));
		visited[i] = (char*) malloc(meshes[i]->GetNumberOfPoints() * sizeof(char));
		memset(visited[i], 0, meshes[i]->GetNumberOfPoints() * sizeof(char));
		rank[i] = (int*) malloc(meshes[i]->GetNumberOfPoints() * sizeof(int));
	}
}

void ComputeRanks(int meshid, int ranklimit)
{
	vtkPolyData* vtkMesh = meshes[meshid];
	for (vtkIdType j = 0; j < vtkMesh->GetNumberOfPoints(); j++)
		rank[meshid][j] = ranklimit;
	
	// find all cells boundary edges
	vtkCellArray* vtkCells = vtkMesh->GetPolys();
	vtkCells->InitTraversal();
	vtkIdType npts;
	vtkIdType* pts;
	multimap<vtkIdType,vtkIdType> mymultimap;
	int count = 0;
	map<int2, int> edges;
	while (vtkCells->GetNextCell(npts, pts) != 0)
	{
		int2 edge1 = make_int2(min(pts[0], pts[1]), max(pts[0], pts[1]));
		int2 edge2 = make_int2(min(pts[0], pts[2]), max(pts[0], pts[2]));
		int2 edge3 = make_int2(min(pts[1], pts[2]), max(pts[1], pts[2]));
		
		map<int2, int>::iterator it;
		
		if ((it = edges.find(edge1)) != edges.end())
		{
			int scr = it->second + 1;
			edges.erase(it);
			edges.insert(pair<int2, int>(edge1, scr));
		}
		else
			edges.insert(pair<int2, int>(edge1, 1));
			
		if ((it = edges.find(edge2)) != edges.end())
		{	
			int scr = it->second + 1;
			edges.erase(it);
			edges.insert(pair<int2, int>(edge2, scr));
		}
		else
			edges.insert(pair<int2, int>(edge2, 1));
			
		if ((it = edges.find(edge3)) != edges.end())
		{	
			int scr = it->second + 1;
			edges.erase(it);
			edges.insert(pair<int2, int>(edge3, scr));
		}
		else
			edges.insert(pair<int2, int>(edge3, 1));
			
	}
	
	// boundary edges are not shared by 2 triangles
	for (map<int2, int>::iterator it = edges.begin(); it != edges.end(); it++ )
	{
		if (it->second != 2)
		{
			rank[meshid][(it->first).x] = 0;
			rank[meshid][(it->first).y] = 0;
		}
	}
	
	// update the ranks for neighbors
	for (int i = 0; i < ranklimit; i++)
	{
		vtkCells->InitTraversal();
		while (vtkCells->GetNextCell(npts, pts) != 0)
		{
			int2 edge1 = make_int2(min(pts[0], pts[1]), max(pts[0], pts[1]));
			int2 edge2 = make_int2(min(pts[0], pts[2]), max(pts[0], pts[2]));
			int2 edge3 = make_int2(min(pts[1], pts[2]), max(pts[1], pts[2]));
			
			rank[meshid][edge1.x] = min(rank[meshid][edge1.x], 1 + rank[meshid][edge1.y]); 
			rank[meshid][edge1.y] = min(rank[meshid][edge1.y], 1 + rank[meshid][edge1.x]);
			
			rank[meshid][edge2.x] = min(rank[meshid][edge2.x], 1 + rank[meshid][edge2.y]); 
			rank[meshid][edge2.y] = min(rank[meshid][edge2.y], 1 + rank[meshid][edge2.x]);
			
			rank[meshid][edge3.x] = min(rank[meshid][edge3.x], 1 + rank[meshid][edge3.y]); 
			rank[meshid][edge3.y] = min(rank[meshid][edge3.y], 1 + rank[meshid][edge3.x]);
		}
	}
	printf("Rank computation for mesh %d is complete.\n", meshid);
	
	// find the best rank
	/*int rr[100];
	memset(rr, 0, 100*sizeof(int));
	for (vtkIdType p = 0; p < vtkMesh->GetNumberOfPoints(); p++)
	{
		int count = 0;
		vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
		vtkMesh->GetPointCells(p, cells);
		for (int k = 0; k < cells->GetNumberOfIds(); k++)
		{
			vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
			vtkMesh->GetCellPoints(cells->GetId(k), ptIds);
			
			// get cell information
			bool onelt = false;
			bool onege = false;
			for (int j = 0; j < ptIds->GetNumberOfIds(); j++)
			{
				int id = ptIds->GetId(j);
				
				if (rank[meshid][p] > rank[meshid][id])
					onelt = true;
					
				if (rank[meshid][p] <= rank[meshid][id])
					onege = true;
			}	
			
			if (onelt && onege)
			{
				count++;
			}
		}
		if (count > 1) 
		{
			for (int z = 99; z > rank[meshid][p]; z --)
				rr[z] ++;
		}
	}
	for (int i = 1; i < 10; i++)
		cout << i << ": " << (rr[i]) << "\n";*/
}

float CheckBestOfMesh(int i, float distlimit)
{
	vtkPointData* pointsData = meshes[i]->GetPointData();
	vtkFloatArray* arr = (vtkFloatArray*) pointsData->GetArray("STRTBL");
		
	while (!sortedvertices[i].empty())
	{
		vtkIdType pid = sortedvertices[i].back();
		double* point = meshes[i]->GetPoints()->GetPoint(pid);
		if ((visited[i][pid]) || (distfield->GetValueAt(point[0], point[1], point[2]) < distlimit))
		{
			visited[i][pid] = 1;
			sortedvertices[i].pop_back();
			continue;
		}
		
		vtkIdType id = sortedvertices[i].back();
		return -arr->GetValue(id);
	}
	
	return numeric_limits<float>::min();
}

void FindConnectedComponent(int meshid, vtkIdType i, float distlimit, int sizelimit, int ranklimit, vector<vtkIdType>& comp)
{
	vtkPolyData* vtkMesh = meshes[meshid];
	vtkPointData* pointsData = vtkMesh->GetPointData();
	vtkFloatArray* arr = (vtkFloatArray*) pointsData->GetArray("STRTBL");
	
	//find component
	// read the complete component
	float maxv = numeric_limits<float>::min();
	queue<vtkIdType> toprocess;
	toprocess.push(i);
	comp.push_back(i);
	visited[meshid][i] = 1;
	while (!toprocess.empty())
	{
		// get point
		vtkIdType p = toprocess.front();
		toprocess.pop();
			
		// get point data
		maxv = max(maxv, -arr->GetValue(p));
			
		// add non visited neihbors
		vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
		vtkMesh->GetPointCells(p, cells);
		for (int k = 0; k < cells->GetNumberOfIds(); k++)
		{
			vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
			vtkMesh->GetCellPoints(cells->GetId(k), ptIds);
			for (int j = 0; j < ptIds->GetNumberOfIds(); j++)
			{
				vtkIdType id = ptIds->GetId(j);
				
				double* point = vtkMesh->GetPoints()->GetPoint(id);
				if (distfield->GetValueAt(point[0], point[1], point[2]) < distlimit)
					visited[meshid][id] = 1;
					
				if (!visited[meshid][id])
				{
					// if it has a rank less than ranklimit then it can only add lower ranks
					if ((rank[meshid][p] >= ranklimit) || (rank[meshid][p] > rank[meshid][id]))
					{
						toprocess.push(id);
						comp.push_back(id);
						visited[meshid][id] = 1;
					}
				}
			}
		}
	}
		
	// evaluate conditions to delete component
	if (comp.size() < sizelimit)
	{
		comp.clear();
	}
}

void UpdateDistancefield(int meshid, vector<vtkIdType>& comp)
{
	vtkPolyData* vtkMesh = meshes[meshid];
	char* bs = (char*) malloc(vtkMesh->GetNumberOfPoints());
	memset(bs, 0, vtkMesh->GetNumberOfPoints());
	for (vector<vtkIdType>::iterator it = comp.begin(); it != comp.end(); it++)
		bs[*it] = 1;
	
	// loop on cells and update distance at neighbors
	vtkCellArray* vtkCells = vtkMesh->GetPolys();
	vtkCells->InitTraversal();
	vtkIdType npts;
	vtkIdType* pts;
	while (vtkCells->GetNextCell(npts, pts) != 0)
	{
		// check the cell has vertices in comp
		bool ignore = true;
		for (int i = 0; i < npts; i++)
		{
			if (bs[pts[i]] == 1)
			{
				ignore = false;
				break;
			}
		}
		if (ignore)
			continue;
			
		// process cell	
		double point[3];
		float3 minb = make_float3(numeric_limits<float>::max(),numeric_limits<float>::max(),numeric_limits<float>::max());
		float3 maxb = make_float3(numeric_limits<float>::min(),numeric_limits<float>::min(),numeric_limits<float>::min());
		float3 pt[3];
	
		// find the cell bounding box
		for (int i = 0; i < npts; i++)
		{
			vtkMesh->GetPoint(pts[i], point);
			pt[i] = make_float3((float)point[0], (float)point[1], (float)point[2]);
			
			//printf("%f %f %f\n", pt[i].x, pt[i].y, pt[i].z);
			
			minb.x = min(minb.x, pt[i].x);
			minb.y = min(minb.y, pt[i].y);
			minb.z = min(minb.z, pt[i].z);
			
			maxb.x = max(maxb.x, pt[i].x);
			maxb.y = max(maxb.y, pt[i].y);
			maxb.z = max(maxb.z, pt[i].z);
		}
	
		// check cell is a triangle
		if (npts != 3) 
			printf("Cell is not a triangle and has been skipped!\n");
	
		// find grid bounding box
		minb = distfield->Space2Grid(minb);
		minb.x = floor(minb.x);
		minb.y = floor(minb.y);
		minb.z = floor(minb.z);		
		maxb = distfield->Space2Grid(maxb);
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
					if (!distfield->IsValid(x, y, z)) continue;
					
					float3 spt = distfield->Grid2Space(make_float3(x,y,z));
					int idx = distfield->Coord2Addr(x,y,z);
					float dist = sqrt(DistPoint3Triangle3Squared(pt, spt));
					//#pragma omp critical
					{
						distfield->data[idx] = min(distfield->data[idx], dist);
					}
				}
			}
		}
	}
	free(bs);
	
	// grow distance field
	ComputeDistanceField(distfield, distfield->data);	
}

void WriteField(int i)
{
	// write the output file
	vector< double > spacing;
	spacing.push_back(nrrdfile->axis[0].spacing);
	spacing.push_back(nrrdfile->axis[1].spacing);
	spacing.push_back(nrrdfile->axis[2].spacing);

	vector< size_t > dims;
	dims.push_back(nrrdfile->axis[0].size);
	dims.push_back(nrrdfile->axis[1].size);
	dims.push_back(nrrdfile->axis[2].size);

	char buf[1000];
	sprintf(buf, "%s%d%s", "dfstep", i, ".nrrd");
	string file_string(buf);
	printf("Write '%s'\n", file_string.c_str());
	writeNrrd(distfield->data, file_string, nrrdTypeFloat, dims, spacing);
}

void CreateDistField(int sizelimit, float distlimit, int* ranklimit)
{

	// initialize the distance field
	for (int i = 0; i < distfield->width * distfield->height * distfield->depth; i++)
	{
		distfield->data[i] = distlimit;
	}

	// loop to find the different componenets
	printf("Start search loop for components\n");
	int ccount = 0;
	while(true)
	{
		
		// find the mesh for the best component
		int bestmesh = -1;
		float beststr = numeric_limits<float>::min();
		for (int i = 0; i < mesh_count; i++)
		{
			float besti = CheckBestOfMesh(i, distlimit);
			if (besti > beststr)
			{
				beststr = besti;
				bestmesh = i;
			}
		}
		if (beststr == numeric_limits<float>::min())
			break;
		
		// find the component
		vector<vtkIdType>* comp = new vector<vtkIdType>;
		FindConnectedComponent(bestmesh, sortedvertices[bestmesh].back(), distlimit, sizelimit, ranklimit[bestmesh], *comp);
		
		// update the distance field
		if (comp->size() > 0)
		{
			printf("Mesh %d has maximum strength %f\n", bestmesh, beststr);
			printf("Component of size %d found\n", comp->size());
			
			ccount++;
			UpdateDistancefield(bestmesh, *comp);
			//WriteField(ccount);
		}
		
		delete comp;
		
	}
	
	printf("Found %d components.\n", ccount);
}

int main (int argc, char *argv[])
{
	if (argc < 8)
	{
		printf("Example: vtkmergestructures 3 mesh1.vtk mesh2.vtk mesh3.vtk in.nrrd out.nrrd distl sizel rankl\n");
		return 0;
	}
	
	cout << "==================\n";
	cout << "vtkmergestructures\n";
	cout << "==================\n";
	
	// read meshes
	mesh_count = atoi(argv[1]);
	printf("The number of meshes is %d\n", mesh_count);
	meshes_files = (char**) malloc(mesh_count * sizeof(char*));
	for (int i = 0; i < mesh_count; i++)
		meshes_files[i] = argv[i+2];
	LoadMeshesAndLists();
	
	// compute the rank 
	int* rankl = (int*) malloc(mesh_count * sizeof(int));
	for (int i = 0; i < mesh_count; i++)
	{
		rankl[i] = atoi(argv[mesh_count+6+i]);
		cout << "Rank limit for mesh " << i << " is " << rankl[i] << "\n";
	}
	for (int i = 0; i < mesh_count; i++)
		ComputeRanks(i, rankl[i]);

	// read nrrd file and create grid
	nrrdfile = readNrrd(argv[mesh_count+2]);
	distfield = new RegularGrid(nrrdfile, 
							   nrrdfile->axis[0].size, nrrdfile->axis[1].size, nrrdfile->axis[2].size, 
							   nrrdfile->axis[0].spacing, nrrdfile->axis[1].spacing, nrrdfile->axis[2].spacing);

										   
	// find the structures
	int avepn = 0;
	for (int i = 0; i < mesh_count; i++)
		avepn += meshes[0]->GetNumberOfPoints();
	avepn /= mesh_count;
	float distl = atof(argv[mesh_count+4]);
	int sizel = avepn * atof(argv[mesh_count+5]);
	printf("Distance limit is %f and size limit is %d\n", distl, sizel);
	CreateDistField(sizel, distl, rankl);
											
	// write the output file
	vector< double > spacing;
	spacing.push_back(nrrdfile->axis[0].spacing);
	spacing.push_back(nrrdfile->axis[1].spacing);
	spacing.push_back(nrrdfile->axis[2].spacing);

	vector< size_t > dims;
	dims.push_back(nrrdfile->axis[0].size);
	dims.push_back(nrrdfile->axis[1].size);
	dims.push_back(nrrdfile->axis[2].size);

	string file_string(argv[mesh_count+3]);
	writeNrrd(distfield->data, file_string, nrrdTypeFloat, dims, spacing);

	printf("Write '%s'\n", argv[mesh_count+3]);
	
	return 0;
}