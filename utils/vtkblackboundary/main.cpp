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
#include <algorithm>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <set>
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
#include <vtkCellArray.h>
#include <vtkPolyDataReader.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkPolyDataWriter.h>
#include <vtkLine.h>

using namespace std;

class Edge 
{
public:
	int x;
	int y;
	double length;

	Edge() {}
	Edge(int _x, int _y) : x(_x), y(_y) {}

	bool operator<(const Edge& e2) const
	{
		if (x < e2.x)
		{
			return true;
		}
		else if (x > e2.x)
		{
			return false;
		}
		else
		{
			if (y < e2.y)
			{
				return true;
			}
			else if (y > e2.y)
			{
				return false;
			}
			else
			{
				return false;
			}
		}
	}

	bool operator==(const Edge e1)
	{
		if ((e1.x == x) && (e1.y == y))
			return true;
		return false;
	}

};

double GetEdgeLength(vtkPolyData* vtkMesh, int id1, int id2)
{
	double pt1[3];
	vtkMesh->GetPoint(id1, pt1);

	double pt2[3];
	vtkMesh->GetPoint(id2, pt2);

	pt1[0] -= pt2[0];
	pt1[1] -= pt2[1];
	pt1[2] -= pt2[2];

	return sqrt(pt1[0] * pt1[0] + pt1[1] * pt1[1] + pt1[2] * pt1[2]);
}

double GetEdgeLength(vtkPolyData* vtkMesh, Edge edge)
{
	double pt1[3];
	vtkMesh->GetPoint(edge.x, pt1);

	double pt2[3];
	vtkMesh->GetPoint(edge.y, pt2);

	pt1[0] -= pt2[0];
	pt1[1] -= pt2[1];
	pt1[2] -= pt2[2];

	return sqrt(pt1[0] * pt1[0] + pt1[1] * pt1[1] + pt1[2] * pt1[2]);
}

int main (int argc, char *argv[])
{
	if (argc != 3)
	{
		printf("Example: vtkblackboundary mesh.vtk out.vtk\n");
		return 0;
	}

	cout << "================\n";
	cout << "vtkblackboundary\n";
	cout << "================\n";
	
	printf("Mesh file is %s\n", argv[1]);
	
	// read the mesh file
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader->SetFileName(argv[1]);
    vtkMeshReader->Update();
	vtkPolyData* vtkMesh = vtkMeshReader->GetOutput();
	
	// get data array 
	vtkPointData* pointsData = vtkMesh->GetPointData();	
	
	// find all cells boundary edges
	vtkCellArray* vtkcells = vtkMesh->GetPolys();
	vtkcells->InitTraversal();
	vtkIdType npts;
	vtkIdType* pts;
	int count = 0;
	map<Edge, int> edges;
	double avelength = 0.0;
	int avecount = 0;
	while (vtkcells->GetNextCell(npts, pts) != 0)
	{
		Edge edge1(min(pts[0], pts[1]), max(pts[0], pts[1]));
		Edge edge2(min(pts[0], pts[2]), max(pts[0], pts[2]));
		Edge edge3(min(pts[1], pts[2]), max(pts[1], pts[2]));
		
		map<Edge, int>::iterator it;
		
		if ((it = edges.find(edge1)) != edges.end())
		{
			int scr = it->second + 1;
			edges.erase(it);
			edges.insert(pair<Edge, int>(edge1, scr));
		}
		else
		{
			edge1.length = GetEdgeLength(vtkMesh, edge1);
			edges.insert(pair<Edge, int>(edge1, 1));
		}
			
		if ((it = edges.find(edge2)) != edges.end())
		{	
			int scr = it->second + 1;
			edges.erase(it);
			edges.insert(pair<Edge, int>(edge2, scr));
		}
		else
		{
			edge2.length = GetEdgeLength(vtkMesh, edge2);
			edges.insert(pair<Edge, int>(edge2, 1));
		}
			
		if ((it = edges.find(edge3)) != edges.end())
		{	
			int scr = it->second + 1;
			edges.erase(it);
			edges.insert(pair<Edge, int>(edge3, scr));
		}
		else
		{
			edge3.length = GetEdgeLength(vtkMesh, edge3);
			edges.insert(pair<Edge, int>(edge3, 1));
		}
		
		avelength += edge1.length;
		avelength += edge2.length;
		avelength += edge3.length;
		avecount += 3;
	}
	avelength /= avecount;
	cout << "Average edge length is " << avelength << "\n";
	
	// get color array
	vtkUnsignedCharArray* colors = (vtkUnsignedCharArray*) pointsData->GetScalars();
	if (colors == NULL)
	{
		printf("No color array found!\n");
		return 0;
	}
	
	// duplicate boundary points
	vtkPoints* points = vtkMesh->GetPoints();
	map<int,int> old2new;
	for (map<Edge, int>::iterator it = edges.begin(); it != edges.end(); it++ )
	{
		if (it->second != 2)
		{
			vtkIdType newid1;
			vtkIdType newid2;
			
			if (old2new.find((*it).first.x) == old2new.end())
			{
				double point[3];
				points->GetPoint((*it).first.x, point);
				newid1 = points->InsertNextPoint(point[0]+2 * avelength, point[1], point[2]);
				old2new[(*it).first.x] = newid1;
				
				double old_color[4];
				colors->GetTuple((*it).first.x, old_color);
				old_color[0] = 0;
				old_color[1] = 0;
				old_color[2] = 0;
				colors->InsertTuple(newid1, old_color);
				
				colors->InsertTuple((*it).first.x, old_color);
			}
			else
				newid1 = old2new[(*it).first.x];
			
			if (old2new.find((*it).first.y) == old2new.end())
			{
				double point[3];
				points->GetPoint((*it).first.y, point);
				newid2 = points->InsertNextPoint(point[0]+2 * avelength, point[1], point[2]);
				old2new[(*it).first.y] = newid2;
				
				double old_color[4];
				colors->GetTuple((*it).first.y, old_color);
				old_color[0] = 0;
				old_color[1] = 0;
				old_color[2] = 0;
				colors->InsertTuple(newid2, old_color);
				
				colors->InsertTuple((*it).first.y, old_color);
			}
			else
				newid2 = old2new[(*it).first.y];
				
			// insert line cell 
			/*vtkIdList* newline = vtkIdList::New();
			newline->Reset();
			newline->SetNumberOfIds(2);
			newline->SetId(0,newid1);
			newline->SetId(1,newid2);*/
			vtkIdType pts[3];
			pts[0] = newid1;
			pts[1] = newid2;
			pts[2] = (*it).first.x;
			vtkMesh->GetPolys()->InsertNextCell(3, pts);
			pts[2] = (*it).first.y;
			vtkMesh->GetPolys()->InsertNextCell(3, pts);
		}
	}
	
	vtkSmartPointer<vtkPolyData> vtkMeshNew = vtkSmartPointer<vtkPolyData>::New();
	vtkMeshNew->SetPoints(vtkMesh->GetPoints());
	vtkMeshNew->SetPolys(vtkMesh->GetPolys());
	vtkMeshNew->GetPointData()->SetScalars(colors);
	
	
	//vtkMesh->SetPoints(points);
	
	// collect boundary points
	/*double* old_color;
	double new_color[4] = {0.0, 0.0, 0.0, 255.0};
	for (map<Edge, int>::iterator it = edges.begin(); it != edges.end(); it++ )
	{
		if (it->second != 2)
		{
			old_color = colors->GetTuple((*it).first.x);
			new_color[3] = old_color[3];
			colors->SetTuple((*it).first.x, new_color);
			
			old_color = colors->GetTuple((*it).first.y);
			new_color[3] = old_color[3];
			colors->SetTuple((*it).first.y, new_color);
		}
	}*/

	// write file
	vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkMeshWriter->SetFileName(argv[2]);
	vtkMeshWriter->SetInput(vtkMeshNew);
	vtkMeshWriter->Write();
	
	return 0;
}