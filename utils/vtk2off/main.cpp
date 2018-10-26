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
#include <vtkExtractEdges.h>

using namespace std;

char* mrk;
int3* cs;

bool ShareTwoVertices(const int3& tri1, const int3& tri2)
{
	int c = 0;
	if ((tri1.x == tri2.x) || (tri1.x == tri2.y) || (tri1.x == tri2.z))
		c++;
		
	if ((tri1.y == tri2.x) || (tri1.y == tri2.y) || (tri1.y == tri2.z))
		c++;
		
	if ((tri1.z == tri2.x) || (tri1.z == tri2.y) || (tri1.z == tri2.z))
		c++;
		
	if (c == 2) return true;
	else return false;
}

bool AreConsistent(int3 tri1, int3 tri2)
{
	
}

/*void AddNeighbors(vtkPolyData* vtkMesh, queue<int2>* mq, int3 tri)
{
	vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New();
	vtkMesh->GetPointCells(tri.x, neighbors);
	for (int k = 0; k < neighbors->GetNumberOfIds(); k++)
	{
		int neighbor = neighbors->GetId(k);
		int3 tric = cs[neighbor];
		
		// check is an edge neighbor
		if (ShareTwoVertices(tri, tric) == false) continue;
		
		// they share an edge
		if (mrk[neighbor] == 0)
			mq->push(neighbor, item.x);
		else
		{
			//if (!ConsistentOrientation(tri, tric))
			{
				printf("Inconsistent orientation\n");
			}
		}
	}	
}*/

int main (int argc, char *argv[])
{
	if (argc != 3)
	{
		printf("Example: vtk2off mesh.vtk mesh.off\n");
		return 0;
	}

	cout << "=======\n";
	cout << "vtk2off\n";
	cout << "=======\n";
	
	printf("Mesh file is %s\n", argv[1]);
	
	// read the mesh file
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader->SetFileName(argv[1]);
    vtkMeshReader->Update();
	vtkPolyData* vtkMesh = vtkMeshReader->GetOutput();
	vtkPoints* points = vtkMesh->GetPoints();
	vtkCellArray* cells = vtkMesh->GetPolys();
	
	// extract edges for edge count
	vtkSmartPointer<vtkExtractEdges> extractEdges = vtkSmartPointer<vtkExtractEdges>::New();
	extractEdges->SetInput(vtkMesh);
	extractEdges->Update();
		
	// write the header information
	fstream ifile;
    ifile.open (argv[2], fstream::out);
	ifile << "OFF\n";
	ifile << "# from vtk file " << argv[1] << "\n";
	ifile << points->GetNumberOfPoints() << " " << cells->GetNumberOfCells() << " " << extractEdges->GetOutput()->GetLines()->GetNumberOfCells() << "\n";
	
	// write vertices to file
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	{
		double point[3];
		points->GetPoint(i, point);
		ifile << point[0] << " " << point[1] << " " << point[2] << "\n";
	}

	// write the cells
	mrk = (char*) malloc(cells->GetNumberOfCells() * sizeof(char));
	memset(mrk, 0, cells->GetNumberOfCells() * sizeof(char));
	cs = (int3*) malloc(cells->GetNumberOfCells() * sizeof(int3));
	vtkIdType npts;
	vtkIdType* pts;
	
	// find all cells
	int j = 0;
	cells->InitTraversal();
	while (cells->GetNextCell(npts, pts) != 0)
	{
		if (npts != 3)
		{
			printf("Not a triangular mesh\n");
			getchar();
		}
		
		cs[j].x = pts[0]; 
		cs[j].y = pts[1]; 
		cs[j].z = pts[2]; 
		
		j++;
	}
	
	// loop fixing direction
	for (int i = 0; i < cells->GetNumberOfCells(); i++)
	{
		// check cell is not done
		if (mrk[i] == 1) continue;
		
		// create queue and grow
		queue<int2> mq;
		mq.push(make_int2(i, -1));
		while (mq.empty() == false)
		{
			// pop item from queue
			int2 item = mq.pop();
			int tri_id1 = item.x;
			int tri_id2 = item.y;
			
			// get triangles that need to match
			int3 tri1 = cs[tri_id1];
			int3 tri2 = cs[tri_id2];
			
			// if already checked skip
			if (mrk[tri_id1]) 
			{
				if (AreConsistent(tri1, tri2) == false)
					printf("Non-orientable surface\n");
				continue;
			}
			mrk[tri_id1] = 1;
						
			// adjust orientation
			AdjustOrientation(tri_id1, tri_id2);
			
			// add neighbors
			AddNeighbors(tri_id1, mq);
		}
	}
	
	
	// write the cells
	/*vtkIdType npts;
	vtkIdType* pts;
	cells->InitTraversal();
	while (cells->GetNextCell(npts, pts))
	{
		if (npts != 3)
		{
			printf("Not a triangular mesh\n");
			getchar();
		}
		ifile << npts << " ";
		for (int i  = 0; i < npts; i++)
			ifile << pts[i] << " ";
		ifile << "\n";
	}
	
	// close file
	ifile.close();*/
	
	return 0;
}