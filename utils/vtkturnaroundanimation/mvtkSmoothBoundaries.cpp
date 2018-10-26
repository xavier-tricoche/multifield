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
#include "mvtkSmoothBoundaries.h"

struct int2_CompareEdge {
    bool operator() (const int2& e1, const int2& e2) const
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
};

void smoothcycle(vtkPoints* points, vector<vtkIdType>* cycle, vector<int3>* tris, int iterations, float weight)
{
	vtkIdType* ptids = (vtkIdType*) malloc(cycle->size() * sizeof(vtkIdType));
	
	int i = 0;
	for (vector<vtkIdType>::iterator it = cycle->begin(); it != cycle->end(); it++)
	{
		ptids[i++] = *it;
	}
	
	for (int k = 0; k < iterations; k++)
	{
		double point1[3];
		double point2[3];
		double point3[3];
		
		for (int i = 0; i < cycle->size(); i++)
		{
		
			points->GetPoint(ptids[i], point1);
			points->GetPoint(ptids[(i + 1) % cycle->size()], point2);
			points->GetPoint(ptids[(i + 2) % cycle->size()], point3);
			
			//printf("before\t%f\t%f\t%f\n", point2[0], point2[1], point2[2]);
			point2[0] = (1 - weight) * point2[0] + weight * (point1[0] + point3[0]) * 0.5;
			point2[1] = (1 - weight) * point2[1] + weight * (point1[1] + point3[1]) * 0.5;
			point2[2] = (1 - weight) * point2[2] + weight * (point1[2] + point3[2]) * 0.5;
			//printf("after\t%f\t%f\t%f\n", point2[0], point2[1], point2[2]);
			
			points->SetPoint(ptids[(i + 1) % cycle->size()], point2);
		}
	}
	
	free(ptids);
}

int mvtkSmoothBoundaries(vtkPolyData* vtkMesh, int iterations, float weight)
{
	printf("===================\n");
	printf("vtksmoothboundaries\n");
	printf("===================\n");
	
	
	// find all cells boundary edges
	vtkPoints* points = vtkMesh->GetPoints();
	vtkCellArray* vtkCells = vtkMesh->GetPolys();
	vtkCells->InitTraversal();
	vtkIdType npts;
	vtkIdType* pts;
	multimap<vtkIdType,vtkIdType> mymultimap;
	int count = 0;
	map<int2, int, int2_CompareEdge> edges;
	while (vtkCells->GetNextCell(npts, pts) != 0)
	{
		int2 edge1 = make_int2(min(pts[0], pts[1]), max(pts[0], pts[1]));
		int2 edge2 = make_int2(min(pts[0], pts[2]), max(pts[0], pts[2]));
		int2 edge3 = make_int2(min(pts[1], pts[2]), max(pts[1], pts[2]));
		
		map<int2, int, int2_CompareEdge>::iterator it;
		
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
	
	for (map<int2, int, int2_CompareEdge>::iterator it = edges.begin(); it != edges.end(); it++ )
	{
		if (it->second != 2)
		{
			count++;
			mymultimap.insert(pair<vtkIdType,vtkIdType>((it->first).x, (it->first).y));
			mymultimap.insert(pair<vtkIdType,vtkIdType>((it->first).y, (it->first).x));
			
			//if (it->second != 1)
				//cout << "edge " << (it->first).x << "," << (it->first).y << " has " << it->second << " triangles!\n";
		}
	}
	
	cout << "Number of edges is " << count << "\n";
	
	// find all closed rings
	int ncycles = 0;
	multimap<vtkIdType,vtkIdType>::iterator it;
	vector<int3> tris;
	while (!mymultimap.empty())
	{
		it = mymultimap.begin();
		vector<vtkIdType> cycle;
		
		// find chain
		vtkIdType f = (*it).first;	
		vtkIdType s = (*it).second;
		cycle.push_back(f);	
		cycle.push_back(s);	
		mymultimap.erase(it);		
		while ((it = mymultimap.find(s)) != mymultimap.end())
		{
			if (((*it).second) == f)
			{
				mymultimap.erase(it);
				continue;
			}
				
			f = (*it).first;			
			s = (*it).second;
			cycle.push_back(s);
						
			mymultimap.erase(it);
		}
		
		// check cycle
		if ((cycle.size() > 0) && (cycle.front() == cycle.back()))
		{
			ncycles++;
			
			smoothcycle(points, &cycle, &tris, iterations, weight);
		}
	}
	cout << "Number of cycles is " << ncycles << "\n";
	
	return 0;
}