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
#include <vtkPolyDataReader.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkPolyDataWriter.h>
#include <vtkIdList.h>
#include <vtkColorTransferFunction.h>
#include <vtkClipPolyData.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkCellArray.h>
#include "vtkAppendPolyData.h"
#include <vtkPlane.h>

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

//  public-domain function by Darel Rex Finley
	//
	//  The passed-in RGB values can be on any desired scale, such as 0 to
	//  to 1, or 0 to 255.  (But use the same scale for all three!)
	//
	//  The "change" parameter works like this:
	//    0.0 creates a black-and-white image.
	//    0.5 reduces the color saturation by half.
	//    1.0 causes no change.
	//    2.0 doubles the color saturation.
	//  Note:  A "change" value greater than 1.0 may project your RGB values
	//  beyond their normal range, in which case you probably should truncate
	//  them to the desired range before trying to use them in an image.
	
	#define  sPr  .299
	#define  sPg  .587
	#define  sPb  .114
	
bool ReduceSaturation(float4& c)
{
	double change = 0.35;
	double P=sqrt(c.x*c.x*sPr + c.y*c.y*sPg + c.z*c.z*sPb ) ;

	c.x=P+(c.x-P)*change;
	c.y=P+(c.y-P)*change;
	c.z=P+(c.z-P)*change;
}
	
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

void BlackBoundary(vtkPolyData* vtkMesh, vtkUnsignedCharArray* colors)
{
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
	
	// duplicate boundary points
	vtkPoints* points = vtkMesh->GetPoints();
	map<int,int> old2new;
	for (map<Edge, int>::iterator it = edges.begin(); it != edges.end(); it++ )
	{
		if (it->second != 2)
		{
			//vtkIdType newid1;
			//vtkIdType newid2;
			
			//if (old2new.find((*it).first.x) == old2new.end())
			{
				//double point[3];
				//points->GetPoint((*it).first.x, point);
				//newid1 = points->InsertNextPoint(point[0]+2 * avelength, point[1], point[2]);
				//old2new[(*it).first.x] = newid1;
				
				double old_color[4];
				colors->GetTuple((*it).first.x, old_color);
				old_color[0] = 0;
				old_color[1] = 0;
				old_color[2] = 0;
				//colors->InsertTuple(newid1, old_color);
				
				colors->InsertTuple((*it).first.x, old_color);
			}
			//else
			//	newid1 = old2new[(*it).first.x];
			
			//if (old2new.find((*it).first.y) == old2new.end())
			{
				//double point[3];
				//points->GetPoint((*it).first.y, point);
				//newid2 = points->InsertNextPoint(point[0]+2 * avelength, point[1], point[2]);
				//old2new[(*it).first.y] = newid2;
				
				double old_color[4];
				colors->GetTuple((*it).first.y, old_color);
				old_color[0] = 0;
				old_color[1] = 0;
				old_color[2] = 0;
				//colors->InsertTuple(newid2, old_color);
				
				colors->InsertTuple((*it).first.y, old_color);
			}
			//else
			//	newid2 = old2new[(*it).first.y];
		}
	}
}

float3 ValueToColor(double val)
{
	vtkColorTransferFunction* colorMap = vtkColorTransferFunction::New();
	
	double vs[9] = {5.00E-01,4.38E-01,3.75E-01,3.13E-01,2.50E-01,1.88E-01,1.25E-01,6.25E-02,0.00E+00};
	
	// color map at http://www.triptropnyc.com/
	/*colorMap->AddRGBPoint(vs[0], 247/ 255.0, 131/ 255.0, 125/ 255.0);
	colorMap->AddRGBPoint(vs[1], 254/ 255.0, 177/ 255.0, 126/ 255.0);
	colorMap->AddRGBPoint(vs[2], 252/ 255.0, 224/ 255.0, 126/ 255.0);
	colorMap->AddRGBPoint(vs[3], 213/ 255.0, 243/ 255.0, 130/ 255.0);
	colorMap->AddRGBPoint(vs[4], 172/ 255.0, 242/ 255.0, 167/ 255.0);
	colorMap->AddRGBPoint(vs[5], 141/ 255.0, 255/ 255.0, 227/ 255.0);
	colorMap->AddRGBPoint(vs[6], 128/ 255.0, 212/ 255.0, 255/ 255.0);
	colorMap->AddRGBPoint(vs[7], 131/ 255.0, 157/ 255.0, 255/ 255.0);
	colorMap->AddRGBPoint(vs[8], 132/ 255.0, 118/ 255.0, 248/ 255.0);
	colorMap->AddRGBPoint(vs[9], 129/ 255.0, 122/ 255.0, 189/ 255.0);*/
	
	
	// color map from http://colorbrewer2.com/ yellow to blue
	colorMap->AddRGBPoint(vs[8], 255/ 255.0, 255/ 255.0, 217/ 255.0);
	colorMap->AddRGBPoint(vs[7], 237/ 255.0, 248/ 255.0, 217/ 255.0);
	colorMap->AddRGBPoint(vs[6], 199/ 255.0, 233/ 255.0, 180/ 255.0);
	colorMap->AddRGBPoint(vs[5], 127/ 255.0, 205/ 255.0, 187/ 255.0);
	colorMap->AddRGBPoint(vs[4], 65/ 255.0, 182/ 255.0, 196/ 255.0);
	colorMap->AddRGBPoint(vs[3], 29/ 255.0, 145/ 255.0, 192/ 255.0);
	colorMap->AddRGBPoint(vs[2], 34/ 255.0, 94/ 255.0, 168/ 255.0);
	colorMap->AddRGBPoint(vs[1], 37/ 255.0, 52/ 255.0, 148/ 255.0);
	colorMap->AddRGBPoint(vs[0], 8/ 255.0, 29/ 255.0, 88/ 255.0);
	
	
	// color map from http://colorbrewer2.com/ magneta
	/*colorMap->AddRGBPoint(vs[8], 255/ 255.0, 247/ 255.0, 243/ 255.0);
	colorMap->AddRGBPoint(vs[7], 253/ 255.0, 224/ 255.0, 221/ 255.0);
	colorMap->AddRGBPoint(vs[6], 252/ 255.0, 197/ 255.0, 192/ 255.0);
	colorMap->AddRGBPoint(vs[5], 250/ 255.0, 159/ 255.0, 181/ 255.0);
	colorMap->AddRGBPoint(vs[4], 247/ 255.0, 104/ 255.0, 161/ 255.0);
	colorMap->AddRGBPoint(vs[3], 221/ 255.0, 52/ 255.0, 151/ 255.0);
	colorMap->AddRGBPoint(vs[2], 174/ 255.0, 1/ 255.0, 126/ 255.0);
	colorMap->AddRGBPoint(vs[1], 122/ 255.0, 1/ 255.0, 119/ 255.0);
	colorMap->AddRGBPoint(vs[0], 73/ 255.0, 0/ 255.0, 106/ 255.0);*/
	
	// color map from http://colorbrewer2.com/ heat
	/*colorMap->AddRGBPoint(vs[8], 255/ 255.0, 255/ 255.0, 229/ 255.0);
	colorMap->AddRGBPoint(vs[7], 255/ 255.0, 247/ 255.0, 188/ 255.0);
	colorMap->AddRGBPoint(vs[6], 254/ 255.0, 227/ 255.0, 145/ 255.0);
	colorMap->AddRGBPoint(vs[5], 254/ 255.0, 196/ 255.0, 79/ 255.0);
	colorMap->AddRGBPoint(vs[4], 254/ 255.0, 153/ 255.0, 41/ 255.0);
	colorMap->AddRGBPoint(vs[3], 236/ 255.0, 112/ 255.0, 20/ 255.0);
	colorMap->AddRGBPoint(vs[2], 204/ 255.0, 76/ 255.0, 2/ 255.0);
	colorMap->AddRGBPoint(vs[1], 153/ 255.0, 52/ 255.0, 4/ 255.0);
	colorMap->AddRGBPoint(vs[0], 102/ 255.0, 37/ 255.0, 6/ 255.0);*/
	
	
	// center map 10 from http://jpgraph.net/download/manuals/chunkhtml/ch22s08.html
	/*colorMap->AddRGBPoint(vs[], 179/ 255.0, 88/ 255.0, 6/ 255.0);
	colorMap->AddRGBPoint(vs[], 224/ 255.0, 130/ 255.0, 20/ 255.0);
	colorMap->AddRGBPoint(vs[], 253/ 255.0, 184/ 255.0, 99/ 255.0);
	colorMap->AddRGBPoint(vs[], 254/ 255.0, 224/ 255.0, 182/ 255.0);
	colorMap->AddRGBPoint(vs[], 255/ 255.0, 255/ 255.0, 255/ 255.0);
	colorMap->AddRGBPoint(vs[], 216/ 255.0, 218/ 255.0, 235/ 255.0);
	colorMap->AddRGBPoint(vs[], 178/ 255.0, 171/ 255.0, 210/ 255.0);
	colorMap->AddRGBPoint(vs[], 128/ 255.0, 115/ 255.0, 172/ 255.0);
	colorMap->AddRGBPoint(vs[], 84./ 255.0, 39./ 255.0, 136/ 255.0);*/
	
	// heat from http://princeofslides.blogspot.com/2011/10/maximizing-sabermetric-visual-content.html
	/*colorMap->AddRGBPoint(vs[0], 165/ 255.0, 0  / 255.0, 38 / 255.0);
	colorMap->AddRGBPoint(vs[1], 215/ 255.0, 49 / 255.0, 39 / 255.0);
	colorMap->AddRGBPoint(vs[2], 244/ 255.0, 111/ 255.0, 68 / 255.0);
	colorMap->AddRGBPoint(vs[3], 253/ 255.0, 177/ 255.0, 99 / 255.0);
	colorMap->AddRGBPoint(vs[4], 254/ 255.0, 226/ 255.0, 147/ 255.0);
	colorMap->AddRGBPoint(vs[5], 251/ 255.0, 253/ 255.0, 196/ 255.0);
	colorMap->AddRGBPoint(vs[6], 217/ 255.0, 239/ 255.0, 246/ 255.0);
	colorMap->AddRGBPoint(vs[7], 163/ 255.0, 210/ 255.0, 229/ 255.0);
	colorMap->AddRGBPoint(vs[8], 108/ 255.0, 163/ 255.0, 204/ 255.0);
	colorMap->AddRGBPoint(vs[9], 65 / 255.0, 105/ 255.0, 174/ 255.0);
	colorMap->AddRGBPoint(vs[10], 49 / 255.0, 54 / 255.0, 149/ 255.0);*/
	
	// color map from http://www.kgs.ku.edu/General/elevatMap.html
	/*colorMap->AddRGBPoint(vs[], 70/ 255.0, 43 / 255.0, 18/ 255.0);
	colorMap->AddRGBPoint(vs[], 131/ 255.0, 82 / 255.0, 55/ 255.0);
	colorMap->AddRGBPoint(vs[], 187/ 255.0, 113 / 255.0, 45/ 255.0);
	colorMap->AddRGBPoint(vs[], 228/ 255.0, 157 / 255.0, 76/ 255.0);
	colorMap->AddRGBPoint(vs[], 236/ 255.0, 185 / 255.0, 109/ 255.0);
	colorMap->AddRGBPoint(vs[], 252/ 255.0, 226 / 255.0, 156/ 255.0);
	colorMap->AddRGBPoint(vs[], 245/ 255.0, 247 / 255.0, 159/ 255.0);
	colorMap->AddRGBPoint(vs[], 200/ 255.0, 211 / 255.0, 116/ 255.0);
	colorMap->AddRGBPoint(vs[], 152/ 255.0, 182 / 255.0, 24/ 255.0);
	colorMap->AddRGBPoint(vs[], 108/ 255.0, 140 / 255.0, 41/ 255.0);
	colorMap->AddRGBPoint(vs[], 96/ 255.0, 116 / 255.0, 50/ 255.0);
	colorMap->AddRGBPoint(vs[], 3/ 255.0, 94 / 255.0, 1/ 255.0);
	colorMap->AddRGBPoint(vs[], 0/ 255.0, 73 / 255.0, 1/ 255.0);
	colorMap->AddRGBPoint(vs[], 22/ 255.0, 88 / 255.0, 52/ 255.0);
	colorMap->AddRGBPoint(vs[], 64/ 255.0, 122 / 255.0, 90/ 255.0);
	colorMap->AddRGBPoint(vs[], 30/ 255.0, 163 / 255.0, 107/ 255.0);
	colorMap->AddRGBPoint(vs[], 95/ 255.0, 201 / 255.0, 165/ 255.0);
	colorMap->AddRGBPoint(vs[], 117/ 255.0, 233 / 255.0, 204/ 255.0);
	colorMap->AddRGBPoint(vs[], 74/ 255.0, 188 / 255.0, 179/ 255.0);
	colorMap->AddRGBPoint(vs[], 22/ 255.0, 139 / 255.0, 154/ 255.0);
	colorMap->AddRGBPoint(vs[], 15/ 255.0, 90 / 255.0, 142/ 255.0);
	colorMap->AddRGBPoint(vs[], 0/ 255.0, 6 / 255.0, 126/ 255.0);
	colorMap->AddRGBPoint(vs[], 17/ 255.0, 3 / 255.0, 89/ 255.0);*/
	
	
	double color[3];
	colorMap->GetColor(val, color);
	return make_float3(color[0], color[1], color[2]);
}

void CutMeshBasedOnValue(vtkPolyData*& vtkMesh, char* table, double minv, double maxv)
{
	vtkMesh->GetPointData()->SetActiveScalars(table);
	
	vtkClipPolyData* clipper3 = vtkClipPolyData::New();
	clipper3->SetInput(vtkMesh);
	clipper3->SetValue(minv);
	clipper3->InsideOutOff();
	clipper3->Update();
	
	vtkClipPolyData* clipper4 = vtkClipPolyData::New();
	clipper4->SetInput(clipper3->GetOutput());
	clipper4->SetValue(maxv);
	clipper4->InsideOutOn();
	clipper4->Update();
	
	vtkMesh = clipper4->GetOutput();
}

set<int>* GetNeighbors(vtkPolyData* vtkMesh, int id)
{
	set<int>* nei = new set<int>;
	vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
	vtkMesh->GetPointCells(id, cells);
	
	for (int i = 0; i < cells->GetNumberOfIds(); i++)
	{
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
		vtkMesh->GetCellPoints(cells->GetId(i), ptIds);
		for (int j = 0; j < ptIds->GetNumberOfIds(); j++)
		{
			int nid = ptIds->GetId(j);
			nei->insert(nid);
		}
	}
	
	nei->erase(id);
	
	return nei;
}

void SmoothColor(vtkPolyData* vtkMesh, vtkUnsignedCharArray* colors, int iterations, float weight)
{
	vtkUnsignedCharArray* newcolors = vtkUnsignedCharArray::New();
	newcolors->SetNumberOfComponents(4);
	newcolors->SetNumberOfTuples(vtkMesh->GetNumberOfPoints());
	
	for (int k = 0; k < iterations; k++)
	{
		// interpolate with neighbors
		for (int i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
		{
			// get all neighbors
			set<int>* nei = GetNeighbors(vtkMesh, i);
			float sum=0.0;
			double newcolor[4];
			newcolor[0] = 0.0;
			newcolor[1] = 0.0;
			newcolor[2] = 0.0;
			newcolor[3] = 0.0;
			for (set<int>::iterator it = nei->begin(); it != nei->end(); it++)
			{
				double point[3];
				vtkMesh->GetPoint(*it, point);
				float3 ptn = make_float3(point[0], point[1], point[2]);
				vtkMesh->GetPoint(i, point);
				float3 pt = make_float3(point[0], point[1], point[2]);
				sum += length(pt - ptn);
				
				double color[4];
				colors->GetTuple(*it, color);
				newcolor[0] += color[0] * length(pt - ptn);
				newcolor[1] += color[1] * length(pt - ptn);
				newcolor[2] += color[2] * length(pt - ptn);
				newcolor[3] += color[3] * length(pt - ptn);
			}
			newcolor[0] /= sum;
			newcolor[1] /= sum;
			newcolor[2] /= sum;
			newcolor[3] /= sum;
			double color[4];
			colors->GetTuple(i, color);
			newcolor[0] = (1.0 - weight) * color[0] + weight * newcolor[0];
			newcolor[1] = (1.0 - weight) * color[1] + weight * newcolor[1];
			newcolor[2] = (1.0 - weight) * color[2] + weight * newcolor[2];
			newcolor[3] = (1.0 - weight) * color[3] + weight * newcolor[3];
			newcolors->SetTuple(i, newcolor);
			
			delete nei;
		}
		
		
		// copy newcolors to colors
		for (int i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
		{
			double color[4];
			newcolors->GetTuple(i, color);
			colors->SetTuple(i, color);
		}
	}
}






enum COLORMAPS {SATURATED_RED_TO_MAGNETA, 
				YELLOW_TO_BLUE, 
				WHITE_TO_MAGNETA, 
				HEAT_WHITE_TO_DARKRED, 
				CENTER_ORANGE_TO_MAGNETA, 
				RED_TO_BLUE, 
				MARON_TO_BLACK,
				MARON_TO_GREEN,
				RED_TO_GREEN};
vtkColorTransferFunction* GetValueColorTable(int map, bool flipmap, float minv, float maxv)
{
	vtkColorTransferFunction* colorMap = vtkColorTransferFunction::New();
	
	vector<float3> mapcol;

	switch (map)
	{
	case SATURATED_RED_TO_MAGNETA:
		// color map at http://www.triptropnyc.com/
		mapcol.push_back(make_float3(248,131,127));
		mapcol.push_back(make_float3(255,179,127));
		mapcol.push_back(make_float3(255,255,127));
		mapcol.push_back(make_float3(220,243,130));
		mapcol.push_back(make_float3(180,243,171));
		mapcol.push_back(make_float3(148,255,230));
		mapcol.push_back(make_float3(127,217,255));
		mapcol.push_back(make_float3(127,166,255));
		mapcol.push_back(make_float3(127,127,250));
		mapcol.push_back(make_float3(127,127,191));
	break;
	
	case YELLOW_TO_BLUE:
		// color map from http://colorbrewer2.com/ yellow to blue
		//mapcol.push_back(make_float3(255,255,217));
		//mapcol.push_back(make_float3(237,248,217));
		mapcol.push_back(make_float3(199,233,180));
		mapcol.push_back(make_float3(127,205,187));
		mapcol.push_back(make_float3(65,182,196));
		mapcol.push_back(make_float3(29,145,192));
		mapcol.push_back(make_float3(34,94,168));
		mapcol.push_back(make_float3(37,52,148));
		mapcol.push_back(make_float3(8,29,88));
	break;
	
	case WHITE_TO_MAGNETA:
		// color map from http://colorbrewer2.com/ magneta
		mapcol.push_back(make_float3(255, 247, 243));
		mapcol.push_back(make_float3(253, 224, 221));
		mapcol.push_back(make_float3(252, 197, 192));
		mapcol.push_back(make_float3(250, 159, 181));
		mapcol.push_back(make_float3(247, 104, 161));
		mapcol.push_back(make_float3(221, 52, 151));
		mapcol.push_back(make_float3(174, 1, 126));
		mapcol.push_back(make_float3(122, 1, 119));
		mapcol.push_back(make_float3(73, 0, 106));
	break;

	case HEAT_WHITE_TO_DARKRED:	
		// color map from http://colorbrewer2.com/ heat
		//mapcol.push_back(make_float3(255, 255, 229));
		//mapcol.push_back(make_float3(255, 247, 188));
		mapcol.push_back(make_float3(254, 227, 145));
		mapcol.push_back(make_float3(254, 196, 79));
		mapcol.push_back(make_float3(254, 153, 41));
		mapcol.push_back(make_float3(236, 112, 20));
		mapcol.push_back(make_float3(204, 76, 2));
		mapcol.push_back(make_float3(153, 52, 4));
		mapcol.push_back(make_float3(102, 37, 6));
	break;
	
	case CENTER_ORANGE_TO_MAGNETA:
		// center map 10 from http://jpgraph.net/download/manuals/chunkhtml/ch22s08.html
		mapcol.push_back(make_float3(179, 88, 6));
		mapcol.push_back(make_float3(224, 130, 20));
		mapcol.push_back(make_float3(253, 184, 99));
		mapcol.push_back(make_float3(254, 224, 182));
		mapcol.push_back(make_float3(255, 255, 255));
		mapcol.push_back(make_float3(216, 218, 235));
		mapcol.push_back(make_float3(178, 171, 210));
		mapcol.push_back(make_float3(128, 115, 172));
		mapcol.push_back(make_float3(84, 39, 136));
	break;
	
	case RED_TO_BLUE:
		// heat from http://princeofslides.blogspot.com/2011/10/maximizing-sabermetric-visual-content.html
		mapcol.push_back(make_float3(165, 0  , 38 ));
		mapcol.push_back(make_float3(215, 49 , 39 ));
		mapcol.push_back(make_float3(244, 111, 68 ));
		mapcol.push_back(make_float3(253, 177, 99 ));
		mapcol.push_back(make_float3(254, 226, 147));
		mapcol.push_back(make_float3(251, 253, 196));
		mapcol.push_back(make_float3(217, 239, 246));
		mapcol.push_back(make_float3(163, 210, 229));
		mapcol.push_back(make_float3(108, 163, 204));
		mapcol.push_back(make_float3(65 , 105, 174));
		mapcol.push_back(make_float3(49 , 54 , 149));
	break;
	
	case MARON_TO_BLACK:
		// color map from http://www.kgs.ku.edu/General/elevatMap.html
		mapcol.push_back(make_float3(70, 43 , 18));
		mapcol.push_back(make_float3(131, 82 , 55));
		mapcol.push_back(make_float3(187, 113 , 45));
		mapcol.push_back(make_float3(228, 157 , 76));
		mapcol.push_back(make_float3(236, 185 , 109));
		mapcol.push_back(make_float3(252, 226 , 156));
		mapcol.push_back(make_float3(245, 247 , 159));
		mapcol.push_back(make_float3(200, 211 , 116));
		mapcol.push_back(make_float3(152, 182 , 24));
		mapcol.push_back(make_float3(108, 140 , 41));
		mapcol.push_back(make_float3(96, 116 , 50));
		mapcol.push_back(make_float3(3, 94 , 1));
		mapcol.push_back(make_float3(0, 73 , 1));
		mapcol.push_back(make_float3(22, 88 , 52));
		mapcol.push_back(make_float3(64, 122 , 90));
		mapcol.push_back(make_float3(30, 163 , 107));
		mapcol.push_back(make_float3(95, 201 , 165));
		mapcol.push_back(make_float3(117, 233 , 204));
		mapcol.push_back(make_float3(74, 188 , 179));
		mapcol.push_back(make_float3(22, 139 , 154));
		mapcol.push_back(make_float3(15, 90 , 142));
		mapcol.push_back(make_float3(0, 6 , 126));
		mapcol.push_back(make_float3(17, 3 , 89));
	break;
	
	case MARON_TO_GREEN:
		// color map from http://www.kgs.ku.edu/General/elevatMap.html
		mapcol.push_back(make_float3(70, 43 , 18));
		mapcol.push_back(make_float3(131, 82 , 55));
		mapcol.push_back(make_float3(187, 113 , 45));
		mapcol.push_back(make_float3(228, 157 , 76));
		mapcol.push_back(make_float3(236, 185 , 109));
		mapcol.push_back(make_float3(252, 226 , 156));
		mapcol.push_back(make_float3(245, 247 , 159));
		mapcol.push_back(make_float3(200, 211 , 116));
		mapcol.push_back(make_float3(152, 182 , 24));
		mapcol.push_back(make_float3(108, 140 , 41));
		mapcol.push_back(make_float3(96, 116 , 50));
		mapcol.push_back(make_float3(3, 94 , 1));
		mapcol.push_back(make_float3(0, 73 , 1));
	break;
	
	case RED_TO_GREEN:
		// color map from http://colorbrewer2.com/ yellow to blue
		mapcol.push_back(make_float3(215,48,39));
		mapcol.push_back(make_float3(244,109,67));
		mapcol.push_back(make_float3(253,174,97));
		mapcol.push_back(make_float3(254,224,139));
		mapcol.push_back(make_float3(255,255,191));
		mapcol.push_back(make_float3(217,239,139));
		mapcol.push_back(make_float3(166,217,106));
		mapcol.push_back(make_float3(102,189,99));
		mapcol.push_back(make_float3(26,152,80));
	break;
	
	}
	
	if (flipmap) 
		reverse(mapcol.begin(),mapcol.end()); 


	for (int i = 0; i < mapcol.size(); i++)
	{
		float3 color = mapcol[i] / 255.0;
		colorMap->AddRGBPoint(minv + (maxv - minv) * i / (mapcol.size() - 1), color.x, color.y, color.z);
	}

	return colorMap;
}

float3 GetValueColor(vtkColorTransferFunction* colorMap, double val)
{
	double color[3];
	colorMap->GetColor(val, color);
	return make_float3(color[0], color[1], color[2]);
}








void MakeQueryTableFromTags(vtkPolyData* vtkMesh)
{
	vtkSmartPointer<vtkFloatArray> newScalars = vtkSmartPointer<vtkFloatArray>::New();
	newScalars->SetName("QUERYTBL");
	
	vtkIntArray* tagarr = (vtkIntArray*) vtkMesh->GetPointData()->GetArray("TAGTBL");
	for (vtkIdType i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
	{
		float value;
		int tag = tagarr->GetValue(i);
		if (tag & 12)
		{
			value = 1.0;
		}
		else
		{
			value = 0.0;
		}
					
		newScalars->InsertValue(i, value);
	}
	vtkMesh->GetPointData()->AddArray(newScalars);
}




void AssignPerClusterColor(vtkPolyData* vtkMesh, vtkUnsignedCharArray* colors)
{
	vtkFloatArray* valarr = (vtkFloatArray*) vtkMesh->GetPointData()->GetArray("VALTBL");
	
	
	vtkIntArray* clsarr = (vtkIntArray*) vtkMesh->GetPointData()->GetArray("CLUSTBL");
	vtkIntArray* tagarr = (vtkIntArray*) vtkMesh->GetPointData()->GetArray("TAGTBL");
	
	// find all tags
	set<int> tags;
	for (vtkIdType i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
		tags.insert(clsarr->GetValue(i));
	printf("Number of distinct clusters is %d\n", tags.size());	
		
	// initialize counts to zero
	vector<int> count;
	for (set<int>::iterator it = tags.begin(); it != tags.end(); it++)
		count.push_back(0);
		
	// now assign colors
	for (vtkIdType i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
	{
		int tag = clsarr->GetValue(i);
		
		float3 color = make_float3(255,255,255);
		
		// for multifield
		switch(tag % 8)
		{
			case 0:
				color = make_float3(144,109,183);
				break;
				
			case 1:
				color = make_float3(146,183,108);
				break;
				
			case 2:
				color = make_float3(107,114,183);
				break;
				
			case 3:
				color = make_float3(183,107,107);
				break;
				
			case 4:
				color = make_float3(183,183,109);
				break;
				
			case 5:
				color = make_float3(108,183,112);
				break;
				
			case 6:
				color = make_float3(183,108,181);
				break;
				
			case 7:
				color = make_float3(108,183,179);
				break;
		};
		
		
		
		//for combustion
	/*	switch(tag % 6)
		{
			case 0:
				color = make_float3(183,107,107);
				break;
				
			case 1:
				color = make_float3(183,177,108);
				break;
				
			case 2:
				color = make_float3(107,183,114);
				break;
				
			case 3:
				color = make_float3(107,183,183);
				break;
				
			case 4:
				color = make_float3(107,111,183);
				break;
				
			case 5:
				color = make_float3(183,108,183);
				break;
		};*/
		
		//for delta
		/*switch(tag % 6)
		{
			case 0:
				color = make_float3(183,107,107);
				break;
				
			case 1:
				color = make_float3(183,107,181);
				break;
				
			case 2:
				color = make_float3(108,114,183);
				break;
				
			case 3:
				color = make_float3(108,181,183);
				break;
				
			case 4:
				color = make_float3(29,234,144);
				break;
		};*/
		
		
		/*switch (tag % 12)
		{
		case 0:
			color = make_float3(183,107,107);
			break;
		case 1:
			color = make_float3(168,160,183);
			break;
		case 2:
			color = make_float3(155,109,183);
			break;
		case 3:
			color = make_float3(137,183,158);
			//color = make_float3(0,0,0);
			break;
		case 4:
			color = make_float3(108,161,120);
			break;
		case 5:
			color = make_float3(183,108,140);
			break;
		case 6:
			color = make_float3(108,183,183);
			break;
		case 7:
			color = make_float3(205,107,77);
			break;
		case 8:
			color = make_float3(166,150,121);
			break;
		case 9:
			color = make_float3(120,109,183);
			//printf("%d ",(tagarr->GetValue(i) >> (6 * 2)) & 3);
			//if ((tagarr->GetValue(i) >> (6 * 2)) & 3)
			//	color = make_float3(0,0,0);
			break;
		case 10:
			color = make_float3(127,183,108);
			break;
		case 11:
			color = make_float3(183,183,109);
			break;
		};*/
		
		//if (tagarr->GetValue(i) == 0)
		//	color *= 0.0;

		float blank[4] = {  color.x,
							color.y,
							color.z,
							255};
		
		// next is for combustion only with mixfrac
		/*blank[3] = fabs(valarr->GetValue(i) - 0.42);
		if (blank[3] < 0.1)
			blank[3] = 255;
		else
			blank[3] = 255 * (1 - (blank[3] - 0.1));*/
		
		colors->SetTuple(i, blank);
	}
}

void AssignPerValueColor(vtkPolyData*& vtkMesh, vtkUnsignedCharArray* colors)
{
	char* colortable = "VALTBL";
	char* cuttable = "QUERYTBL";
	int fadedopacity = 200;
	//MakeQueryTableFromTags(vtkMesh);
	vtkColorTransferFunction* colormap = GetValueColorTable(HEAT_WHITE_TO_DARKRED, false, 0, 25700.0);

	// cut the mesh to get what you want
	bool docut = false;
	vtkPolyData* cutMesh = vtkMesh;
	if (docut) CutMeshBasedOnValue(cutMesh, cuttable, 0.5, 2.0);
	printf("Size of mesh after cutting is %d\n", cutMesh->GetNumberOfPoints());
	
	// assign colors to the cut mesh
	if (docut) 
	{
		vtkFloatArray* valarrCut = (vtkFloatArray*) cutMesh->GetPointData()->GetArray(colortable);
		vtkUnsignedCharArray* colorsCut = vtkUnsignedCharArray::New();
		colorsCut->SetNumberOfComponents(4);
		colorsCut->SetNumberOfTuples(cutMesh->GetNumberOfPoints());
		for (vtkIdType i = 0; i < cutMesh->GetNumberOfPoints(); i++)
		{
			double* point = cutMesh->GetPoint(i);
			float val = valarrCut->GetValue(i);

			float blank[4];	
			
			float3 fc = GetValueColor(colormap, val);
			blank[0] = 255 * fc.x;
			blank[1] = 255 * fc.y;
			blank[2] = 255 * fc.z;
			blank[3] = 255;
			
			colorsCut->SetTuple(i, blank);
		}
		cutMesh->GetPointData()->SetScalars(colorsCut);
	}
	
	// assign colors to the original mesh
	vtkFloatArray* valarrOri = (vtkFloatArray*) vtkMesh->GetPointData()->GetArray(colortable);
	vtkUnsignedCharArray* colorsOri = vtkUnsignedCharArray::New();
	colorsOri->SetNumberOfComponents(4);
	colorsOri->SetNumberOfTuples(vtkMesh->GetNumberOfPoints());
	for (vtkIdType i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
	{
		double* point = vtkMesh->GetPoint(i);
		float val = valarrOri->GetValue(i);

		float blank[4];	
		
		float3 fc = GetValueColor(colormap, val);
		blank[0] = 255 * fc.x;
		blank[1] = 255 * fc.y;
		blank[2] = 255 * fc.z;
		/*if (docut) 
			blank[3] = fadedopacity;
		else
			blank[3] = 255;*/
		
		//if (((vtkFloatArray*) vtkMesh->GetPointData()->GetArray(cuttable))->GetValue(i) > 0.5)
		{
			blank[3] = 255;
		}
		/*else
		{
			blank[3] = fadedopacity;
		}*/
		
		colorsOri->SetTuple(i, blank);
	}
	vtkMesh->GetPointData()->SetScalars(colorsOri);
	
	// merge the two meshes
	if (docut) 
	{
		vtkAppendPolyData *app = vtkAppendPolyData::New();
		app->AddInput(vtkMesh);
		app->AddInput(cutMesh);
		app->Update();
		vtkPolyData* oriMesh = vtkMesh;
		vtkMesh = app->GetOutput();

		//oriMesh->Delete();
		//cutMesh->Delete();
		//app->Delete();
	}
	
	// refill input color array with proper values
	colors->Reset();
	colors->SetNumberOfTuples(vtkMesh->GetNumberOfPoints());	
	for (vtkIdType i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
	{
		double old_color[4];
		vtkMesh->GetPointData()->GetScalars()->GetTuple(i, old_color);
		
		colors->SetTuple(i, old_color);
	}

	colormap->Delete();
}


void SaturateColorsAndDarkBoundaries(vtkPolyData* vtkMesh, vtkUnsignedCharArray* colors, bool makedark)
{
	// cells to make dark
	set<int> todark;
	vtkCellArray* vtkcells = vtkMesh->GetPolys();
	vtkcells->InitTraversal();
	vtkIdType npts;
	vtkIdType* pts;
	if (makedark)
	{
		while (vtkcells->GetNextCell(npts, pts) != 0) 
		{
			int3 cls;
			cls.x = ((vtkIntArray*) vtkMesh->GetPointData()->GetArray("CLUSTBL"))->GetValue(pts[0]);
			cls.y = ((vtkIntArray*) vtkMesh->GetPointData()->GetArray("CLUSTBL"))->GetValue(pts[1]);
			cls.z = ((vtkIntArray*) vtkMesh->GetPointData()->GetArray("CLUSTBL"))->GetValue(pts[2]);
			bool dark = (cls.x != cls.y) || (cls.x != cls.z) || (cls.y != cls.z);
			if (dark)
			{
				todark.insert(pts[0]);
				todark.insert(pts[1]);
				todark.insert(pts[2]);
			}
		}
	}
	
	// now iterations to saturate colors
	for (vtkIdType i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
	{
		double color[4];
		colors->GetTuple(i, color);
		float4 c = make_float4(color[0], color[1], color[2], color[3]);
		ReduceSaturation(c);
		color[0] = c.x;
		color[1] = c.y;
		color[2] = c.z;
		colors->SetTuple(i, color);
	}
	
	// now iteration to darken the boundary
	for(set<int>::iterator it = todark.begin(); it != todark.end(); it++)
	{
		double color[4];
		colors->GetTuple(*it, color);
		float4 c = make_float4(color[0], color[1], color[2], color[3]);
		c *= 0.9;
		color[0] = c.x;
		color[1] = c.y;
		color[2] = c.z;
		colors->SetTuple(*it, color);
	}
}

int main (int argc, char *argv[])
{
	if (argc != 3)
	{
		printf("Example: vtkassigncolor mesh.vtk out.vtk\n");
		return 0;
	}

	cout << "==============\n";
	cout << "vtkassigncolor\n";
	cout << "==============\n";
	
	printf("Mesh file is %s\n", argv[1]);
	
	/*Nrrd* nrrdfile = readNrrd(argv[3]);
	RegularGrid* dataset = new RegularGrid(nrrdfile, 
										   nrrdfile->axis[0].size, nrrdfile->axis[1].size, nrrdfile->axis[2].size, 
										   nrrdfile->axis[0].spacing, nrrdfile->axis[1].spacing, nrrdfile->axis[2].spacing);
	*/
	
	// read the mesh file
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader->SetFileName(argv[1]);
    vtkMeshReader->Update();
	vtkPolyData* vtkMesh = vtkMeshReader->GetOutput();
	
	// get the points
	vtkPoints* points = vtkMesh->GetPoints();
	cout << "Number of points is " << points->GetNumberOfPoints() << "\n";
	vtkPointData* pointsData = vtkMesh->GetPointData();
	
	// create the array for the colors
	double bounds[6];
	vtkMesh->GetBounds(bounds);
	vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
	colors->SetNumberOfComponents(4);
	colors->SetNumberOfTuples(vtkMesh->GetNumberOfPoints());
	
	
	
	
	// rescale to camera view
	/*float3 center;
	center.x = (bounds[1] + bounds[0]) / 2.0;
	center.y = (bounds[3] + bounds[2]) / 2.0;
	center.z = (bounds[5] + bounds[4]) / 2.0;
	cout << "Center at " << center.x << " " << center.y << " " << center.z << "\n";
	vtkTransform* transf = vtkTransform::New();
	transf->Scale(0.001,0.001, 0.001);
	vtkTransformPolyDataFilter* tf = vtkTransformPolyDataFilter::New();
	tf->SetInput(vtkMesh);
	tf->SetTransform(transf);
	tf->Update();
	vtkMesh = tf->GetOutput();*/
	
	
	
	vtkPlane* clipplane = vtkPlane::New();
	//clipplane->SetOrigin(0.75,0,0);
	//clipplane->SetOrigin(0.6,0,0);
	//clipplane->SetNormal(1,0,0);
	clipplane->SetNormal(0,-1,-1);
	clipplane->SetOrigin(0.3,0.124,0.124);
	vtkClipPolyData* clipper2 = vtkClipPolyData::New();
	clipper2->SetInput(vtkMesh);
	clipper2->SetClipFunction(clipplane);
	clipper2->Update();
	vtkMesh = clipper2->GetOutput();
	
	
	// find the color 
	//AssignPerValueColor(vtkMesh, colors);
	AssignPerClusterColor(vtkMesh, colors);
	BlackBoundary(vtkMesh, colors);
	SmoothColor(vtkMesh, colors, 15, 0.1);
	SaturateColorsAndDarkBoundaries(vtkMesh, colors, false);
	
	
	// add color array
	vtkMesh->GetPointData()->SetScalars(colors);

	// write file
	vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkMeshWriter->SetFileName(argv[2]);
	vtkMeshWriter->SetInput(vtkMesh);
	vtkMeshWriter->Write();
	
	//nrrdNuke(nrrdfile);
	
	return 0;
}