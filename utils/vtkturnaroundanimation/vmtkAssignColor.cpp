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
#include "vmtkAssignColor.h"

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

	clipper3->Delete();
	clipper4->Delete();
}

void MakeBoundaryBlack(vtkPolyData* vtkMesh, vtkUnsignedCharArray* colors)
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
			double old_color[4];

			colors->GetTuple((*it).first.x, old_color);
			old_color[0] = 0;
			old_color[1] = 0;
			old_color[2] = 0;	
			colors->InsertTuple((*it).first.x, old_color);

			colors->GetTuple((*it).first.y, old_color);
			old_color[0] = 0;
			old_color[1] = 0;
			old_color[2] = 0;	
			colors->InsertTuple((*it).first.y, old_color);
		}
	}
}

void MakeColorsSmooth(vtkPolyData* vtkMesh, vtkUnsignedCharArray* colors, int iterations, float weight)
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

	newcolors->Delete();
}

enum COLORMAPS {SATURATED_RED_TO_MAGNETA, 
				YELLOW_TO_BLUE, 
				WHITE_TO_MAGNETA, 
				HEAT_WHITE_TO_DARKRED, 
				CENTER_ORANGE_TO_MAGNETA, 
				RED_TO_BLUE, 
				MARON_TO_BLACK };
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
		mapcol.push_back(make_float3(255,255,217));
		mapcol.push_back(make_float3(237,248,217));
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
		mapcol.push_back(make_float3(255, 255, 229));
		mapcol.push_back(make_float3(255, 247, 188));
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
	// reads TAGTBL and computes the table QUERYTBL of zeros and ones

	vtkSmartPointer<vtkFloatArray> newScalars = vtkSmartPointer<vtkFloatArray>::New();
	newScalars->SetName("QUERYTBL");
	
	vtkIntArray* tagarr = (vtkIntArray*) vtkMesh->GetPointData()->GetArray("TAGTBL");
	for (vtkIdType i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
	{
		float value;
		int tag = tagarr->GetValue(i);
		if (tag & 768)
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

void AssignPerValueColor(vtkPolyData*& vtkMesh, vtkUnsignedCharArray* colors, vtkScalarBarActor*& scalarBar)
{
	char* colortable = "VALTBL";
	char* cuttable = "QUERYTBL";
	float minv = 0.28;
	float maxv = 0.306;
	bool docut = false;
	int fadedopacity = 50;
	//MakeQueryTableFromTags(vtkMesh);
	vtkColorTransferFunction* colormap = GetValueColorTable(WHITE_TO_MAGNETA, false, minv, maxv);

	// create the legend object
	scalarBar = vtkScalarBarActor::New();
	scalarBar->SetLookupTable(colormap);
	scalarBar->SetTitle("Density");
	vtkTextProperty* prpt = vtkTextProperty::New();
	prpt->SetFontSize(10);
	prpt->SetColor(0.0,0.0,0.0);
	prpt->SetFontFamilyToCourier();
	scalarBar->SetTitleTextProperty(prpt);
	scalarBar->SetNumberOfLabels(2);
	scalarBar->SetOrientationToHorizontal();
	scalarBar->SetTextureGridWidth(0.2);
	scalarBar->SetHeight(0.1);
	scalarBar->GetLabelTextProperty()->SetFontSize(8);
	scalarBar->GetLabelTextProperty()->SetColor(0,0,0);
	scalarBar->GetLabelTextProperty()->SetFontFamilyToCourier();
	scalarBar->SetPosition(0.78, 0.9);
	scalarBar->SetLabelFormat("%.2g");

	// cut the mesh to get what you want
	vtkPolyData* cutMesh = vtkMesh;
	if (docut) CutMeshBasedOnValue(cutMesh, cuttable, minv, maxv);
	
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
		if (docut) 
			blank[3] = fadedopacity;
		else
			blank[3] = 255;
		
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

		oriMesh->Delete();
		cutMesh->Delete();
		app->Delete();
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

void AssignPerClusterColor(vtkPolyData* vtkMesh, vtkUnsignedCharArray* colors)
{
	// cluster array, tag array
	vtkIntArray* clsarr = (vtkIntArray*) vtkMesh->GetPointData()->GetArray("CLUSTBL");
	vtkIntArray* tagarr = (vtkIntArray*) vtkMesh->GetPointData()->GetArray("TAGTBL");
	
	// find all tags
	set<int> tags;
	for (vtkIdType i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
		tags.insert(clsarr->GetValue(i));
	printf("Number of distinct clusters is %d\n", tags.size());	
		
	// now assign colors
	for (vtkIdType i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
	{
		int tag = clsarr->GetValue(i);
		
		float3 color = make_float3(255,255,255);
		
		// for multifield
		/*switch(tag % 8)
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

			default:
				printf("Error: unknown cluster\n");
				getchar();
				break;
		};*/
		
		
		//for combustion
		switch(tag % 6)
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

			default:
				printf("Error: unknown cluster\n");
				getchar();
				break;
		};
		
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

			default:
				printf("Error: unknown cluster\n");
				getchar();
				break;
		};*/
		
		// any different colors
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
				break;
			case 10:
				color = make_float3(127,183,108);
				break;
			case 11:
				color = make_float3(183,183,109);
				break;

			default:
				printf("Error: unknown cluster\n");
				getchar();
				break;
		};*/

		float blank[4] = {  color.x,
							color.y,
							color.z,
							255};
		
		colors->SetTuple(i, blank);
	}
}




void mvtkAssignColor(vtkPolyData*& vtkMesh, vtkScalarBarActor*& scalarBar)
{

	cout << "==============\n";
	cout << "vtkassigncolor\n";
	cout << "==============\n";

	
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
		

	// find the color 
	AssignPerValueColor(vtkMesh, colors, scalarBar);
	//AssignPerClusterColor(vtkMesh, colors);
	
	
	// enhance colors
	//MakeBoundaryBlack(vtkMesh, colors);
	//MakeColorsSmooth(vtkMesh, colors, 5, 0.1);
	
	// add color array
	vtkMesh->GetPointData()->SetScalars(colors);
}