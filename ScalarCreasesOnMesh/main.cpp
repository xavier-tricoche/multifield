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
#include <map>
#include <queue>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <omp.h>
#include <vector_types.h>
#include <vector_functions.h>
#include <cutil_inline.h>
#include <cutil_math.h>
#include <math/fixed_vector.hpp>
#include <math/bezier_spline.hpp>
#include <math/dopri5.hpp>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <teem/nrrd.h>
#include <teem/ten.h>
#include <teem/seek.h>
#include <teem/air.h>
#include <teem/limn.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkLODActor.h>
#include <vtkConeSource.h>
#include <vtkGlyph3D.h>
#include <vtkCellPicker.h>
#include <vtkTextMapper.h>
#include <vtkActor2D.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTextProperty.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkLookupTable.h>
#include <vtkObjectFactory.h>
#include <vtkProperty.h>
#include <vtkLine.h>
#include <vtkCellArray.h>
#include <vtkTubeFilter.h>
#include <vtkPNGReader.h>
#include <vtkImageActor.h>
#include <vtkPointPicker.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPolyDataReader.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>
#include <vtkExtractEdges.h>
#include <vtkPolyDataNormals.h>

#include "MyMath.h"
#include "MyTeem.h"
#include "MyGeometry.h"

using namespace std;


////////////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////////////

vtkRenderer *ren;
vtkRenderWindow *renWin;
vtkRenderWindowInteractor *iren;
vtkPolyData *linespolydata;
vtkPolyDataMapper* linespolydatamapper;
vtkActor* linespolydataactor;
vtkTubeFilter* triangleEdgeTubes;
vtkPolyDataMapper* triangleEdgeMapper;
vtkExtractEdges* triangleEdges;
vtkTubeFilter* ridgeEdgeTubes;
vtkPolyDataMapper* ridgeEdgeMapper;
vtkExtractEdges* ridgeEdges;
vtkPointPicker *picker;
vtkActor2D *textActor;
vtkTextMapper *textMapper;
int MouseMotion;
vtkActor* meshActor;
vtkPolyData* vtkMesh;

double GetDistance(vtkPolyData* mesh, vtkIdType id1, vtkIdType id2)
{
	double pt[3];
	mesh->GetPoint(id1, pt);
	float3 pt1 = make_float3(pt[0], pt[1], pt[2]);
	mesh->GetPoint(id2, pt);
	float3 pt2 = make_float3(pt[0], pt[1], pt[2]);

	return length(pt2 - pt1);
}

double GetValueAtDist(vtkPolyData* mesh, vtkIdType id1, vtkIdType id2, double dist)
{
	//vtkFloatArray* valarr = (vtkFloatArray*) mesh->GetPointData()->GetArray("VALTBL");
	vtkDoubleArray* valarr = (vtkDoubleArray*) mesh->GetPointData()->GetArray("stripes");

	double pt[3];
	mesh->GetPoint(id1, pt);
	float3 pt1 = make_float3(pt[0], pt[1], pt[2]);
	mesh->GetPoint(id2, pt);
	float3 pt2 = make_float3(pt[0], pt[1], pt[2]);

	double val1 = valarr->GetValue(id1);
	double val2 = valarr->GetValue(id2);

	double r = dist / length(pt2 - pt1);
	return val1 * (1 - r) + val2 * r;
	//return val2;
}

void GetNeighbors(vtkPolyData* mesh, vtkIdType idx, set<int>* n)
{
	vtkIdList* cells = vtkIdList::New();
	mesh->GetPointCells(idx, cells);
	for (int i = 0; i < cells->GetNumberOfIds(); i++)
	{
		int cellid = cells->GetId(i);
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
		mesh->GetCellPoints(cellid, ptIds);
		
		for (int j = 0; j < ptIds->GetNumberOfIds(); j++)
		{
			int id = ptIds->GetId(j);
			if (id != idx)
				n->insert(id);
		}
	}
}

int IsCreasePoint(vtkPolyData* mesh, vtkIdType idx)
{
	vtkIdList* cells = vtkIdList::New();
	mesh->GetPointCells(idx, cells);
	vector<int2> sides;
	double dist = numeric_limits<double>::max();
	for (int i = 0; i < cells->GetNumberOfIds(); i++)
	{
		int cellid = cells->GetId(i);
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
		mesh->GetCellPoints(cellid, ptIds);
		
		int2 side;
		int k = 0;
		for (int j = 0; j < ptIds->GetNumberOfIds(); j++)
		{
			int id = ptIds->GetId(j);
			if (id == idx)
				continue;

			dist = min(dist, GetDistance(mesh, idx, id));
			if (k == 0)
				side.x = id;
			else
				side.y = id;
			k++;
		}

		sides.push_back(side);
	}

	for (int i = 0; i < sides.size(); i++)
	{
		int k = -1;
		for (int j = i + 1; j < sides.size(); j++)
		{
			if (sides[j].x == sides[i].y) 
			{
				k = j;
				break;
			}

			if (sides[j].y == sides[i].y)
			{
				int tmp = sides[j].y;
				sides[j].y = sides[j].x;
				sides[j].x = tmp;
				k = j;
				break;
			}
		}
		if ((k == -1) && (i != sides.size() - 1))
			return false;
		if (k == -1)
			continue;

		int2 side = sides[k];
		sides.erase(sides.begin() + k);
		sides.insert(sides.begin() + i + 1, side);
	}

	vtkDoubleArray* valarr = (vtkDoubleArray*) mesh->GetPointData()->GetArray("stripes");
	set<int> maxpts;
	vector<int> bla;
	double vi = valarr->GetValue(idx);
	double mxvl = vi;
	for (int i = 0; i < sides.size(); i++)
	{
		int id1 = sides[i].x;
		int id2 = sides[i].y;
		int id3 = sides[(i + 1) % sides.size()].y;
		bla.push_back(id2);

		double v1 = GetValueAtDist(mesh, idx, id1, dist);
		double v2 = GetValueAtDist(mesh, idx, id2, dist);
		double v3 = GetValueAtDist(mesh, idx, id3, dist);
		
		mxvl = max(mxvl, v2);

		if ((v2 >= v1) && (v2 >= v3))
			maxpts.insert(sides[i].y);
	}

	if (mxvl == vi)
		return 2;

	if (sides[0].x != sides[sides.size()-1].y)
		return false;
	
	if (maxpts.size() < 2)
		return false;

	while (maxpts.find(bla.front()) == maxpts.end())
	{
		bla.push_back(bla.front());
		bla.erase(bla.begin());
	}



	bool c1 = true;
	bool c2 = true;
	int inlowc = 0;
	bool inlow = false;
	bool hslow = false;
	double av = 0.0;
	int avc = 0;
 	for (int i = 0; i < bla.size(); i++)
	{
		if (maxpts.find(bla[i]) == maxpts.end())
		{
			if (!inlow)
				inlowc++;
			inlow = true;
			double v1 = GetValueAtDist(mesh, idx, bla[i], dist);
			av += v1;
			avc++;
			if (v1 < vi)
				hslow = true;
		}
		else
		{
			if ((inlow) && (!hslow))
				c1 = false;
			if ((inlow) && ((av/avc) > vi))
				c2 = false;
			inlow = false;
			hslow = false;
			av = 0.0;
			avc = 0;
		}
	}
	if ((inlow) && (!hslow))
		c1 = false;
	if ((inlow) && ((av/avc) > vi))
		c2 = false;

	if (inlowc < 2)
		return false;

	if ((!c1) && (!c2))
		return false;

	if (c2)
		return 2;

	if (c1)
		return 1;
}

void FindRidge(vtkPolyData* mesh)
{
	map<int, int> pt2d;
	set<vtkIdType> pts;
	vtkPoints* points = mesh->GetPoints();
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	{
		pt2d[i] = IsCreasePoint(mesh, i);
		if (pt2d[i])
			pts.insert(i);
	}

	vector<int2> edges;
	vector<int> pt2ln(points->GetNumberOfPoints());
	map<int, set<int>> neigbors;
	for (set<int>::iterator it = pts.begin(); it != pts.end(); it++)
	{
		set<int> n;
		GetNeighbors(mesh, *it, &n);
		for (set<int>::iterator nit = n.begin(); nit != n.end(); nit++)
		{
			if ((*it < *nit) && (pts.find(*nit) != pts.end()))
			{
				int2 nedge = make_int2(*it, *nit);
				
				if ((pt2ln[nedge.x] >= 2) && (pt2ln[nedge.y] >= 1))
					continue;

				if ((pt2ln[nedge.y] >= 2) && (pt2ln[nedge.x] >= 1))
					continue;

				if ((pt2d[nedge.x] != 2) && (pt2d[nedge.y] != 2))
					continue;


				pt2ln[nedge.x]++;
				pt2ln[nedge.y]++;

				neigbors[nedge.x].insert(nedge.y);
				neigbors[nedge.y].insert(nedge.x);

				edges.push_back(nedge);
			}
		}
	}

	// now do some smoothing
	for (int z = 0; z < 10; z++)
	{
		for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
		{
			if (pt2ln[i] <= 1)
				continue;

			double ws = 0.0;

			double pt1[3];
			points->GetPoint(i, pt1);
			float3 f1 = make_float3(pt1[0],pt1[1], pt1[2]);

			// find total weights
			for (set<int>::iterator it = neigbors[i].begin(); it != neigbors[i].end(); it++)
			{
				double pt2[3];
				points->GetPoint(*it, pt2);
				float3 f2 = make_float3(pt2[0],pt2[1], pt2[2]);
			
				ws += length(f1 - f2);
			}

			// compute weighted average
			float3 wap = make_float3(0.0);
			for (set<int>::iterator it = neigbors[i].begin(); it != neigbors[i].end(); it++)
			{
				double pt2[3];
				points->GetPoint(*it, pt2);
				float3 f2 = make_float3(pt2[0],pt2[1], pt2[2]);
			
				wap += length(f1 - f2) * f2 / ws;
			}

			pt1[0] = 0.8 * pt1[0] + 0.2 * wap.x;
			pt1[1] = 0.8 * pt1[1] + 0.2 * wap.y;
			pt1[2] = 0.8 * pt1[2] + 0.2 * wap.z;

			points->SetPoint(i, pt1);
		}
	}

	// reset the elements
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
		pt2ln[i] = 0;
	linespolydata->Delete();
	linespolydata = vtkPolyData::New();	
	vtkPoints* lpts = vtkPoints::New();
	vtkCellArray* lines = vtkCellArray::New();
	for (int i = 0; i < edges.size(); i++)
	{
		double pt1[3];
		points->GetPoint(edges[i].x, pt1);

		double pt2[3];
		points->GetPoint(edges[i].y, pt2);

		vtkIdType connectivity[2];
		connectivity[0] = lpts->InsertNextPoint(pt1);
		connectivity[1] = lpts->InsertNextPoint(pt2);

		vtkLine* line = vtkLine::New();
		line->GetPointIds()->SetId(0,connectivity[0]);
		line->GetPointIds()->SetId(1,connectivity[1]);
		lines->InsertNextCell(line);
	}

	// add geometry
	linespolydata->SetPoints(lpts);
	linespolydata->SetLines(lines);

	// update the mapper
	linespolydata->Update();
	linespolydata->Modified();
	linespolydatamapper->SetInput(linespolydata);
	linespolydatamapper->Update();

	// write the output
	vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
	writer->SetInput(linespolydata);
	writer->SetFileName("C:\\Users\\sbarakat\\Desktop\\test_data_samer\\test_data_samer\\duds_surface_largebox_smoothed_ftle_out.vtk");
	writer->Write();

	ridgeEdgeTubes->SetInput(linespolydata);
}

////////////////////////////////////////////////////////////////////////////////
// picker variables
////////////////////////////////////////////////////////////////////////////////
class PickCommand : public vtkCommand
{
public:

    static PickCommand *New() { return new PickCommand; }
    void Delete() { delete this; }
 
    virtual void Execute(vtkObject *caller, unsigned long l, void *callData)
    {
        if (picker->GetPointId() < 0 )
        {
            textActor->VisibilityOff();
        }
        else
        {
            double selpt[3];
            picker->GetSelectionPoint(selpt);
            double x = selpt[0];
            double y = selpt[1];
            double pickPos[3];
            picker->GetPickPosition(pickPos);
            double xp = pickPos[0];
            double yp = pickPos[1];
            double zp = pickPos[2];
       
			//vtkFloatArray* valarr = (vtkFloatArray*) vtkMesh->GetPointData()->GetArray("VALTBL");
			vtkDoubleArray* valarr = (vtkDoubleArray*) vtkMesh->GetPointData()->GetArray("stripes");
			double pt[3];
			vtkMesh->GetPoint(picker->GetPointId(), pt);

			if (picker->GetActor() != meshActor)
				return;

            char text[120];
			sprintf( text, "%d, %5.5lf %5.5lf %5.5lf, %5.5lf", picker->GetPointId(), pt[0], pt[1], pt[2], valarr->GetValue(picker->GetPointId()));
            textMapper->SetInput( text );
            textActor->SetPosition(x, y);
            textActor->VisibilityOn();
        }
   
        renWin->Render();
    }
};

void PickerInteractionCallback( vtkObject* vtkNotUsed(object), unsigned long event, void* clientdata, void* vtkNotUsed(calldata))
{
	vtkInteractorStyleTrackballCamera * style = 
	(vtkInteractorStyleTrackballCamera*)clientdata;
	switch( event )
	{
		case vtkCommand::LeftButtonPressEvent:
			MouseMotion = 0;
			style->OnLeftButtonDown();
			break;
		case vtkCommand::LeftButtonReleaseEvent:
			if (MouseMotion == 0)
			{
				int *pick = iren->GetEventPosition();
				//if (picker->GetActor() == meshActor)
				{
					picker->Pick((double)pick[0], (double)pick[1], 0.0, ren);
					//printf("User picked %d\n", picker->GetPointId());
				}
			}
			style->OnLeftButtonUp();
			break;
		case vtkCommand::MouseMoveEvent:
			MouseMotion = 1;
			style->OnMouseMove();
			break;
	}
}

class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
  public:
    static KeyPressInteractorStyle* New();
	vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);
 

	//MeshRegistration* meshreg;

	virtual void OnKeyPress() 
	{
		// Get the keypress
		vtkRenderWindowInteractor *rwi = this->Interactor;
		std::string key = rwi->GetKeySym();
 
		// Output the key that was pressed
		std::cout << "Pressed " << key << std::endl;
 
		if (key == "1")
		{

		}
		
		// Forward events
		vtkInteractorStyleTrackballCamera::OnKeyPress();
    } 
};

vtkStandardNewMacro(KeyPressInteractorStyle);

////////////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////////////
int main (int argc, char *argv[])
{
	MouseMotion = 0;

	// Create a cell picker.
    PickCommand* pickObserver = PickCommand::New();
    picker = vtkPointPicker::New();
    picker->AddObserver( vtkCommand::EndPickEvent, pickObserver );

	// Create the Renderer, RenderWindow, and RenderWindowInteractor
    KeyPressInteractorStyle* style = KeyPressInteractorStyle::New();
    vtkCallbackCommand* pickerCommand = vtkCallbackCommand::New();
    pickerCommand->SetClientData(style);
    pickerCommand->SetCallback(PickerInteractionCallback);
    style->AddObserver(vtkCommand::LeftButtonPressEvent, pickerCommand);
    style->AddObserver(vtkCommand::MouseMoveEvent, pickerCommand);
    style->AddObserver(vtkCommand::LeftButtonReleaseEvent, pickerCommand);

	// Create a text mapper and actor to display the results of picking.
    textMapper = vtkTextMapper::New();
    vtkTextProperty *tprop = textMapper->GetTextProperty();
    tprop->SetFontFamilyToArial();
    tprop->SetFontSize(12);
    tprop->BoldOn();
    tprop->SetColor(1, 0, 0);
    textActor = vtkActor2D::New();
    textActor->VisibilityOff();
    textActor->SetMapper(textMapper);

	// read the mesh file
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	//vtkMeshReader->SetFileName("E:\\Chaos2Rost\\Projects\\ScalarCreasesOnMesh\\Ds_0060.vtk");
	//vtkMeshReader->SetFileName("C:\\Users\\sbarakat\\Desktop\\test_data_segmentation\\test_data_segmentation\\test_sin_segmentation.vtk");
	//vtkMeshReader->SetFileName("C:\\Users\\sbarakat\\Desktop\\test_data_samer\\test_data_samer\\duds_surface_largebox_smoothed_ftle.vtk");
	vtkMeshReader->SetFileName("C:\\Users\\sbarakat\\Desktop\\test_sphere.vtk");
    vtkMeshReader->Update();
	vtkMesh = vtkMeshReader->GetOutput();
	vtkPoints* points = vtkMesh->GetPoints();

	// normals fixing object
	vtkPolyDataNormals* normals = vtkPolyDataNormals::New();
	normals->SetInput(vtkMesh);

	// mesh mapper and actor
	vtkPolyDataMapper* meshMapper = vtkPolyDataMapper::New(); 
	meshMapper->SetInputConnection(normals->GetOutputPort());
	meshActor = vtkActor::New();
	meshActor->SetMapper(meshMapper);

	// mesh edges extracted from msh
	triangleEdges = vtkExtractEdges::New();
	triangleEdges->SetInput(vtkMesh);
	triangleEdges->Update();
	triangleEdgeTubes = vtkTubeFilter::New();
	triangleEdgeTubes->SetInput(triangleEdges->GetOutput());
	triangleEdgeTubes->SetRadius(.001);
	triangleEdgeTubes->SetNumberOfSides(6);
	triangleEdgeTubes->UseDefaultNormalOn();
	triangleEdgeTubes->SetDefaultNormal(.577, .577, .577);
	triangleEdgeMapper = vtkPolyDataMapper::New();
	triangleEdgeMapper->SetInput(triangleEdgeTubes->GetOutput());
	triangleEdgeMapper->ScalarVisibilityOff();
	vtkActor* triangleEdgeActor = vtkActor::New();
	triangleEdgeActor->SetMapper(triangleEdgeMapper);
	triangleEdgeActor->GetProperty()->SetDiffuseColor(0, 0, 0);
	triangleEdgeActor->GetProperty()->SetSpecular(.4);
	triangleEdgeActor->GetProperty()->SetSpecularPower(10);

	// tubes for the ridge
	ridgeEdgeTubes = vtkTubeFilter::New();
	ridgeEdgeTubes->SetRadius(.03);
	ridgeEdgeTubes->SetNumberOfSides(6);
	ridgeEdgeTubes->UseDefaultNormalOn();
	ridgeEdgeTubes->SetDefaultNormal(.577, .577, .577);
	ridgeEdgeMapper = vtkPolyDataMapper::New();
	ridgeEdgeMapper->SetInput(ridgeEdgeTubes->GetOutput());
	ridgeEdgeMapper->ScalarVisibilityOff();
	vtkActor* ridgeEdgeActor = vtkActor::New();
	ridgeEdgeActor->SetMapper(ridgeEdgeMapper);
	ridgeEdgeActor->GetProperty()->SetDiffuseColor(0, 0, 0);
	ridgeEdgeActor->GetProperty()->SetSpecular(.4);
	ridgeEdgeActor->GetProperty()->SetSpecularPower(10);

	// color the mesh with scalar
	meshMapper->SetScalarModeToDefault();
	meshMapper->SetScalarVisibility(1);
	meshMapper->SetScalarModeToUsePointFieldData();
	meshMapper->SetScalarRange(0,0.05);
	//meshMapper->ColorByArrayComponent("stripes", 0);
	meshMapper->Update();

	// lines to connect matches
	linespolydata = vtkPolyData::New();
	linespolydatamapper = vtkPolyDataMapper::New();
	linespolydatamapper->SetInput(linespolydata);
	linespolydataactor = vtkActor::New();
	linespolydataactor->SetMapper(linespolydatamapper);
	linespolydataactor->GetProperty()->SetColor(0,0,0);

	// render window
	ren = vtkRenderer::New();
	renWin = vtkRenderWindow::New();
    renWin->AddRenderer(ren);
    iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);
    iren->SetInteractorStyle(style);
    iren->SetPicker(picker);

	// Add the actors to the renderer, set the background and size
	//ren->AddActor(linespolydataactor);
	ren->AddActor2D(textActor);
    ren->AddActor(meshActor);
    //ren->AddActor(triangleEdgeActor);
	ren->AddActor(ridgeEdgeActor);
    ren->SetBackground(1, 1, 1);
	ren->TwoSidedLightingOn();
    renWin->SetSize(768,768);

	FindRidge(vtkMesh);
    iren->Initialize();
    iren->Start();

	return 0;
}