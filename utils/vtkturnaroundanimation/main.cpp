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
// Includes
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <vector_types.h>
#include <vector_functions.h>
#include <cutil_inline.h>
#include <cutil_math.h>
#include <omp.h>
#include <teem/nrrd.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataNormals.h>
#include <vtkPlane.h>
#include <vtkClipPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkFeatureEdges.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkWindowToImageFilter.h>
#include <vtkTIFFWriter.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkLODActor.h>
#include <vtkConeSource.h>
#include <vtkGlyph3D.h>
#include <vtkCellPicker.h>
#include <vtkTextMapper.h>
#include <vtkActor2D.h>
#include <vtkTextProperty.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkLookupTable.h>
#include <vtkObjectFactory.h>
#include <vtkLine.h>
#include <vtkScalarBarActor.h>
#include <vtkLegendBoxActor.h>
#include <vtkCameraWidget.h>
#include <vtkCameraRepresentation.h>
#include <vtkCameraInterpolator.h>
#include <vtkInteractorEventRecorder.h>
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>

#include "mvtkSmoothBoundaries.h"
#include "mvtkAssignNrrd.h"
#include "vmtkAssignColor.h"

using namespace std;

void LoadMeshToScene(string meshfile, string nrrdfile, int idx);
void SaveViewPort();

//*****************************************************************************
// Some general variables
//*****************************************************************************

string meshesdirectory;
string meshesfilter;
vector<string> meshesfiles;

string nrrdsdirectory;
string nrrdsfilter;
vector<string> nrrdsfiles;

string imgsprefix;
string camerapathfile;

int imgscount = 0;
int meshindex = 0;

bool firstmesh = true;
float3 focal;

//*****************************************************************************
// Camera path reader
//*****************************************************************************
FILE* ppc_file = NULL;

int ppc_on = 0;
int ppc_index = 0;
int ppc_count = 10000;
int ppc_speed = 8;

float3 PPCP_read()
{
	float3 ret;

	char str[128];
	fscanf (ppc_file, "%s\n", str);
	fscanf (ppc_file, "%s\n", str); // x
	str[strlen(str) - 1] = '\0';
	ret.x = atof(str+3);
	fscanf (ppc_file, "%s\n", str); // y
	str[strlen(str) - 1] = '\0';
	ret.y = atof(str+3);
	fscanf (ppc_file, "%s\n", str); // z
	str[strlen(str) - 1] = '\0';
	ret.z = atof(str+3);
	fscanf (ppc_file, "%s\n", str);
	fscanf (ppc_file, "%s\n", str);

	return ret;
}

float3 PPCP_next()
{	
	if ((ppc_index >= ppc_count) /*|| (ppc_index == 2008)*/)
	{
		ppc_on = 0;
		ppc_index = 0;

		fclose (ppc_file);

		return make_float3(0.0);
	}

	float3 ret;
	for (int i = 0 ; (i < ppc_speed) && (ppc_index < ppc_count); i++)
	{
		ret = PPCP_read();
		ppc_index++;
	}

	if ((ppc_index % 100) == 0)
		printf("Current Position %f %f %f at %d out of %d\n", ret.x, ret.y, ret.z, ppc_index, ppc_count);

	return ret;
}

void PPCP_start(const char* filename)
{
	printf("Start camera path \n");

	char str [128];
	
	ppc_file = fopen (filename,"r+");	

	for(int i = 0; i < 20; i++)
		fscanf (ppc_file, "%s\n", str);

	
	fscanf (ppc_file, "%s\n", str); // number of points
	str[strlen(str) - 1] = '\0';
	ppc_count = atoi(str+12);
	fscanf (ppc_file, "%s\n", str);

	ppc_on = 1;
	ppc_index = 0;	

	ppc_speed = floorf(ppc_count / meshesfiles.size());

	PPCP_next();
}

//*****************************************************************************
// Objects defenitions
//*****************************************************************************
vtkPolyDataMapper* edgesmapper;
vtkPolyDataMapper* meshmapper;
vtkActor* edgesactor;
vtkActor* meshactor;
vtkScalarBarActor* scalarBar = NULL;
vtkRenderWindowInteractor* iren;
vtkRenderWindow* renwin;
vtkRenderer* ren;
vtkCameraRepresentation* rep;
vtkCameraWidget* widget;
vtkCameraInterpolator* interpolator;
vtkInteractorEventRecorder* recorder;

bool white = false;

int gci = 0;
class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
  public:
    static KeyPressInteractorStyle* New();
	vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);
 
	virtual void OnKeyPress() 
	{
		// Get the keypress
		std::string key = iren->GetKeySym();
		printf("Key %s pressed\n", key.c_str()); 

		if (key == "P")
		{
			printf("Play all meshes\n");
			for (int i = 0; i < meshesfiles.size(); i++)
			{
				meshindex = i;

				// update the camera
				if (ppc_file != NULL)
				{
					float3 camerapos = PPCP_next();
					ren->GetActiveCamera()->SetPosition(camerapos.x, camerapos.y, camerapos.z);
				}

				// direction
				/*double position[3];
				ren->GetActiveCamera()->GetPosition(position);
				ren->GetActiveCamera()->SetFocalPoint(focal.x, focal.y, focal.z);*/

				// load the new file
				printf("%d: Mesh File %s\n", meshindex, meshesfiles[meshindex].c_str());
				LoadMeshToScene(meshesfiles[meshindex], nrrdsfiles[meshindex], meshindex);
				SaveViewPort();
			}
			meshindex = 0;
		}

		if (key == "p")
		{
			printf("Play single meshe\n");
			//while (meshindex < meshesfiles.size())
			{
				cout << "Enter an index between 0 and " << (meshesfiles.size() - 1) << ":\n";
				cin >> meshindex;

				// update the camera
				if (ppc_file != NULL)
				{
					float3 camerapos = PPCP_next();
					ren->GetActiveCamera()->SetPosition(camerapos.x, camerapos.y, camerapos.z);
				}

				// direction
				/*double position[3];
				ren->GetActiveCamera()->GetPosition(position);
				ren->GetActiveCamera()->SetFocalPoint(focal.x, focal.y, focal.z);*/

				// load the new file
				printf("Mesh File %s\n", meshesfiles[meshindex].c_str());
				LoadMeshToScene(meshesfiles[meshindex], nrrdsfiles[meshindex], meshindex);
				SaveViewPort();
				meshindex++;
			}
			meshindex = 0;
		}

		if (key == "c")
		{
			printf("Capture picture\n");
			SaveViewPort();
		}

		if (key == "s")
		{
			imgscount = 0;
			printf("Reset image file count to %d\n", imgscount);
		}


		if (key == "f")
		{
			if (white)
			{
				ren->SetBackground(0, 0, 0);
				white = false;
			}
			else
			{
				ren->SetBackground(1, 1, 1);
				white = true;
			}
		}
		if (key == "z")
		{
			rep->AnimatePath(iren);
		}

		if (key == "X")
		{
			gci = 0;
		}
		if (key == "x")
		{
			double tmin = interpolator->GetMinimumT();
			double tmax = interpolator->GetMaximumT();
			int nsteps = 400;
			printf("Animation from time %lf to time %lf\n", tmin, tmax);
			//for (int i = 0; i < nsteps; i++)
			{
				double t = tmin + gci * (tmax - tmin) / (nsteps - 1);
				//printf("%lf\n", t);

				vtkCamera* cc = vtkCamera::New();
				interpolator->InterpolateCamera(t, cc);
				//cc->SetClippingRange(0.0, 100);

				double cr[2];
				cc->GetClippingRange(cr);
				printf("clipping range is %lf %lf \n", cr[0], cr[1]);

				
				ren->GetActiveCamera()->SetPosition(cc->GetPosition());
				ren->GetActiveCamera()->SetFocalPoint(cc->GetFocalPoint());
				ren->GetActiveCamera()->SetViewUp(cc->GetViewUp());
				ren->GetActiveCamera()->SetViewAngle(cc->GetViewAngle());

				gci++;
			}
			SaveViewPort();
		}
	}
};

vtkStandardNewMacro(KeyPressInteractorStyle);

//*****************************************************************************
// Initialize scene main objects
//*****************************************************************************
void InitializeScene()
{
	// create actors and mappers
	edgesmapper = vtkPolyDataMapper::New();
	meshmapper = vtkPolyDataMapper::New();
	edgesactor = vtkActor::New();
	meshactor = vtkActor::New();

	// Create the Renderer, RenderWindow, and RenderWindowInteractor
    KeyPressInteractorStyle* style = KeyPressInteractorStyle::New();

	// create a window to render into
	renwin = vtkRenderWindow::New();
	ren = vtkRenderer::New();
	renwin->AddRenderer(ren);
	ren->AddActor(edgesactor);
	ren->AddActor(meshactor);
	ren->SetBackground(0, 0, 0);
	ren->TwoSidedLightingOn();
	renwin->SetSize(960, 720);

	// create an interactor
	iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renwin);
	iren->SetInteractorStyle(style);

	// Create the widget
	rep = vtkCameraRepresentation::New();
	rep->SetNumberOfFrames(400);
	interpolator = vtkCameraInterpolator::New();
	interpolator->SetInterpolationTypeToSpline();
	rep->SetInterpolator(interpolator);
	widget = vtkCameraWidget::New();
	widget->SetInteractor(iren);
	widget->SetRepresentation(rep);
	rep->SetCamera(ren->GetActiveCamera());
	widget->On();

	// mesh for the delta wing
	vtkPolyDataReader *reader = vtkPolyDataReader::New();
	reader->SetFileName("Z:\\nabla\\Samer\\Datasets\\Multifield\\TEMP\\edeltamesh.vtk");
	reader->Update();
	vtkPolyData* vtkMesh = reader->GetOutput();
	vtkPolyDataMapper* wingmapper = vtkPolyDataMapper::New();
	wingmapper->SetInput(vtkMesh); 
	wingmapper->SetScalarModeToDefault();
	wingmapper->SetScalarVisibility(1);
	wingmapper->Update();
	vtkActor* wingactor = vtkActor::New();
	wingactor->SetMapper(wingmapper);
	wingactor->GetProperty()->SetColor(0.7, 0.66, 0.5);
	wingactor->GetProperty()->SetOpacity(1.0);
	ren->AddActor(wingactor);


	ren->ResetCamera();
	ren->GetActiveCamera()->SetClippingRange(0.0, 10000);
	rep->GetCamera()->SetClippingRange(0.0, 10000);

	double cr[2];
	rep->GetCamera()->GetClippingRange(cr);
	printf("clipping range is %lf %lf \n", cr[0], cr[1]);
}

//*****************************************************************************
// Get a certain mesh ready to render
//*****************************************************************************
void LoadMeshToScene(string meshfile, string nrrdfile, int idx)
{
	// first read the mesh file
	vtkPolyDataReader *reader = vtkPolyDataReader::New();
	reader->SetFileName(meshfile.c_str());
	reader->Update();
	vtkPolyData* vtkMesh = reader->GetOutput();

	// focal point
	if (firstmesh)
	{
		firstmesh = false;
		double bounds[6];
		vtkMesh->GetBounds(bounds);
		focal.x = (bounds[0] + bounds[1]) / 2.0;
		focal.y = (bounds[2] + bounds[3]) / 2.0;
		focal.z = (bounds[4] + bounds[5]) / 2.0;

		printf("Focal Point %f %f %f\n", focal.x, focal.y, focal.z);
	}

	// do nrrd assignmed, coloring, and boundary smoothing
	//mvtkAssignNrrd(vtkMesh, string("VALTBL"), nrrdfile);
	if (scalarBar != NULL)
	{
		scalarBar->Delete();
		scalarBar = NULL;
	}
	mvtkAssignColor(vtkMesh, scalarBar);
	//mvtkSmoothBoundaries(vtkMesh, 5.0, 0.5);

	// move mesh down a bit
	vtkTransform* trans = vtkTransform::New();
	trans->Translate(-0.001*idx, 0, 0);
	vtkTransformPolyDataFilter* tf = vtkTransformPolyDataFilter::New();
	tf->SetInput(vtkMesh);
	tf->SetTransform(trans);
	tf->Update();
	vtkMesh = tf->GetOutput();


	// add the legend
	if (scalarBar != NULL)
	{
		ren->AddActor2D(scalarBar);
	}

	// compute the mesh normals
	vtkPolyDataNormals* normals = vtkPolyDataNormals::New();
	normals->SetInput(vtkMesh);
	normals->NonManifoldTraversalOn();
	normals->ComputePointNormalsOn();
	normals->ConsistencyOn();
	//normals->Update();

	// do the clipping
	vtkPlane* clipplane = vtkPlane::New();
	//clipplane->SetNormal(0,-1,-1);
	//clipplane->SetOrigin(0.3,0.124,0.124);
	clipplane->SetOrigin(0.75,0,0);
	clipplane->SetNormal(1,0,0);
	vtkClipPolyData* clipper = vtkClipPolyData::New();
	clipper->SetInput(vtkMesh);
	clipper->SetClipFunction(clipplane);

	// compute the normals for the clipped object
	vtkPolyDataNormals* normalscliped = vtkPolyDataNormals::New();
	normalscliped->SetInput(clipper->GetOutput());
	normalscliped->NonManifoldTraversalOn();
	normalscliped->ComputePointNormalsOn();
	normalscliped->ConsistencyOn();
	//normalscliped->Update();

	// map to graphics library (input is normals or normalscliped) (vtkPolyDataNormals use splitting to fix it)
	meshmapper->SetInput(normalscliped->GetOutput()); 
	meshmapper->SetScalarModeToDefault();
	meshmapper->SetScalarVisibility(1);
	meshmapper->Update();

	// The edge extractor for the mesh (input is reader or clipper)
	vtkFeatureEdges* edgesextractor = vtkFeatureEdges::New();
	edgesextractor->SetInput(clipper->GetOutput());
	edgesextractor->ColoringOff();
	edgesextractor->BoundaryEdgesOn();
	edgesextractor->ManifoldEdgesOff();
	edgesextractor->NonManifoldEdgesOff();
	edgesextractor->FeatureEdgesOff();
	edgesmapper->SetInput(edgesextractor->GetOutput());
	edgesmapper->SetResolveCoincidentTopologyToPolygonOffset();
	edgesmapper->ScalarVisibilityOff();
	edgesmapper->Update();
	edgesactor->SetMapper(edgesmapper);
	edgesactor->GetProperty()->SetLineWidth(3);
	edgesactor->GetProperty()->SetColor(0,0,0);
	edgesactor->VisibilityOn();

	// actor for the object
    meshactor->SetMapper(meshmapper);
	meshactor->GetProperty()->SetColor(0.7, 0.66, 0.5);
	meshactor->GetProperty()->SetOpacity(1.0);

	// delete objects
	/*reader->Delete();
	normals->Delete();
	clipplane->Delete();
	clipper->Delete();
	normalscliped->Delete();
	edgesextractor->Delete();
	trans->Delete();
	tf->Delete();*/

	printf("Done\n");
}

//*****************************************************************************
// Save picture
//*****************************************************************************
void SaveViewPort()
{
	vtkWindowToImageFilter *capture = vtkWindowToImageFilter::New();
	capture->SetInput(renwin);

	vtkTIFFWriter *writer = vtkTIFFWriter::New();
	writer->SetInputConnection(capture->GetOutputPort());

	char buffer[16];
	sprintf(buffer, "%d", imgscount++);
	string ids = std::string( 5 - strlen(buffer), '0').append(string(buffer));

	std::string basename(imgsprefix+string(ids)+string(".tif"));

	writer->SetFileName(basename.c_str());
	writer->Write();

	capture->Delete();
	writer->Delete();
}

//*****************************************************************************
// All files in a directory satisfaying a filter
//*****************************************************************************
void GetFilesFiltered(string director, string filter, vector<string>& files)
{
	// get all meshes names
	const std::string target_path( director.c_str() );
	const boost::regex my_filter( filter.c_str() );
	boost::filesystem::directory_iterator end_itr; // Default ctor yields past-the-end
	for( boost::filesystem::directory_iterator i( target_path ); i != end_itr; ++i )
	{
		// Skip if not a file
		if( !boost::filesystem::is_regular_file( i->status() ) ) continue;

		// Skip if no match
		boost::smatch what;
		if( !boost::regex_match( i->path().leaf().string(), what, my_filter ) ) continue;

		// File matches, store it
		files.push_back(  i->path().string() );
	}
	sort(files.begin(), files.end());
	//for (int i = 0; i < files.size(); i++)
		//printf("File %s\n", files[i].c_str());
	printf("Found %d files\n", files.size());
}

//*****************************************************************************
// Main program
//*****************************************************************************
int main(int argc, char** argv) 
{
	// input information
	meshesdirectory = string("Z:\\nabla\\Samer\\Datasets\\Multifield\\delta65_mf\\");
	meshesfilter = string("Ds_0.*\.vtk");
	nrrdsdirectory = string("Z:\\nabla\\Samer\\Datasets\\Multifield\\delta65\\");
	nrrdsfilter = string("delta65_5_0.*\.nrrd");
	imgsprefix = string("D:\\Media\\Movies\\IEVIS12\\Clip13\\");
	camerapathfile = string("");

	// get all meshes names
	GetFilesFiltered(meshesdirectory, meshesfilter, meshesfiles);

	// get all meshes names
	GetFilesFiltered(nrrdsdirectory, nrrdsfilter, nrrdsfiles);

	// check every mesh has a nrrd
	if (nrrdsfiles.size() != meshesfiles.size())
	{
		printf("Wrong number of nrrds\n");
		getchar();
		return 0;
	}

	// path file
	if (camerapathfile.length() > 0)
	{
		PPCP_start(camerapathfile.c_str());
	}
	
	// start visualization
	InitializeScene();
	iren->Initialize();
	renwin->Render();
	iren->Start();


	return 0;
}
