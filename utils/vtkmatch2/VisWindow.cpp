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
#include "VisWindow.h"

int MouseMotion;
vtkRenderer *ren1;
vtkRenderWindow *renWin;
vtkRenderWindowInteractor *iren;
vtkCellPicker *picker;
vtkActor2D *textActor;
vtkTextMapper *textMapper;
vtkActor *mesh1Actor;
vtkActor *mesh2Actor;
int selectedMesh;
vtkIdType selectedPoint;
double transMesh2[3];
vtkPolyDataMapper *mesh1Mapper;
vtkPolyDataMapper *mesh2Mapper;

vtkPolyData *linespolydata;
vtkPolyDataMapper* linespolydatamapper;
vtkActor* linespolydataactor;

vtkPolyData *treelinespolydata;
vtkPolyDataMapper* treelinespolydatamapper;
vtkActor* treelinespolydataactor;

class PickCommand : public vtkCommand
{
public:

    static PickCommand *New() { return new PickCommand; }
    void Delete() { delete this; }
 
    virtual void Execute(vtkObject *caller, unsigned long l, void *callData)
    {
        if (picker->GetCellId() < 0 )
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
            picker->GetPickPosition( pickPos );
            double xp = pickPos[0];
            double yp = pickPos[1];
            double zp = pickPos[2];
       
            char text[120];
            sprintf( text, "(%5.5f,  %5.5f,  %5.5f)", xp, yp, zp );
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
				picker->Pick((double)pick[0], (double)pick[1], 0.0, ren1);
				printf("User picked %d\n", picker->GetPointId());
				if (picker->GetActor() == mesh1Actor)
				{
					selectedPoint = picker->GetPointId();
					selectedMesh = 1;
					printf("Mesh1\n");
				}
				else
				{
					selectedPoint = picker->GetPointId();
					selectedMesh = 2;
					printf("Mesh2\n");
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
 

	MeshRegistration* meshreg;

	virtual void OnKeyPress() 
	{
		// Get the keypress
		vtkRenderWindowInteractor *rwi = this->Interactor;
		std::string key = rwi->GetKeySym();
 
		// Output the key that was pressed
		std::cout << "Pressed " << key << std::endl;
 
		// print selected info
		if (key == "i")
		{
			this->meshreg->PrintCandidates(selectedPoint, 1);
		}

		// Handle an arrow key
		if(key == "Up")
		{
			std::cout << "The up arrow was pressed." << std::endl;
		}

		// clear all to restart learning
		if (key == "c")
		{
			for (int id = 0; id < this->meshreg->vtkMesh1->GetNumberOfPoints(); id++)
			{
				this->meshreg->candidates[id].clear();
			}
			this->meshreg->updown = 0;
			printf("Learning reset\n");
		}

		// final fix and output 
		if (key == "f")
		{
			this->meshreg->FinalFix();

			double bounds[6];
			this->meshreg->vtkMesh1->GetBounds(bounds);

			// mesh 2 colors
			vtkUnsignedCharArray* colors2 = vtkUnsignedCharArray::New();
			colors2->SetNumberOfComponents(4);
			colors2->SetNumberOfTuples(this->meshreg->vtkMesh2->GetNumberOfPoints());
			for (vtkIdType i = 0; i < this->meshreg->vtkMesh2->GetNumberOfPoints(); i++)
			{
				double point[3];
				this->meshreg->vtkMesh1->GetPoint(this->meshreg->m2m1[i], point);

				float blank[4] = {  (point[0] - bounds[0]) / (bounds[1] - bounds[0]) * 255.0,
									(point[1] - bounds[2]) / (bounds[3] - bounds[2]) * 255.0,
									(point[2] - bounds[5]) / (bounds[5] - bounds[4]) * 255.0,
									255.0 };

				colors2->SetTuple(i, blank);
			}
			this->meshreg->vtkMesh2->GetPointData()->SetScalars(colors2); 
			mesh2Mapper->SetScalarModeToDefault();
			mesh2Mapper->SetScalarVisibility(1);
			mesh2Mapper->Update();
			mesh2Mapper->Modified();

			renWin->Render();

		}

		// do next step
		if (key == "n")
		{
			for (int i = 0; i < 4; i++)
			{
				this->meshreg->DoRegistrationStep();
				

				double bounds[6];
				this->meshreg->vtkMesh1->GetBounds(bounds);

				// mesh 1 colors
				vtkUnsignedCharArray* colors2 = vtkUnsignedCharArray::New();
				colors2->SetNumberOfComponents(4);
				colors2->SetNumberOfTuples(this->meshreg->vtkMesh2->GetNumberOfPoints());
				for (vtkIdType i = 0; i < this->meshreg->vtkMesh2->GetNumberOfPoints(); i++)
				{
					double point[3];
					this->meshreg->vtkMesh2->GetPoint(i, point);

					float blank[4] = {  255.0,255.0,255.0,255.0 };
					colors2->SetTuple(i, blank);
				}
				this->meshreg->vtkMesh2->GetPointData()->SetScalars(colors2); 
				mesh2Mapper->SetScalarModeToDefault();
				mesh2Mapper->SetScalarVisibility(1);
				mesh2Mapper->Update();
				mesh2Mapper->Modified();


				// mesh 1 colors
				vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
				colors->SetNumberOfComponents(4);
				colors->SetNumberOfTuples(this->meshreg->vtkMesh1->GetNumberOfPoints());
				for (vtkIdType i = 0; i < this->meshreg->vtkMesh1->GetNumberOfPoints(); i++)
				{
					double point[3];
					this->meshreg->vtkMesh1->GetPoint(i, point);

					float blank[4] = {  (point[0] - bounds[0]) / (bounds[1] - bounds[0]) * 255.0,
										(point[1] - bounds[2]) / (bounds[3] - bounds[2]) * 255.0,
										(point[2] - bounds[5]) / (bounds[5] - bounds[4]) * 255.0,
										255.0
									};
					colors->SetTuple(i, blank);


					if (this->meshreg->candidates[i].size() != 0)
					{
						colors2->SetTuple(this->meshreg->candidates[i][0].id, blank);
					}
				}
				this->meshreg->vtkMesh1->GetPointData()->SetScalars(colors); 
				mesh1Mapper->SetScalarModeToDefault();
				mesh1Mapper->SetScalarVisibility(1);
				mesh1Mapper->Update();
				mesh1Mapper->Modified();

				renWin->Render();
			}
		}

		// build the tree
		if (key == "b")
		{
			this->meshreg->BuildTree();
		}

		// set selected id
		if (key == "s")
		{
			char str[256];
			cin.get (str,256);
			selectedPoint = atoi(str);
		}

		// show one line
		if (key == "l")
		{
			// reset the elements
			linespolydata->Delete();
			linespolydata = vtkPolyData::New();
			
			vtkPoints* pts = vtkPoints::New();
			vtkCellArray* lines = vtkCellArray::New();

			//if (this->meshreg->candidates[selectedPoint].size() != 0)
			{
				int pt1id, pt2id;
				double pt1[3];
				double pt2[3];

				if (selectedMesh == 1)
				{
					pt1id = selectedPoint;
					pt2id = this->meshreg->candidates[selectedPoint][0].id;
					this->meshreg->vtkMesh1->GetPoint(pt1id, pt1);
					this->meshreg->vtkMesh2->GetPoint(pt2id, pt2);

					char text[120];
					sprintf( text, "%d (%5.5f,  %5.5f,  %5.5f)", pt1id, pt1[0], pt1[1], pt1[2] );
					textMapper->SetInput( text );
				}
				else
				{
					pt1id = this->meshreg->m2m1[selectedPoint];
					pt2id = selectedPoint;
					this->meshreg->vtkMesh1->GetPoint(pt1id, pt1);
					this->meshreg->vtkMesh2->GetPoint(pt2id, pt2);

					char text[120];
					sprintf( text, "%d (%5.5f,  %5.5f,  %5.5f)", pt2id, pt2[0], pt2[1], pt2[2] );
					textMapper->SetInput( text );
				}

				// change color up to the root
				//this->meshreg->vtkMesh1->GetPointData()->SetScalars(
				TNode* node = this->meshreg->id2node[pt1id];
				while (node->level > 0)
				{
					printf("%d -> %d\n", node->id, this->meshreg->candidates[node->id][0].id);
					float blank[4] = { 0, 0, 0, 255};
					this->meshreg->vtkMesh1->GetPointData()->GetScalars()->SetTuple(node->id, blank);
					this->meshreg->vtkMesh2->GetPointData()->GetScalars()->SetTuple(this->meshreg->candidates[node->id][0].id, blank);
					node = node->parent;
				}


				this->meshreg->PrintCandidates(pt1id, 1000);
				

				pt2[0] += transMesh2[0];
				pt2[1] += transMesh2[1];
				pt2[2] += transMesh2[2];

				vtkIdType connectivity[2];
				connectivity[0] = pts->InsertNextPoint(pt1);
				connectivity[1] = pts->InsertNextPoint(pt2);

				vtkLine* line = vtkLine::New();
				line->GetPointIds()->SetId(0,connectivity[0]);
				line->GetPointIds()->SetId(1,connectivity[1]);
				lines->InsertNextCell(line);

				// add geometry
				linespolydata->SetPoints(pts);
				linespolydata->SetLines(lines);

				// update the mapper
				linespolydata->Update();
				linespolydata->Modified();
				linespolydatamapper->SetInput(linespolydata);
				linespolydatamapper->Update();

				// set the color to black
				float blank[4] = { 0, 0, 0, 255};
				this->meshreg->vtkMesh1->GetPointData()->GetScalars()->SetTuple(pt1id, blank);
				this->meshreg->vtkMesh2->GetPointData()->GetScalars()->SetTuple(pt2id, blank);

				renWin->Render();
			}
		}

		// show all lines
		if (key == "L")
		{
			// reset the elements
			linespolydata->Delete();
			linespolydata = vtkPolyData::New();
			
			vtkPoints* pts = vtkPoints::New();
			vtkCellArray* lines = vtkCellArray::New();

			// add the lines
			for (int id = 0; id < this->meshreg->vtkMesh1->GetNumberOfPoints(); id++)
			{
				int matchid = (this->meshreg->candidates[id].size() > 0) ? this->meshreg->candidates[id][0].id : -1;
				if (matchid != -1)
				{
					double pt1[3];
					this->meshreg->vtkMesh1->GetPoint(id, pt1);

					double pt2[3];
					this->meshreg->vtkMesh2->GetPoint(matchid, pt2);

					pt2[0] += transMesh2[0];
					pt2[1] += transMesh2[1];
					pt2[2] += transMesh2[2];

					vtkIdType connectivity[2];
					connectivity[0] = pts->InsertNextPoint(pt1);
					connectivity[1] = pts->InsertNextPoint(pt2);

					vtkLine* line = vtkLine::New();
					line->GetPointIds()->SetId(0,connectivity[0]);
					line->GetPointIds()->SetId(1,connectivity[1]);
					lines->InsertNextCell(line);

				}
			}

			// add geometry
			linespolydata->SetPoints(pts);
			linespolydata->SetLines(lines);

			// update the mapper
			linespolydata->Update();
			linespolydata->Modified();
			linespolydatamapper->SetInput(linespolydata);
			linespolydatamapper->Update();

			renWin->Render();
		}

		// show all lines
		if (key == "T")
		{
			// reset the elements
			treelinespolydata->Delete();
			treelinespolydata = vtkPolyData::New();
			
			vtkPoints* pts = vtkPoints::New();
			vtkCellArray* lines = vtkCellArray::New();

			// add the lines
			queue< TNode* > mq;
			mq.push(this->meshreg->rootnode);
			while (!mq.empty())
			{
				TNode* node = mq.front();
				mq.pop();

				if ((node->parent != NULL) && (node->parent->id != -1))
				{
					// add line
					double pt1[3];
					this->meshreg->vtkMesh1->GetPoint(node->id, pt1);

					double pt2[3];
					this->meshreg->vtkMesh1->GetPoint(node->parent->id, pt2);

					vtkIdType connectivity[2];
					connectivity[0] = pts->InsertNextPoint(pt1);
					connectivity[1] = pts->InsertNextPoint(pt2);

					vtkLine* line = vtkLine::New();
					line->GetPointIds()->SetId(0,connectivity[0]);
					line->GetPointIds()->SetId(1,connectivity[1]);
					lines->InsertNextCell(line);
				}

				// add children
				for (int i = 0; i < node->children.size(); i++)
					mq.push(node->children[i]);
			}

			// add geometry
			treelinespolydata->SetPoints(pts);
			treelinespolydata->SetLines(lines);

			// update the mapper
			treelinespolydata->Update();
			treelinespolydata->Modified();
			treelinespolydatamapper->SetInput(treelinespolydata);
			treelinespolydatamapper->Update();

			renWin->Render();
		}

		// color by 3D parametrization
		if (key == "u")
		{
			double bounds[6];
			this->meshreg->vtkMesh1->GetBounds(bounds);

			// mesh 2 colors
			vtkUnsignedCharArray* colors2 = vtkUnsignedCharArray::New();
			colors2->SetNumberOfComponents(4);
			colors2->SetNumberOfTuples(this->meshreg->vtkMesh2->GetNumberOfPoints());
			for (vtkIdType i = 0; i < this->meshreg->vtkMesh2->GetNumberOfPoints(); i++)
			{
				double point[3];
				this->meshreg->vtkMesh2->GetPoint(i, point);

				float blank[4] = {  255.0,255.0,255.0,255.0 };
				colors2->SetTuple(i, blank);
			}
			this->meshreg->vtkMesh2->GetPointData()->SetScalars(colors2); 
			mesh2Mapper->SetScalarModeToDefault();
			mesh2Mapper->SetScalarVisibility(1);
			mesh2Mapper->Update();
			mesh2Mapper->Modified();


			// mesh 1 colors
			vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
			colors->SetNumberOfComponents(4);
			colors->SetNumberOfTuples(this->meshreg->vtkMesh1->GetNumberOfPoints());
			for (vtkIdType i = 0; i < this->meshreg->vtkMesh1->GetNumberOfPoints(); i++)
			{
				double point[3];
				this->meshreg->vtkMesh1->GetPoint(i, point);

				float blank[4] = {  (point[0] - bounds[0]) / (bounds[1] - bounds[0]) * 255.0,
									(point[1] - bounds[2]) / (bounds[3] - bounds[2]) * 255.0,
									(point[2] - bounds[5]) / (bounds[5] - bounds[4]) * 255.0,
									255.0
								};
				colors->SetTuple(i, blank);


				if (this->meshreg->candidates[i].size() != 0)
				{
					colors2->SetTuple(this->meshreg->candidates[i][0].id, blank);

					if (this->meshreg->candidates[i][0].id == i)
					{
						blank[0] = 255;
						blank[1] = 0;
						blank[2] = 0;
						colors->SetTuple(i, blank);
					}
				}

				
			}
			this->meshreg->vtkMesh1->GetPointData()->SetScalars(colors); 
			mesh1Mapper->SetScalarModeToDefault();
			mesh1Mapper->SetScalarVisibility(1);
			mesh1Mapper->Update();
			mesh1Mapper->Modified();

			mesh1Actor->GetProperty()->SetInterpolationToFlat();
			mesh2Actor->GetProperty()->SetInterpolationToFlat();

			renWin->Render();
		}

		// color by connected component
		if (key == "o")
		{
			// mesh 1 colors
			vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
			colors->SetNumberOfComponents(4);
			colors->SetNumberOfTuples(this->meshreg->vtkMesh1->GetNumberOfPoints());
			for (vtkIdType i = 0; i < this->meshreg->vtkMesh1->GetNumberOfPoints(); i++)
			{

				int cc = this->meshreg->id2cc[i];

				float3 color = make_float3(255,255,255);
				switch (cc % 8)
				{
				case 0:
					color = make_float3(0, 255, 255);
					break;
				case 1:
					color = make_float3(0, 0, 255);
					break;
				case 2:
					color = make_float3(255, 0, 255);
					break;
				case 3:
					color = make_float3(255, 0, 0);
					break;
				case 4:
					color = make_float3(255, 255, 0);
					break;
				case 5:
					color = make_float3(0, 255, 0);
					break;
				case 6:
					color = make_float3(0, 127, 255);
					break;
				case 7:
					color = make_float3(255, 127, 0);
					break;
				};

				float blank[4] = {  color.x,
									color.y,
									color.z,
									255.0
								};
				colors->SetTuple(i, blank);
			}
			this->meshreg->vtkMesh1->GetPointData()->SetScalars(colors); 
			mesh1Mapper->SetScalarModeToDefault();
			mesh1Mapper->SetScalarVisibility(1);
			mesh1Mapper->Update();
			mesh1Mapper->Modified();

			renWin->Render();
		}

		// color by boundary points
		if (key == "B")
		{
			// mesh 1 colors
			vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
			colors->SetNumberOfComponents(4);
			colors->SetNumberOfTuples(this->meshreg->vtkMesh1->GetNumberOfPoints());
			for (vtkIdType i = 0; i < this->meshreg->vtkMesh1->GetNumberOfPoints(); i++)
			{

				int cc = this->meshreg->id2cc[i];

				float3 color = make_float3(255,255,255);
				if (this->meshreg->bpts->find(i) != this->meshreg->bpts->end())
					color = make_float3(255,0,0);

				float blank[4] = {  color.x,
									color.y,
									color.z,
									255.0
								};
				colors->SetTuple(i, blank);
			}
			this->meshreg->vtkMesh1->GetPointData()->SetScalars(colors); 
			mesh1Mapper->SetScalarModeToDefault();
			mesh1Mapper->SetScalarVisibility(1);
			mesh1Mapper->Update();
			mesh1Mapper->Modified();

			renWin->Render();
		}

		// Forward events
		vtkInteractorStyleTrackballCamera::OnKeyPress();
    } 
};

vtkStandardNewMacro(KeyPressInteractorStyle);

void VisWindow::InitWindow()
{
    MouseMotion = 0;

	// color lookup table
	vtkLookupTable* lut = vtkLookupTable::New();
    lut->SetNumberOfColors(256);
    lut->SetHueRange(0.15, 1.0);
    lut->SetSaturationRange(1.0, 1.0);
    lut->SetValueRange(1.0, 1.0);
    lut->SetAlphaRange(1.0, 1.0);
    lut->SetRange(-1,1);

	// mesh 1 graphics primitives
    mesh1Mapper = vtkPolyDataMapper::New();
	mesh1Mapper->SetInput(meshreg->vtkMesh1);
    mesh1Mapper->GlobalImmediateModeRenderingOn();
	mesh1Mapper->SetScalarModeToUsePointFieldData();
	mesh1Mapper->SetUseLookupTableScalarRange(1);
	mesh1Mapper->ColorByArrayComponent("MTCHMSK", 0);
    mesh1Actor = vtkActor::New();
    mesh1Actor->SetMapper(mesh1Mapper);

	// mesh 2 graphics primitives
    mesh2Mapper = vtkPolyDataMapper::New();
	mesh2Mapper->SetInput(meshreg->vtkMesh2);
    mesh2Mapper->GlobalImmediateModeRenderingOn();
	mesh2Mapper->SetScalarModeToUsePointFieldData();
	mesh2Mapper->SetUseLookupTableScalarRange(1);
	mesh2Mapper->ColorByArrayComponent("MTCHMSK", 0);
    mesh2Actor = vtkActor::New();
    mesh2Actor->SetMapper(mesh2Mapper);
	double bounds[6];
	meshreg->vtkMesh2->GetBounds(bounds);
	transMesh2[0] = 0.0;
	transMesh2[1] = -1.2 * (bounds[3] - bounds[2]);
	transMesh2[2] = 0.0;
	mesh2Actor->SetPosition(transMesh2);

	// Create a cell picker.
    PickCommand* pickObserver = PickCommand::New();
    picker = vtkCellPicker::New();
    picker->AddObserver( vtkCommand::EndPickEvent, pickObserver );

	// Create a text mapper and actor to display the results of picking.
    textMapper = vtkTextMapper::New();
    vtkTextProperty *tprop = textMapper->GetTextProperty();
    tprop->SetFontFamilyToArial();
    tprop->SetFontSize(12);
    tprop->BoldOn();
	//tprop->ShadowOn();
    tprop->SetColor(1, 0, 0);
    textActor = vtkActor2D::New();
    textActor->VisibilityOff();
    textActor->SetMapper(textMapper);

	// lines to connect matches
	linespolydata = vtkPolyData::New();
	linespolydatamapper = vtkPolyDataMapper::New();
	linespolydatamapper->SetInput(linespolydata);
	linespolydataactor = vtkActor::New();
	linespolydataactor->SetMapper(linespolydatamapper);
	linespolydataactor->GetProperty()->SetColor(0,0,0);

	// lines for the tree
	treelinespolydata = vtkPolyData::New();
	treelinespolydatamapper = vtkPolyDataMapper::New();
	treelinespolydatamapper->SetInput(treelinespolydata);
	treelinespolydataactor = vtkActor::New();
	treelinespolydataactor->SetMapper(treelinespolydatamapper);
	treelinespolydataactor->GetProperty()->SetColor(0,0,0);

	// Create the Renderer, RenderWindow, and RenderWindowInteractor
    KeyPressInteractorStyle* style = KeyPressInteractorStyle::New();
	style->meshreg = this->meshreg;
    vtkCallbackCommand* pickerCommand = vtkCallbackCommand::New();
    pickerCommand->SetClientData(style);
    pickerCommand->SetCallback(PickerInteractionCallback);
    style->AddObserver(vtkCommand::LeftButtonPressEvent, pickerCommand);
    style->AddObserver(vtkCommand::MouseMoveEvent, pickerCommand);
    style->AddObserver(vtkCommand::LeftButtonReleaseEvent, pickerCommand);
    ren1 = vtkRenderer::New();

	// render window
	renWin = vtkRenderWindow::New();
    renWin->AddRenderer(ren1);
    iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);
    iren->SetInteractorStyle(style);
    iren->SetPicker(picker);

	// Add the actors to the renderer, set the background and size
    ren1->AddActor2D(textActor);
    ren1->AddActor(mesh1Actor);
    ren1->AddActor(mesh2Actor);
	ren1->AddActor(linespolydataactor);
	ren1->AddActor(treelinespolydataactor);
    ren1->SetBackground(1, 1, 1);
	ren1->TwoSidedLightingOn();
    renWin->SetSize(512,864);

	// Get the camera and zoom in closer to the image.
    //vtkCamera *cam1 = ren1->GetActiveCamera();
    //cam1->Zoom(1.4);

    iren->Initialize();
    iren->Start();

    /*picker->RemoveObserver( pickObserver );
    sphere->Delete();
    sphereMapper->Delete();
    sphereActor->Delete();
    cone->Delete();
    glyph->Delete();
    spikeMapper->Delete();
    spikeActor->Delete();
    picker->Delete();
    textMapper->Delete();
    textActor->Delete();
    pickerCommand->Delete();
    style->Delete();
    ren1->Delete();
    renWin->Delete();
    pickObserver->Delete();
	//iren->Delete();*/
}