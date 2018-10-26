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
#include <queue>

#include "MeshRegistration.h"

class VisWindow
{
public:

	MeshRegistration* meshreg;

	void InitWindow();
};



