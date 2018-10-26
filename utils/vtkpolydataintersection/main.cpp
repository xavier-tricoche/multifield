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
#include <vtkPolyDataReader.h>
#include <vtkIntersectionPolyDataFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkPlane.h>
#include <vtkClipPolyData.h>

using namespace std;


////////////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Main entry point
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv )
{
	printf("Find intersection between two surfaces.\n");
	vtkSmartPointer<vtkPolyDataReader> reader1 = vtkSmartPointer<vtkPolyDataReader>::New();
	reader1->SetFileName("E:\\Chaos2Rost\\Tests\\Bernd\\testdata_largebox_ftle\\ridges_c.vtk");
	reader1->Update();
	vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New();
	reader2->SetFileName("E:\\Chaos2Rost\\Tests\\Bernd\\testdata_largebox_ftle\\isocontour.vtk");
	reader2->Update();

	vtkSmartPointer<vtkIntersectionPolyDataFilter> intersectionPolyDataFilter =	vtkSmartPointer<vtkIntersectionPolyDataFilter>::New();
	intersectionPolyDataFilter->SetInputConnection( 0, reader1->GetOutputPort() );
	intersectionPolyDataFilter->SetInputConnection( 1, reader2->GetOutputPort() );
	intersectionPolyDataFilter->Update();

	vtkPolyDataWriter* w = vtkPolyDataWriter::New();
	w->SetInputConnection(intersectionPolyDataFilter->GetOutputPort(0));
	w->SetFileName("E:\\Chaos2Rost\\Tests\\Bernd\\testdata_largebox_ftle\\intersections.vtk");
	w->Write();

	return 0;
}