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
/*
    FTLE computation

    This program loads a file describing a flow field and performs an integration 
	for the computation of FTLE.
*/

// OpenGL Graphics Includes
#include <GL/glew.h>
#if defined(__APPLE__) || defined(__MACOSX)
    #include <OpenGL/OpenGL.h>
    #include <GLUT/glut.h>
#else
	#include <windows.h>
    #include <GL/freeglut.h>
    #ifdef UNIX
       #include <GL/glx.h>
    #endif
#endif

// Includes
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <vector_types.h>
#include <vector_functions.h>
#include <omp.h>
#include <teem/nrrd.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

// Utilities, OpenCL and system includes
#include <oclUtils.h>
#include <shrQATest.h>

#if defined (__APPLE__) || defined(MACOSX)
   #define GL_SHARING_EXTENSION "cl_APPLE_gl_sharing"
#else
   #define GL_SHARING_EXTENSION "cl_khr_gl_sharing"
#endif

#include "MeshRegistration.h"
#include "VisWindow.h"

using namespace std;

// Main program
//*****************************************************************************
int main(int argc, char** argv) 
{

	///////////////////////////////////////////////////////////////////////////////
	//MyCGALSimplifyMesh("Y:\\sbarakat\\Projects\\cgalreducemesh\\D_0701.off", "Y:\\sbarakat\\Projects\\cgalreducemesh\\smaller.off", 0.5);

	//MeshRegistration* meshreg = new MeshRegistration("Z:\\nabla\\Samer\\Datasets\\Multifield\\bonzew3d.vtk", "Z:\\nabla\\Samer\\Datasets\\Multifield\\bonzew3d_t.vtk");
	//MeshRegistration* meshreg = new MeshRegistration("Z:\\nabla\\Samer\\Datasets\\Multifield\\Ds_0691_r2.vtk", "Z:\\nabla\\Samer\\Datasets\\Multifield\\Ds_0701_r2.vtk");
	//MeshRegistration* meshreg = new MeshRegistration("Z:\\nabla\\Samer\\Datasets\\Multifield\\comb_mf\\Ds_0090.vtk", "Z:\\nabla\\Samer\\Datasets\\Multifield\\comb_mf\\Ds_0091.vtk");
	
	//MeshRegistration* meshreg = new MeshRegistration("Z:\\nabla\\Samer\\Datasets\\Multifield\\vorts_mf\\case2\\mesh1.vtk", "Z:\\nabla\\Samer\\Datasets\\Multifield\\vorts_mf\\case2\\mesh2_t.vtk");


	MeshRegistration* meshreg = new MeshRegistration("Z:\\nabla\\Samer\\Datasets\\Multifield\\multi_mf\\Ds_0090.vtk", "Z:\\nabla\\Samer\\Datasets\\Multifield\\multi_mf\\Ds_0091.vtk");


	VisWindow* viswin = new VisWindow();
	viswin->meshreg = meshreg;
	viswin->InitWindow();

}
