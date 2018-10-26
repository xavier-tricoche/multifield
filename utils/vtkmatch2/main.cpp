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
// OpenGL Graphics Includes
#include <GL/glew.h>
#if defined(__APPLE__) || defined(__MACOSX)
    #include <OpenGL/OpenGL.h>
    #include <GLUT/glut.h>
#else
	//#include <windows.h>
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

#include "MeshRegistration.h"

using namespace std;


// Main program
//*****************************************************************************
int main(int argc, char** argv) 
{
	if (argc != 3)
	{
		printf("Example: vtkmatch2 mesh1.vtk mesh2.vtk\n");
		return 0;
	}

	cout << "=========\n";
	cout << "vtkmatch2\n";
	cout << "=========\n";

	printf("%s -> %s\n", argv[1], argv[2]);
	
	MeshRegistration* meshreg = new MeshRegistration(argv[1], argv[2]);
	meshreg->BuildTree();
	meshreg->DoRegistrationStep();
	meshreg->FinalFix();
	
	return 0;
}
