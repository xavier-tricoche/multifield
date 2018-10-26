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
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector_types.h>
#include <vector_functions.h>
#include <cutil_inline.h>
#include <cutil_math.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPNGWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkActor.h>
#include <vtkOBJExporter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include "EasyBMP.h"
#include "GenerateTexture.h"

using namespace std;

struct int4_compare {
    bool operator() (const int4& lhs, const int4& rhs) const{
        if (lhs.x < rhs.x)
			return true;
		else if (lhs.x > rhs.x)
			return false;
		else
		{
			if (lhs.y < rhs.y)
				return true;
			else if (lhs.y > rhs.y)
				return false;
			else
			{
				 if (lhs.z < rhs.z)
					return true;
				else if (lhs.z > rhs.z)
					return false;
				else return lhs.w < rhs.w;
			}
		}
    }
};

int cint(double x);

void replaceAll(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}

int main (int argc, char *argv[])
{
	if (argc != 3)
	{
		printf("Example: vtk2objmtlexporter mesh.vtk mesh.obj\n");
		return 0;
	}

	cout << "==================\n";
	cout << "vtk2objmtlexporter\n";
	cout << "==================\n";
	
	printf("Mesh file is %s\n", argv[1]);
	
	// read the mesh file
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader->SetFileName(argv[1]);
    vtkMeshReader->Update();
	vtkPolyData* vtkMesh = vtkMeshReader->GetOutput();

	// get the points
	vtkPoints* points = vtkMesh->GetPoints();
	cout << "Number of points is " << points->GetNumberOfPoints() << "\n";
	vtkPointData* pointsData = vtkMesh->GetPointData();
	vtkUnsignedCharArray* colors = (vtkUnsignedCharArray*) vtkMesh->GetPointData()->GetScalars();	
	
	// create the obj file
	FILE * pFile;
	pFile = fopen (argv[2],"w");
	string mtlfile(argv[2]);
	replaceAll(mtlfile, string(".obj"), string(".mtl"));
	fprintf(pFile, "mtllib %s\n", mtlfile.c_str());
	
	// write the vertices
	for (int i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
	{
		double point[3];
		vtkMesh->GetPoint(i, point);
		fprintf(pFile, "v %lf %lf %lf\n", point[0], point[1], point[2]);
	}

	// do something different
	string prefix(argv[2]);
	replaceAll(prefix, string(".obj"), string(""));
	
	// check if there is color or not
	if (colors != NULL)
	{
		printf("We found colors\n");
		GenerateColorTexture(vtkMesh, prefix, pFile);
	}
	else
	{
		printf("Mesh does not have colors\n");
		vtkCellArray* vtkcells = vtkMesh->GetPolys();
		vtkcells->InitTraversal();
		vtkIdType npts;
		vtkIdType* pts;
		while (vtkcells->GetNextCell(npts, pts) != 0)
		{
			fprintf(pFile, "f %d %d %d\n", pts[0] + 1, pts[1] + 1, pts[2] + 1);
		}
	}
	
	
	// write the other mesh here
	/*vtkSmartPointer<vtkPolyDataReader> wingreader = vtkSmartPointer<vtkPolyDataReader>::New();
	wingreader->SetFileName("/home/sbarakat/nabla/Samer/Datasets/Multifield/TEMP/edeltamesh.vtk");
    wingreader->Update();
	vtkPolyData* wingmesh = wingreader->GetOutput();
	fprintf(pFile, "usemtl %sg\n", prefix.c_str());
	for (int i = 0; i < wingmesh->GetNumberOfPoints(); i++)
	{
		double point[3];
		wingmesh->GetPoint(i, point);
		fprintf(pFile, "v %lf %lf %lf\n", point[0], point[1], point[2]);
	}
	vtkCellArray* wingcells = wingmesh->GetPolys();
	wingcells->InitTraversal();
	vtkIdType wnpts;
	vtkIdType* wpts;
	while (wingcells->GetNextCell(wnpts, wpts) != 0)
	{
		fprintf(pFile, "f %d %d %d\n", wpts[0] + vtkMesh->GetNumberOfPoints() + 1, wpts[1] + vtkMesh->GetNumberOfPoints() + 1, wpts[2] + vtkMesh->GetNumberOfPoints() + 1);
	}*/
	
	
	// close the obj file
	fclose (pFile);
	
	// now output the bounding box
	double bounds[6];
	vtkMesh->GetBounds(bounds);
	pFile = fopen (string(prefix + string("B.obj")).c_str(),"w");
	fprintf(pFile, "# cube.obj\n");
	fprintf(pFile, "g cube\n");
	fprintf(pFile, "v  %lf %lf %lf\n", bounds[0], bounds[2], bounds[4]);
	fprintf(pFile, "v  %lf %lf %lf\n", bounds[0], bounds[2], bounds[5]);
	fprintf(pFile, "v  %lf %lf %lf\n", bounds[0], bounds[3], bounds[4]);
	fprintf(pFile, "v  %lf %lf %lf\n", bounds[0], bounds[3], bounds[5]);
	fprintf(pFile, "v  %lf %lf %lf\n", bounds[1], bounds[2], bounds[4]);
	fprintf(pFile, "v  %lf %lf %lf\n", bounds[1], bounds[2], bounds[5]);
	fprintf(pFile, "v  %lf %lf %lf\n", bounds[1], bounds[3], bounds[4]);
	fprintf(pFile, "v  %lf %lf %lf\n", bounds[1], bounds[3], bounds[5]);
	fprintf(pFile, "vn  0.0  0.0  1.0\n");
	fprintf(pFile, "vn  0.0  0.0 -1.0\n");
	fprintf(pFile, "vn  0.0  1.0  0.0\n");
	fprintf(pFile, "vn  0.0 -1.0  0.0\n");
	fprintf(pFile, "vn  1.0  0.0  0.0\n");
	fprintf(pFile, "vn -1.0  0.0  0.0\n");
	fprintf(pFile, "f  1//2  7//2  5//2\n");
	fprintf(pFile, "f  1//2  3//2  7//2\n");
	fprintf(pFile, "f  1//6  4//6  3//6\n"); 
	fprintf(pFile, "f  1//6  2//6  4//6\n");
	fprintf(pFile, "f  3//3  8//3  7//3\n"); 
	fprintf(pFile, "f  3//3  4//3  8//3\n");
	fprintf(pFile, "f  5//5  7//5  8//5\n");
	fprintf(pFile, "f  5//5  8//5  6//5\n");
	fprintf(pFile, "f  1//4  5//4  6//4\n");
	fprintf(pFile, "f  1//4  6//4  2//4\n");
	fprintf(pFile, "f  2//1  6//1  8//1\n");
	fprintf(pFile, "f  2//1  8//1  4//1\n");
	fclose (pFile);
	
	return 0;
}