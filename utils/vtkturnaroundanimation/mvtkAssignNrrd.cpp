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
#include "mvtkAssignNrrd.h"

Nrrd *readNrrd(const char* filename)
{
    Nrrd *nin = nrrdNew();
    if (nrrdLoad(nin, filename, NULL)) {
        printf("Error reading %s\n",filename);
                return NULL;
    }
    return nin;
}

void mvtkAssignNrrd(vtkPolyData* vtkMesh, string tblname, string nrrdfilename)
{
	cout << "=============\n";
	cout << "vtkassignnrrd\n";
	cout << "=============\n";
	
	Nrrd* nrrdfile = readNrrd(nrrdfilename.c_str());
	RegularGrid* dataset = new RegularGrid(nrrdfile, 
										   nrrdfile->axis[0].size, nrrdfile->axis[1].size, nrrdfile->axis[2].size, 
										   nrrdfile->axis[0].spacing, nrrdfile->axis[1].spacing, nrrdfile->axis[2].spacing);
	
	// get the points
	vtkPoints* points = vtkMesh->GetPoints();
	cout << "Number of points is " << points->GetNumberOfPoints() << "\n";
	vtkPointData* pointsData = vtkMesh->GetPointData();
	vtkSmartPointer<vtkFloatArray> newScalars = vtkSmartPointer<vtkFloatArray>::New();
	newScalars->SetName(tblname.c_str());
	float minv = numeric_limits<float>::max();
	float maxv = numeric_limits<float>::min();
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
	{
		double* point = points->GetPoint(i);
		float value = (float) (dataset->GetValueAt(point[0], point[1], point[2]));
		minv = min(minv, value);
		maxv = max(maxv, value);
		
		newScalars->InsertValue(i, value);
	}
	pointsData->AddArray(newScalars);
	printf("min is %f and max is %f\n", minv, maxv);
	
	// free nrrd
	delete dataset;
}