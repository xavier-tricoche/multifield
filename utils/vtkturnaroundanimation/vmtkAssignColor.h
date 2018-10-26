#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <algorithm>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <teem/nrrd.h>
#include <vector_types.h>
#include <vector_functions.h>
#include <cutil_inline.h>
#include <cutil_math.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPolyDataReader.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkPolyDataWriter.h>
#include <vtkIdList.h>
#include <vtkColorTransferFunction.h>
#include <vtkClipPolyData.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkCellArray.h>
#include "vtkAppendPolyData.h"
#include <vtkPlane.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include <vtkLegendBoxActor.h>

#include "RegularGrid.h"

using namespace std;

void mvtkAssignColor(vtkPolyData*& vtkMesh, vtkScalarBarActor*& scalarBar);