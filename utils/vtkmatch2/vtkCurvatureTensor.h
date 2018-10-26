

#ifndef VTKCURVATURETENSOR
#define VTKCURVATURETENSOR

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>

#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include "Enriched_polyhedron.h"
#include "parser_obj.h"
#include "curvature_estimator.h"
#include "curvature.h"

#include "vtkCellArray.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolygon.h"
#include "vtkTensor.h"
#include "vtkTriangle.h"
#include "vtkDataArray.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkPolyDataMapper.h"
#include "vtkOBJExporter.h"

#include <vector_types.h>
#include <vector_functions.h>
#include <driver_functions.h>
#include <cutil_inline.h>
#include <cutil_math.h>

using namespace std;

typedef double number_type;
typedef CGAL::Simple_cartesian<number_type> Enriched_kernel;
typedef Enriched_polyhedron<Enriched_kernel,Enriched_items> Mesh;
typedef Mesh::Vertex_iterator Vertex_iterator;
typedef CCurvature<Enriched_kernel> Curvature;

class vtkCurvatureTensor
{
public:
	double Tolerance;
	Mesh *m_pMesh;
	int ErrorCount;

	vtkCurvatureTensor();
	void GetTensor(vtkPolyData *mesh);
	void GetTensor2(vtkPolyData *mesh);
	int RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector);
	void PrintSelf(ostream& os, vtkIndent indent);
};
#endif
