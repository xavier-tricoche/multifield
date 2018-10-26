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
#include "vtkCurvatureTensor.h"


void dbprint(double** arr)
{
  for (int n = 0; n < 3; n++)
    {
    for (int m = 0; m < 3; m++)
      {cout << arr[n][m] << " ";}
    cout << endl;
    }
  cout << "**************" << endl;
}

vtkCurvatureTensor::vtkCurvatureTensor()
{
  this -> Tolerance = 1.0e-3;
}

/********************************************************
 Taubin tensor:

  M_v_i = sum_neighbours_v_j_of_v_i w_ij*k_ij*T_ij*T_ij^t

  where
    w_ij = is the sum over j of the surface areas of the cells containing
           vertices i and j, that should be twice the surface area of the cells
           containing i

    k_ij = an approximation of the  directional curvature in direction v_i - v_j
                   (v_i - v_j)
         = 2*N_i^t*------------  N_i normal at v_i
                   |v_i - v_j|

    T_ij = the projection orthogonal to N_i of v_i - v_j, in matrix terms
           (Id - N_i*N_^t)*(v_i - v_j)
         = ----------------------------
           |Id - N_i*N_^t)*(v_i - v_j)|

    here the vertex at which we compute is noted v_i, as it is a triangulated
    surface, the two others are called v_j, v_k
*******************************************************/
void vtkCurvatureTensor::GetTensor(vtkPolyData *mesh)
{
	ErrorCount = 0;
  
  //printf("Start vtkCurvatureTensor::GetTensor");

  double TOL = this->Tolerance; // will depend on the quality of vtkMath::Jacobi;
  
  if (mesh->GetNumberOfPolys()==0 || mesh->GetNumberOfPoints()==0)
      {
      printf("No points/cells to operate on");
      return;
      }
  
  vtkDataArray *normals=mesh->GetPointData()->GetNormals();
  int numPts = mesh->GetNumberOfPoints();
  
  vtkIdList* vertices, *vertices_n, *neighbours;
 
  //     create-allocate
  neighbours = vtkIdList::New();

  // vertices, tangents(edges) and normals
  double v_i[3], v_j[3], v_k[3], n_i[3], e_ij[3], en, T_ij[3];
  double T_1[3], T_2[3];
  double** M_i, *pc;
  double S_i[3][3];
  double** T_i; // the Taubin tensor
  M_i   = new double*[3];
  T_i   = new double*[3];
  //S_i   = new double*[3];
  M_i[0] = new double[3]; M_i[1] = new double[3]; M_i[2] = new double[3];
  T_i[0] = new double[3]; T_i[1] = new double[3]; T_i[2] = new double[3];
  //S_i[0] = new double[3]; S_i[1] = new double[3]; S_i[2] = new double[3];
  pc         = new double[3];

  // Data arrays
  vtkFloatArray* newScalars = vtkFloatArray::New();
  newScalars->SetName("Mean Curvature");

  vtkDoubleArray* curvatureTensor = vtkDoubleArray::New();
  curvatureTensor->SetName("Curvature Tensor");
  curvatureTensor->SetNumberOfComponents(9);
  curvatureTensor->SetNumberOfTuples(numPts);

  vtkDoubleArray* curvatureVector1 = vtkDoubleArray::New();
  curvatureVector1->SetName("Maximum Curvature Vector");
  curvatureVector1->SetNumberOfComponents(3);
  curvatureVector1->SetNumberOfTuples(numPts);

  vtkDoubleArray* curvatureVector2 = vtkDoubleArray::New();
  curvatureVector2->SetName("Minimum Curvature Vector");
  curvatureVector2->SetNumberOfComponents(3);
  curvatureVector2->SetNumberOfTuples(numPts);

  vtkDoubleArray* principalCurvatures = vtkDoubleArray::New();
  principalCurvatures->SetName("Principal Curvatures");
  principalCurvatures->SetNumberOfComponents(2);
  principalCurvatures->SetNumberOfTuples(numPts);

  // loop over vertices v_i (or loop over faces, increment?)
  //      get neighbours v_j to v_i
  unsigned short ncells;
  vtkIdType* cells, *pts, npts;
  
  mesh->BuildLinks();

  // temporary storage.
  double w_ij = 0.0, A_i = 0.0, a_i = 0.0, l_ij = 0.0, kappa_ij = 0.0,\
    k_1 = 0.0, k_2 = 0.0, m_1 = 0.0, m_2 = 0.0;
  double d0[3], d1[3], d2[3], c0[3], c1[3], c2[3];
  for (int i = 0; i < numPts; i++)
  {
    //printf(<< "vertex " << i);
    
    mesh->GetPointCells(i,ncells,cells);
    mesh->GetPoint(i,v_i);
    normals->GetTuple(i,n_i);
    
    // resetting stuff for next vertex
    neighbours->Reset(); A_i = 0.0;
    for (int n = 0; n < 3; n++) for (int m = 0; m < 3; m++) M_i[n][m] = 0.0;
    
    // explicit computation, on neighbours, etc
    for (int c = 0; c < ncells; c++)
    {
      mesh->GetCellPoints(cells[c],npts,pts);
      for (int j = 0; j < npts; j++)
      {
	if (pts[j] != i)
	{
	  neighbours->InsertNextId(pts[j]);
	  //printf(<< "\tadding " << pts[j] << " to neighbours");
	}
      } // loop over vertices of cell
    } // loop over cells: neighbours should have all been found
      
    for (int j = 0; j < neighbours -> GetNumberOfIds(); j++)
    {
      mesh->GetPoint(neighbours->GetId(j),  v_j);
      mesh->GetPoint(neighbours->GetId((j+1) % neighbours -> GetNumberOfIds() ),v_k);  // use mod? or not necessary
      w_ij = vtkTriangle::TriangleArea(v_i,v_j,v_k);
      A_i += w_ij;
      
      e_ij[0] = v_i[0] - v_j[0];
      e_ij[1] = v_i[1] - v_j[1];
      e_ij[2] = v_i[2] - v_j[2];
      
      
      /*--------------------------------------------
	T = (I - n_i*n_i^t)* e_ij / norm
	= e_ij - (e_ij*n_i)*n_i =: e_ij - en*n_i;
	--------------------------------------------*/
      
      en = e_ij[0]*n_i[0] + e_ij[1]*n_i[1] + e_ij[2]*n_i[2];
      T_ij[0] = (e_ij[0] - en*n_i[0]);
      T_ij[1] = (e_ij[1] - en*n_i[1]);
      T_ij[2] = (e_ij[2] - en*n_i[2]);
      vtkMath::Normalize(T_ij);

      // normal curvature in direction i to j: Eq.(6) in Taubin
      kappa_ij = -2.0*en / vtkMath::Norm(e_ij);
      
      // contribution from neighbour j: Eq.(2) in Taubin
      M_i[0][0] += w_ij*kappa_ij*T_ij[0]*T_ij[0];
      M_i[0][1] += w_ij*kappa_ij*T_ij[0]*T_ij[1];
      M_i[1][0] += w_ij*kappa_ij*T_ij[1]*T_ij[0];
      M_i[1][1] += w_ij*kappa_ij*T_ij[1]*T_ij[1];
      M_i[0][2] += w_ij*kappa_ij*T_ij[0]*T_ij[2];
      M_i[2][0] += w_ij*kappa_ij*T_ij[2]*T_ij[0];
      M_i[2][2] += w_ij*kappa_ij*T_ij[2]*T_ij[2];
      M_i[1][2] += w_ij*kappa_ij*T_ij[1]*T_ij[2];
      M_i[2][1] += w_ij*kappa_ij*T_ij[2]*T_ij[1];
       
    } // loop over neighbours of v
    
      // Normalise
    if (fabs(A_i) > 2*VTK_DOUBLE_MIN) 
    {
      a_i = 1.0 / A_i;
    }
    else
    {
      printf("total surface area of facets less than tol");
    }
      
    for (int n = 0; n < 3; n++)
    for (int m = 0; m < 3; m++)
    {
      M_i[n][m] *= a_i;
    }

    // Diagonalisation
    if (!vtkMath::Jacobi(M_i,pc,T_i))
	{
		// update data
		curvatureTensor->SetTuple9(i,0,0,0,0,0,0,0,0,0);
		curvatureVector1->SetTuple3(i,0,0,0);
		curvatureVector2->SetTuple3(i,0,0,0);
		principalCurvatures->SetTuple2(i,0,0);

		newScalars->InsertValue(i, 0.0);

		ErrorCount++;
		continue;
	}
    
    d0[0] = T_i[0][0]; d0[1] = T_i[1][0]; d0[2] = T_i[2][0];
    d1[0] = T_i[0][1]; d1[1] = T_i[1][1]; d1[2] = T_i[2][1]; 
    d2[0] = T_i[0][2]; d2[1] = T_i[1][2]; d2[2] = T_i[2][2];
    vtkMath::Cross(d0,n_i,c0);
    vtkMath::Cross(d1,n_i,c1);
    vtkMath::Cross(d2,n_i,c2);
    if (vtkMath::Norm(c0) < TOL)
    {
      m_1 = pc[1]; m_2 = pc[2];
      T_1[0] = d1[0]; T_1[1] = d1[1]; T_1[2] = d1[2]; 
      T_2[0] = d2[0]; T_2[1] = d2[1]; T_2[2] = d2[2]; 
    }
    else if (vtkMath::Norm(c1) < TOL)
    {
      m_1 = pc[0]; m_2 = pc[2];
      T_1[0] = d0[0]; T_1[1] = d0[1]; T_1[2] = d0[2]; 
      T_2[0] = d2[0]; T_2[1] = d2[1]; T_2[2] = d2[2];
    }
    else if (vtkMath::Norm(c2) < TOL)
    {
      m_1 = pc[0]; m_2 = pc[1];
      T_1[0] = d0[0]; T_1[1] = d0[1]; T_1[2] = d0[2]; 
      T_2[0] = d1[0]; T_2[1] = d1[1]; T_2[2] = d1[2];
    }
    else
    {
      printf("No eigenvector parallel to normal, keeping previous vertex's");
	  ErrorCount++;
    }
    
    // Equation (5) in Taubin:
    // NB: m_1 >= m_2 according to vtkMath::Jacobi ==> k_1 >= k_2
    k_1 = 3.0*m_1 - m_2;
    k_2 = 3.0*m_2 - m_1;
    
    // the curvature tensor is the symmetric tensor with
    // eigenvectors: {N   T_1  T_2}					\
    //                              ==> S = k_1*T_1*T_1^t + k_2*T_2*T_2^t + 0*N*N^t;
    // eigenvalues:  {0   k_1  k_2}/
    S_i[0][0] = k_1*T_1[0]*T_1[0] + k_2*T_2[0]*T_2[0];
    S_i[1][1] = k_1*T_1[1]*T_1[1] + k_2*T_2[1]*T_2[1];
    S_i[2][2] = k_1*T_1[2]*T_1[2] + k_2*T_2[2]*T_2[2];
    S_i[0][1] = k_1*T_1[0]*T_1[1] + k_2*T_2[0]*T_2[1]; S_i[1][0] = S_i[0][1];
    S_i[0][2] = k_1*T_1[0]*T_1[2] + k_2*T_2[0]*T_2[2]; S_i[2][0] = S_i[0][2];
    S_i[1][2] = k_1*T_1[1]*T_1[2] + k_2*T_2[1]*T_2[2]; S_i[2][1] = S_i[1][2];
    
    // update data
    curvatureTensor->SetTuple9(i,S_i[0][0],S_i[0][1],S_i[0][2],		\
			         S_i[1][0],S_i[1][1],S_i[1][2],		\
			         S_i[2][0],S_i[2][1],S_i[2][2]);
    curvatureVector1->SetTuple3(i,T_1[0],T_1[1],T_1[2]);
    curvatureVector2->SetTuple3(i,T_2[0],T_2[1],T_2[2]);
    principalCurvatures->SetTuple2(i,k_1,k_2);

	newScalars->InsertValue(i, 0.5 * abs(k_1 + k_2));
    
  } // loop over vertices v_i
  
  mesh->GetPointData()->AddArray(newScalars);
  mesh->GetPointData()->AddArray(curvatureTensor);
  mesh->GetPointData()->AddArray(curvatureVector1);
  mesh->GetPointData()->AddArray(curvatureVector2);
  mesh->GetPointData()->AddArray(principalCurvatures);

  // cleaning:
  //if (principalCurvatures) { principalCurvatures->Delete(); }
  //if (curvatureVector1)    { curvatureVector1->Delete(); }
  //if (curvatureVector2)    { curvatureVector2->Delete(); }
  //if (curvatureTensor)     { curvatureTensor->Delete(); }
  //if (neighbours)          { neighbours->Delete(); }

  if (pc) 
  {
    delete [] pc; pc = 0;
  }
  if (M_i)
  {
    for (int n = 0; n < 3; n++) { delete [] M_i[n]; M_i[n] = 0; }
    delete [] M_i; M_i = 0;
  }
  if (T_i)
  {
    for (int n = 0; n < 3; n++) { delete [] T_i[n]; T_i[n] = 0; }
    delete [] T_i; T_i = 0;
  }

  printf("Total errors is %d\n", ErrorCount);
}

void vtkMesh2Obj(vtkPolyData* mesh, char* file)
{
	vtkRenderer* ren = vtkRenderer::New();
	vtkRenderWindow* renWin = vtkRenderWindow::New(); 
	renWin->AddRenderer(ren);
	vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
	mapper->SetInput(mesh);
	vtkActor* ma = vtkActor::New();
	ma->SetMapper(mapper);
	ren->AddActor(ma);

	vtkOBJExporter* obj = vtkOBJExporter::New();
	obj->SetInput(renWin);
	obj->SetFilePrefix(file);
	obj->Write();
}

void vtkCurvatureTensor::GetTensor2(vtkPolyData *mesh)
{
	double* bounds = mesh->GetBounds();
	vtkMesh2Obj(mesh, "tmpobj");

	// allocate a new mesh
	m_pMesh = new Enriched_polyhedron<Enriched_kernel,Enriched_items>;
	CGAL_assertion(m_pMesh != NULL);

	// obj extension
	Parser_obj<Enriched_kernel,Enriched_items> parser;
	parser.read("tmpobj.obj",m_pMesh);

	// update mesh properties in the status bar 
	// and compute normals
	m_pMesh->compute_type();
	m_pMesh->compute_normals();

	// estimate curvature tensors
	typedef CCurvature_estimator<Mesh,Enriched_kernel> Estimator;
	Estimator estimator(m_pMesh);
	estimator.run(0);

	// Data arrays
	vtkFloatArray* newScalars = vtkFloatArray::New();
	newScalars->SetName("Mean Curvature");
	int numPts = mesh->GetNumberOfPoints();

	vtkDoubleArray* curvatureVector1 = vtkDoubleArray::New();
	curvatureVector1->SetName("Maximum Curvature Vector");
	curvatureVector1->SetNumberOfComponents(3);
	curvatureVector1->SetNumberOfTuples(numPts);

	vtkDoubleArray* curvatureVector2 = vtkDoubleArray::New();
	curvatureVector2->SetName("Minimum Curvature Vector");
	curvatureVector2->SetNumberOfComponents(3);
	curvatureVector2->SetNumberOfTuples(numPts);

	vtkDoubleArray* principalCurvatures = vtkDoubleArray::New();
	principalCurvatures->SetName("Principal Curvatures");
	principalCurvatures->SetNumberOfComponents(2);
	principalCurvatures->SetNumberOfTuples(numPts);

	// range info
	double mincur = numeric_limits<double>::max();
	double maxcur = numeric_limits<double>::min();

	// fill data
	Vertex_iterator hVertex = NULL;
	int i = 0;
	for(hVertex = m_pMesh->vertices_begin(); hVertex != m_pMesh->vertices_end(); hVertex++)
	{
		Curvature& curvature = hVertex->curvature();
		
		// max vector
		float3 dir = make_float3(curvature.m_direction_kmax.x(), curvature.m_direction_kmax.y(), curvature.m_direction_kmax.z());
		dir = normalize(dir);
		curvatureVector1->SetTuple3(i,dir.x, dir.y, dir.z);

		// min vector
		dir = make_float3(curvature.m_direction_kmin.x(), curvature.m_direction_kmin.y(), curvature.m_direction_kmin.z());
		dir = normalize(dir);
		curvatureVector2->SetTuple3(i,dir.x, dir.y, dir.z);
		
		// max and min values
		principalCurvatures->SetTuple2(i, curvature.kmax(), curvature.kmin());

		newScalars->InsertValue(i, 0.5 * abs(curvature.kmax() + curvature.kmin()));

		mincur = min(mincur, curvature.kmin());
		mincur = min(mincur, curvature.kmax());
		maxcur = max(maxcur, curvature.kmin());
		maxcur = max(maxcur, curvature.kmax());

		i++;
	}

	printf("Curvatures [%lf , %lf]\n", mincur, maxcur);

	mesh->GetPointData()->AddArray(newScalars);
	mesh->GetPointData()->AddArray(curvatureVector1);
	mesh->GetPointData()->AddArray(curvatureVector2);
	mesh->GetPointData()->AddArray(principalCurvatures);
}

int vtkCurvatureTensor::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and ouptut
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Null input check
  if (!input)
    {
    return 0;
    }
  output->CopyStructure(input);
  output->GetPointData()->PassData(input->GetPointData());
  output->GetFieldData()->PassData(input->GetFieldData());

  this->GetTensor(output);
  
  return 1;
}

void vtkCurvatureTensor::PrintSelf(ostream& os, vtkIndent indent)
{
	//this->Superclass::PrintSelf(os,indent);
	//os << indent << "CurvatureType: " << this->CurvatureType << "\n";
	//os << indent << "InvertMeanCurvature: " << this->InvertMeanCurvature << "\n";
}