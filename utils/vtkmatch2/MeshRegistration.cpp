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
#include "MeshRegistration.h"

bool myisnan(double var)
{
    volatile double temp = var;
    return temp != temp;
}

bool myisinf(double var)
{
	return (numeric_limits<double>::infinity() == var) || (numeric_limits<double>::infinity() == -var);
}

bool myisnan(float var)
{
    volatile float temp = var;
    return temp != temp;
}

bool myisinf(float var)
{
	return (numeric_limits<float>::infinity() == var) || (numeric_limits<float>::infinity() == -var);
}

void quicksortswap(int *x,int *y)
{
   int temp;
   temp = *x;
   *x = *y;
   *y = temp;
}
	 
int quicksortswapchoose_pivot(int i,int j )
{
   return((i+j) /2);
}
	 
void quicksort(BoundaryDistSorter& bs, vector<int>& list,int m,int n)
{
   int key,i,j,k;
   if( m < n)
   {
      k = quicksortswapchoose_pivot(m,n);
      quicksortswap(&list[m],&list[k]);
	  key = bs.vals[m];
      i = m+1;
      j = n;
      while(i <= j)
      {
         while((i <= n) && (bs.vals[i] <= key))
                i++;
         while((j >= m) && (bs.vals[j] > key))
                j--;
         if( i < j)
                quicksortswap(&list[i],&list[j]);
      }
      // swap two elements
      quicksortswap(&list[m],&list[j]);
      // recursively sort the lesser list
      quicksort(bs, list,m,j-1);
      quicksort(bs, list,j+1,n);
   }
}

float3 GetFloat3Point(vtkPolyData* vtkMesh, const int& id)
{
	double pt[3];
	vtkMesh->GetPoint(id, pt);

	return make_float3(pt[0], pt[1], pt[2]);
}

float3 GetFloat3Normal(vtkPolyData* vtkMesh, const int& id)
{
	double pt[3];
	vtkMesh->GetPointData()->GetNormals()->GetTuple(id, pt);

	return normalize(make_float3(pt[0], pt[1], pt[2]));
}

double GetEdgeLength(vtkPolyData* vtkMesh, int id1, int id2)
{
	double pt1[3];
	vtkMesh->GetPoint(id1, pt1);

	double pt2[3];
	vtkMesh->GetPoint(id2, pt2);

	pt1[0] -= pt2[0];
	pt1[1] -= pt2[1];
	pt1[2] -= pt2[2];

	return sqrt(pt1[0] * pt1[0] + pt1[1] * pt1[1] + pt1[2] * pt1[2]);
}

double GetEdgeLength(vtkPolyData* vtkMesh, Edge edge)
{
	double pt1[3];
	vtkMesh->GetPoint(edge.x, pt1);

	double pt2[3];
	vtkMesh->GetPoint(edge.y, pt2);

	pt1[0] -= pt2[0];
	pt1[1] -= pt2[1];
	pt1[2] -= pt2[2];

	return sqrt(pt1[0] * pt1[0] + pt1[1] * pt1[1] + pt1[2] * pt1[2]);
}

int ClosestNeighbot(vtkPolyData* vtkMesh, int from, int to1, int to2)
{
	if (to1 == -1)
		return to2;
	if (to2 == -1)
		return to1;

	if (GetEdgeLength(vtkMesh, Edge(from, to1)) < GetEdgeLength(vtkMesh, Edge(from, to2)))
		return to1;
	else
		return to2;
}

float3 ClosestPointOnSegment(float3 P, float3 S0, float3 S1)
{
    float3 v = S1 - S0;
    float3 w = P - S0;

    double c1 = dot(w,v);
    if ( c1 <= 0 )
        return S0;

    double c2 = dot(v,v);
    if ( c2 <= c1 )
        return S1;

    double b = c1 / c2;
    float3 Pb = S0 + b * v;
    return Pb;
}

float3 ClosestPointOnTriangle(const float3& mPoint, const float3& V0, const float3& V1, const float3& V2)
{
	float l1 = length(V0 - V1);
	float l2 = length(V0 - V2);
	float l3 = length(V1 - V2);
	if ((l1 < 0.00001) || (l2 < 0.00001) || (l3 < 0.00001))
	{
		if ((l1 < 0.00001) && (l2 < 0.00001) && (l3 < 0.00001))
		{
			return V0;
		}

		if (l1 < 0.00001)
		{
			return ClosestPointOnSegment(mPoint, V0, V2);
		}
		else if (l2 < 0.00001)
		{
			return ClosestPointOnSegment(mPoint, V0, V1);
		}
		else if (l3 < 0.00001)
		{
			return ClosestPointOnSegment(mPoint, V0, V1);
		}
	}

    float3 diff = V0 - mPoint;
    float3 edge0 = V1 - V0;
    float3 edge1 = V2 - V0;
    double a00 = dot(edge0,edge0);
    double a01 = dot(edge0,edge1);
    double a11 = dot(edge1,edge1);
    double b0 = dot(diff,edge0);
    double b1 = dot(diff,edge1);
    double c = dot(diff,diff);
    double det = abs(a00*a11 - a01*a01);
    double s = a01*b1 - a11*b0;
    double t = a01*b0 - a00*b1;
    double sqrDistance;
	
	if (abs(dot(normalize(edge0), normalize(edge1))) > 0.999)
	{
		return ClosestPointOnSegment(mPoint, V0, V1);
	}
	
    if (s + t <= det)
    {
        if (s < (double)0)
        {
            if (t < (double)0)  // region 4
            {
                if (b0 < (double)0)
                {
                    t = (double)0;
                    if (-b0 >= a00)
                    {
                        s = (double)1;
                        sqrDistance = a00 + ((double)2)*b0 + c;
                    }
                    else
                    {
                        s = -b0/a00;
                        sqrDistance = b0*s + c;
                    }
                }
                else
                {
                    s = (double)0;
                    if (b1 >= (double)0)
                    {
                        t = (double)0;
                        sqrDistance = c;
                    }
                    else if (-b1 >= a11)
                    {
                        t = (double)1;
                        sqrDistance = a11 + ((double)2)*b1 + c;
                    }
                    else
                    {
                        t = -b1/a11;
                        sqrDistance = b1*t + c;
                    }
                }
            }
            else  // region 3
            {
                s = (double)0;
                if (b1 >= (double)0)
                {
                    t = (double)0;
                    sqrDistance = c;
                }
                else if (-b1 >= a11)
                {
                    t = (double)1;
                    sqrDistance = a11 + ((double)2)*b1 + c;
                }
                else
                {
                    t = -b1/a11;
                    sqrDistance = b1*t + c;
                }
            }
        }
        else if (t < (double)0)  // region 5
        {
            t = (double)0;
            if (b0 >= (double)0)
            {
                s = (double)0;
                sqrDistance = c;
            }
            else if (-b0 >= a00)
            {
                s = (double)1;
                sqrDistance = a00 + ((double)2)*b0 + c;
            }
            else
            {
                s = -b0/a00;
                sqrDistance = b0*s + c;
            }
        }
        else  // region 0
        {
            // minimum at interior point
            double invDet = ((double)1)/det;
            s *= invDet;
            t *= invDet;
            sqrDistance = s*(a00*s + a01*t + ((double)2)*b0) +
                t*(a01*s + a11*t + ((double)2)*b1) + c;
        }
    }
    else
    {
        double tmp0, tmp1, numer, denom;

        if (s < (double)0)  // region 2
        {
            tmp0 = a01 + b0;
            tmp1 = a11 + b1;
            if (tmp1 > tmp0)
            {
                numer = tmp1 - tmp0;
                denom = a00 - ((double)2)*a01 + a11;
                if (numer >= denom)
                {
                    s = (double)1;
                    t = (double)0;
                    sqrDistance = a00 + ((double)2)*b0 + c;
                }
                else
                {
                    s = numer/denom;
                    t = (double)1 - s;
                    sqrDistance = s*(a00*s + a01*t + ((double)2)*b0) +
                        t*(a01*s + a11*t + ((double)2)*b1) + c;
                }
            }
            else
            {
                s = (double)0;
                if (tmp1 <= (double)0)
                {
                    t = (double)1;
                    sqrDistance = a11 + ((double)2)*b1 + c;
                }
                else if (b1 >= (double)0)
                {
                    t = (double)0;
                    sqrDistance = c;
                }
                else
                {
                    t = -b1/a11;
                    sqrDistance = b1*t + c;
                }
            }
        }
        else if (t < (double)0)  // region 6
        {
            tmp0 = a01 + b1;
            tmp1 = a00 + b0;
            if (tmp1 > tmp0)
            {
                numer = tmp1 - tmp0;
                denom = a00 - ((double)2)*a01 + a11;
                if (numer >= denom)
                {
                    t = (double)1;
                    s = (double)0;
                    sqrDistance = a11 + ((double)2)*b1 + c;
                }
                else
                {
                    t = numer/denom;
                    s = (double)1 - t;
                    sqrDistance = s*(a00*s + a01*t + ((double)2)*b0) +
                        t*(a01*s + a11*t + ((double)2)*b1) + c;
                }
            }
            else
            {
                t = (double)0;
                if (tmp1 <= (double)0)
                {
                    s = (double)1;
                    sqrDistance = a00 + ((double)2)*b0 + c;
                }
                else if (b0 >= (double)0)
                {
                    s = (double)0;
                    sqrDistance = c;
                }
                else
                {
                    s = -b0/a00;
                    sqrDistance = b0*s + c;
                }
            }
        }
        else  // region 1
        {
            numer = a11 + b1 - a01 - b0;
            if (numer <= (double)0)
            {
                s = (double)0;
                t = (double)1;
                sqrDistance = a11 + ((double)2)*b1 + c;
            }
            else
            {
                denom = a00 - ((double)2)*a01 + a11;
                if (numer >= denom)
                {
                    s = (double)1;
                    t = (double)0;
                    sqrDistance = a00 + ((double)2)*b0 + c;
                }
                else
                {
                    s = numer/denom;
                    t = (double)1 - s;
                    sqrDistance = s*(a00*s + a01*t + ((double)2)*b0) +
                        t*(a01*s + a11*t + ((double)2)*b1) + c;
                }
            }
        }
    }

    // Account for numerical round-off error.
    if (sqrDistance < (double)0)
    {
        sqrDistance = (double)0;
    }

    float3 mClosestPoint0 = mPoint;
    float3 mClosestPoint1 = V0 + s*edge0 + t*edge1;

	float3 mTriangleBary;
    mTriangleBary.y = s;
    mTriangleBary.z = t;
    mTriangleBary.x = (double)1 - s - t;

    return mClosestPoint1;
}

float3 ApplyTranformation(const float3& pt, float Tg[4][4], float Tl[4][4])
{
	float3 ret;
	float3 pt2;

	// apply global transformation
	pt2.x = pt.x * Tg[0][0] + pt.y * Tg[1][0] + pt.z * Tg[2][0] + Tg[3][0];
	pt2.y = pt.x * Tg[0][1] + pt.y * Tg[1][1] + pt.z * Tg[2][1] + Tg[3][1];
	pt2.z = pt.x * Tg[0][2] + pt.y * Tg[1][2] + pt.z * Tg[2][2] + Tg[3][2];
	pt2 /= (pt.x * Tg[0][3] + pt.y * Tg[1][3] + pt.z * Tg[2][3] + Tg[3][3]);

	// apply local transformation
	ret.x = pt2.x * Tl[0][0] + pt2.y * Tl[1][0] + pt2.z * Tl[2][0] + Tl[3][0];
	ret.y = pt2.x * Tl[0][1] + pt2.y * Tl[1][1] + pt2.z * Tl[2][1] + Tl[3][1];
	ret.z = pt2.x * Tl[0][2] + pt2.y * Tl[1][2] + pt2.z * Tl[2][2] + Tl[3][2];
	ret /= (pt2.x * Tl[0][3] + pt2.y * Tl[1][3] + pt2.z * Tl[2][3] + Tl[3][3]);

	return ret;
}

float3 ApplyTranformation(const float3& pt, float T[4][4])
{
	float3 ret;

	// apply global transformation
	ret.x = pt.x * T[0][0] + pt.y * T[1][0] + pt.z * T[2][0] + T[3][0];
	ret.y = pt.x * T[0][1] + pt.y * T[1][1] + pt.z * T[2][1] + T[3][1];
	ret.z = pt.x * T[0][2] + pt.y * T[1][2] + pt.z * T[2][2] + T[3][2];
	ret /= (pt.x * T[0][3] + pt.y * T[1][3] + pt.z * T[2][3] + T[3][3]);

	return ret;
}

int InverseTransform(const Candidate& ca, float invT[4][4])
{
	float T[4][4];
	int i,j,k;
	for (i=0; i<=3; i++)
	{
		for (j=0; j<=3; j++)
		{
			float sum = 0;
			for (k=0; k<=3; k++)
				sum = sum + ca.Tg[i][k] * ca.Tl[k][j];
			T[i][j] = sum;
		}
	}


	// allocate data
	gsl_matrix* A = gsl_matrix_alloc(4, 4);
	gsl_matrix* invA = gsl_matrix_alloc(4, 4);
	gsl_permutation *p = gsl_permutation_alloc(4);

	// set matrix
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			gsl_matrix_set(A, i, j, T[i][j]);

	// compute the inverse
	int signum;
	int err = gsl_linalg_LU_decomp(A, p, &signum);
	err |= gsl_linalg_LU_invert(A, p, invA);
	if (err > 0)
	{
		printf("Failed to compute the inverse!\n");
	}

	// get result
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			invT[i][j] = gsl_matrix_get(invA, i, j);

	// free memory
	gsl_matrix_free(A);
	gsl_matrix_free(invA);
	gsl_permutation_free(p);

	return err;
}

void matrixFromAxisAngle(float3 axis, float angle, float T[4][4]) 
{
	//http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToMatrix/index.htm

    double c = cos(angle);
    double s = sin(angle);
    double t = 1.0 - c;
	//  if axis is not already normalised then uncomment this
	// double magnitude = Math.sqrt(a1.x*a1.x + a1.y*a1.y + a1.z*a1.z);
	// if (magnitude==0) throw error;
	// a1.x /= magnitude;
	// a1.y /= magnitude;
	// a1.z /= magnitude;

    T[0][0] = c + axis.x*axis.x*t;
    T[1][1] = c + axis.y*axis.y*t;
    T[2][2] = c + axis.z*axis.z*t;


    double tmp1 = axis.x*axis.y*t;
    double tmp2 = axis.z*s;
    T[1][0] = tmp1 + tmp2;
    T[0][1] = tmp1 - tmp2;
    tmp1 = axis.x*axis.z*t;
    tmp2 = axis.y*s;
    T[2][0] = tmp1 - tmp2;
    T[0][2] = tmp1 + tmp2;    
	tmp1 = axis.y*axis.z*t;
    tmp2 = axis.x*s;
    T[2][1] = tmp1 + tmp2;
    T[1][2] = tmp1 - tmp2;


	T[0][3] = 0.0;
	T[1][3] = 0.0;
	T[2][3] = 0.0;
}

void multiplyMatrix(float T1[4][4], float T2[4][4], float T[4][4])
{
	int i,j,k;
	for (i=0; i<=3; i++)
	{
		for (j=0; j<=3; j++)
		{
			float sum = 0;
			for (k=0; k<=3; k++)
				sum = sum + T1[i][k] * T2[k][j];
			T[i][j] = sum;
		}
	}
}

void FindPrincipalDirections(vector<float3>* pts, vector<float3>* directions)
{
	int n = pts->size();
	int m = 3;

	// create matrix for the principal components
	gsl_matrix* pc = gsl_matrix_alloc(m, m);

	// create matrix for the scores
	gsl_matrix* score = gsl_matrix_alloc(m, n);

	// create matrix for the data
	gsl_matrix* data = gsl_matrix_alloc(n, m);

	// put data in the matrix
	for (int i = 0; i < n; i ++)
	{
		gsl_matrix_set(data, i, 0, pts->at(i).x);
		gsl_matrix_set(data, i, 1, pts->at(i).y);
		gsl_matrix_set(data, i, 2, pts->at(i).z);
	}

	// call the PCA function
	PCA::FindPCA(pc, score, data);

	// get results
	for (int  i = 0; i < 3; i++)
	{
		float3 dir;
		dir.x = gsl_matrix_get(pc, 0, i);
		dir.y = gsl_matrix_get(pc, 1, i);
		dir.z = gsl_matrix_get(pc, 2, i);
		
		directions->push_back(dir);
	}

	// free memory
	gsl_matrix_free(pc);
	gsl_matrix_free(score);
	gsl_matrix_free(data);
}

void SolveLinearSystem(int n, double a_data[], double b_data[], double res[])
{
	gsl_matrix_view m = gsl_matrix_view_array (a_data, n, n);
     
	gsl_vector_view b = gsl_vector_view_array (b_data, n);
     
	gsl_vector *x = gsl_vector_alloc (n);
       
	int s;
     
	gsl_permutation * p = gsl_permutation_alloc (n);
     
	gsl_linalg_LU_decomp (&m.matrix, p, &s);
     
	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
     
	for (int i = 0; i < n; i++)
		res[i] = gsl_vector_get(x, i);
     
	gsl_permutation_free (p);
	gsl_vector_free (x);
}

bool FindRotationAnglesPCA(vector<float3>* pts1, vector<float3>* pts2, float T[4][4])
{
	vector<float3> dir1;
	FindPrincipalDirections(pts1, &dir1);

	vector<float3> dir2;
	FindPrincipalDirections(pts2, &dir2);

	// first normalize
	for (int i = 0; i < 3; i++)
	{
		dir1[i] = normalize(dir1[i]);
		dir2[i] = normalize(dir2[i]);
	}

	// make sure dir1 directions are correct
	if (dot(dir1[2], cross(dir1[0], dir1[1])) < 0.0)
	{
		dir1[2] = -dir1[2];
	}

	// make sure dir2 directions are correct
	if (dot(dir2[2], cross(dir2[0], dir2[1])) < 0.0)
	{
		dir2[2] = -dir2[2];
	}

	// create 3x3 distance matrix
	float mat[3][3];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			mat[i][j] = abs(dot(dir1[i], dir2[j]));
		}
	}

	// find best correspondences
	vector<int2> mapping;
	for (int i = 0; i < 3; i++)
	{
		// loop on array to find highest item 
		float maxv = 0.0;
		int2 coord;
		for (int r = 0; r < 3; r++)
		{
			for (int s = 0; s < 3; s++)
			{
				if (mat[r][s] > maxv)
				{
					maxv = mat[r][s];
					coord.x = r;
					coord.y = s;
				}
			}
		}
		mapping.push_back(coord);

		// set the row and colomn of coord to zero
		for (int r = 0; r < 3; r++)
			mat[coord.x][r] = 0.0;
		for (int r = 0; r < 3; r++)
			mat[r][coord.y] = 0.0;
	}

	// set dir2 correctly
	vector<float3> dir3 = dir2;
	for (int i = 0; i < 3; i++)
	{
		int2 coord = mapping[i];
		dir2[coord.x] = dir3[coord.y];
	}

	// check for the confusion case
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (i == j) continue;
			if (abs(dot(dir1[i], dir1[j])) > 0.000001)
			{
				printf("h %f\n", dot(dir1[i], dir1[j]));
			}
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (i == j) continue;
			if (abs(dot(dir2[i], dir2[j])) > 0.000001)
			{
				printf("h %f\n", dot(dir2[i], dir2[j]));
			}
		}
	}

	// make sure dir1 directions are correct
	if (dot(dir1[2], cross(dir1[0], dir1[1])) < 0.0)
	{
		dir1[2] = -dir1[2];
	}

	// now adjust dir2 according to dir 1
	float ang[3];
	for (int i = 0; i < 3; i++)
	{
		if (dot(dir1[i], dir2[i]) < 0.0)
		{
			dir2[i] = -dir2[i];
		}

		double an = acos(dot(dir1[i], dir2[i])) * 180.0 / 3.14159265;
		printf("An %f\n", an);
		ang[i] = an;
	}
	float minang = min(ang[0], min(ang[1], ang[2]));
	int minangi = (ang[0] < ang[1])? 0 : 1;
	minangi = (ang[minangi] < ang[2])? minangi : 2;

	/*for (int i = 0; i < 3; i++)
	{
		if (i == minangi) continue;
		float3 vec = cross(dir1[i], dir2[minangi]);
		vec = cross(dir2[minangi], vec);
		dir2[i] = vec;
	}*/

	for (int i = 0; i < 3; i++)
	{
		if (dot(dir1[i], dir2[i]) < 0.0)
		{
			dir2[i] = -dir2[i];
		}

		double an = acos(dot(dir1[i], dir2[i])) * 180.0 / 3.14159265;
		printf("An %f\n", an);
		ang[i] = an;
	}


	// now check the angle
	double angle = acos(dot(dir1[0], dir2[0])) * 180.0 / 3.14159265;
	printf("Angle is %f\n", angle);

	//if ((angle > 30.0) || (angle < 0.0))
	//	return false;

	// gsl solve linear system to find rotation
	double a_data[9];
	double b_data[3];
	double x[3];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			a_data[j * 3 + 0] = dir1[j].x;
			a_data[j * 3 + 1] = dir1[j].y;
			a_data[j * 3 + 2] = dir1[j].z;
		}

		if (i == 0)
		{
			b_data[0] = dir2[0].x;
			b_data[1] = dir2[1].x;
			b_data[2] = dir2[2].x;
		}
		if (i == 1)
		{
			b_data[0] = dir2[0].y;
			b_data[1] = dir2[1].y;
			b_data[2] = dir2[2].y;
		}
		if (i == 2)
		{
			b_data[0] = dir2[0].z;
			b_data[1] = dir2[1].z;
			b_data[2] = dir2[2].z;
		}

		SolveLinearSystem(3, a_data, b_data, x);

		T[0][i] = x[0];
		T[1][i] = x[1];
		T[2][i] = x[2];
	}

	T[0][3] = 0.0;
	T[1][3] = 0.0;
	T[2][3] = 0.0;

	return true;
}

MeshRegistration::MeshRegistration(char* file1, char* file2)
{
	m2m1 = NULL;

	meshfile1 = string(file1);
	meshfile2 = string(file2);

	vtkObject::SetGlobalWarningDisplay(0);

	// read the mesh file
	vtkMeshReader1 = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader1->SetFileName(file1);
    vtkMeshReader1->Update();
	vtkMesh1 = vtkMeshReader1->GetOutput();

	// read the mesh file
	vtkMeshReader2 = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader2->SetFileName(file2);
    vtkMeshReader2->Update();
	vtkMesh2 = vtkMeshReader2->GetOutput();

	// find curvatures for meshes
	vtkCurvatureTensor* vtkcur = new vtkCurvatureTensor();
	vtkcur->GetTensor(vtkMesh1);
	vtkcur->GetTensor(vtkMesh2);

	// add vertices to the sorted list
	for (int i = 0; i < vtkMesh1->GetNumberOfPoints(); i++)
	{
		int2 val = make_int2(i, 1);
		idsSorted.push_back(val);
	}
	for (int i = 0; i < vtkMesh2->GetNumberOfPoints(); i++)
	{
		int2 val = make_int2(i, 2);
		idsSorted.push_back(val);
	}
	CurvatureSorter csorter(vtkMesh1, vtkMesh2);
	std::sort(idsSorted.begin(), idsSorted.end(), csorter);
	idsSortedIndex = 0;
	for (int i = 0; i <idsSorted.size(); i++)
	{
		int2 val = idsSorted[i];
		if (val.y == 2) continue;
		meshid2curvid[val.x] = i;
	}

	// initialize min spanning tree
	mstreeSize = 0;
	maxlevel = 0;

	// create the point locator for the tree
	vtkPoints* mstreePoints = vtkPoints::New();
	mstreeLocator = vtkIncrementalOctreePointLocator::New();
	mstreeLocator->FreeSearchStructure();
	mstreeLocator->InitPointInsertion(mstreePoints, vtkMesh1->GetBounds(), vtkMesh1->GetNumberOfPoints());

	// cadidate search locator
	mesh2Locator = (vtkIncrementalOctreePointLocator**) malloc(THREADNUMBER * sizeof(vtkIncrementalOctreePointLocator*));
	for (int i = 0; i < THREADNUMBER; i++)
	{
		vtkPoints* mesh2LocatorPoints = vtkPoints::New();
		mesh2Locator[i] = vtkIncrementalOctreePointLocator::New();
		mesh2Locator[i]->FreeSearchStructure();
		mesh2Locator[i]->InitPointInsertion(mesh2LocatorPoints, vtkMesh2->GetBounds(), vtkMesh2->GetNumberOfPoints());

		for (int j = 0; j < vtkMesh2->GetNumberOfPoints(); j++)
		{
			//double tup;
			//vtkMesh2->GetPointData()->GetArray("MTCHMSK")->GetTuple(j, &tup);
			//if (tup != 1)
			//	continue;

			double point[3];
			vtkMesh2->GetPoint(j, point);
			mesh2Locator[i]->InsertPoint(j, point);
		}
	}
	printf("Locator for mesh 2 has %d points\n", mesh2Locator[0]->GetNumberOfPoints());
	omp_set_num_threads(THREADNUMBER);

	// cadidate search locator
	vtkPoints* mesh1LocatorPoints = vtkPoints::New();
	mesh1Locator = vtkIncrementalOctreePointLocator::New();
	mesh1Locator->FreeSearchStructure();
	mesh1Locator->InitPointInsertion(mesh1LocatorPoints, vtkMesh1->GetBounds(), vtkMesh1->GetNumberOfPoints());
	for (int j = 0; j < vtkMesh1->GetNumberOfPoints(); j++)
	{
		//double tup;
		//vtkMesh1->GetPointData()->GetArray("MTCHMSK")->GetTuple(j, &tup);
		//if (tup != 1)
		//	continue;

		double point[3];
		vtkMesh1->GetPoint(j, point);
		mesh1Locator->InsertPoint(j, point);
	}
	printf("Locator for mesh 1 has %d points\n", mesh1Locator->GetNumberOfPoints());


	// all vertices are new
	newvertex = (bool*) malloc(vtkMesh1->GetNumberOfPoints() * sizeof(bool));
	for (int i = 0; i < vtkMesh1->GetNumberOfPoints(); i++)
		newvertex[i] = true;

	// candidates
	candidates = new vector<Candidate>[vtkMesh1->GetNumberOfPoints()];

	// registration status
	add_count = 0.05 * vtkMesh1->GetNumberOfPoints();
	done_reg = false;
	updown = 0;
}

MeshRegistration::~MeshRegistration(void)
{
}

void MeshRegistration::PrintCandidates(int id, int count)
{
	#pragma omp critical 
	{
		fflush(stdout);
		printf("Candidates for %d at level %d\n", id, id2node[id]->level);
		for (int r = 0; ((r < candidates[id].size()) && (r < count)); r++)
		{
			Candidate c = candidates[id][r];

			printf("%d -> b = %f, phi = %f, conf = %f\n", c.id, c.belief, c.phi, c.confidence);
			
			double point[3];
			vtkMesh1->GetPoint(id, point);
			printf("%lf %lf %lf -> %f %f %f\n", point[0], point[1], point[2], c.point.x, c.point.y, c.point.z);
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					printf("%lf\t", c.Tl[i][j]);
				}
				printf("\n");
			}
		}
		fflush(stdout);
	}
}

void MeshRegistration::BuildTree()
{
	vtkPolyData* vtkMesh = this->vtkMesh1;

	// bounds
	double bounds[6];
	vtkMesh1->GetBounds(bounds);
	diag = make_float3(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]);
	printf("Mesh1 Bounds %f %f %f %f %f %f\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

	double bounds2[6];
	vtkMesh2->GetBounds(bounds2);
	printf("Mesh2 Bounds %f %f %f %f %f %f\n", bounds2[0], bounds2[1], bounds2[2], bounds2[3], bounds2[4], bounds2[5]);

	// get data array 
	vtkPointData* pointsData = vtkMesh->GetPointData();	
	
	// find all cells boundary edges
	vtkCellArray* vtkCells = vtkMesh->GetPolys();
	vtkCells->InitTraversal();
	vtkIdType npts;
	vtkIdType* pts;
	int count = 0;
	map<Edge, int> edges;
	while (vtkCells->GetNextCell(npts, pts) != 0)
	{
		Edge edge1(min(pts[0], pts[1]), max(pts[0], pts[1]));
		Edge edge2(min(pts[0], pts[2]), max(pts[0], pts[2]));
		Edge edge3(min(pts[1], pts[2]), max(pts[1], pts[2]));
		
		map<Edge, int>::iterator it;
		
		if ((it = edges.find(edge1)) != edges.end())
		{
			int scr = it->second + 1;
			edges.erase(it);
			edges.insert(pair<Edge, int>(edge1, scr));
		}
		else
		{
			edge1.length = GetEdgeLength(vtkMesh, edge1);
			edges.insert(pair<Edge, int>(edge1, 1));
		}
			
		if ((it = edges.find(edge2)) != edges.end())
		{	
			int scr = it->second + 1;
			edges.erase(it);
			edges.insert(pair<Edge, int>(edge2, scr));
		}
		else
		{
			edge2.length = GetEdgeLength(vtkMesh, edge2);
			edges.insert(pair<Edge, int>(edge2, 1));
		}
			
		if ((it = edges.find(edge3)) != edges.end())
		{	
			int scr = it->second + 1;
			edges.erase(it);
			edges.insert(pair<Edge, int>(edge3, scr));
		}
		else
		{
			edge3.length = GetEdgeLength(vtkMesh, edge3);
			edges.insert(pair<Edge, int>(edge3, 1));
		}
	}
	
	// collect boundary points
	int aveedgelength_count = 0;
	aveedgelength = 0.0;
	bpts = new set<vtkIdType>();
	vector<Edge>* alledges = new vector<Edge>();
	for (map<Edge, int>::iterator it = edges.begin(); it != edges.end(); it++ )
	{
		aveedgelength += it->first.length;
		aveedgelength_count++;

		if (it->second != 2)
		{
			count++;
			bpts->insert((*it).first.x);
			bpts->insert((*it).first.y);
		}

		alledges->push_back((*it).first);
	}
	aveedgelength /= aveedgelength_count;

	// sort all edges based on the length
	EdgeLengthSorter csorter;
	std::sort(alledges->begin(), alledges->end(), csorter);

	// initialize distance from boundary
	double* vals = (double*) malloc(vtkMesh->GetNumberOfPoints() * sizeof(double));
	for (int i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
	{
		vals[i] = numeric_limits<double>::max();
	}
	for (set<vtkIdType>::iterator r = bpts->begin(); r != bpts->end(); r++)
	{
		vals[*r] = 0.0;
	}

	// loop on the points minimizing the distance from boundary
	bool changed = true;
	while (changed)
	{
		changed = false;
		for (int i = 0; i < alledges->size(); i++)
		{
			Edge edge = alledges->at(i);

			if (vals[edge.y] + edge.length < vals[edge.x])
			{
				vals[edge.x] = vals[edge.y] + edge.length;
				changed = true;
			}
				
			if (vals[edge.x] + edge.length < vals[edge.y])
			{
				vals[edge.y] = vals[edge.x] + edge.length;
				changed = true;
			}
		}
	}

	// create the neighbor list
	neighborlist.resize(vtkMesh->GetNumberOfPoints());
	for (int i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
	{
		vector<int> ps;
		neighborlist.push_back(ps);
	}
	for (int i = 0; i < alledges->size(); i++)
	{
		Edge edge = alledges->at(i);

		neighborlist[edge.x].push_back(edge.y);
		neighborlist[edge.y].push_back(edge.x);
	}

	// sort all vertice based on distance from the boundary
	vector<int> bsorted;
	for (int i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
	{
		bsorted.push_back(i);
	}
	BoundaryDistSorter bdsorter(vals);
	std::sort(bsorted.begin(), bsorted.end(), bdsorter);

	// initialize for tree building
	char* marks = (char*) malloc(vtkMesh->GetNumberOfPoints() * sizeof(char));
	memset(marks, 0, vtkMesh->GetNumberOfPoints() * sizeof(char));
	rootnode = new TNode(NULL, -1);
	levelnodes[0].push_back(rootnode);

	// build the tree
	int nodescount = 0;
	int cc = 0;
	int discarded = 0;
	int notdiscarded = 0;
	double bdislimit = 3 * length(diag) * 0.007;
	while (nodescount < vtkMesh->GetNumberOfPoints()) 
	{
		map<int, vector<TNode*> > tmplevelnodes;

		// pop one vertex
		int id = bsorted.back();
		bsorted.pop_back();
		if (marks[id] == 1)
			continue;

		// create node for root
		vector<int> toprocess;
		toprocess.push_back(id);
		marks[id] = 1;
		nodescount++;

		// create node
		TNode* newnode = new TNode(rootnode, id);
		id2node[id] = newnode;
		maxlevel = max(maxlevel, newnode->level);
		tmplevelnodes[newnode->level].push_back(newnode);
		id2cc[id] = cc;

		int csize = 1;
		while (toprocess.empty() == false)
		{
			// get vertex from list to process
			int parentid = toprocess.front();
			toprocess.erase(toprocess.begin());

			// loop on closest
			for (int s = 0; s < neighborlist[parentid].size(); s++)
			{
				int childid = neighborlist[parentid][s];
				if (marks[childid] == 1)
					continue;

				// condition on growth
				if ((vals[parentid] < bdislimit) && (vals[childid] > vals[parentid]))
					continue;

				// update size
				csize ++;

				// add node
				toprocess.push_back(childid);
				marks[childid] = 1;
				nodescount++;

				// create the node 
				TNode* parentn = id2node[parentid];
				TNode* newnode = new TNode(parentn, childid);
				parentn->children.push_back(newnode);
				id2node[childid] = newnode;
				maxlevel = max(maxlevel, newnode->level);
				tmplevelnodes[newnode->level].push_back(newnode);
				id2cc[childid] = cc;
			}
		}

		//if ((csize < 5) && (csize < 0.0001 * vtkMesh->GetNumberOfPoints()))
		//{
		//	discarded += csize;
		//}
		//else
		{
			notdiscarded += csize;

			rootnode->children.push_back(newnode);

			for (int l = 0; l <= maxlevel; l++)
			{
				levelnodes[l].insert(levelnodes[l].end(), tmplevelnodes[l].begin(), tmplevelnodes[l].end());
			}
			
			//printf("Tree size is %d\n", csize);
		}

		cc++;
	}

	printf("Number of tree is %d\n", rootnode->children.size());
	printf("Total discarded is %d and included is %d\n", discarded, notdiscarded);

	// free memory
	free(vals);
	free(marks);
}

double MeshRegistration::SquarePointDist(int mesh, int id1, int id2)
{
	if (mesh == 1)
	{
		double pt1[3];
		vtkMesh1->GetPoint(id1, pt1);
		double pt2[3];
		vtkMesh1->GetPoint(id2, pt2);
		float3 diff = make_float3(pt1[0] - pt2[0], pt1[1] - pt2[1], pt1[2] - pt2[2]);
		return dot(diff, diff);
	}
	else if (mesh == 2)
	{
		double pt1[3];
		vtkMesh2->GetPoint(id1, pt1);
		double pt2[3];
		vtkMesh2->GetPoint(id2, pt2);
		float3 diff = make_float3(pt1[0] - pt2[0], pt1[1] - pt2[1], pt1[2] - pt2[2]);
		return dot(diff, diff);
	}
	else
	{
		printf("ERROR: Unknown mesh\n");
		return 0;
	}
}

void MeshRegistration::FindClosestOnSurfacePoint(int& id, float3& pt)
{
	// find the closest id to point 
	int thid = omp_get_thread_num();
	id = mesh2Locator[thid]->FindClosestPoint(pt.x, pt.y, pt.z);

	// search for closest location
	float3 ret;
	double mindist = numeric_limits<double>::max();
	vtkIdList* cells = vtkIdList::New();
	vtkMesh2->GetPointCells(id, cells);
	vtkGenericCell* cell = vtkGenericCell::New();

	for (int i = 0; i < cells->GetNumberOfIds(); i++)
	{
		vtkMesh2->GetCell(cells->GetId(i), cell);
		if (cell->GetCellType() != VTK_TRIANGLE) 
			continue;
		
		double point[3];
		vtkMesh2->GetPoint(cell->GetPointId(0), point);
		float3 V0 = make_float3(point[0], point[1], point[2]);
		
		vtkMesh2->GetPoint(cell->GetPointId(1), point);
		float3 V1 = make_float3(point[0], point[1], point[2]);
		
		vtkMesh2->GetPoint(cell->GetPointId(2), point);
		float3 V2 = make_float3(point[0], point[1], point[2]);
		
		float3 c = ClosestPointOnTriangle(pt, V0, V1, V2);

		// check for errors
		if (myisnan(c.x) || myisnan(c.y) || myisnan(c.z))
		{
			printf("\nError %s %d\n", __FILE__, __LINE__);
			char chr = getchar();
		}
		
		if (length(pt - c) < mindist)
		{
			mindist = length(pt - c);
			ret = c;
		}
	}
	
	// update id again
	//id = mesh2Locator[thid]->FindClosestPoint(ret.x, ret.y, ret.z);
	pt = ret;

	cell->Delete();
	cells->Delete();
}

void MeshRegistration::CheckCandidateCase(int Xi, Candidate& ca, double a, double b, double c, double d, float3 va, float3 vb, float3 vc, float3 vd, vector<Candidate>* ret, double curvlimit)
{
	CUERR("CheckCandidateCase");

	double c1 = max(abs(a), abs(b));
	double c2 = max(abs(c), abs(d));
	c1 = abs(c1 - c2);

	ca.phi = 1.0;//exp(-c1);
	ret->push_back(ca);
}

void MeshRegistration::CheckCandidate(TNode* child, Candidate& ca, vector<Candidate>* ret, double curvlimit)
{
	CUERR("CheckCandidate");

	double curXi[2];
	vtkMesh1->GetPointData()->GetArray("Principal Curvatures")->GetTuple(child->id, curXi);
	double vec1_Xi[3];
	vtkMesh1->GetPointData()->GetArray("Maximum Curvature Vector")->GetTuple(child->id, vec1_Xi);
	double vec2_Xi[3];
	vtkMesh1->GetPointData()->GetArray("Minimum Curvature Vector")->GetTuple(child->id, vec2_Xi);

	double curxi[2];
	vtkMesh2->GetPointData()->GetArray("Principal Curvatures")->GetTuple(ca.id, curxi);
	double vec1_xi[3];
	vtkMesh2->GetPointData()->GetArray("Maximum Curvature Vector")->GetTuple(ca.id, vec1_xi);
	double vec2_xi[3];
	vtkMesh2->GetPointData()->GetArray("Minimum Curvature Vector")->GetTuple(ca.id, vec2_xi);

	float3 v1_Xi = make_float3(vec1_Xi[0], vec1_Xi[1], vec1_Xi[2]);
	float3 v2_Xi = make_float3(vec2_Xi[0], vec2_Xi[1], vec2_Xi[2]);
	float3 v1_xi = make_float3(vec1_xi[0], vec1_xi[1], vec1_xi[2]);
	float3 v2_xi = make_float3(vec2_xi[0], vec2_xi[1], vec2_xi[2]);

	CheckCandidateCase(child->id, ca, curXi[0], curXi[1], curxi[0], curxi[1], v1_Xi, v2_Xi, v1_xi, v2_xi, ret, curvlimit);
}

set<int>* MeshRegistration::FindNNeighborsInMesh2(int n, int id)
{
	set<int>* neighbors = new set<int>();

	vector<int> q;
	q.push_back(id);
	vtkIdList* cellIds = vtkIdList::New();
	vtkIdList* ptIds = vtkIdList::New();
	while ((neighbors->size() < n) && (!q.empty()))
	{
		int pid = q.front();
		q.erase(q.begin());

		// get cells of pid
		vtkMesh2->GetPointCells(pid, cellIds);
		for (int i = 0; i < cellIds->GetNumberOfIds(); i++)
		{
			vtkMesh2->GetCellPoints(cellIds->GetId(i), ptIds);
			for (int j = 0; j < ptIds->GetNumberOfIds(); j++)
			{
				int ptId = ptIds->GetId(j);
				if (ptId == pid) continue;
				if (ptId == id) continue;

				if (neighbors->find(ptId) == neighbors->end())
				{
					neighbors->insert(ptId);
					q.push_back(ptId);
				}
			}
		}
	}

	cellIds->Delete();
	ptIds->Delete();

	return neighbors;
}

void MeshRegistration::ComputeLocalTransform(int Xi, Candidate& ca)
{
	float3 ptXi = GetFloat3Point(vtkMesh1, Xi);
	//float3 nrXi = GetFloat3Normal(vtkMesh1, Xi);
	ptXi = ApplyTranformation(ptXi, ca.Tg);

	float3 ptXj = ca.point;
	//float3 nrXj = GetFloat3Normal(vtkMesh2, ca.id);

	//http://www.chriskugler.com/xna/xna-get-rotation-matrix-between-two-vectors/

	/*float dotp = dot(nrXi, nrXj);
	if (dotp < 0.0)
	{
		nrXj = -nrXj;
		dotp = -dotp;
	}*/

	//if (!myisnan(dotp))
	//{
		//float angle = acos(dotp);
		//if (!myisnan(angle))
		//{
			//float3 crossp = normalize(cross(nrXi, nrXj));

			//if (updown < 8)
			//angle = 0.0;
			//matrixFromAxisAngle(crossp, angle, ca.Tl);
			//printf("%f,",angle);
			//ca.rotationangle = angle;
			//ca.rotationaxis = crossp;

			ca.Tl[3][0] = ptXj.x - ptXi.x;
			ca.Tl[3][1] = ptXj.y - ptXi.y;
			ca.Tl[3][2] = ptXj.z - ptXi.z;
			ca.Tl[3][3] = 1.0;
		//}
		//else CUERRP("angle is not a number")
	//}
	//else 
	//	CUERRP("dot product is not a number")

}

vector<Candidate>* MeshRegistration::FindVertexCandidates(int id, double curvlimit, bool goingup)
{
	CUERR("FindVertexCandidates");
	// For each match in parent find the closest match at the gradparent. Get
	// transformation from parent, and scaling from grandparent and parent.
	// Find the expected center and set the radius to be twice that of the 
	// distance to the match of parent. Search the sorted list for the matching
	// items.

	int thid = omp_get_thread_num();

	TNode* child = id2node[id];
	TNode* parent = child->parent;
	vector<Candidate>* ret = new vector<Candidate>();
	
	// point to transform
	float3 ipt = GetFloat3Point(vtkMesh1, child->id);

	// list to carry
	vector<Candidate> pc;
	map<int, float3> iden;

	// my own candidates neighborhood
	if (goingup == true)
	{
		// children candidates neighborhood
		for (vector<TNode*>::iterator it = child->children.begin(); it != child->children.end(); it++)
		{
			for (int r = 0; r < candidates[(*it)->id].size(); r++)
			{
				Candidate c = candidates[(*it)->id][r];

				// find center (at the old candidate)
				float3 centerf = ApplyTranformation(ipt, c.Tg, c.Tl);

				// find the closest point
				int cid = mesh2Locator[thid]->FindClosestPoint(centerf.x, centerf.y, centerf.z);
				Candidate ca(c);
				ca.id = cid;
				ca.belief = 0.0;
				ca.phi = 1.0;
				ca.point = centerf;
				pc.push_back(ca);
			}
		}
	}
	else
	{
		// parent candidates neighborhood
		if ((parent != NULL) && (parent->id != -1))
		{
			for (int r = 0; r < candidates[parent->id].size(); r++)
			{
				Candidate c = candidates[parent->id][r];

				// find center (at the old candidate)
				float3 centerf = ApplyTranformation(ipt, c.Tg, c.Tl);

				// find the closest point
				int cid = mesh2Locator[thid]->FindClosestPoint(centerf.x, centerf.y, centerf.z);
				Candidate ca(c);
				ca.id = cid;
				ca.belief = 0.0;
				ca.phi = 1.0;
				ca.point = centerf;
				pc.push_back(ca);
			}
		}
	}

	// add many surrounding candidates
	TNode* n = id2node[id];

	// if leaf and going up use previous candidates
	if (((goingup == true) && (n->children.size() == 0)) || (updown >= 2))
	{
		for (int i = 0; i < candidates[id].size(); i++)
		{
			pc.push_back(candidates[id][i]);
		}
	}

	// if more than two steps bring surrounding region
	if ((updown >= 2) && (candidates[id].size() > 0))
	{
		// just add a bunch of surrounding vertices
		vtkIdList* result = vtkIdList::New();
		double point[3];
		point[0] = candidates[id][0].point.x;
		point[1] = candidates[id][0].point.y;
		point[2] = candidates[id][0].point.z;
		mesh2Locator[thid]->FindClosestNPoints(CANDIDATECOUNTMAX*2, point, result);
		for (int i = 0; i < result->GetNumberOfIds(); i++)
		{
			int cid = result->GetId(i);
			float3 centerf = GetFloat3Point(vtkMesh2, cid);

			// find the closest point
			Candidate ca(candidates[id][0]);
			ca.id = cid;
			ca.belief = 0.0;
			ca.phi = 1.0;
			ca.point = centerf;
			ComputeLocalTransform(id, ca);
			pc.push_back(ca);
		}
		result->Delete();
	}

	// add a single point near the original position
	if (pc.size() == 0)
	{
		// shift ipt by the bounding box info
		//double bounds1[6];
		//vtkMesh1->GetBounds(bounds1);
		//double bounds2[6];
		//vtkMesh2->GetBounds(bounds2);
		float3 iptt = ipt;
		//iptt.x = (iptt.x - bounds1[0]) / (bounds1[1] - bounds1[0]);
		//iptt.x = bounds2[0] + iptt.x * (bounds2[1] - bounds2[0]);
		//iptt.y = (iptt.y - bounds1[2]) / (bounds1[3] - bounds1[2]);
		//iptt.y = bounds2[2] + iptt.y * (bounds2[3] - bounds2[2]);
		//iptt.z = (iptt.z - bounds1[4]) / (bounds1[5] - bounds1[4]);
		//iptt.z = bounds2[4] + iptt.z * (bounds2[5] - bounds2[4]);

		// do the search
		int cid = mesh2Locator[thid]->FindClosestPoint(iptt.x, iptt.y, iptt.z);
		float3 centerf = GetFloat3Point(vtkMesh2, cid);

		// find the closest point
		Candidate ca(cid, 0.0, 1.0, centerf);
		ComputeLocalTransform(id, ca);
		pc.push_back(ca);
	}

	// search for the points
	for (int i = 0; i < pc.size(); i++)
	{
		map<int,float3>::iterator it = iden.find(pc[i].id);
		if (it != iden.end()) 
		{
			if (length(pc[i].point - it->second) < 0.01)
				continue;
		}
		else
		{
			iden[pc[i].id] = pc[i].point;
		}
		CheckCandidate(child, pc[i], ret, curvlimit);
	}

	return ret;
}

double MeshRegistration::Psi(const int& Xi, Candidate& xi, const int& Xj, Candidate& xj)
{
	CUERR("Psi");

	double ret = 0.0;


	double point[3];

	if (updown >= 2)
	{
		// get point Xj
		float3 Xjp = GetFloat3Point(vtkMesh1, Xj);

		// get point xi
		float3 ti_xj = ApplyTranformation(Xjp, xi.Tg, xi.Tl);

		ret = -length(ti_xj - xj.point) / length(Xjp - xj.point);
		ret = exp(ret);
	}
	else
	{
		// difference in local transformation
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				ret += abs(xi.Tl[i][j] - xj.Tl[i][j]);
			}
		}

		ret = exp(-ret);
	}

	return ret;
}

bool sortmappings (int3 i,int3 j) { return (i.z > j.z); }

void MeshRegistration::ComputeTransform(int Xi, Candidate& c, vector<int3>* mappings, bool goingup)
{
	CUERR("ComputeTransform");
	
	float3 old = c.point;

	// adjust location
	FindClosestOnSurfacePoint(c.id, c.point);

	if (updown >= 2)
	{
		if ((updown % 2 == 0)/* || (length(c.point - old) < aveedgelength)*/)
		{
			// compute new transfomration
			ComputeLocalTransform(Xi, c);
		}
	}
}

void MeshRegistration::AdjustNode(TNode* n, bool goingup, bool lastiter)
{
	CUERR("AdjustNode");

	int Xi = n->id;

	// update the visiting count
	#pragma omp critical 
	{
		visited++;
	}

	// find the candidate list for the point
	vector<Candidate>* cl = NULL;
	double curvlimit = CURVEPS;
	cl = FindVertexCandidates(n->id, curvlimit, goingup);

	// delete old candidates
	int oldbest = (candidates[Xi].size() > 0) ? candidates[Xi][0].id : -1;
	candidates[Xi].clear();

	// for each candidate evaluate its belief
	if (cl->size() == 0)
	{
		printf("Point %d has no candidates!!! \n", Xi);
	}

	// check every candidate
	for (int k = 0; k < cl->size(); k++)
	{
		Candidate ca = (*cl)[k];

		// candidate id
		int xi = ca.id;

		// loop on neighbors of Xi 
		bool skip = false;
		double Pi_j = 1.0;
		vector<int3> mappings;

		// compute the transofrm here for psi to be correct
		ComputeTransform(Xi, ca, &mappings, goingup);

		////////// Children as neighbors
		if ((goingup) || (lastiter))
		{
			for (vector<TNode*>::iterator childit = n->children.begin(); childit != n->children.end(); childit++)
			{
				int Xj = (*childit)->id;

				// loop on candidates of Xj 
				int3 best;
				double sum_xj = 0.0;
				for (int r = 0; r < candidates[Xj].size(); r++)
				{
					double psi = Psi(Xi, ca, Xj, candidates[Xj][r]);
					double s_xj = psi * candidates[Xj][r].belief;
					if (s_xj > sum_xj)
					{
						sum_xj = s_xj;
						best.x = Xj;
						best.y = r;
						best.z = psi * 100000;
					}
				}

				//mappings.push_back(best);
				Pi_j *= sum_xj;
			}
		}

		////////// parent as neighbor
		if ((!goingup) || (lastiter))
		{
			int Xj = n->parent->id;

			// loop on candidates of Xj 
			int3 best;
			double sum_xj = 0.0;
			for (int r = 0; r < candidates[Xj].size(); r++)
			{
				double psi = Psi(Xi, ca, Xj, candidates[Xj][r]);
				double s_xj = psi * candidates[Xj][r].belief;
				if (s_xj > sum_xj)
				{
					sum_xj = s_xj;
					best.x = Xj;
					best.y = r;
					best.z = psi * 100000;
				}
			}

			//mappings.push_back(best);
			Pi_j *= sum_xj;
		}

		// check for not important candidates
		if (skip) continue;
		if (Pi_j == 0.0) continue;

		// finally the belief and insert in the list of candidates
		Candidate newcd(ca);
		newcd.belief = newcd.phi * Pi_j;
		bool fi = true;
		int r = 0;
		while (((r < CANDIDATECOUNTMIN) || (newcd.belief > 0.99 * candidates[Xi][0].belief)) && (r < CANDIDATECOUNTMAX))
		{
			if (candidates[Xi].size() == r)
			{
				candidates[Xi].push_back(newcd);
				break;
			}
			if (candidates[Xi][r].belief < newcd.belief)
			{
				Candidate tmp(candidates[Xi][r]);
				candidates[Xi][r].Copy(newcd);
				newcd.Copy(tmp);			
			}

			r ++;
		}
	}

	// normalize beliefs
	double den = 0.0;
	for (int r = 0; r < candidates[Xi].size(); r++)
	{
		den += candidates[Xi][r].belief;
		candidates[Xi][r].confidence = ParentAgreement(Xi, r);
	}
	if (den == 0.0)
	{
		// check the cl
		if (cl->size() == 0)
		{
			printf("At node %d the candidate list was empty\n", Xi);
			char chr = getchar();
		}

		// add all items in cl
		vector<int3> mappings;
		for (int k = 0; k < cl->size(); k++)
		{
			Candidate ca = (*cl)[k];
			Candidate newcd(ca);
			newcd.belief = newcd.phi * 1.0;
			ComputeTransform(Xi, newcd, &mappings, goingup);
			candidates[Xi].push_back(newcd);
		}

		// give equal priorities
		den = 1.0 / candidates[Xi].size();
		for (int r = 0; r < candidates[Xi].size(); r++)
			candidates[Xi][r].belief = den;
	}
	for (int r = 0; r < candidates[Xi].size(); r++)
	{
		candidates[Xi][r].belief /= den;
	}

	delete cl;


	// record change if happened
	int newbest = (candidates[Xi].size() > 0) ? candidates[Xi][0].id : -1;
	if ((oldbest != newbest) && (oldbest != -1))
	{
		#pragma omp critical 
		{
			meshshift++;
		}
	}
	if (newbest == Xi)
	{
		#pragma omp critical 
		{
			correct++;
		}
	}

	// check the new distance
	if (newbest != -1)
	{
		/*double point[3];
		vtkMesh1->GetPoint(Xi, point);
		float3 pt1 = make_float3(point[0], point[1], point[2]);
		vtkMesh2->GetPoint(newbest, point);
		float3 pt2 = make_float3(point[0], point[1], point[2]);

		if (length(pt1 - pt2) > 100.0)
		{
			TNode* cur = id2node[Xi];
			while (cur->level > 0)
			{
				PrintCandidates(cur->id, 5);
				cur = cur->parent;
			}
		}*/
	}
	
}

void MeshRegistration::PropagateUp()
{
	meshshift = 0;
	correct = 0;
	visited = 0;

	for (int level = maxlevel; level >= 1; level--)
	{
		if ((level < 25) || (level % 100 == 0)) printf("%d ",level); fflush(stdout);
		levelmeshshift[level] = 0;
		levelavegrdist[level] = 0.0;

		#pragma omp parallel for schedule(dynamic,32)
		for (int i = 0; i < levelnodes[level].size(); i++)
		{
			AdjustNode(levelnodes[level][i], true, false);
		}
	}
}

void MeshRegistration::PropagateDown()
{
	meshshift = 0;
	correct = 0;
	visited = 0;

	for (int level = 2; level <= maxlevel; level++)
	{
		if ((level < 25) || (level % 100 == 0)) printf("%d ",level); fflush(stdout);
		levelmeshshift[level] = 0;
		levelavegrdist[level] = 0.0;

		#pragma omp parallel for schedule(dynamic,32)
		for (int i = 0; i < levelnodes[level].size(); i++)
		{
			AdjustNode(levelnodes[level][i], false, false);
		}
	}
}

double MeshRegistration::ComputeTransformationVariance()
{
	double Pr1 = 0.0;
	double Pr2 = 0.0;
	double Pr3 = 0.0;

	//ofstream myfile;
	//myfile.open ("Z:\\nabla\\Samer\\Datasets\\Multifield\\comb_mf\\reslog2.txt");
  
	for (int Xi = 0; Xi < vtkMesh1->GetNumberOfPoints(); Xi++)
	{
		if (candidates[Xi].size() == 0)
			continue;

		TNode* n = id2node[Xi];

		double Pi = 1.0;

		// children
		for (vector<TNode*>::iterator childit = n->children.begin(); childit != n->children.end(); childit++)
		{
			int Xj = (*childit)->id;

			Pi *= Psi(Xi, candidates[Xi][0], Xj, candidates[Xj][0]);
		}

		// parent
		if ((n->parent != NULL) && (n->parent != rootnode))
		{
			int Xj = n->parent->id;

			Pi *= Psi(Xi, candidates[Xi][0], Xj, candidates[Xj][0]);
		}

		//if (n->level >= 118)
		//myfile << "level:" << n->level << " id:" << Xi << " " << Pi << " " << candidates[Xi][0].id << " " << candidates[Xi][0].belief << 
		//	" " << candidates[Xi][0].point.x << " " << candidates[Xi][0].point.y << " " << candidates[Xi][0].point.z << 
		//	" " << candidates[Xi][0].T[3][0] << " " << candidates[Xi][0].T[3][2] << " " << candidates[Xi][0].T[3][2] << "\n";

		//myfile << "level:" << n->level << " ";
		//for (vector<TNode*>::iterator childit = n->children.begin(); childit != n->children.end(); childit++)
		//	myfile << (*childit)->id << " ";
		//myfile << "\n";
		
		Pr1 += Pi;
		Pr2 += candidates[Xi][0].belief;

		// compute distance
		double point[3];
		vtkMesh1->GetPoint(Xi, point);
		float3 pt1 = make_float3(point[0], point[1], point[2]);
		vtkMesh2->GetPoint(candidates[Xi][0].id, point);
		float3 pt2 = make_float3(point[0], point[1], point[2]);
		Pr3 += length(pt1 - pt2);
	}

	Pr3 /= vtkMesh1->GetNumberOfPoints();

	if (myisnan(Pr1) || myisinf(Pr1))
		printf("Error here at %d\n", __LINE__);

	//myfile.close();
	printf("Scores are %lf %lf %lf \n\n", Pr1, Pr2, Pr3);

	return Pr3;
}

bool MeshRegistration::AdjustGlobalTransform()
{
	float3 translation;
	float3 axis;
	double angle;

	bool ret = false;

	// for each root 
	for (int root = 0; root < rootnode->children.size(); root++)
	{
		// initialize variables
		translation = make_float3(0.0);
		axis = make_float3(0.0);
		angle = 0.0;
		int rootid = rootnode->children[root]->id;

		// collect vertices in the tree
		set<int> allids;
		vector<TNode*> q;
		q.push_back(rootnode->children[root]);
		while (!q.empty())
		{
			TNode* n = q.back();
			q.pop_back();
			allids.insert(n->id);

			// add children
			for (int i = 0; i < n->children.size(); i++)
			{
				if (allids.find(n->children[i]->id) == allids.end())
				{
					q.push_back(n->children[i]);
				}
			}
		}
		if (allids.size() < 0.01 * vtkMesh1->GetNumberOfPoints())
			continue;
		
		////// MATRIX OLD
		// copy old global matrix
		float T1[4][4];
		for (int r = 0; r < 4; r++)
		{
			for (int s = 0; s < 4; s++)
			{
				T1[r][s] = candidates[rootid][0].Tg[r][s];
			}
		}




		////// MATRIX OLD ROTATION OLNY
		float onlyRotaion[4][4];
		memcpy(onlyRotaion, T1, 16 * sizeof(float));
		onlyRotaion[3][0] = onlyRotaion[3][1] = onlyRotaion[3][2] = 0.0;




		////// MATRIX TRANSLATION
		// The translation matrix computation
		// compute the translation matrix
		for (set<int>::iterator it = allids.begin(); it != allids.end(); it++)
		{
			int id = *it;

			// apply old rotation than rediscover trasnlation
			float3 pt = GetFloat3Point(vtkMesh1, id);
			pt = ApplyTranformation(pt, onlyRotaion);
			translation += candidates[id][0].point - pt;
		}
		translation /= allids.size();
		float T2[4][4];
		for (int r = 0; r < 4; r ++)
		{
			for (int s = 0; s < 4; s++)
			{
				if (r == s)
					T2[r][s] = 1.0;
				else
					T2[r][s] = 0.0;
			}
		}
		T2[3][0] = translation.x;
		T2[3][1] = translation.y;
		T2[3][2] = translation.z;
		T2[3][3] = 1.0;
		float3 oldtrans = make_float3(T1[3][0], T1[3][1], T1[3][2]);
		float diff = length(translation - oldtrans) / (0.001 + length(oldtrans));





		////// MATRIX ROTATION
		// compute matrix for new rotation
		vector<float3> pts1;
		vector<float3> pts2;
		for (set<int>::iterator it = allids.begin(); it != allids.end(); it++)
		{
			int id = *it;

			float3 pt1 = GetFloat3Point(vtkMesh1, id);
			float3 pt2 = candidates[id][0].point;
			float3 pt3 = ApplyTranformation(pt1, candidates[id][0].Tg, candidates[id][0].Tl);

			if (length(pt2 - pt3) < 3.0 * aveedgelength)
			{
				pts1.push_back(pt1);
				pts2.push_back(pt2);
			}
		}
		float T3[4][4];
		bool rres = false;//FindRotationAnglesPCA(&pts1,&pts2,T3);
		if ((!rres) /*|| (diff > 0.05)*/ || (pts1.size() < 50))
		{
			// if translation is big rotation remain as is
			memcpy(T3, onlyRotaion, 16 * sizeof(float));
		}
		float anglediff = 0.0;
		float sument = 0.0;
		for (int r = 0; r < 3; r++)
		{
			for (int s = 0; s < 3; s++)
			{
				anglediff += abs(T1[r][s] - T3[r][s]);
				sument += abs(T1[r][s]);
			}
		}
		anglediff /= sument;
		T3[3][0] = 0.0;
		T3[3][1] = 0.0;
		T3[3][2] = 0.0;
		T3[3][3] = 1.0;
		
		// compute new T
		float T[4][4];
		float T4[4][4];
		multiplyMatrix(T3, T2, T);

		
		// change in translation is significant
		if ((diff > 0.05) || (anglediff > 0.05))
		{
			
			printf("================ %d\n", allids.size());
			printf("Translation %f, %f %f %f\n", diff, translation.x, translation.y, translation.z);
			printf("Rotation change %f\n\n", anglediff);
		
			for (set<int>::iterator it = allids.begin(); it != allids.end(); it++)
			{
				int id = *it;

				candidates[id].clear();

				// do the search
				float3 ipt = GetFloat3Point(vtkMesh1, id);
				ipt = ApplyTranformation(ipt, T);
				int cid = mesh2Locator[0]->FindClosestPoint(ipt.x, ipt.y, ipt.z);
				float3 centerf = GetFloat3Point(vtkMesh2, cid);

				// find the closest point
				Candidate ca(cid, 1.0, 1.0, centerf);
				memcpy(ca.Tg, T, 16 * sizeof(float));
				ComputeLocalTransform(id, ca);
				candidates[id].push_back(ca);

				/*for (int k = 0; k < candidates[id].size(); k++)
				{
					memcpy(candidates[id][k].Tg, T, 16 * sizeof(float));
					ComputeLocalTransform(id, candidates[id][k]);
				}*/
			}

			ret = true;
		}
	}

	return ret;
}

void MeshRegistration::DoRegistrationStep()
{
	// Execute step
	float old_scr = 0.0;
	float new_scr = 0.0;

	//while ((updown < 4) || ((abs(new_scr - old_scr) / old_scr) > 0.02))
	while(true)
	{
		printf("Step # %d\n", updown);
		
		if (updown % 2 == 0)
		{
			PropagateUp();
			printf("===> Up: %d total visited, %d moved and %d correct\n", visited, meshshift, correct);
			

			ComputeTransformationVariance();
		}
		else
		{
			PropagateDown();
			printf("===> Down: %d total visited, %d moved and %d correct\n", visited, meshshift, correct);


			old_scr = new_scr;
			new_scr = ComputeTransformationVariance();
		}

		//if (updown != 0)
			//AdjustGlobalTransform();

		
		
		if (updown == 1)
			break;
			
		updown ++;
	}
}

double MeshRegistration::ParentAgreement(int Xi, int r)
{
	int Xj = id2node[Xi]->parent->id;
	if (Xj == -1)
		return 1.0;

	// loop on candidates of Xj 
	double psi = 1.0;
	if (candidates[Xj].size() > 0)
		psi = Psi(Xi, candidates[Xi][r], Xj, candidates[Xj][0]);

	return psi;
}

double MeshRegistration::PointDist(int idm1, int idm2)
{
	double p1[3];
	vtkMesh1->GetPoint(idm1, p1);

	double p2[3];
	vtkMesh2->GetPoint(idm2, p2);

	return length(make_float3(p1[0],p1[1],p1[2]) - make_float3(p2[0],p2[1],p2[2]));
}

void MeshRegistration::FinalFix()
{
	int thid = omp_get_thread_num();

	////// STEP 1:
	// This loop role is to find the best match from mesh2 to mesh1
	// for points in mesh2 that has a pointing to vertices from 
	// mesh1. May be the best way to match is to find the vertex in
	// mesh1 that is most confident about its assignment.
	m2m1 = (int*) malloc(vtkMesh2->GetNumberOfPoints() * sizeof(int));
	for (int i = 0; i < vtkMesh2->GetNumberOfPoints(); i++)
	{
		m2m1[i] = -1;
	}
	for (int i = 0; i < vtkMesh1->GetNumberOfPoints(); i++)
	{
		// if too far just ignore it cause it can't be correct
		if ((PointDist(i, candidates[i][0].id) > 0.05 * length(diag)) /*|| (candidates[i][0].confidence < 0.5)*/)
		{
			candidates[i].clear();
			continue;
		}

		// now check if there is a confusion
		if (m2m1[candidates[i][0].id] != -1)
		{
			int idm2 = candidates[i][0].id;
			int idm1_1 = m2m1[idm2];
			int idm1_2 = i;

			double ptm1_1[3];
			vtkMesh1->GetPoint(idm1_1, ptm1_1);
			float3 ptm1_1_f = make_float3(ptm1_1[0], ptm1_1[1], ptm1_1[2]);

			double ptm1_2[3];
			vtkMesh1->GetPoint(idm1_2, ptm1_2);
			float3 ptm1_2_f = make_float3(ptm1_2[0], ptm1_2[1], ptm1_2[2]);

			double ptm2[3];
			vtkMesh2->GetPoint(idm2, ptm2);
			float3 ptm2_f = make_float3(ptm2[0], ptm2[1], ptm2[2]);

			// if not neighbors than it is a confusion
			if (length(ptm1_1_f - ptm1_2_f) > 3.0 * aveedgelength)
			{
				// assign it to the closest
				if (PointDist(idm1_1, idm2) < PointDist(idm1_2, idm2))
					m2m1[idm2] = idm1_1;
				else
					m2m1[idm2] = idm1_2;
			}
		}
		else
			m2m1[candidates[i][0].id] = i;
	}

	////// STEP 2:
	// some vertices in mesh2 will remain without assignments after
	// the last loop because they simply have not vertices in mesh1
	// pointing to them. So they check their neighbors and use the 
	// inverse transform to find their match in mesh1
	for (int i = 0; i < vtkMesh2->GetNumberOfPoints(); i++)
	{
		if (m2m1[i] != -1) continue;

		set<int>* nei = FindNNeighborsInMesh2(25, i);

		// find closest neighbor who has an assignment
		int clid = -1;
		double cldist = numeric_limits<double>::max();
		for (set<int>::iterator it = nei->begin(); it != nei->end(); it++)
		{
			int nid = *it;
			if ((m2m1[nid] == -1) || (candidates[m2m1[nid]].size() == 0))
				continue;

			double niddist = GetEdgeLength(vtkMesh2, i, nid);
			if (niddist > 3.0 * aveedgelength)
				continue;

			if (niddist < cldist)
			{
				clid = nid;
				cldist = niddist;
			}
		}
		delete nei;
		if (clid == -1) 
			continue;

		// get the point in mesh 2 
		double point[3];
		vtkMesh2->GetPoint(i, point);
		float3 pt = make_float3(point[0], point[1], point[2]);

		// apply transform to the point
		float invT[4][4];
		int err = InverseTransform(candidates[m2m1[clid]][0], invT);
		pt = ApplyTranformation(pt, invT);

		// find the closest point to the transformed coordinates
		int cid = mesh1Locator->FindClosestPoint(pt.x, pt.y, pt.z);
		vtkMesh1->GetPoint(cid, point);
		double dist = length(make_float3(point[0], point[1], point[2]) - pt);

		// check it with existing possibilities
		if (dist < 2.0 * aveedgelength)
		{
			m2m1[i] = cid;
		}
	}

	////// STEP 3:
	// some points in mesh1 now has no points pointing to from mesh2 and also
	// these points are point to a vertex in mesh2 that has a far different 
	// assignment in mesh1. These vertices in mesh1 should be with no assignment.
	vtkSmartPointer<vtkIntArray> newScalars1 = vtkSmartPointer<vtkIntArray>::New();
	newScalars1->SetName("FMTBL");
	for (vtkIdType i = 0; i < vtkMesh1->GetNumberOfPoints(); i++)
	{
		if ((candidates[i].size() > 0) && (m2m1[candidates[i][0].id] != i) && (m2m1[candidates[i][0].id] != -1))
		{
			int idm2 = candidates[i][0].id;
			int idm1_1 = m2m1[idm2];
			int idm1_2 = i;

			double ptm1_1[3];
			vtkMesh1->GetPoint(idm1_1, ptm1_1);
			float3 ptm1_1_f = make_float3(ptm1_1[0], ptm1_1[1], ptm1_1[2]);

			double ptm1_2[3];
			vtkMesh1->GetPoint(idm1_2, ptm1_2);
			float3 ptm1_2_f = make_float3(ptm1_2[0], ptm1_2[1], ptm1_2[2]);

			double ptm2[3];
			vtkMesh2->GetPoint(idm2, ptm2);
			float3 ptm2_f = make_float3(ptm2[0], ptm2[1], ptm2[2]);

			if (abs(length(ptm1_1_f - ptm2_f) - length(ptm1_2_f - ptm2_f)) > 3.0 * aveedgelength)
			{
				candidates[idm1_2].clear();
			}
		}

		if (candidates[i].size() > 0)
			newScalars1->InsertValue(i, candidates[i][0].id);
		else
			newScalars1->InsertValue(i, -1);
	}
	vtkMesh1->GetPointData()->AddArray(newScalars1);

	// create matching array for mesh2
	vtkSmartPointer<vtkIntArray> newScalars2 = vtkSmartPointer<vtkIntArray>::New();
	newScalars2->SetName("BMTBL");
	for (vtkIdType i = 0; i < vtkMesh2->GetNumberOfPoints(); i++)
	{
		newScalars2->InsertValue(i, m2m1[i]);
	}
	vtkMesh2->GetPointData()->AddArray(newScalars2);

	// write files
	printf("Writing mesh %s\n", meshfile1.c_str());
	vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter1 = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkMeshWriter1->SetFileName(meshfile1.c_str());
	vtkMeshWriter1->SetInput(vtkMesh1);
	vtkMeshWriter1->Write();
	printf("Writing mesh %s\n", meshfile2.c_str());
	vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter2 = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkMeshWriter2->SetFileName(meshfile2.c_str());
	vtkMeshWriter2->SetInput(vtkMesh2);
	vtkMeshWriter2->Write();
	printf("Done\n");

}


