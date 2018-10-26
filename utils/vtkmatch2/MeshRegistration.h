#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <limits>
#include <exception>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPolyDataReader.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>
#include <vtkPointSet.h>
#include <vtkIncrementalOctreePointLocator.h>
#include <vtkGenericCell.h>
#include <vector_types.h>
#include <vector_functions.h>
#include <driver_functions.h>
#include <cutil_inline.h>
#include <cutil_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <omp.h>

#include "vtkCurvatureTensor.h"
#include "SolveDALM.h"
#include "PCA.h"

using namespace std;

#define CANDIDATECOUNTMIN 10
#define CANDIDATECOUNTMAX 20
#define MYPI 3.14159265
#define CURVEPS 2.0
#define THREADNUMBER 64

#define CUERRP(x) {printf("%s in file %s at line %d\n", x, __FILE__, __LINE__); fflush(stdout);}
#define CUERR(x) {}

class Edge 
{
public:
	int x;
	int y;
	double length;

	Edge() {}
	Edge(int _x, int _y) : x(_x), y(_y) {}

	bool operator<(const Edge& e2) const
	{
		if (x < e2.x)
		{
			return true;
		}
		else if (x > e2.x)
		{
			return false;
		}
		else
		{
			if (y < e2.y)
			{
				return true;
			}
			else if (y > e2.y)
			{
				return false;
			}
			else
			{
				return false;
			}
		}
	}

	bool operator==(const Edge e1)
	{
		if ((e1.x == x) && (e1.y == y))
			return true;
		return false;
	}

};

class HeapEdge 
{
public:
	int x;
	int y;
	double length;

	HeapEdge() {}
	HeapEdge(int _x, int _y) : x(_x), y(_y) {}

	bool operator<(const HeapEdge& e2) const
	{
		return length > e2.length;
	}
};

class EdgeLengthSorter {
public:
	

    EdgeLengthSorter() {	}

    bool operator()(Edge o1, Edge o2) {
		if (o1.length < o2.length)
			return true;

		if (o1.length > o2.length)
			return false;

		if (o1.length == o2.length)
		{
			return o1 < o2;
		}
    }
};

class BoundaryDistSorter {
public:
	double* vals;

	BoundaryDistSorter(double* _vals) : vals(_vals) {	}

	bool operator()(int o1, int o2) {
		if (vals[o1] < vals[o2])
			return true;
		if (vals[o1] > vals[o2])
			return false;

		if (vals[o1] == vals[o2])
		{
			return o1 < o2;
		}
    }
};

class CurvatureSorter {
public:
	
	vtkPolyData* vtkMesh1;
	vtkPolyData* vtkMesh2;

    CurvatureSorter(vtkPolyData* _vtkMesh1, vtkPolyData* _vtkMesh2) : vtkMesh1(_vtkMesh1),vtkMesh2(_vtkMesh2) {	}

    bool operator()(int2 o1, int2 o2) {
		vtkPolyData* mo1 = (o1.y == 1)? vtkMesh1 : vtkMesh2;
		vtkPolyData* mo2 = (o2.y == 1)? vtkMesh1 : vtkMesh2;

		double cur1[2];
		mo1->GetPointData()->GetScalars("Principal Curvatures")->GetTuple(o1.x, cur1);
		double v1 = max(abs(cur1[0]), abs(cur1[1]));

		double cur2[2];
		mo2->GetPointData()->GetScalars("Principal Curvatures")->GetTuple(o2.x, cur2);
		double v2 = max(abs(cur2[0]), abs(cur2[1]));

        return v1 > v2;
    }
};

class TNode {
public:
	TNode* parent;
	vector<TNode*> children;
	int level;
	int id;

	TNode(TNode* _parent, int _id) { 
		parent = _parent; 
		if (parent == NULL)
			level = 0;
		else
			level = parent->level + 1; 
		id = _id;
	}
};

class Candidate {
public:
	int id;
	float belief;
	float phi;
	float3 point;
	float confidence;

	// transformation info
	float3 rotationaxis;
	float rotationangle;
	float Tg[4][4];
	float Tl[4][4];

	Candidate(int _id, float _belief, float _phi, float3 _point)
	{
		id = _id;
		belief = _belief;
		phi = _phi;
		point = _point;
		confidence = 0.0;

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				if (i != j)
				{
					Tg[i][j] = 0.0;
					Tl[i][j] = 0.0;
				}
				else
				{
					Tg[i][j] = 1.0;
					Tl[i][j] = 1.0;
				}
			}
		}
	}

	Candidate(const Candidate& c)
	{
		id = c.id;
		belief = c.belief;
		phi = c.phi;
		point = c.point;
		confidence = c.confidence;

		rotationaxis = c.rotationaxis;
		rotationangle = c.rotationangle;
		memcpy(Tl, c.Tl, 16 * sizeof(float));
		memcpy(Tg, c.Tg, 16 * sizeof(float));
	}

	void Copy(const Candidate& c)
	{
		id = c.id;
		belief = c.belief;
		phi = c.phi;
		point = c.point;
		confidence = c.confidence;

		rotationaxis = c.rotationaxis;
		rotationangle = c.rotationangle;
		memcpy(Tl, c.Tl, 16 * sizeof(float));
		memcpy(Tg, c.Tg, 16 * sizeof(float));
	}
};

class MeshRegistration
{
public:

	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader1;
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader2;

	// the meshes to register
	string meshfile1;
	string meshfile2;
	vtkPolyData* vtkMesh1;	
	vtkPolyData* vtkMesh2;

	// list of vertex neighbors
	vector< vector<int> > neighborlist;
	double aveedgelength;
	float3 diag;

	// minimum spanning tree for inference
	map<int, TNode*> id2node;
	map<int, int> id2cc;
	int mstreeSize;
	vtkIncrementalOctreePointLocator* mstreeLocator;
	int maxlevel;
	map<int, vector<TNode*> > levelnodes;
	TNode* rootnode;

	// boundary points
	set<vtkIdType>* bpts;

	// sorted list for curvature
	bool* newvertex;
	int idsSortedIndex;
	vector<int2> idsSorted;
	map<int, int> meshid2curvid;

	// candidates information array of vectors
	vector<Candidate>* candidates;

	// locator to find points close to transformation
	vtkIncrementalOctreePointLocator** mesh2Locator;
	vtkIncrementalOctreePointLocator* mesh1Locator;

	// mapping from the second to the first mesh
	int* m2m1;

	// registration status 
	int add_count;
	bool done_reg;

	// statistics information
	int meshshift;
	int correct;
	int visited;
	int levelmeshshift[10000];
	double levelavegrdist[10000];

	// counter for the iterations
	int updown;

	MeshRegistration(char* file1, char* file2);
	~MeshRegistration(void);

	void BuildTree();
	void PrintCandidates(int id, int count);
	double SquarePointDist(int mesh, int id1, int id2);
	void FindClosestOnSurfacePoint(int& id, float3& pt);


	set<int>* FindNNeighborsInMesh2(int n, int id);
	vector<Candidate>* FindVertexCandidates(int id, double curvlimit, bool goingup);
	void CheckCandidate(TNode* child, Candidate& ca, vector<Candidate>* ret, double curvlimit);
	void CheckCandidateCase(int Xi, Candidate& ca, double a, double b, double c, double d, float3 va, float3 vb, float3 vc, float3 vd, vector<Candidate>* ret, double curvlimit);

	double Psi(const int& Xi, Candidate& xi, const int& Xj, Candidate& xj);
	void ComputeTransform(int Xi, Candidate& c, vector<int3>* mappings, bool goingup);
	void AdjustNode(TNode* n, bool goingup, bool lastiter);
	void PropagateUp();
	void PropagateDown();
	bool AdjustGlobalTransform();
	double ComputeTransformationVariance();
	void ComputeLocalTransform(int Xi, Candidate& ca);
	void DoRegistrationStep();

	double PointDist(int idm1, int idm2);
	double ParentAgreement(int Xi, int r);
	void FinalFix();
};

