#ifndef PNMIS
#define PNMIS

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
#include "BulletCollision/CollisionDispatch/btCompoundCollisionAlgorithm.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletCollision/CollisionShapes/btCompoundShape.h"
#include "BulletCollision/BroadphaseCollision/btDbvt.h"
#include "LinearMath/btIDebugDraw.h"
#include "LinearMath/btAabbUtil2.h"
#include <BulletCollision/BroadphaseCollision/btBroadphaseInterface.h>
#include <BulletCollision/BroadphaseCollision/btDbvtBroadphase.h>

#include "MyGeometry.h"
#include "Vertex_descriptor.h"

using namespace std;



struct	BroadphaseAabbTester : public btDbvt::ICollide
{
	btBroadphaseAabbCallback& m_aabbCallback;
	BroadphaseAabbTester(btBroadphaseAabbCallback& orgCallback)
		:m_aabbCallback(orgCallback)
	{
	}
	void Process(const btDbvtNode* leaf)
	{
		btDbvtProxy*	proxy=(btDbvtProxy*)leaf->data;
		m_aabbCallback.process(proxy);
	}
};	

class MyAabbCallback : public btBroadphaseAabbCallback
{
public:
	set<int> querylines;

	bool process(const btBroadphaseProxy  *proxy)
	{
		int i = int(proxy->m_clientObject);
		querylines.insert(i);
		return true;
	}
};

class NMIS_Vertex
{
public:
	double V;
	double U;
	TargetPoint data;

	bool IsConflicting(NMIS_Vertex& v, double distance, double limit)
	{
		if ((data.origin != v.data.origin) &&
			//(intersect2D_Segments(Line(data.origin, data.target), Line(v.data.origin, v.data.target)) > 0))
			(dist3D_Segment_to_Segment(Line(data.origin, data.target), Line(v.data.origin, v.data.target)) < distance))
		{
			//if (abs(dot(normalize(data.origin - data.target), normalize(v.data.origin - v.data.target))) < limit)
				return true;
		}
		return false;
	}
};

double hFunc(double a)
{
	if (a == 0.0)
		return 1.0;
	return 0.0;
}

vector<NMIS_Vertex>* FindNearMaximumIndependentSet(vector<NMIS_Vertex>* vertices, int iterations)
{
	// create the collision tree
	btDbvt* tree = new btDbvt();
	tree->optimizeIncremental(5);
	for (int i = 0; i < vertices->size(); i++)
	{
		float2 pt0 = vertices->at(i).data.origin;
		float2 pt1 = vertices->at(i).data.target;
		btVector3 minv(min(pt0.x, pt1.x), min(pt0.y, pt1.y), 0.0);
		btVector3 maxv(max(pt0.x, pt1.x), max(pt0.y, pt1.y), 1.0);
		btDbvtVolume volume( btDbvtAabbMm::FromMM(minv, maxv));		

		btDbvtNode* node = tree->insert(volume, new int(i));
	}

	// now build the adjacency list
	double distance = 0.008;
	double limit = 0.94;
	MyAabbCallback callback;
	vector<set<int>*> adjacentlist;
	for (int i = 0; i < vertices->size(); i++)
	{
		adjacentlist.push_back(new set<int>());

		float2 pt0 = vertices->at(i).data.origin;
		float2 pt1 = vertices->at(i).data.target;

		btVector3 minv(min(pt0.x, pt1.x), min(pt0.y, pt1.y), 0.0);
		minv.setX(minv.x() - distance);
		minv.setY(minv.y() - distance);

		btVector3 maxv(max(pt0.x, pt1.x), max(pt0.y, pt1.y), 1.0);
		maxv.setX(maxv.x() + distance);
		maxv.setY(maxv.y() + distance);

		btDbvtVolume volume(btDbvtAabbMm::FromMM(minv, maxv));		

		callback.querylines.clear();
		tree->collideTV(tree->m_root, volume, (BroadphaseAabbTester)callback);

		for (set<int>::iterator it = callback.querylines.begin(); it != callback.querylines.end(); it++)
		{
			if (vertices->at(i).IsConflicting(vertices->at(*it), distance, limit))
				adjacentlist[i]->insert(*it);
		}
	}


	// Initial values for U_i
	for (int i = 0; i < vertices->size(); i++)
	{
		vertices->at(i).U = -0.1;
	}

	// interpolator function
	double ipx[5] = {-1.0, 0.0, 0.94, 1.0, 2.0};
	double ipy[5] = {1.0, 1.0, 0.0, -1.0, -1.0};
	gsl_akima_interpolator* ip = new gsl_akima_interpolator(5, ipx, ipy);
	for (float i = 0.0; i < 1.0; i+=0.01)
	{
		printf("%f\t%f\n", i, ip->GetValue(i));
	}

	// start learning loop
	double stopc = 0.0;
	for (int iteration = 0; iteration < iterations; iteration++)
	{
		// evaluate values for V_i
		int included = 0;
		for (int i = 0; i < vertices->size(); i++)
		{
			if (vertices->at(i).U > 0.0)
				vertices->at(i).V = 1.0;
			else
			{
				included++;
				vertices->at(i).V = 0.0;
			}
		}
		printf("%d included.\n", included);
		stopc = 0.95 * stopc + 0.05 * included;
		if (abs(stopc - included) < 1.0)
			break;

		// compute the derivative
		for (int i = 0; i < vertices->size(); i++)
		{
			double dU = 0.0;
			double sigma = 0.0;
			for (set<int>::iterator it = adjacentlist[i]->begin(); it != adjacentlist[i]->end(); it++)
			{
				double ag = abs(dot(normalize(vertices->at(i).data.origin - vertices->at(i).data.target), normalize(vertices->at(*it).data.origin - vertices->at(*it).data.target)));

				//sigma += (1.0 - vertices->at(*it).V) * ip->GetValue(ag);
				sigma += (1.0 - vertices->at(*it).V) * ((ag < limit)? 1.0 : 0.0);
			}

			dU =  sigma * (1.0 - vertices->at(i).V)		// is included and others want it out
				- hFunc(sigma) * vertices->at(i).V;		// is not included and no one disagrees

			vertices->at(i).U += 0.2 * dU;
		}
	}

	vector<NMIS_Vertex>* output = new vector<NMIS_Vertex>();
	for (int i = 0; i < vertices->size(); i++)
	{
		if (vertices->at(i).U > 0.0)
			vertices->at(i).V = 1.0;
		else
			vertices->at(i).V = 0.0;

		if (vertices->at(i).V == 0.0)
		{
			output->push_back(vertices->at(i));
		}
	}

	return output;

}
#endif