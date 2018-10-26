#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <string.h>
#include <math.h>
#include <time.h>
#include <map>
#include <vector_types.h>
#include <vector_functions.h>
#include <cutil_inline.h>
#include <cutil_math.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

using namespace std;

class VariableLink 
{
public:
	int x;
	int y;
	double length;

	VariableLink() {}
	VariableLink(int _x, int _y)  
	{
		x = _x;
		y = _y;
	}

	bool operator<(const VariableLink& e2) const
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

	bool operator==(const VariableLink e1)
	{
		if ((e1.x == x) && (e1.y == y))
			return true;
		return false;
	}

};

class FieldTypeNode
{ 
public:
	
	int field;
	int type;

	FieldTypeNode()
	{
	}

	FieldTypeNode(int _field, int _type) 
	{
		field = _field;
		type = _type;
	}
	
	bool operator<(const FieldTypeNode& node) const
	{
		if (field < node.field)
			return true;

		if ((field == node.field) && (type < node.type))
			return true;

		return false;
	};

	bool operator==(const FieldTypeNode& node) const
	{
		if ((field == node.field) && (type == node.type))
			return true;
		
		return false;
	}
};




////// (1) undirected graph 
typedef boost::adjacency_list
  <
    //Store all edges as a std::vector
    boost::setS,
    //Store all vertices in a std::vector
    boost::setS,
    //Relations are both ways (in this example)
    //(note: but you can freely change it to boost::directedS)
    boost::undirectedS,
    //All vertices are person names of type std::string
    FieldTypeNode,
    //All edges are weights equal to the encounter frequencies
    boost::property<boost::edge_weight_t,double>,
    //Graph itself has a std::string name
    boost::property<boost::graph_name_t,std::string>
  > FieldTypeGraph_t;
typedef boost::graph_traits<FieldTypeGraph_t>::vertex_descriptor fieldtypevertex_t;
typedef boost::graph_traits<FieldTypeGraph_t>::edge_descriptor fieldtypeedge_t;

////// (2) graph to be written in file
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, 
	  FieldTypeNode,
	  boost::property<boost::edge_weight_t, double>,
	  boost::property<boost::graph_name_t,std::string>
> FieldTypeGraph_o;
typedef boost::graph_traits<FieldTypeGraph_o>::vertex_descriptor fieldtypevertex_o;
typedef boost::graph_traits<FieldTypeGraph_o>::edge_descriptor fieldtypeedge_o;

struct FieldTypeGraph_writer {

	// my_node_writer() {}
	FieldTypeGraph_writer(FieldTypeGraph_o* g_) : g (g_) {};

	void operator()(std::ostream& out, const fieldtypevertex_o& v) {
		string strt;
		if ((*g)[v].type == 0)
			strt = string("r");
		else
			strt = string("v");
		out << "[label=\"" << (*g)[v].field << "/" << strt << "\"]";
	};

	void operator()(std::ostream& out, const fieldtypeedge_o& e) {
		boost::property_map<FieldTypeGraph_o, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, *g);
		out << "[weight=" << (get(weightmap, e) / 10000.0) << "]";
	};

	// bleah.  I can't get references right...
	// according to http://www.knowledgesearch.org/doc/examples.html
	// it should be a reference here, not the type itself.
	// but g++ either barfs, or the program segfaults.
	FieldTypeGraph_o* g;
};

class FieldTypeGraph
{
public:

	int nvar;

	map<VariableLink, float4> correlationmatrix;

	FieldTypeGraph_t graph;
	map<FieldTypeNode, fieldtypevertex_t> nodes;

	FieldTypeGraph(int _nvar);
	~FieldTypeGraph(void);

	void UpdateGraphWithInput(vtkPolyData* mesh);
	void UpdateGraphWithInput2(vtkPolyData* mesh);
	void BuildGraph();
	void WriteFile(string filename);
};

