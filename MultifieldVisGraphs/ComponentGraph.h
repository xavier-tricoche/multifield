#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <string.h>
#include <math.h>
#include <time.h>
#include <map>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

using namespace std;

class ComponentNode
{ 
public:
	
	int timestep;
	int component;

	ComponentNode()
	{
	}

	ComponentNode(int _timestep, int _component) 
	{
		timestep = _timestep;
		component = _component;
	}
	
	bool operator<(const ComponentNode& node) const
	{
		if (timestep < node.timestep)
			return true;

		if ((timestep == node.timestep) && (component < node.component))
			return true;

		return false;
	};

	bool operator==(const ComponentNode& node) const
	{
		if ((timestep == node.timestep) && (component == node.component))
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
    ComponentNode,
    //All edges are weights equal to the encounter frequencies
    boost::property<boost::edge_weight_t,double>,
    //Graph itself has a std::string name
    boost::property<boost::graph_name_t,std::string>
  > componentgraph_t;
typedef boost::graph_traits<componentgraph_t>::vertex_descriptor componentvertex_t;
typedef boost::graph_traits<componentgraph_t>::edge_descriptor componentedge_t;

////// (2) graph to be written in file
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, 
	  ComponentNode,
	  boost::property< boost::edge_weight_t, double >,
	  boost::property<boost::graph_name_t,std::string>
> componentgraph_o;
typedef boost::graph_traits<componentgraph_o>::vertex_descriptor componentvertex_o;
typedef boost::graph_traits<componentgraph_o>::edge_descriptor componentedge_o;
struct componentgraph_writer {
	// my_node_writer() {}
	componentgraph_writer(componentgraph_o* g_) : g (g_) {};
	void operator()(std::ostream& out, const componentvertex_o& v) {
		out << "[label=\"" << (*g)[v].timestep << "/" << (*g)[v].component << "\"]";
	};
	// bleah.  I can't get references right...
	// according to http://www.knowledgesearch.org/doc/examples.html
	// it should be a reference here, not the type itself.
	// but g++ either barfs, or the program segfaults.
	componentgraph_o* g;
};




// main class definition
class ComponentGraph
{
public:

	componentgraph_t graph;
	map<ComponentNode, componentvertex_t> nodes;

	ComponentGraph(void);
	~ComponentGraph(void);

	void UpdateGraphWithInput(int ts1, vtkPolyData* mesh1, int ts2, vtkPolyData* mesh2);
	void ComponentGraph::WriteFile(string filename);
};

