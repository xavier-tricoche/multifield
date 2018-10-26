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
#include "ComponentGraph.h"


ComponentGraph::ComponentGraph(void)
{
}


ComponentGraph::~ComponentGraph(void)
{
}


void ComponentGraph::UpdateGraphWithInput(int ts1, vtkPolyData* mesh1, int ts2, vtkPolyData* mesh2)
{
	// add nodes for the first mesh
	vtkIntArray* mesh1components = (vtkIntArray*) mesh1->GetPointData()->GetArray("CCTBL");
	for (int i = 0; i < mesh1->GetNumberOfPoints(); i++)
	{
		ComponentNode cn(ts1, mesh1components->GetValue(i));

		if (nodes.find(cn) == nodes.end())
		{
			componentvertex_t v = boost::add_vertex(graph);
			graph[v].timestep = ts1;
			graph[v].component = mesh1components->GetValue(i);
			nodes[cn] = v;
		}
	}

	// add nodes for the second mesh
	vtkIntArray* mesh2components = (vtkIntArray*) mesh2->GetPointData()->GetArray("CCTBL");
	for (int i = 0; i < mesh2->GetNumberOfPoints(); i++)
	{
		ComponentNode cn(ts2, mesh2components->GetValue(i));

		if (nodes.find(cn) == nodes.end())
		{
			componentvertex_t v = boost::add_vertex(graph);
			graph[v].timestep = ts2;
			graph[v].component = mesh2components->GetValue(i);
			nodes[cn] = v;
		}
	}

	// get links weights map
	boost::property_map<componentgraph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, graph);

	// add links
	vtkIntArray* m12 = (vtkIntArray*) mesh1->GetPointData()->GetArray("FMTBL");
	for (int id1 = 0; id1 < mesh1->GetNumberOfPoints(); id1++)
	{
		int id2 = m12->GetValue(id1);
		int comp1 = mesh1components->GetValue(id1);
		int comp2 = mesh2components->GetValue(id2);

		componentvertex_t v1 = nodes[ComponentNode(ts1, comp1)];
		componentvertex_t v2 = nodes[ComponentNode(ts2, comp2)];

		bool bRet;
		componentedge_t e;
		tie(e,bRet) = boost::edge(v1, v2, graph);

		if (bRet == true)
		{
			put(weightmap, e, get(weightmap, e) + 1.0);
		}
		else
		{
			boost::add_edge(v1, v2, 1.0, graph);
		}
	}
}

void ComponentGraph::WriteFile(string filename)
{
	// initialize Graph
	componentgraph_o ng;

	// copy the vertices
	map<ComponentNode, componentvertex_o> nodes_o;
	boost::graph_traits<componentgraph_t>::vertex_iterator vi, vend;
	for (boost::tie(vi, vend) = vertices(graph); vi != vend; ++vi) 
	{
		componentvertex_o v = boost::add_vertex(ng);
		ng[v].timestep = graph[*vi].timestep;
		ng[v].component = graph[*vi].component;
		nodes_o[ng[v]] = v;
	}
	
	// copy the edges
	boost::property_map<componentgraph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, graph);
	boost::graph_traits<componentgraph_t>::edge_iterator ei, ei_end;
	for (boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei) 
	{
		componentvertex_t gv1 = source(*ei, graph);
		componentvertex_t gv2 = target(*ei, graph);
		componentvertex_o v1 = nodes_o[graph[gv1]];
		componentvertex_o v2 = nodes_o[graph[gv2]];
		boost::add_edge(v1, v2, get(weightmap, *ei), ng);
	}
		
	// write the graph to a file
	std::ofstream output(filename.c_str(), std::ios::out);
	boost::write_graphviz(output, ng, componentgraph_writer(&ng));
	output.close();
}