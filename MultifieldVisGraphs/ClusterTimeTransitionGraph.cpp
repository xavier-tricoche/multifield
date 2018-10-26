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
#include "ClusterTimeTransitionGraph.h"

ClusterTimeTransitionGraph::ClusterTimeTransitionGraph(void)
{
}


ClusterTimeTransitionGraph::~ClusterTimeTransitionGraph(void)
{
}

void ClusterTimeTransitionGraph::CreateGraphNodes(vector<ClusterNode>* clusters, int ts_str, int ts_end)
{
	for (int t = ts_str; t <= ts_end; t++)
	{
		for (int i = 0; i < clusters->size(); i++)
		{
			clustertimetransitionvertex_t v = boost::add_vertex(graph);
			graph[v].cluster = clusters->at(i).cluster;
			graph[v].size = clusters->at(i).size;
			graph[v].nvar = clusters->at(i).nvar;
			graph[v].entropy = clusters->at(i).entropy;
			graph[v].mpclusterint = clusters->at(i).mpclusterint;
			graph[v].mpclusterstr = clusters->at(i).mpclusterstr;
			graph[v].ts = t;
			nodes[graph[v]] = v;
		}
	}
}

void ClusterTimeTransitionGraph::UpdateGraphWithInput(int ts1, vtkPolyData* mesh1, int ts2, vtkPolyData* mesh2, vector<ClusterNode>* clusters)
{
	// get the arrays
	vtkIntArray* mesh1clustertimetransitions = (vtkIntArray*) mesh1->GetPointData()->GetArray("CLUSTBL");
	vtkIntArray* mesh2clustertimetransitions = (vtkIntArray*) mesh2->GetPointData()->GetArray("CLUSTBL");

	// get links weights map
	boost::property_map<ClusterTimeTransitionGraph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, graph);

	// add links
	vtkIntArray* m12 = (vtkIntArray*) mesh1->GetPointData()->GetArray("FMTBL");
	for (int id1 = 0; id1 < mesh1->GetNumberOfPoints(); id1++)
	{
		int id2 = m12->GetValue(id1);
		if (id2 == -1) continue;

		int comp1 = mesh1clustertimetransitions->GetValue(id1);
		int comp2 = mesh2clustertimetransitions->GetValue(id2);

		clustertimetransitionvertex_t v1 = nodes[ClusterTimeNode(clusters->at(comp1), ts1)];
		clustertimetransitionvertex_t v2 = nodes[ClusterTimeNode(clusters->at(comp2), ts2)];

		bool bRet;
		clustertimetransitionedge_t e;
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

void ClusterTimeTransitionGraph::WriteFile(string filename)
{
	// initialize Graph
	ClusterTimeTransitionGraph_o ng;

	// copy the vertices
	map<ClusterTimeNode, clustertimetransitionvertex_o> nodes_o;
	boost::graph_traits<ClusterTimeTransitionGraph_t>::vertex_iterator vi, vend;
	for (boost::tie(vi, vend) = vertices(graph); vi != vend; ++vi)
	{
		clustertimetransitionvertex_o v = boost::add_vertex(ng);
		ng[v].cluster = graph[*vi].cluster;
		ng[v].size = graph[*vi].size;
		ng[v].nvar = graph[*vi].nvar;
		ng[v].entropy = graph[*vi].entropy;
		ng[v].mpclusterint = graph[*vi].mpclusterint;
		ng[v].mpclusterstr = graph[*vi].mpclusterstr;
		ng[v].ts = graph[*vi].ts;
		nodes_o[ng[v]] = v;
	}

	// copy the edges
	boost::property_map<ClusterTimeTransitionGraph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, graph);
	boost::graph_traits<ClusterTimeTransitionGraph_t>::edge_iterator ei, ei_end;
	for (boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei)
	{
		clustertimetransitionvertex_t gv1 = source(*ei, graph);
		clustertimetransitionvertex_t gv2 = target(*ei, graph);
		clustertimetransitionvertex_o v1 = nodes_o[graph[gv1]];
		clustertimetransitionvertex_o v2 = nodes_o[graph[gv2]];
		boost::add_edge(v1, v2, get(weightmap, *ei), ng);
	}

	// write the graph to a file
	std::ofstream output(filename.c_str(), std::ios::out);
	boost::write_graphviz(output, ng, ClusterTimeTransitionGraph_writer(&ng), ClusterTimeTransitionGraph_writer(&ng));
	output.close();
}
