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
#include "FieldTypeGraph.h"


FieldTypeGraph::FieldTypeGraph(int _nvar)
{
	nvar = _nvar;

	for (int i = 0; i < nvar; i++)
	{
		fieldtypevertex_t v1 = boost::add_vertex(graph);
		graph[v1].field = i;
		graph[v1].type = 0;
		nodes[graph[v1]] = v1;

		fieldtypevertex_t v2 = boost::add_vertex(graph);
		graph[v2].field = i;
		graph[v2].type = 1;
		nodes[graph[v2]] = v2;
	}

	for (int b1 = 0; b1 < 2 * nvar; b1++)
	{
		for (int b2 = b1 + 1; b2 < 2 * nvar; b2++)
		{
			correlationmatrix[VariableLink(b1,b2)] = make_float4(0.0);
		}
	}
}


FieldTypeGraph::~FieldTypeGraph(void)
{
}

void FieldTypeGraph::UpdateGraphWithInput(vtkPolyData* mesh)
{
	// get the tag table
	vtkIntArray* tagsarr = (vtkIntArray*) mesh->GetPointData()->GetArray("TAGTBL");

	// get links weights map
	boost::property_map<FieldTypeGraph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, graph);

	// loop on vertices creating links
	for (int i = 0; i < mesh->GetNumberOfPoints(); i++)
	{
		int tag = tagsarr->GetValue(i);

		for (int b1 = 0; b1 < 2 * nvar; b1++)
		{
			int bit1 = (tag >> b1) & 1;
			if (!bit1) continue;

			for (int b2 = 0; b2 < 2 * nvar; b2++)
			{
				if (b1 == b2) continue;

				int bit2 = (tag >> b2) & 1;
				if (!bit2) continue;

				fieldtypevertex_t v1 = nodes[FieldTypeNode(b1/2, b1 % 2)];
				fieldtypevertex_t v2 = nodes[FieldTypeNode(b2/2, b2 % 2)];

				// both structure types coexist
				bool bRet;
				fieldtypeedge_t e;
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
	}
}

void FieldTypeGraph::UpdateGraphWithInput2(vtkPolyData* mesh)
{
	// get the tag table
	vtkIntArray* tagsarr = (vtkIntArray*) mesh->GetPointData()->GetArray("TAGTBL");

	// get links weights map
	boost::property_map<FieldTypeGraph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, graph);

	// loop on vertices creating links
	for (int i = 0; i < mesh->GetNumberOfPoints(); i++)
	{
		int tag = tagsarr->GetValue(i);

		for (int b1 = 0; b1 < 2 * nvar; b1++)
		{
			int bit1 = (tag >> b1) & 1;

			for (int b2 = b1 + 1; b2 < 2 * nvar; b2++)
			{
				int bit2 = (tag >> b2) & 1;

				
				if ((bit1 == 0) && (bit2 == 0))
					correlationmatrix[VariableLink(b1,b2)].x += 1.0;
				if ((bit1 == 0) && (bit2 == 1))
					correlationmatrix[VariableLink(b1,b2)].y += 1.0;
				if ((bit1 == 1) && (bit2 == 0))
					correlationmatrix[VariableLink(b1,b2)].z += 1.0;
				if ((bit1 == 1) && (bit2 == 1))
					correlationmatrix[VariableLink(b1,b2)].w += 1.0;
			}
		}
	}
}

void FieldTypeGraph::BuildGraph()
{
	// get links weights map
	boost::property_map<FieldTypeGraph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, graph);

	// on corelation matrices adding links
	for (int b1 = 0; b1 < 2 * nvar; b1++)
	{
		for (int b2 = b1 + 1; b2 < 2 * nvar; b2++)
		{
			float4 mat = correlationmatrix[VariableLink(b1,b2)];
			float num = mat.x * mat.w - mat.y * mat.z;
			float den = (mat.x + mat.y) * (mat.z + mat.w) * (mat.x + mat.z) * (mat.y + mat.w);
			den = sqrtf(den);

			fieldtypevertex_t v1 = nodes[FieldTypeNode(b1/2, b1 % 2)];
			fieldtypevertex_t v2 = nodes[FieldTypeNode(b2/2, b2 % 2)];

			// both structure types coexist
			bool bRet;
			fieldtypeedge_t e;
			tie(e,bRet) = boost::edge(v1, v2, graph);

			if (bRet == true)
			{
				put(weightmap, e, get(weightmap, e) + 1.0);
				printf("How come!\n");
			}
			else
			{
				boost::add_edge(v1, v2, num / den, graph);
			}

		}
	}
}

void FieldTypeGraph::WriteFile(string filename)
{
	// initialize Graph
	FieldTypeGraph_o ng;

	// copy the vertices
	map<FieldTypeNode, fieldtypevertex_o> nodes_o;
	boost::graph_traits<FieldTypeGraph_t>::vertex_iterator vi, vend;
	for (boost::tie(vi, vend) = vertices(graph); vi != vend; ++vi) 
	{
		fieldtypevertex_o v = boost::add_vertex(ng);
		ng[v].field = graph[*vi].field;
		ng[v].type = graph[*vi].type;
		nodes_o[ng[v]] = v;
	}
	
	// copy the edges
	boost::property_map<FieldTypeGraph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, graph);
	boost::graph_traits<FieldTypeGraph_t>::edge_iterator ei, ei_end;
	for (boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei) 
	{
		fieldtypevertex_t gv1 = source(*ei, graph);
		fieldtypevertex_t gv2 = target(*ei, graph);
		fieldtypevertex_o v1 = nodes_o[graph[gv1]];
		fieldtypevertex_o v2 = nodes_o[graph[gv2]];
		boost::add_edge(v1, v2, get(weightmap, *ei), ng);
	}
		
	// write the graph to a file
	std::ofstream output(filename.c_str(), std::ios::out);
	boost::write_graphviz(output, ng, FieldTypeGraph_writer(&ng), FieldTypeGraph_writer(&ng));
	output.close();
}