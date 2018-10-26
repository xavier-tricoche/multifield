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
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <algorithm>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <teem/nrrd.h>
#include <vector_types.h>
#include <vector_functions.h>
#include <cutil_inline.h>
#include <cutil_math.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPolyDataReader.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkPolyDataWriter.h>

using namespace std;

int nvar = 0;
int ncluster = 0;
int ts_count = 0;
int ts_str = 0;
int ts_end = 0;
int fs_count = 0;
char** meshes_files;

class Cluster {
public:
	string pattern;
	int size;
	
	Cluster()
	{
		for (int i = 0; i < 2 * nvar; i++)
			pattern += "x";
			
		size = 0;
	}
	
	Cluster(Cluster c, int bit, int val)
	{
		pattern = string(c.pattern);
		if (val == 0)
			pattern[bit] = '0';
		else
			pattern[bit] = '1';
	}
	
	bool TagMatch(int tag)
	{
		bool match = true;
		for (int i = 0; i < 2 * nvar; i++)
		{
			int bit = ((tag >> i) & 1);
			
			if (((bit == 1) && (pattern[i] == '0')) || 
				((bit == 0) && (pattern[i] == '1')))
			{
				match = false;
				break;
			}
		}
		return match;
	}
	
	void ComputeSize(int* tagcnt)
	{
		size = 0;
		
		int ntags = pow(2, 2 * nvar);
		for (int i = 0; i < ntags; i++)
		{
			if (TagMatch(i))
				size += tagcnt[i];
		}
	}
	
	int FindHighestEntropyBit(int* tagcnt)
	{
		double mindiff = numeric_limits<double>::max();
		int bestbit = -1;
		int ntags = pow(2, 2 * nvar);
		
		// loop on bits
		for (int bit = 0; bit < 2 * nvar; bit++)
		{
			if (pattern[bit] != 'x')
				continue;
				
			int sum0 = 0;
			int sum1 = 0;
			
			// loop on all matching tags
			for (int tag = 0; tag < ntags; tag++)
			{
				if (!TagMatch(tag)) continue;
				
				if (((tag >> bit) & 1) == 0)
					sum0 += tagcnt[tag];
				else
					sum1 += tagcnt[tag];
			}
			
			//cout << "(" << bit << ") " << abs(sum0 - sum1) << "\n";
			
			if (abs(sum0 - sum1) < mindiff)
			{
				bestbit = bit;
				mindiff = abs(sum0 - sum1);
			}
		}
		
		return bestbit;
	}
	
	bool operator<(const Cluster& e2) const
	{
		if (size < e2.size)
			return true;
		else
			return false;
	}
	
	void Print()
	{
		cout << string(pattern.rbegin(), pattern.rend()) << " " << size << "\n";
	}
};

int main (int argc, char *argv[])
{
	if (argc != 7)
	{
		printf("Example: vtkhierarchicalclustering list.txt 5 10 25 8 cluster.txt\n");
		return 0;
	}

	cout << "=========================\n";
	cout << "vtkhierarchicalclustering\n";
	cout << "=========================\n";
	
	nvar = atoi(argv[2]);
	ts_str = atoi(argv[3]);
	ts_end = atoi(argv[4]);
	ncluster = atoi(argv[5]);
	
	cout << "Number of variables is " << nvar << "\n";
	
	// read the file for information
	fstream input(argv[1], std::ios::in);
	input >> ts_count;
	input >> fs_count;
	meshes_files = (char**) malloc(ts_count * sizeof(char*));
	for (int i = 0; i < ts_count; i++)
	{
		string buffer;
		input >> buffer;
		meshes_files[i] = new char[buffer.size()+1];
		strcpy(meshes_files[i], buffer.c_str());
	}
	input.close();

	// create an array for all patterns counts
	int ntags = pow(2, 2 * nvar);
	printf("Number of possible tags is %d\n", ntags);
	int* tagcnt = (int*) malloc(ntags * sizeof(int));
	memset(tagcnt, 0, ntags * sizeof(int));
	
	// loop loading meshes and counting
	vtkPolyData* mesh;
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader;
	for (int i = ts_str; i <= ts_end; i++)
	{
		cout << "time step " << i << "\n";
		
		// read mesh
		cout << "Loading " << meshes_files[i] << "\n";
		vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
		vtkMeshReader->SetFileName(meshes_files[i]);
		vtkMeshReader->Update();
		mesh = vtkMeshReader->GetOutput();
		
		// get table for the tags
		vtkIntArray* tagarr = (vtkIntArray*) mesh->GetPointData()->GetArray("TAGTBL");
		
		// set count
		for (vtkIdType j = 0; j < mesh->GetNumberOfPoints(); j++)
		{
			int tag = tagarr->GetValue(j);
			tagcnt[tag]++;
		}
	}
	cout << "Done counting\n";
	
	// now do the clustering 
	printf("Start the clustering\n");
	vector<Cluster> clusters;
	Cluster fc;
	fc.ComputeSize(tagcnt);
	clusters.push_back(fc);
	make_heap (clusters.begin(),clusters.end());
	while (clusters.size() < ncluster)
	{
		// find biggest cluster
		pop_heap(clusters.begin(),clusters.end()); 
		Cluster c = clusters.back();
		clusters.pop_back();
		
		//cout << "split "; c.Print();
		
		// find the bit with highest entropy 
		int bit = c.FindHighestEntropyBit(tagcnt);
		
		// do split
		Cluster c1(c, bit, 0);
		c1.ComputeSize(tagcnt);
		clusters.push_back(c1); 
		push_heap(clusters.begin(), clusters.end());
		//cout << "new 0: "; c1.Print();
		
		Cluster c2(c, bit, 1);
		c2.ComputeSize(tagcnt);
		clusters.push_back(c2); 
		push_heap(clusters.begin(), clusters.end());
		//cout << "new 1: "; c2.Print();
	}
	
	// print the clusters
	for (vector<Cluster>:: iterator it = clusters.begin(); it != clusters.end(); it++)
	{
		(*it).Print();
	}
	
	// write the result to the files
	for (int i = ts_str; i <= ts_end; i++)
	{
		cout << "time step " << i << "\n";
		
		// read mesh
		cout << "Loading " << meshes_files[i] << "\n";
		vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
		vtkMeshReader->SetFileName(meshes_files[i]);
		vtkMeshReader->Update();
		mesh = vtkMeshReader->GetOutput();
		
		// get table for the tags
		vtkIntArray* tagarr = (vtkIntArray*) mesh->GetPointData()->GetArray("TAGTBL");
		
		// set count
		vtkPointData* pointsData = mesh->GetPointData();
		vtkSmartPointer<vtkIntArray> newScalars = vtkSmartPointer<vtkIntArray>::New();
		newScalars->SetName("CLUSTBL");
		for (vtkIdType j = 0; j < mesh->GetNumberOfPoints(); j++)
		{
			int tag = tagarr->GetValue(j);
			
			for (int k = 0; k < clusters.size(); k++)
			{
				if (clusters[k].TagMatch(tag))
				{
					newScalars->InsertValue(j, k);
				}
			}
		}
		pointsData->AddArray(newScalars);
		
		// write file
		vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
		vtkMeshWriter->SetFileName(meshes_files[i]);
		vtkMeshWriter->SetInput(mesh);
		vtkMeshWriter->Write();
	}
	
	// free memory
	free(tagcnt);
	
	return 0;
}