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
#include "GenerateTexture.h"

int cint(double x){

	double fractpart, intpart;
	fractpart = modf (x , &intpart);

	if (fractpart>=.5)
		return x>=0?ceil(x):floor(x);
	else
		return x<0?ceil(x):floor(x);
}

double dist_Point_to_Line(float3 P, float3 LP0, float3 LP1)
{
    float3 v = LP1 - LP0;
    float3 w = P - LP0;

    double c1 = dot(w,v);
    double c2 = dot(v,v);
    double b = c1 / c2;

    float3 Pb = LP0 + b * v;
    return length(P - Pb);
}

typedef struct {
   char id_len;                 // ID Field (Number of bytes - max 255)
   char map_type;               // Colormap Field (0 or 1)
   char img_type;               // Image Type (7 options - color vs. compression)
   int  map_first;              // Color Map stuff - first entry index
   int  map_len;                // Color Map stuff - total entries in file
   char map_entry_size;         // Color Map stuff - number of bits per entry
   int  x;                      // X-coordinate of origin 
   int  y;                      // Y-coordinate of origin
   int  width;                  // Width in Pixels
   int  height;                 // Height in Pixels
   char bpp;                    // Number of bits per pixel
   char misc;                   // Other stuff - scan origin and alpha bits
} targa_header;

void writeheader(targa_header h, FILE *tga) 
  {
   fputc(h.id_len, tga);          // Write chars for ID, map, and image type
   fputc(h.map_type, tga);
   fputc(h.img_type, tga);
   fputc(h.map_first % 256, tga); // Write integer, low order byte first
   fputc(h.map_first / 256, tga); // Write second byte of integer, high order
   fputc(h.map_len % 256, tga);   // Another integer 
   fputc(h.map_len / 256, tga);
   fputc(h.map_entry_size, tga);  // Write a char - only one byte
   fputc(h.x % 256, tga);         // More integers
   fputc(h.x / 256, tga);
   fputc(h.y % 256, tga);
   fputc(h.y / 256, tga);
   fputc(h.width % 256, tga);     // Even more integers
   fputc(h.width / 256, tga);
   fputc(h.height % 256, tga);
   fputc(h.height / 256, tga);
   fputc(h.bpp, tga);             // Write two chars
   fputc(h.misc, tga);
}

int write_truecolor_tga(string filename, RGB_t* data, unsigned width, unsigned height )
{
	printf("Color image dimensions %d %d -> %s\n", width, height, filename.c_str());
   FILE *tga;               // Pointer to a FILE
   targa_header header;     // Variable of targa_header type
   int x, y;

/* First, set all the fields in the header to appropriate values */
   header.id_len = 0;          /* no ID field */
   header.map_type = 0;        /* no colormap */
   header.img_type = 2;        /* trust me */
   header.map_first = 0;       /* not used */
   header.map_len = 0;         /* not used */
   header.map_entry_size = 0;  /* not used */
   header.x = 0;               /* image starts at (0,0) */
   header.y = 0;
   header.width = width;         /* image is 200 x 100 */
   header.height = height;
   header.bpp = 24;            /* 24 bits per pixel */
   header.misc = 0x20;         /* scan from upper left corner */

/* Open a file for writing targa data.  Call the file "test.tga" and
      write in binary mode (wb) so that nothing is lost as characters
      are written to the file */

   tga = fopen(filename.c_str(), "wb"); /* Write the header information  */

   writeheader(header, tga);  

   /* Write the data for a graphic that measures 100 by 200 pixels. */

   for(y = 0; y < height; y++)      // Create 100 Rows of Pixels 
      for(x = 0; x < width; x++)   // Create 200 Pixels in each Row

       { /* For each pixel, write a character representing the RGB color.
            Notice that the order that is written to the file is B-G-R.
            This sequence just cycles through the colors in some pattern
            but all char values must be integers between 0 and 255. */

		int yf = height - 1 - y;
        fputc(data[x + width * yf].blue, tga);               // Write char for BLUE            
        fputc(data[x + width * yf].green, tga);   // Write char for GREEN
        fputc(data[x + width * yf].red, tga);               // Write char for RED
       }

/* close the file */
   fclose(tga);

/* that was easy, right? */
   return 0;
}

void WriteBMPImage(string filename, RGB_t* data, int w, int h)
{
	BMP Output;
	Output.SetSize( w , h );
	Output.SetBitDepth( 24 );

	for (int x = 0; x < w; x++)
	{
		for (int y = 0; y < h; y++)
		{
			Output(x,h - 1 - y)->Red = data[x + w * y].alpha;
			Output(x,h - 1 - y)->Blue = data[x + w * y].alpha;
			Output(x,h - 1 - y)->Green = data[x + w * y].alpha;
		}
	}
	
	Output.SetBitDepth( 8 );

	if( Output.TellBitDepth() < 16 )
		{ CreateGrayscaleColorTable( Output ); }

	cout << "writing 8bpp ... " << endl;					
	Output.WriteToFile( filename.c_str() );

}

void GenerateColorTexture(vtkPolyData* vtkMesh, string texprefix, FILE* objfile)
{
	// create the set of trinagles T
	vector<Triangle> T;
	vtkCellArray* vtkcells = vtkMesh->GetPolys();
	vtkcells->InitTraversal();
	vtkIdType npts;
	vtkIdType* pts;
	double averageheight = 0.0;
	while (vtkcells->GetNextCell(npts, pts) != 0)
	{
		T.push_back(Triangle(vtkMesh, pts));
		if (isnan(T.back().GetHeight()))
		{
			printf("OMG!\n");
			T.pop_back();
			continue;
		}
		averageheight += T.back().GetHeight();
	}
	averageheight /= T.size();
	make_heap(T.begin(), T.end(), CompareHeight());
	double scale = 3.0 / averageheight;
	printf("Average Height %lf and scale is %lf\n", averageheight, scale);
	
	// mtl file
	string mtlfile = texprefix + string(".mtl");
	FILE* pFile = fopen (mtlfile.c_str(),"w");
		
	// create the textures	
	int vtcount = 0;
	int s = 1024;
	vector<QuadraticTexture*> textures;
	int height;
	int length;
	while (!T.empty())
	{
		QuadraticTexture* texture = new QuadraticTexture(textures.size(), texprefix, s, s);
		textures.push_back(texture);
		fprintf(objfile, "usemtl %s\n", texture->filename.c_str());
		
		height = 0;
		while (!T.empty())
		{
			Triangle t = T.front();
			height = height + t.GetPixelHeight(scale);
			if (height > s)
				break;
			pop_heap(T.begin(), T.end(), CompareHeight()); T.pop_back();
			texture->InsertTriangleInNewRow(t, scale, objfile, vtcount);
			length = t.GetPixelBaseLength(scale);
			while (!T.empty())
			{
				Triangle t = T.front();
				length = length + t.GetPixelBaseLength(scale);
				if (length > s)
					break;
				pop_heap(T.begin(), T.end(), CompareHeight()); T.pop_back();
				texture->InsertTriangleInSameRow(t, scale, objfile, vtcount);
			}
		}
		
		// write files for textures
		string tgafile = texture->filename + string(".tga");
		string bmpfile = texture->filename + string(".bmp");
		write_truecolor_tga(tgafile, texture->data, texture->width, texture->height);
		WriteBMPImage(bmpfile, texture->data, texture->width, texture->height);
		texture->Destruct();
		
		// write the material file
		fprintf(pFile, "newmtl %s\n", texture->filename.c_str());
		fprintf(pFile, "Ka 1.000 1.000 1.000\n");
		fprintf(pFile, "Kd 1.000 1.000 1.000\n");
		fprintf(pFile, "Ks 0.000 0.000 0.000\n");
		fprintf(pFile, "Ns 0.0\n");
		fprintf(pFile, "d 1.0\n");
		fprintf(pFile, "illum 9\n");
		fprintf(pFile, "map_Ka %s\n", tgafile.c_str());      
		fprintf(pFile, "map_Kd %s\n", tgafile.c_str());      		  
		fprintf(pFile, "map_Ks %s\n", tgafile.c_str());      
		fprintf(pFile, "map_d %s\n", bmpfile.c_str());
	}
	
	// mesh of only gray color
	fprintf(pFile, "newmtl %sg\n", texprefix.c_str());
	fprintf(pFile, "Ka 0.100 0.100 0.100\n");
	fprintf(pFile, "Kd 0.100 0.100 0.100\n");
	fprintf(pFile, "Ks 0.000 0.000 0.000\n");
	fprintf(pFile, "Ns 0\n");
	fprintf(pFile, "d 1.0\n");
	fprintf(pFile, "illum 9\n");
		
	fclose (pFile);
}
