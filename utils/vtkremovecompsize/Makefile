all:
	gcc-4.4 -I/usr/include/atlas \
		-I/media/SDA2/scratch/sbarakat/Libraries/cuda/include \
		-I/media/SDA2/scratch/sbarakat/Libraries/NVIDIA_GPU_Computing_SDK/C/common/inc \
		-I/media/SDA2/scratch/sbarakat/Libraries/nvis/nvis/ \
		-I/media/SDA2/scratch/sbarakat/Libraries/vtk-5.8.0/vtk-binaries/include/vtk-5.8/ \
		-I/media/SDA2/scratch/sbarakat/Libraries/CGAL-4.1/include \
		-I/usr/include -I/usr/local/include -I/opt/local/include \
		-L/media/SDA2/scratch/sbarakat/Libraries/vtk-5.8.0/vtk-binaries/lib/vtk-5.8/ \
		-L/media/SDA2/scratch/sbarakat/Libraries/CGAL-4.1/lib \
		-L/usr/lib -L/usr/local/lib -L/opt/local/lib \
		-lMapReduceMPI -lmpistubs -lvtkalglib -lvtkCharts -lvtkCommon -lvtkDICOMParser -lvtkmetaio -lvtkNetCDF -lvtkNetCDF_cxx -lvtkpng -lvtkproj4 -lvtkRendering -lvtksqlite -lvtkexpat -lvtkFiltering -lvtkfreetype -lvtkftgl -lvtkGenericFiltering -lvtksys -lvtktiff -lvtkverdict -lvtkViews -lvtkVolumeRendering -lvtkWidgets -lvtkexoIIc -lvtklibxml2 -lvtkzlib -lvtkGeovis -lvtkGraphics -lvtkhdf5 -lvtkHybrid -lvtkImaging -lvtkInfovis -lvtkIO -lvtkjpeg \
		-lgsl -lgslcblas -lteem -fopenmp -lCGAL -latlas -llapack -lcblas \
		-frounding-math \
		-o vtkremovecompsize main.cpp
clean: 
	rm vtkremovecompsize;