all:
	g++ -I/usr/include/atlas \
		-I/media/SDA2/scratch/sbarakat/Libraries/cuda/include \
		-I/media/SDA2/scratch/sbarakat/Libraries/NVIDIA_GPU_Computing_SDK/C/common/inc \
		-I/media/SDA2/scratch/sbarakat/Libraries/nvis/nvis/ \
		-I/media/SDA2/scratch/sbarakat/Libraries/vtk-5.8.0/vtk-binaries/include/vtk-5.8/ \
		-I/media/SDA2/scratch/sbarakat/Libraries/valette-ACVD-6f93093/Common \
		-I/media/SDA2/scratch/sbarakat/Libraries/valette-ACVD-6f93093/DiscreteRemeshing \
		-I/usr/include -I/usr/local/include -I/opt/local/include \
		-L/media/SDA2/scratch/sbarakat/Libraries/vtk-5.8.0/vtk-binaries/lib/vtk-5.8/ \
		-L/media/SDA2/scratch/sbarakat/Libraries/valette-ACVD-6f93093/Common \
		-L/media/SDA2/scratch/sbarakat/Libraries/valette-ACVD-6f93093/DiscreteRemeshing \
		-L/usr/lib -L/usr/local/lib -L/opt/local/lib \
		-lMapReduceMPI -lmpistubs -lvtkalglib -lvtkCharts -lvtkCommon -lvtkDICOMParser -lvtkmetaio -lvtkNetCDF -lvtkNetCDF_cxx -lvtkpng -lvtkproj4 -lvtkRendering -lvtksqlite -lvtkexpat -lvtkFiltering -lvtkfreetype -lvtkftgl -lvtkGenericFiltering -lvtksys -lvtktiff -lvtkverdict -lvtkViews -lvtkVolumeRendering -lvtkWidgets -lvtkexoIIc -lvtklibxml2 -lvtkzlib -lvtkGeovis -lvtkGraphics -lvtkhdf5 -lvtkHybrid -lvtkImaging -lvtkInfovis -lvtkIO -lvtkjpeg \
		-lgsl -lgslcblas -lteem -fopenmp -lCGAL -latlas -llapack -lcblas -lvtkSurface -lvtkDiscreteRemeshing \
		-frounding-math \
		-O2 -o vmtkremeshing main.cpp vtkvmtkDataSetItem.cxx vtkvmtkItem.cxx vtkvmtkMath.cxx vtkvmtkNeighborhood.cxx vtkvmtkPolyDataBoundaryExtractor.cxx vtkvmtkPolyDataManifoldExtendedNeighborhood.cxx \
		vtkvmtkPolyDataManifoldNeighborhood.cxx vtkvmtkPolyDataManifoldStencil.cxx vtkvmtkPolyDataSurfaceRemeshing.cxx vtkvmtkPolyDataUmbrellaStencil.cxx vtkvmtkStencil.cxx
clean: 
	rm vmtkremeshing;