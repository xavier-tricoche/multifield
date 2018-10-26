#!/bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/media/SDA2/scratch/sbarakat/Libraries/vtk-5.8.0/vtk-binaries/lib/vtk-5.8/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/media/SDA2/scratch/sbarakat/Libraries/valette-ACVD-6f93093/Common/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/media/SDA2/scratch/sbarakat/Libraries/valette-ACVD-6f93093/DiscreteRemeshing
export PATH=$PATH:/opt/local/bin/
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/mtrv
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/scalenrrd2unit
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/vtkadjustunitbbox
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/vtksimplifymesh
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/vtksimplifymeshquadric
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/vtkremovecompsize
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/vtkfillholes
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/vtkfindbbox
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/vtkassignnrrd
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/vtksmoothmesh
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/vtktag
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/vtksegment
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/vtkdistancefield
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/vtkmergestructures
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/nrrdblurinside
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/vtk2vtktri
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/vmtkremeshing
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/vtkcomputenormals
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/vtkcomponents
export PATH=$PATH:/media/SDA2/scratch/sbarakat/Projects/vtkremovemask

for f in /home/sbarakat/nabla/Samer/Datasets/Multifield/multi_mf/Ds_*.vtk
do
filename=$(basename $f)
filename=${filename%.*}
filenamelng=${#filename}
ts=${filename:3:4}
ts="/home/sbarakat/nabla/Samer/Datasets/Multifield/multi/multifield_1_$ts.nrrd"
vtkassignnrrd $f VALTBL $ts
echo $f
echo $ts


done
