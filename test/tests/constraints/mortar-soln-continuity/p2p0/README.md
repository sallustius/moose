In order to be able to run the Jupyter notebook in this directory, you must
first create the mesh generation binary (which is just `a.out` in the
notebook). The mesh binary is created using the `second-order.cpp` script that
employs the gmsh API. You must have a recent gmsh version (> 4.0) and gmsh must
be built with the cmake option `-DENABLE_BUILD_DYNAMIC`. You can then compile
`second-order.cpp` with the command:
`clang++ -std=c++11 -I/PATH/TO/GMSH/INCLUDE -L/PATH/TO/GMSH/LIB second-order.cpp -lgmsh`
to generate the meshing binary.
