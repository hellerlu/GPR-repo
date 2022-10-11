A C++11 program to generate saturated RSA configurations of 2D polygons and 1-4 dimensional spheres.

Usage:
	Compile the program with flag -std=c++0x (which enables C++ 11 features), run the program, and input parameters as directed.

	The first three parameters affects the program's efficiency but not its correctness as long as they are positive integers. A reasonable choice would be:
		d-dimensional spheres: 3, 10^d, and 5*10^d
		2D polygons: 3, 50000, and 5000000.

	The fourth parameter is the random seed. The fifth parameter is a time limit in seconds. If the program has been running for longer than its limit, it will save its progress (warning: save files can be huge!). Calculation will resume when you run it again.

	One then choose a shape. For spheres, one enters its dimension and radius (assuming that the side length of the simulation box is one). For polygons, one enters a scaling factor, and then enter its vertices in polar coordinates.

Multithreading:
	To improve running speed, compile with OpenMP enabled, and then enter a negative first parameter. For example, for 2D spheres, one would enter -3 100 500 for the first three parameters, rather than 3 100 500.

Extending the program to other shapes/dimensions:

	To extend the program to other shapes, one need to inherit class Shape, and implement its pure virtual functions. Depending on the shape, the implementation may be nontrivial. 
	If the voxel space's dimension is higher than 4, one should additionally change line 19. If the voxel space's dimension is higher than 8, one should additionally change line 560: voxeltype should be an integer type with more than d bits.
