# Kirkpatrick_Triangulation
The aim of this project is to create a point location structure for traingulations of point set with 
integer coordinates. The given triangulation has a triangular outer face, which is the first triangle in the
triangle list. The structure supports the following operations:

ploc_t* create_ploc(int* points, int* triangles, int n, int m) creates the point location structure.
Here n is the number of points, m is the number of triangles, point is a pointer to the n X 2 array 
of point coordinates and triangles is the pointer to m X 3 array of point numbers that make the 
triangles. So each triangle is given bythree elements of points array.

int query_ploc(ploc_t* pl, int x, int y) returns the number of triangle which contains the point (x, y) or -1 if 
the point is not contained in any triangle.

