// Gmsh project created on Wed Nov 22 16:01:23 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {5, 0, 0, 1.0};
Point(2) = {0, 5, 0, 1.0};
Point(3) = {0., 0., 0, 1.0};
Point(4) = {2., 0., 0, 1.0};
Point(5) = {0., 2., 0, 1.0};
//+
Circle(1) = {1, 3, 2};
Circle(2) = {4, 3, 5};
//+
Line(3) = {1, 4};
Line(4) = {2, 5};
//+
Curve Loop(1) = {3, 2, -4, -1};
//+
Plane Surface(1) = {1};
//+
Physical Surface("domain", 1) = {1};
//+
Physical Curve("bc", 2) = {3, 1, 4, 2};

