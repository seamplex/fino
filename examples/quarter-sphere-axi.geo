//
SetFactory("OpenCASCADE");
a = 1;
b = 2;
lc = (b-a)/20;

Point(1) = {0, 0, 0, lc};

Point(2) = {+a, 0, 0, lc};
Point(3) = {0, +a, 0, lc};

Point(4) = {+b, 0, 0, lc};
Point(5) = {0, +b, 0, lc};

Circle(8) = {2,1,3};
Circle(9) = {4,1,5};

Mesh.CharacteristicLengthMin = 0.25 * lc;
Mesh.CharacteristicLengthMax = 1.25 * lc;

/*
Mesh.ElementOrder = 2;
Mesh.CharacteristicLengthFromPoints = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.Lloyd = 1;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen =1;
Mesh.HighOrderOptimize = 1;
Mesh.SecondOrderLinear = 0;
Mesh.Algorithm = 2;
Mesh.Algorithm3D = 1;
*/

Mesh.Algorithm = 8;
Mesh.RecombineAll = 1;

//+
Line(10) = {3, 5};
//+
Line(11) = {4, 2};
//+
Line Loop(1) = {8, 10, -9, 11};
//+
Plane Surface(1) = {1};
//+
Physical Surface("bulk") = {1};
//+
Physical Line("inner") = {8};
//+
Physical Line("outer") = {9};
//+
Physical Line("vertical") = {10};
//+
Physical Line("horizontal") = {11};
