//
SetFactory("OpenCASCADE");
a = 1;
b = 2;
h = 1;
lc = (b-a)/20;

Rectangle(1) = {a, 0, 0, (b-a), h};

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
Physical Surface("bulk") = {1};
//+
Physical Line("inner") = {4};
//+
Physical Line("outer") = {2};
//+
Physical Line("infinite") = {1, 3};
