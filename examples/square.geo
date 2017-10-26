//
SetFactory("OpenCASCADE");
a = 1;
lc = a/20;

Point(1) = {-a, -a, 0, lc};
Point(2) = {+a, -a, 0, lc};
Point(3) = {+a, +a, 0, lc};
Point(4) = {-a, +a, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

// p() = PointsOf{ Volume{a()}; };
// Characteristic Length { p() } = 181.806;
Mesh.CharacteristicLengthMin = 0.25 * lc;
Mesh.CharacteristicLengthMax = 1.25 * lc;
Mesh.ElementOrder = 1;
Mesh.CharacteristicLengthFromPoints = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.Lloyd = 0;
Mesh.Optimize = 0;
Mesh.OptimizeNetgen = 0;
Mesh.HighOrderOptimize = 0;
Mesh.SecondOrderLinear = 0;
Mesh.Algorithm = 2;
Mesh.Algorithm3D = 1;

// Mesh.Algorithm = 8;
// Mesh.RecombineAll = 1;

//+
Physical Line("left") = {4};
//+
Physical Line("right") = {2};
//+
Physical Line("bottom") = {1};
//+
Physical Line("top") = {3};

Physical Surface("bulk") = {6};

