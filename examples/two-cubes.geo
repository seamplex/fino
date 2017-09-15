//
SetFactory("OpenCASCADE");
a() = ShapeFromFile("two-cubes.brep");
p() = PointsOf{ Volume{a()}; };
Characteristic Length { p() } = 1.113;
Mesh.CharacteristicLengthMin = 0.5 * 0.711;
Mesh.CharacteristicLengthMax = 1.25 * 0.711;
Mesh.ElementOrder = 1;
Mesh.CharacteristicLengthFromPoints = 1;
Mesh.CharacteristicLengthFromCurvature = 1;
Mesh.CharacteristicLengthExtendFromBoundary = 1;
Mesh.Lloyd = 1;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.HighOrderOptimize = 1;
Mesh.SecondOrderLinear = 0;
Mesh.Algorithm = 8;
Mesh.Algorithm3D = 4;

Physical Surface("fixed", 1) = {6};
Physical Surface("load", 2) = {7};
Physical Surface(10) = {1,2,3,4,5,8,9,10,11};

Physical Volume("solid1") = {1};
Physical Volume("solid2") = {2};

