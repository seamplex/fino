//
SetFactory("OpenCASCADE");
l = 150;
d = 2;
Cylinder(1) = {0, 0, 0, 0, l, 0, d/2};

Mesh.OptimizeNetgen = 1;
Mesh.HighOrderOptimize = 2;
Mesh.Algorithm = 6;
Mesh.ElementOrder = 2;
Mesh.CharacteristicLengthMax = d/4;

Physical Volume ("bulk") = {1};
Physical Surface("fixed") = {3};
