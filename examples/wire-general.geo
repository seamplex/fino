//
SetFactory("OpenCASCADE");
l = 150;
d = 2;
Cylinder(1) = {0, 0, 0, 0, l, 0, d/2};

Mesh.Optimize = 1;
Mesh.Algorithm = 6;
Mesh.ElementOrder = 2;
Mesh.CharacteristicLengthMax = d/3;

Physical Volume ("bulk") = {1};
Physical Surface("fixed") = {3};
