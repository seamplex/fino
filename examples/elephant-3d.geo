// an solid cylinder of radius r and length (height) z
r = 0.5;
l = 2;
lc = 1e-1;

SetFactory("OpenCASCADE");

Cylinder(1) = {0, 0, 0,  0, l, 0, r};

Physical Surface("bottom") = {3};
Physical Surface("top") = {2};
Physical Volume("bulk") = {1};

Mesh.CharacteristicLengthMax = lc;
Mesh.ElementOrder = 2;
Mesh.HighOrderOptimize = 2;
