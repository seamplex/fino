// use gmsh >= 2.16.0
SetFactory("OpenCASCADE");

h=H/2;

Point(1) = {0, 0, 0, lc};
Point(2) = {0, 0, h, lc};
Cylinder(4) = {0,0,0, 0,0,h, r};

// embed both centers in the surfaces
Point{1} In Surface{3};
Point{2} In Surface{2};

Mesh.CharacteristicLengthMin = 0.8*lc;
Mesh.CharacteristicLengthMax = 1.2*lc;

Physical Point("origin") = {1};
Physical Point("up") = {2};

Physical Surface("lat") = {1};
Physical Surface("sup") = {2};
Physical Surface("inf") = {3};

Physical Volume("bulk") = {4};
