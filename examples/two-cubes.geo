SetFactory("OpenCASCADE");
Merge "two-cubes.brep";
Mesh.CharacteristicLengthMax = 1;
Mesh.ElementOrder = 2;

Physical Surface("fixed", 1) = {6};
Physical Surface("load", 2) = {7};

Physical Volume("solid1") = {1};
Physical Volume("solid2") = {2};

