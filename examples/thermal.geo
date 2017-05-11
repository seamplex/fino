Merge "prism.brep";
Mesh.CharacteristicLengthFactor = 0.25;
Mesh.Algorithm = 8;
//Mesh.RecombineAll = 1;

Physical Surface("left",   1) = {1};
Physical Surface("right",  2) = {2};
Physical Surface("front",  3) = {3};
Physical Surface("rear",   4) = {4};
Physical Surface("bottom", 5) = {5};
Physical Surface("top",    6) = {6};

Physical Volume("bulk") = {1};
