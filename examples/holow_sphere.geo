Merge "holow_sphere.stp";
Mesh.CharacteristicLengthFactor = 0.25;
Mesh.Algorithm = 8;
//Mesh.RecombineAll = 1;

Physical Surface("inner",  1) = {6};
Physical Surface("outer",  2) = {1};

Physical Volume("bulk") = {1};
