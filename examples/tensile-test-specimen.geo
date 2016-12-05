Merge "tensile-test-specimen.step";

Mesh.CharacteristicLengthMin = 1;
Mesh.CharacteristicLengthMax = 3;
Mesh.Algorithm = 8;

Physical Surface("left") = {1};
Physical Surface("right") = {7};

Physical Volume("bulk") = {1};
