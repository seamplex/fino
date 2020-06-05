// https://autofem.com/examples/determining_natural_frequencie.html
SetFactory("OpenCASCADE");

Box(1) = {0,0,0, 0.500, 0.050, 0.020};

Physical Surface("fixed") = {1};
Physical Volume("bulk") = {1};

Mesh.ElementOrder = 2;
Mesh.Algorithm = 8;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Mesh.CharacteristicLengthMax = 0.010;

