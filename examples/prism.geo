//+
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 4, 1, 2};

Mesh.CharacteristicLengthMax = 0.2;
Mesh.Algorithm = 5;
Mesh.ElementOrder = 2;

Physical Surface("xy", 1) = {5};
Physical Surface("yz", 2) = {1};
Physical Surface("zx", 3) = {3};
Physical Surface(10) = {2,4,6};

Physical Volume("bulk") = {1};
