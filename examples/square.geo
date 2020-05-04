//
SetFactory("OpenCASCADE");
a = 1;
lc = a/20;

Rectangle(1) = {-a/2, -a/2, 0, +a, +a, 0};

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;
Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;
Mesh.Algorithm3D = 2;
Mesh.RecombineAll = 1;

Physical Line("right") = {2};
Physical Line("bottom") = {1};
Physical Line("top") = {3};
Physical Line("left") = {4};

Physical Surface("bulk") = {1};
