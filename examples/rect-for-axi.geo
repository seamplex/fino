//
SetFactory("OpenCASCADE");
a = 2;
b = 1;
lc = a/20;

Rectangle(1) = {0, 0, 0, a, b, 0};

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;
Mesh.ElementOrder = 1;
Mesh.Algorithm = 1;
Mesh.RecombineAll = 1;

Physical Line("right") = {2};
Physical Line("bottom") = {1};
Physical Line("top") = {3};
Physical Line("left") = {4};

Physical Surface("bulk") = {1};
