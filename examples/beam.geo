//+
SetFactory("OpenCASCADE");
//+
a = 1;
l = 10;
Box(1) = {0, -0.5*a, -0.5*a, l,  a, a};

Characteristic Length {1:8} = a/4;
Mesh.CharacteristicLengthMin = a/4;
Mesh.CharacteristicLengthMax = a/4;

Mesh.Algorithm = 6;
Mesh.ElementOrder = 2;

//+
Physical Surface("fixed") = {1};
//+
Physical Line("uno") = {6};
//+
Physical Line("dos") = {8};
//+
Physical Surface("load") = {2};
//+
Physical Volume("bulk") = {1};
